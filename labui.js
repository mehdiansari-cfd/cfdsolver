/* ============================================================
   labui.js  —  shared UI layer for the Euler Lab pages
   ------------------------------------------------------------
   UI only: knows nothing about fluxes or the Euler equations. The
   pages hand it plain callbacks (init / advance / record) and arrays.

   Plain JavaScript with React.createElement (no JSX), so it loads
   with a normal <script src> — also from file:// — after React.

   Exposes window.LabUI:
     Player            playback engine: computes in time-budgeted
                       slices per animation frame (the page never
                       freezes), records frames, and plays them back
                       at a steady rate with a scrubbable timeline
     drawPlot          canvas line plot (numerical + dashed reference)
     colormap, CMAP    warm cream → gold → rust → near-black LUT
     fmt, niceTicks    number formatting / axis ticks
     useCanvas         HiDPI canvas + ResizeObserver redraw hook
     useSmoothRange    axis range that grows at once, shrinks gently
     components        NumField, SelectField, Section, Segmented,
                       PlotCard, Transport, TopBar, HelpDrawer,
                       Colorbar, Probe
   ============================================================ */
(function () {
  "use strict";
  const h = React.createElement;
  const { useState, useEffect, useRef, useCallback, useMemo } = React;

  // ---------- Theme (mirrors lab.css tokens; canvases can't read CSS vars cheaply) ----------
  const THEME = {
    plotBg:   "#101317",
    grid:     "#1f242b",
    axis:     "#3a414c",
    tick:     "#7d8794",
    series:   "#d9653f",   // numerical solution — one role, one colour, every plot
    data:     "#4a8fe0",   // measured / imported data points (validated vs series on plotBg)
    ref:      "#98a0ab",   // analytical reference — neutral, dashed
    cursor:   "#d7a33c",
    font:     "10px ui-monospace, 'SF Mono', SFMono-Regular, Menlo, Consolas, monospace"
  };

  // ---------- Number formatting ----------
  function fmt(v, digits) {
    if (!Number.isFinite(v)) return "—";
    const a = Math.abs(v);
    if (a === 0) return "0";
    if (a < 1e-3 || a >= 1e4) return v.toExponential(digits == null ? 2 : digits);
    return v.toFixed(digits == null ? (a < 1 ? 4 : a < 100 ? 3 : 1) : digits);
  }

  // "Nice" tick positions covering [lo, hi] (Heckbert's 1-2-5 rule).
  function niceTicks(lo, hi, count) {
    const span = hi - lo;
    if (!(span > 0)) return [lo];
    const raw = span / Math.max(1, count);
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const r = raw / mag;
    const step = (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * mag;
    const out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) {
      out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
    }
    return out;
  }

  // ---------- Warm perceptual colormap as a 256-entry LUT ----------
  const STOPS = [[247, 242, 227], [218, 178, 98], [173, 88, 62], [50, 34, 30]];
  const CMAP = new Uint8ClampedArray(256 * 3);
  for (let k = 0; k < 256; k++) {
    const seg = (k / 255) * 3, i = Math.min(2, Math.floor(seg)), f = seg - i;
    const a = STOPS[i], b = STOPS[i + 1];
    CMAP[3 * k]     = a[0] + (b[0] - a[0]) * f;
    CMAP[3 * k + 1] = a[1] + (b[1] - a[1]) * f;
    CMAP[3 * k + 2] = a[2] + (b[2] - a[2]) * f;
  }
  function colormap(norm) {
    const k = 3 * Math.max(0, Math.min(255, Math.round(norm * 255)));
    return "rgb(" + CMAP[k] + "," + CMAP[k + 1] + "," + CMAP[k + 2] + ")";
  }

  // ============================================================
  // Player — compute + playback engine
  // ------------------------------------------------------------
  // cfg = {
  //   init()               -> sim state at t = 0
  //   advance(sim, tTarget)-> sim state advanced towards tTarget
  //   stableDt(sim)        -> the solver's CFL time step at this state
  //   record(sim)          -> compact snapshot (copied arrays)
  //   check(sim)           -> null, or an error string (NaN, negative p, ...)
  //   tEnd, maxFrames, playSeconds (wall time of a full run at 1x)
  //   onChange()           -> called once per animation frame while active
  // }
  // A frame is recorded each time the solution passes the next multiple of
  // tEnd / maxFrames. Steps are never clipped to land on a frame time (that
  // would add a sliver step per frame); frames keep their actual t and are
  // looked up by binary search.
  // ============================================================
  const SLICE_MS_PLAYING = 8;    // compute budget per animation frame while playing
  const SLICE_MS_SOLVING = 10;   // ... while solving flat out (leaves room to paint at 60 fps)
  const CHUNK_MS = 4;            // target wall time of one advance() call
  const CHUNK_SLACK = 0.02;      // aim a chunk 2% of a step short of k whole steps,
                                 // so its last step is trimmed, not followed by a sliver

  function Player(cfg) {
    this.cfg = cfg;
    this.tEnd = cfg.tEnd;
    this.maxFrames = Math.max(2, cfg.maxFrames | 0);
    this.dtFrame = cfg.tEnd / this.maxFrames;
    this.playSeconds = cfg.playSeconds || 8;
    this.speed = 1;
    this.frames = [];
    this.sim = null;
    this.cursorT = 0;
    this.playing = false;
    this.follow = false;      // cursor tracks the compute frontier (solve mode)
    this.done = false;
    this.error = null;
    this.computeMs = 0;
    this._stepsPerChunk = 1;
    this._nextFrameT = 0;
    this._raf = null;
    this._last = 0;
    this._dead = false;
    this._loop = this._loop.bind(this);
  }

  Player.prototype.start = function (mode) {
    this.frames = []; this.done = false; this.error = null; this.computeMs = 0;
    this._stepsPerChunk = 1; this._nextFrameT = this.dtFrame;
    try {
      this.sim = this.cfg.init();
      this.frames.push({ t: this.sim.t, step: this.sim.step, snap: this.cfg.record(this.sim) });
    } catch (e) {
      this.error = "Initialisation failed: " + e.message;
    }
    this.cursorT = 0;
    this.playing = mode === "play";
    this.follow = mode !== "play";
    this._kick();
  };

  Player.prototype.frontierT = function () {
    return this.frames.length ? this.frames[this.frames.length - 1].t : 0;
  };
  // Index of the recorded frame nearest to time t.
  Player.prototype.indexAt = function (t) {
    const F = this.frames;
    let lo = 0, hi = F.length - 1;
    while (lo < hi) { const mid = (lo + hi + 1) >> 1; if (F[mid].t <= t) lo = mid; else hi = mid - 1; }
    return lo + 1 < F.length && F[lo + 1].t - t < t - F[lo].t ? lo + 1 : lo;
  };
  Player.prototype.frameAt = function (t) {
    return this.frames.length ? this.frames[this.indexAt(t)] : null;
  };
  Player.prototype.current = function () { return this.frameAt(this.cursorT); };
  Player.prototype.progress = function () { return this.frontierT() / this.tEnd; };

  Player.prototype.play = function () {
    if (this.error && !this.frames.length) return;
    if (this.cursorT >= this.tEnd - 0.5 * this.dtFrame) this.cursorT = 0;   // replay from the top
    this.playing = true; this.follow = false; this._kick();
  };
  Player.prototype.pause = function () { this.playing = false; this.follow = false; this._notify(); };
  Player.prototype.toggle = function () { this.playing ? this.pause() : this.play(); };
  Player.prototype.seek = function (t) {
    this.cursorT = Math.max(0, Math.min(this.frontierT(), t));
    this.follow = false; this._notify();
  };
  Player.prototype.stepFrame = function (dir) {
    if (!this.frames.length) return;
    this.playing = false;
    const k = Math.max(0, Math.min(this.frames.length - 1, this.indexAt(this.cursorT) + dir));
    this.seek(this.frames[k].t);
  };
  Player.prototype.setSpeed = function (s) { this.speed = s; this._notify(); };
  Player.prototype.destroy = function () {
    this._dead = true;
    if (this._raf) cancelAnimationFrame(this._raf);
    this._raf = null;
  };
  Player.prototype.status = function () {
    if (this.error) return "error";
    if (this.playing && !this.done && this.cursorT >= this.frontierT() - 1e-12) return "buffering";
    if (this.playing) return "playing";
    if (!this.done) return "computing";
    return "ready";
  };

  Player.prototype._notify = function () { if (!this._dead && this.cfg.onChange) this.cfg.onChange(); };
  Player.prototype._kick = function () {
    if (this._dead) return;
    if (!this._raf) { this._last = performance.now(); this._raf = requestAnimationFrame(this._loop); }
    this._notify();
  };

  // Advance the simulation by roughly CHUNK_MS of work: k whole CFL steps,
  // with k re-tuned from the measured cost per step.
  Player.prototype._computeChunk = function () {
    const cfg = this.cfg, sim = this.sim;
    const dt = cfg.stableDt(sim);
    if (!(dt > 0 && dt < Infinity)) {
      this.error = "The solver stalled (non-finite time step) at t = " + fmt(sim.t) + ".";
      return;
    }
    const target = Math.min(this.tEnd, sim.t + (this._stepsPerChunk - CHUNK_SLACK) * dt);
    const s0 = sim.step, w0 = performance.now();
    const next = cfg.advance(sim, target);
    const wall = performance.now() - w0;
    this.computeMs += wall;
    const steps = next.step - s0;
    if (!Number.isFinite(next.t) || steps === 0) {
      this.error = "The solver stalled at t = " + fmt(sim.t) + ".";
      return;
    }
    this.sim = next;
    const bad = cfg.check ? cfg.check(next) : null;
    if (bad) { this.error = bad + " at t = " + fmt(next.t) + "."; return; }
    this._stepsPerChunk = Math.max(1, Math.floor(CHUNK_MS * steps / Math.max(wall, 1e-3)));

    const finished = next.t >= this.tEnd * (1 - 1e-12);
    if (finished || next.t >= this._nextFrameT) {
      this.frames.push({ t: next.t, step: next.step, snap: cfg.record(next) });
      this._nextFrameT = (Math.floor(next.t / this.dtFrame + 1e-9) + 1) * this.dtFrame;
      if (finished) this.done = true;
    }
  };

  Player.prototype._loop = function (now) {
    this._raf = null;
    if (this._dead) return;
    const dtWall = Math.min(0.1, (now - this._last) / 1000);   // clamp after tab switches
    this._last = now;

    if (!this.done && !this.error) {
      const budget = this.playing ? SLICE_MS_PLAYING : SLICE_MS_SOLVING;
      const w0 = performance.now();
      try {
        while (!this.done && !this.error && performance.now() - w0 < budget) this._computeChunk();
      } catch (e) {
        this.error = "Solver error: " + e.message;
      }
    }

    if (this.playing) {
      this.cursorT += dtWall * this.speed * this.tEnd / this.playSeconds;
      const front = this.frontierT();
      if (this.cursorT >= front) {
        this.cursorT = front;
        if (this.done || this.error) this.playing = false;
      }
    } else if (this.follow) {
      this.cursorT = this.frontierT();
    }

    this._notify();
    if (this.playing || (!this.done && !this.error)) this._raf = requestAnimationFrame(this._loop);
  };

  // React binding: returns [player, start(cfg, mode)]. Each Player update
  // bumps one counter, so the page re-renders once per animation frame.
  function usePlayer() {
    const [, setTick] = useState(0);
    const ref = useRef(null);
    const start = useCallback((cfg, mode) => {
      if (ref.current) ref.current.destroy();
      const p = new Player(Object.assign({}, cfg, { onChange: () => setTick(n => n + 1) }));
      ref.current = p;
      p.start(mode);
      return p;
    }, []);
    useEffect(() => () => { if (ref.current) ref.current.destroy(); }, []);
    return [ref.current, start];
  }

  // ============================================================
  // Canvas helpers
  // ============================================================

  // Size a canvas's backing store to its CSS box × devicePixelRatio.
  function fitCanvas(cv) {
    const dpr = window.devicePixelRatio || 1;
    const w = cv.clientWidth, hh = cv.clientHeight;
    const bw = Math.max(1, Math.round(w * dpr)), bh = Math.max(1, Math.round(hh * dpr));
    if (cv.width !== bw || cv.height !== bh) { cv.width = bw; cv.height = bh; }
    const ctx = cv.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx, W: w, H: hh };
  }

  // Calls draw(ctx, W, H) after every render and whenever the canvas resizes.
  function useCanvas(draw, deps) {
    const ref = useRef(null);
    const drawRef = useRef(draw);
    drawRef.current = draw;
    useEffect(() => {
      const cv = ref.current; if (!cv) return;
      const f = fitCanvas(cv); drawRef.current(f.ctx, f.W, f.H);
    }, deps);
    useEffect(() => {
      const cv = ref.current; if (!cv || typeof ResizeObserver === "undefined") return;
      const ro = new ResizeObserver(() => { const f = fitCanvas(cv); drawRef.current(f.ctx, f.W, f.H); });
      ro.observe(cv);
      return () => ro.disconnect();
    }, []);
    return ref;
  }

  // Axis range with hysteresis: expands immediately so data never clips, and
  // relaxes inward a little each frame so autoscaling doesn't jitter.
  const SHRINK_RATE = 0.12;
  function useSmoothRange(lo, hi, resetKey) {
    const st = useRef({ key: null, lo: 0, hi: 1 });
    const s = st.current;
    if (!Number.isFinite(lo) || !Number.isFinite(hi)) return [s.lo, s.hi];
    if (s.key !== resetKey) { s.key = resetKey; s.lo = lo; s.hi = hi; return [lo, hi]; }
    s.lo = lo < s.lo ? lo : s.lo + (lo - s.lo) * SHRINK_RATE;
    s.hi = hi > s.hi ? hi : s.hi + (hi - s.hi) * SHRINK_RATE;
    return [s.lo, s.hi];
  }

  function finiteRange(arrs) {
    let lo = Infinity, hi = -Infinity;
    for (const a of arrs) {
      if (!a) continue;
      for (let i = 0; i < a.length; i++) {
        const v = a[i];
        if (Number.isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; }
      }
    }
    return [lo, hi];
  }

  // ---------- Line plot ----------
  const MIN_REL_SPAN = 0.02;
  // o = { x, y, ref, lo, hi, logY, hoverX, xLabel }
  function drawPlot(ctx, W, H, o) {
    const pad = { l: 46, r: 10, t: 8, b: 22 };
    const pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;
    ctx.fillStyle = THEME.plotBg; ctx.fillRect(0, 0, W, H);
    if (!o.x || o.x.length < 2 || pw <= 0 || ph <= 0) return null;

    const tr = o.logY ? (v => Math.log10(Math.max(v, 1e-300))) : (v => v);
    const pts = o.points || null;
    let lo = tr(o.lo), hi = tr(o.hi);
    if (!(hi > lo)) { const c = Number.isFinite(lo) ? lo : 0; lo = c - 0.5; hi = c + 0.5; }
    // A nearly constant series (e.g. steady mass flow) spans at least
    // MIN_REL_SPAN of its magnitude, so round-off isn't blown up into noise.
    if (!o.logY) {
      const minSpan = MIN_REL_SPAN * Math.max(Math.abs(lo), Math.abs(hi));
      if (hi - lo < minSpan) { const c = 0.5 * (lo + hi); lo = c - 0.5 * minSpan; hi = c + 0.5 * minSpan; }
    }
    const yPad = (hi - lo) * 0.06; lo -= yPad; hi += yPad;
    const xMin = o.x[0], xMax = o.x[o.x.length - 1];
    const X = v => pad.l + ((v - xMin) / (xMax - xMin)) * pw;
    const Y = v => pad.t + (1 - (tr(v) - lo) / (hi - lo)) * ph;

    // grid + ticks
    ctx.font = THEME.font; ctx.lineWidth = 1;
    ctx.strokeStyle = THEME.grid; ctx.fillStyle = THEME.tick;
    ctx.textAlign = "right"; ctx.textBaseline = "middle";
    const yt = o.logY
      ? niceTicks(Math.ceil(lo), Math.floor(hi), 4).filter(v => v === Math.round(v))
      : niceTicks(lo, hi, 4);
    for (const v of yt) {
      const py = Math.round(pad.t + (1 - (v - lo) / (hi - lo)) * ph) + 0.5;
      ctx.beginPath(); ctx.moveTo(pad.l, py); ctx.lineTo(pad.l + pw, py); ctx.stroke();
      ctx.fillText(o.logY ? "1e" + v : fmt(v, Math.abs(v) >= 1e4 || (Math.abs(v) < 1e-3 && v !== 0) ? 0 : undefined)
                   .replace(/(\.\d*?)0+$/, "$1").replace(/\.$/, ""), pad.l - 6, py);
    }
    ctx.textAlign = "center"; ctx.textBaseline = "top";
    for (const v of niceTicks(xMin, xMax, 5)) {
      const px = Math.round(X(v)) + 0.5;
      ctx.beginPath(); ctx.moveTo(px, pad.t); ctx.lineTo(px, pad.t + ph); ctx.stroke();
      ctx.fillText(String(+v.toFixed(6)), px, pad.t + ph + 6);
    }
    ctx.strokeStyle = THEME.axis;
    ctx.strokeRect(pad.l + 0.5, pad.t + 0.5, pw - 1, ph - 1);

    ctx.save();
    ctx.beginPath(); ctx.rect(pad.l, pad.t, pw, ph); ctx.clip();
    ctx.lineJoin = "round"; ctx.lineCap = "round";
    const line = (arr, color, width, dash) => {
      ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash || []);
      ctx.beginPath();
      let pen = false;
      for (let i = 0; i < o.x.length; i++) {
        const v = arr[i];
        if (!Number.isFinite(v)) { pen = false; continue; }
        const px = X(o.x[i]), py = Y(v);
        if (pen) ctx.lineTo(px, py); else { ctx.moveTo(px, py); pen = true; }
      }
      ctx.stroke(); ctx.setLineDash([]);
    };
    if (o.ref) line(o.ref, THEME.ref, 1.25, [4, 3]);
    line(o.y, THEME.series, 2);
    if (pts) {                                   // measured data: hollow 8 px markers with a surface ring
      for (const q of pts) {
        if (!Number.isFinite(q.y) || q.x < xMin || q.x > xMax) continue;
        const px = X(q.x), py = Y(q.y);
        ctx.beginPath(); ctx.arc(px, py, 4, 0, Math.PI * 2);
        ctx.lineWidth = 4; ctx.strokeStyle = THEME.plotBg; ctx.stroke();
        ctx.lineWidth = 1.75; ctx.strokeStyle = THEME.data; ctx.stroke();
      }
    }

    // linked hover cursor
    let hit = null;
    if (o.hoverX != null && o.hoverX >= xMin && o.hoverX <= xMax) {
      const i = Math.max(0, Math.min(o.x.length - 1,
                Math.round(((o.hoverX - xMin) / (xMax - xMin)) * (o.x.length - 1))));
      const px = Math.round(X(o.x[i])) + 0.5;
      ctx.strokeStyle = THEME.cursor; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(px, pad.t); ctx.lineTo(px, pad.t + ph); ctx.stroke();
      if (Number.isFinite(o.y[i])) {
        ctx.fillStyle = THEME.series; ctx.strokeStyle = THEME.plotBg; ctx.lineWidth = 2;
        ctx.beginPath(); ctx.arc(px, Y(o.y[i]), 4, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
      }
      hit = i;
    }
    ctx.restore();
    return { pad, pw, xMin, xMax, hit };
  }

  // ============================================================
  // Components
  // ============================================================

  // Number input that keeps the user's text while typing ("0.", "-", "1e")
  // and only reports finite numbers upward.
  function NumField({ label, value, onChange, step, title }) {
    const [text, setText] = useState(String(value));
    const last = useRef(value);
    if (value !== last.current) { last.current = value; if (parseFloat(text) !== value) setText(String(value)); }
    return h("label", { className: "field", title },
      h("span", { className: "field-label" }, label),
      h("input", {
        type: "number", value: text, step: step || 0.1,
        onChange: e => {
          setText(e.target.value);
          const v = parseFloat(e.target.value);
          if (Number.isFinite(v)) { last.current = v; onChange(v); }
        }
      }));
  }

  function SelectField({ label, value, onChange, options }) {
    return h("label", { className: "field" },
      h("span", { className: "field-label" }, label),
      h("select", { value, onChange: e => onChange(e.target.value) },
        options.map(o => {
          const v = typeof o === "string" ? o : o.value;
          return h("option", { key: v, value: v }, typeof o === "string" ? o : o.label);
        })));
  }

  function Section({ title, children, open, aside }) {
    return h("details", { className: "rail-section", open: open !== false },
      h("summary", null, h("span", null, title), aside ? h("span", { className: "aside" }, aside) : null),
      h("div", { className: "rail-body" }, children));
  }

  function Segmented({ value, onChange, options, label }) {
    return h("div", { className: "segmented", role: "radiogroup", "aria-label": label },
      options.map(o => {
        const v = typeof o === "string" ? o : o.value;
        return h("button", {
          key: v, type: "button", role: "radio", "aria-checked": v === value,
          className: v === value ? "on" : "", onClick: () => onChange(v)
        }, typeof o === "string" ? o : o.label);
      }));
  }

  function Toggle({ checked, onChange, label }) {
    return h("label", { className: "toggle" },
      h("input", { type: "checkbox", checked, onChange: e => onChange(e.target.checked) }),
      h("span", null, label));
  }

  // Canvas line plot with a linked hover cursor and a live readout.
  // Hovering reports a domain x through onHover so every plot (and the
  // field view) can show the same probe position.
  // yScale converts SI values for display (e.g. 1e-3 for Pa → kPa); points are
  // measured data [{x, y}] in display units, drawn as hollow markers.
  function PlotCard({ title, x, y: yRaw, refY: refRaw, hoverX, onHover, resetKey, logY, refLabel, wide,
                      yScale, points, pointsLabel }) {
    const sc = yScale || 1;
    const y = useMemo(() => (sc === 1 || !yRaw ? yRaw : Array.from(yRaw, v => v * sc)), [yRaw, sc]);
    const refArr = useMemo(() => (sc === 1 || !refRaw ? refRaw : Array.from(refRaw, v => v * sc)), [refRaw, sc]);
    const ptsY = points ? points.map(q => q.y) : null;
    const [rlo, rhi] = finiteRange(logY ? [y] : [y, refArr, ptsY]);
    const [lo, hi] = useSmoothRange(rlo, rhi, resetKey);
    const geo = useRef(null);
    const cvRef = useCanvas((ctx, W, H) => {
      geo.current = drawPlot(ctx, W, H, { x, y, ref: refArr, lo, hi, logY, hoverX, points });
    }, [x, y, refArr, lo, hi, logY, hoverX, points]);

    const onMove = e => {
      const g = geo.current; if (!g || !onHover) return;
      const r = e.currentTarget.getBoundingClientRect();
      const f = (e.clientX - r.left - g.pad.l) / g.pw;
      onHover(f < 0 || f > 1 ? null : g.xMin + f * (g.xMax - g.xMin));
    };
    let readout = null;
    if (hoverX != null && x && x.length) {
      const i = Math.max(0, Math.min(x.length - 1,
                Math.round(((hoverX - x[0]) / (x[x.length - 1] - x[0])) * (x.length - 1))));
      readout = h("span", { className: "readout" },
        h("b", null, fmt(y[i])),
        refArr && Number.isFinite(refArr[i]) ? h("span", { className: "muted" }, " · ref " + fmt(refArr[i])) : null);
    }
    return h("figure", { className: "plot-card" + (wide ? " wide" : "") },
      h("figcaption", null,
        h("span", { className: "plot-title" }, title),
        readout,
        refArr || points ? h("span", { className: "legend" },
          h("i", { className: "sw-num" }), "numerical",
          refArr ? h("i", { className: "sw-ref" }) : null, refArr ? (refLabel || "exact") : null,
          points ? h("i", { className: "sw-data" }) : null, points ? (pointsLabel || "data") : null) : null),
      h("canvas", { ref: cvRef, className: "plot-canvas", onMouseMove: onMove,
                    onMouseLeave: () => onHover && onHover(null) }));
  }

  // ---------- Transport bar: solve / play / scrub / speed ----------
  const SPEEDS = [0.25, 0.5, 1, 2, 4];
  const STATUS_TEXT = {
    idle: "idle", computing: "computing", playing: "playing",
    buffering: "buffering", ready: "ready", error: "error", stale: "params changed"
  };
  // tUnit = { scale, label } shows time in display units (e.g. ms for SI runs).
  function Transport({ player, stale, onSolve, onPlay, info, tUnit }) {
    const ts = tUnit ? tUnit.scale : 1, tl = tUnit ? " " + tUnit.label : "";
    const p = player;
    const status = stale ? "stale" : (p ? p.status() : "idle");
    const tEnd = p ? p.tEnd : 1;
    const cursor = p ? p.cursorT : 0;
    const buffered = p ? p.progress() : 0;
    const playing = p && p.playing;
    return h("div", { className: "transport" },
      h("button", { className: "btn btn-primary", onClick: onSolve, title: "Solve to t_end (Enter)" }, "Solve"),
      h("button", {
        className: "btn icon", onClick: () => (p && !stale ? p.toggle() : onPlay()),
        title: playing ? "Pause (Space)" : "Play (Space)", "aria-label": playing ? "Pause" : "Play"
      }, playing ? "❚❚" : "▶"),
      h("button", { className: "btn icon", disabled: !p, onClick: () => p && p.stepFrame(-1),
                    title: "Previous frame (←)", "aria-label": "Previous frame" }, "‹"),
      h("button", { className: "btn icon", disabled: !p, onClick: () => p && p.stepFrame(1),
                    title: "Next frame (→)", "aria-label": "Next frame" }, "›"),
      h("div", { className: "scrub" },
        h("div", { className: "scrub-buffer", style: { width: (100 * buffered).toFixed(2) + "%" } }),
        h("input", {
          type: "range", min: 0, max: 1000, step: 1, "aria-label": "Time",
          value: Math.round(1000 * cursor / tEnd),
          disabled: !p,
          onChange: e => p && (p.pause(), p.seek(tEnd * e.target.value / 1000))
        })),
      h("select", { className: "speed", value: p ? p.speed : 1, title: "Playback speed",
                    onChange: e => p && p.setSpeed(parseFloat(e.target.value)) },
        SPEEDS.map(s => h("option", { key: s, value: s }, s + "×"))),
      h("div", { className: "clock" },
        h("span", null, h("span", { className: "muted" }, "t "), fmt(cursor * ts, 4)),
        h("span", { className: "muted" }, " / " + fmt(tEnd * ts, 3) + tl)),
      info ? h("div", { className: "meta" }, info) : null,
      h("span", { className: "chip chip-" + status }, STATUS_TEXT[status] || status));
  }

  // ---------- Top bar with page tabs ----------
  function TopBar({ active, scheme, onHelp }) {
    const tabs = [
      { id: "shocktube", label: "Shock tube", href: "index.html#shocktube" },
      { id: "nozzle",    label: "Nozzle",     href: "index.html#nozzle" },
      { id: "facility",  label: "Facility",   href: "index.html#facility" },
      { id: "verify",    label: "Verify",     href: "index.html#verify" },
      { id: "2d",        label: "2D field",   href: "index2d.html" }
    ];
    return h("header", { className: "topbar" },
      h("div", { className: "brand" }, h("span", { className: "brand-mark" }, "∂U/∂t"), "Euler Lab"),
      h("nav", { className: "tabs" },
        tabs.map(t => h("a", { key: t.id, href: t.href, className: "tab" + (t.id === active ? " on" : "") }, t.label))),
      h("div", { className: "topbar-right" },
        scheme ? h("span", { className: "scheme" }, scheme) : null,
        h("button", { className: "btn icon", onClick: onHelp, title: "Help & shortcuts (?)", "aria-label": "Help" }, "?")));
  }

  function HelpDrawer({ open, onClose, children }) {
    if (!open) return null;
    return h("div", { className: "drawer-backdrop", onClick: onClose },
      h("aside", { className: "drawer", onClick: e => e.stopPropagation(), role: "dialog", "aria-label": "Help" },
        h("button", { className: "btn icon drawer-close", onClick: onClose, "aria-label": "Close" }, "×"),
        children,
        h("h3", null, "Shortcuts"),
        h("dl", { className: "keys" },
          h("dt", null, "Enter"), h("dd", null, "solve to t_end"),
          h("dt", null, "Space"), h("dd", null, "play / pause"),
          h("dt", null, "← →"),   h("dd", null, "step one frame"),
          h("dt", null, "Home End"), h("dd", null, "jump to start / end"),
          h("dt", null, "?"),     h("dd", null, "toggle this panel"))));
  }

  // Global keyboard shortcuts (ignored while typing in a field).
  function useShortcuts(handlers) {
    const ref = useRef(handlers); ref.current = handlers;
    useEffect(() => {
      const onKey = e => {
        const tag = (e.target.tagName || "").toLowerCase();
        if (tag === "input" || tag === "select" || tag === "textarea" || e.metaKey || e.ctrlKey || e.altKey) return;
        const fn = ref.current[e.key];
        if (fn) { e.preventDefault(); fn(e); }
      };
      window.addEventListener("keydown", onKey);
      return () => window.removeEventListener("keydown", onKey);
    }, []);
  }

  // ---------- Colour bar for a field view ----------
  function Colorbar({ lo, hi, label }) {
    return h("div", { className: "colorbar" },
      h("span", { className: "muted" }, label),
      h("span", null, fmt(lo, 3)),
      h("i", { className: "colorbar-ramp" }),
      h("span", null, fmt(hi, 3)));
  }

  // ---------- Probe readout (values at the hovered position) ----------
  function Probe({ items, where }) {
    if (!items) return h("div", { className: "probe muted" }, "hover a plot or the field to probe values");
    return h("div", { className: "probe" },
      h("span", { className: "muted" }, where),
      items.map(([k, v]) => h("span", { key: k }, h("span", { className: "muted" }, k + " "), fmt(v))));
  }

  // ---------- Files: download text, pick a local text file ----------
  function downloadText(name, text, mime) {
    const url = URL.createObjectURL(new Blob([text], { type: mime || "text/plain" }));
    const a = document.createElement("a");
    a.href = url; a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }
  function pickTextFile(accept) {
    return new Promise((resolve, reject) => {
      const inp = document.createElement("input");
      inp.type = "file"; inp.accept = accept || "";
      inp.onchange = () => {
        const f = inp.files && inp.files[0];
        if (!f) return reject(new Error("No file chosen"));
        f.text().then(text => resolve({ name: f.name, text }), reject);
      };
      inp.click();
    });
  }
  // Rows of a CSV / whitespace table whose cells are all numbers (headers skipped).
  function numericTable(text) {
    const rows = [];
    for (const line of text.split(/\r?\n/)) {
      const cells = line.trim().split(/[\s,;]+/).filter(Boolean).map(Number);
      if (cells.length >= 2 && cells.every(Number.isFinite)) rows.push(cells);
    }
    return rows;
  }

  // ---------- Drive a generator in time slices between animation frames ----------
  // onYield(value) after each slice; onDone(result) / onError(err). Returns cancel().
  function runSliced(gen, onYield, onDone, onError, sliceMs) {
    let dead = false, raf = null;
    const budget = sliceMs || 12;
    const tick = () => {
      raf = null;
      if (dead) return;
      const t0 = performance.now();
      try {
        let last;
        while (performance.now() - t0 < budget) {
          const r = gen.next();
          if (r.done) { onDone(r.value); return; }
          last = r.value;
        }
        if (last !== undefined) onYield(last);
      } catch (e) { onError(e); return; }
      raf = requestAnimationFrame(tick);
    };
    raf = requestAnimationFrame(tick);
    return () => { dead = true; if (raf) cancelAnimationFrame(raf); };
  }

  window.LabUI = {
    THEME, fmt, niceTicks, CMAP, colormap,
    Player, usePlayer, fitCanvas, useCanvas, useSmoothRange, finiteRange, drawPlot, useShortcuts,
    NumField, SelectField, Section, Segmented, Toggle, PlotCard, Transport, TopBar, HelpDrawer, Colorbar, Probe,
    downloadText, pickTextFile, numericTable, runSliced
  };
})();
