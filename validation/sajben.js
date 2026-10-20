#!/usr/bin/env node
/* ============================================================
   validation/sajben.js — quasi-1D comparison with the Sajben diffuser
   (see validation/sajben/README.md). Node built-ins only.

   Geometry: upper / lower wall coordinates (inches), A(x) ∝ h(x) per
   unit depth. Data: wall static pressure p/p0 against x/H* (H* = throat
   height, x from the throat) on the bottom and top walls.
   Outputs, per case, in results/sajben/: computed and measured CSV, an
   SVG plot, and report.json (RMS Δ(p/p0) per wall, shock positions).
   Skips cleanly (exit 0) when the experimental files are absent.
   ============================================================ */
"use strict";
const fs = require("fs");
const path = require("path");
const Z = require("../nozzle");
const C = require("../case");

const HERE = path.join(__dirname, "sajben");
const OUT = path.join(__dirname, "..", "results", "sajben");

// Rows of a whitespace/comma separated file whose cells are all numbers
// (headers, VARIABLES / ZONE lines and # comments are skipped).
function numericRows(file) {
  const rows = [];
  for (const line of fs.readFileSync(file, "utf8").split(/\r?\n/)) {
    if (/^\s*#/.test(line)) continue;
    const cells = line.trim().split(/[\s,;]+/).filter(Boolean);
    if (!cells.length) continue;
    const nums = cells.map(Number);
    if (nums.every(Number.isFinite)) rows.push(nums);
  }
  return rows;
}

function interp(xs, ys, x) {
  if (x <= xs[0]) return ys[0];
  if (x >= xs[xs.length - 1]) return ys[ys.length - 1];
  let lo = 0, hi = xs.length - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
  return ys[lo] + (x - xs[lo]) / (xs[hi] - xs[lo]) * (ys[hi] - ys[lo]);
}

function wall(file, cols, scale) {
  const rows = numericRows(file).map(r => [r[cols.x] * scale, r[cols.y] * scale]).sort((a, b) => a[0] - b[0]);
  if (rows.length < 2) throw new Error(file + ": fewer than two numeric (x, y) rows — check geometry.columns");
  return { x: rows.map(r => r[0]), y: rows.map(r => r[1]) };
}

// Shock = largest pressure rise between neighbouring points, x at its midpoint.
function shockAt(X, P) {
  let k = 1, best = -Infinity;
  for (let i = 1; i < P.length; i++) if (X[i] > X[i - 1] && P[i] - P[i - 1] > best) { best = P[i] - P[i - 1]; k = i; }
  return 0.5 * (X[k] + X[k - 1]);
}

// ---------- SVG plot: computed vs measured p/p0 along x/H* ----------
const SHOCK_MARGIN = 0.3;      // x/H* excluded either side of the shocks
const RECOVERY_LENGTH = 1.5;   // further x/H* behind the later shock before "recovered flow"
const SVG_W = 760, SVG_H = 420, SVG_PAD = { l: 60, r: 20, t: 40, b: 50 };
function svgPlot(title, comp, exact, walls) {
  const xs = comp.map(c => c.xH).concat(...walls.map(w => w.pts.map(p => p.xH)));
  const ys = comp.map(c => c.p).concat(exact.map(e => e.p), ...walls.map(w => w.pts.map(p => p.p)));
  const x0 = Math.floor(Math.min(...xs)), x1 = Math.ceil(Math.max(...xs));
  const y0 = Math.floor(Math.min(...ys) * 20) / 20, y1 = Math.ceil(Math.max(...ys) * 20) / 20;
  const pw = SVG_W - SVG_PAD.l - SVG_PAD.r, ph = SVG_H - SVG_PAD.t - SVG_PAD.b;
  const X = v => (SVG_PAD.l + (v - x0) / (x1 - x0) * pw).toFixed(1), Y = v => (SVG_PAD.t + (1 - (v - y0) / (y1 - y0)) * ph).toFixed(1);
  const line = (pts, attrs) => `<polyline fill="none" ${attrs} points="${pts.map(p => X(p.xH) + "," + Y(p.p)).join(" ")}"/>`;
  const out = [`<svg xmlns="http://www.w3.org/2000/svg" width="${SVG_W}" height="${SVG_H}" font-family="ui-monospace, Menlo, monospace" font-size="11">`,
    `<rect width="100%" height="100%" fill="#101317"/>`,
    `<text x="${SVG_PAD.l}" y="22" fill="#e7e9ec" font-size="13">${title}</text>`];
  for (let v = x0; v <= x1; v++) out.push(`<line x1="${X(v)}" x2="${X(v)}" y1="${SVG_PAD.t}" y2="${SVG_PAD.t + ph}" stroke="#1f242b"/><text x="${X(v)}" y="${SVG_PAD.t + ph + 16}" fill="#7d8794" text-anchor="middle">${v}</text>`);
  for (let v = y0; v <= y1 + 1e-9; v += 0.05) out.push(`<line x1="${SVG_PAD.l}" x2="${SVG_PAD.l + pw}" y1="${Y(v)}" y2="${Y(v)}" stroke="#1f242b"/><text x="${SVG_PAD.l - 6}" y="${(+Y(v) + 4).toFixed(1)}" fill="#7d8794" text-anchor="end">${v.toFixed(2)}</text>`);
  out.push(`<rect x="${SVG_PAD.l}" y="${SVG_PAD.t}" width="${pw}" height="${ph}" fill="none" stroke="#3a414c"/>`,
    `<text x="${SVG_PAD.l + pw / 2}" y="${SVG_H - 12}" fill="#7d8794" text-anchor="middle">x / H*  (from the throat)</text>`,
    `<text transform="translate(16 ${SVG_PAD.t + ph / 2}) rotate(-90)" fill="#7d8794" text-anchor="middle">p / p0</text>`,
    line(exact, `stroke="#98a0ab" stroke-width="1.25" stroke-dasharray="4 3"`),
    line(comp, `stroke="#d9653f" stroke-width="2"`));
  walls.forEach((w, k) => w.pts.forEach(p => out.push(k === 0
    ? `<circle cx="${X(p.xH)}" cy="${Y(p.p)}" r="4" fill="none" stroke="#4a8fe0" stroke-width="1.75"/>`
    : `<rect x="${(+X(p.xH) - 3.5).toFixed(1)}" y="${(+Y(p.p) - 3.5).toFixed(1)}" width="7" height="7" fill="#4a8fe0"/>`)));
  const lx = SVG_PAD.l + pw - 250, ly = SVG_PAD.t + ph - 72;
  out.push(`<rect x="${lx - 8}" y="${ly - 14}" width="250" height="78" fill="#101317" stroke="#3a414c"/>`,
    `<line x1="${lx}" x2="${lx + 22}" y1="${ly}" y2="${ly}" stroke="#d9653f" stroke-width="2"/><text x="${lx + 30}" y="${ly + 4}" fill="#e7e9ec">quasi-1D Euler (numerical)</text>`,
    `<line x1="${lx}" x2="${lx + 22}" y1="${ly + 18}" y2="${ly + 18}" stroke="#98a0ab" stroke-dasharray="4 3"/><text x="${lx + 30}" y="${ly + 22}" fill="#e7e9ec">quasi-1D exact</text>`,
    `<circle cx="${lx + 11}" cy="${ly + 36}" r="4" fill="none" stroke="#4a8fe0" stroke-width="1.75"/><text x="${lx + 30}" y="${ly + 40}" fill="#e7e9ec">measured, ${walls[0].name} wall</text>`,
    `<rect x="${lx + 7.5}" y="${ly + 50.5}" width="7" height="7" fill="#4a8fe0"/><text x="${lx + 30}" y="${ly + 58}" fill="#e7e9ec">measured, ${walls[1] ? walls[1].name : ""} wall</text>`,
    `</svg>`);
  return out.join("\n") + "\n";
}

// opts: { N, write (default true), quiet } — test.js runs a coarser, silent copy.
function main(opts) {
  const o = Object.assign({ write: true, quiet: false }, opts || {});
  const log = o.quiet ? () => {} : console.log;
  const cfg = JSON.parse(fs.readFileSync(path.join(HERE, "config.json"), "utf8"));
  if (o.N) cfg.numerics.N = o.N;
  const need = [cfg.geometry.upper, cfg.geometry.lower].concat(Object.values(cfg.cases).map(c => c.file));
  const missing = need.filter(f => !fs.existsSync(path.join(HERE, f)));
  if (missing.length) {
    console.log("SKIPPED: Sajben data not present. Missing:\n  " + missing.map(f => path.join("validation/sajben", f)).join("\n  "));
    console.log("See validation/sajben/README.md for where to get them.");
    return null;
  }
  const unset = ["p0", "T0"].filter(k => cfg[k] == null)
    .concat(Object.entries(cfg.cases).filter(([, c]) => c.exitPressureRatio == null).map(([k]) => "cases." + k + ".exitPressureRatio"));
  if (unset.length) throw new Error("Fill these in validation/sajben/config.json from the NPARC case page: " + unset.join(", "));

  // Geometry: h(x) = y_upper - y_lower on the union of the two x grids (SI).
  const scale = cfg.geometry.lengthScale || 1;
  const up = wall(path.join(HERE, cfg.geometry.upper), cfg.geometry.columns, scale);
  const lo = wall(path.join(HERE, cfg.geometry.lower), cfg.geometry.columns, scale);
  const xs = Array.from(new Set(up.x.concat(lo.x))).sort((a, b) => a - b)
    .filter(x => x >= Math.max(up.x[0], lo.x[0]) && x <= Math.min(up.x[up.x.length - 1], lo.x[lo.x.length - 1]));
  const h = xs.map(x => interp(up.x, up.y, x) - interp(lo.x, lo.y, x));
  const xThroat = cfg.geometry.xThroat != null ? cfg.geometry.xThroat * scale : xs[h.indexOf(Math.min(...h))];
  const hThroat = interp(xs, h, xThroat);
  const xH = x => (x - xThroat) / hThroat;
  log("geometry: " + xs.length + " stations, x/H* in [" + xH(xs[0]).toFixed(3) + ", " + xH(xs[xs.length - 1]).toFixed(3) +
              "], H* = " + (hThroat * 1000).toFixed(2) + " mm, A_exit/A* = " + (h[h.length - 1] / hThroat).toFixed(4));

  const gas = C.resolveGas(cfg.gas).gas;
  if (o.write) fs.mkdirSync(OUT, { recursive: true });
  const report = {};
  for (const [name, cs] of Object.entries(cfg.cases)) {
    const N = cfg.numerics.N;
    let s = Z.initNozzleState({ type: "table", areaTable: { x: xs, A: h }, gas, p0: cfg.p0, T0: cfg.T0,
                                pb: cs.exitPressureRatio * cfg.p0, N, cfl: cfg.numerics.cfl, flux: cfg.numerics.flux });
    const transit = s.L / Math.sqrt(1.4 * gas.R * cfg.T0);
    let prev = s.rho.slice(), t = transit, res = 1, k = 0;
    for (; k < cfg.numerics.maxTransits && res > cfg.numerics.steadyTol; k++) {
      s = Z.nozzleStepUntil(s, t);
      let d = 0, m = 0;
      for (let i = 2; i < N + 2; i++) { d += Math.abs(s.rho[i] - prev[i]); m += s.rho[i]; }
      res = d / m; prev = s.rho.slice(); t += transit;
    }
    const sol = Z.extractNozzle(s), exact = Z.quasi1DSolution(s, "auto"), x0 = xs[0];
    const comp = sol.x.map((x, i) => ({ xH: xH(x + x0), p: sol.p[i] / cfg.p0, M: sol.M[i] }));
    const ex = exact.x.map((x, i) => ({ xH: xH(x + x0), p: exact.p[i] / cfg.p0 }));
    const rows = numericRows(path.join(HERE, cs.file));
    const walls = cs.walls.map(w => ({
      name: w.name,
      pts: rows.map(r => ({ xH: cs.xScale === "throat-heights" ? r[w.x] : xH(r[w.x] * scale), p: cs.pIsRatio ? r[w.p] : r[w.p] / cfg.p0 }))
               .filter((q, i, a) => i === 0 || q.xH !== a[i - 1].xH)          // drop repeated outflow point
    }));
    const cx = comp.map(c => c.xH), cp = comp.map(c => c.p);
    const r = {
      exitPressureRatio: cs.exitPressureRatio, transits: k, residual: res, regime: exact.regime,
      shock_xH_computed: shockAt(cx, cp), shock_xH_exactQuasi1D: exact.xShock != null ? xH(exact.xShock + x0) : null,
      walls: {}
    };
    // Error by region: ahead of both shocks (isentropic, the fair test of the
    // model), between the measured and computed shocks (shock-position error),
    // and the recovered flow well behind both.
    for (const w of walls) {
      const xsM = shockAt(w.pts.map(q => q.xH), w.pts.map(q => q.p));
      const lo = Math.min(xsM, r.shock_xH_computed) - SHOCK_MARGIN, hi = Math.max(xsM, r.shock_xH_computed) + SHOCK_MARGIN;
      const acc = { all: [], upstream: [], downstream: [] };
      for (const q of w.pts) {
        const d = interp(cx, cp, q.xH) - q.p;
        acc.all.push(d);
        if (q.xH < lo) acc.upstream.push(d);
        else if (q.xH > hi + RECOVERY_LENGTH) acc.downstream.push(d);
      }
      const rms = a => (a.length ? Math.sqrt(a.reduce((t, d) => t + d * d, 0) / a.length) : null);
      const maxAbs = a => (a.length ? Math.max(...a.map(Math.abs)) : null);
      r.walls[w.name] = {
        points: w.pts.length, shock_xH_measured: xsM, shock_xH_error: r.shock_xH_computed - xsM,
        rmsDeltaP_all: rms(acc.all),
        upstream: { points: acc.upstream.length, rmsDeltaP: rms(acc.upstream), maxDeltaP: maxAbs(acc.upstream) },
        downstream: { points: acc.downstream.length, rmsDeltaP: rms(acc.downstream), maxDeltaP: maxAbs(acc.downstream) }
      };
    }
    report[name] = r;
    if (!o.write) { log(name + ": " + JSON.stringify(r)); continue; }
    // A case file for run.js / the Nozzle tab, and the measurements in that
    // case's own x coordinate (m from the duct inlet) for "Load p/p0 CSV".
    const caseFile = { kind: "nozzle", name: "sajben-" + name, gas: cfg.gas,
      geometry: { type: "table", x: xs.map(v => +(v - x0).toFixed(7)), A: h.map(v => +v.toFixed(7)) },
      p0: cfg.p0, T0: cfg.T0, pb: cs.exitPressureRatio * cfg.p0,
      numerics: { N: cfg.numerics.N, cfl: cfg.numerics.cfl, flux: cfg.numerics.flux }, steady: { tol: 1e-6, tMax: 0.5 } };
    fs.writeFileSync(path.join(OUT, name + ".case.json"), JSON.stringify(caseFile, null, 2) + "\n");
    fs.writeFileSync(path.join(OUT, name + ".measured_duct_x.csv"),
      "# x [m from duct inlet], p/p0 — measured, both walls (load with 'Load p/p0 CSV' in the Nozzle tab)\n" +
      [].concat(...walls.map(w => w.pts.map(q => (xThroat + q.xH * hThroat - x0).toFixed(6) + "," + q.p))).join("\n") + "\n");
    fs.writeFileSync(path.join(OUT, name + ".computed.csv"), C.tableToCSV(comp.map((c, i) => ({ x_over_H: c.xH, p_over_p0: c.p, M: c.M, p_over_p0_exact_quasi1D: ex[i].p }))));
    fs.writeFileSync(path.join(OUT, name + ".measured.csv"), C.tableToCSV([].concat(...walls.map(w => w.pts.map(q => ({ wall: w.name, x_over_H: q.xH, p_over_p0: q.p }))))));
    fs.writeFileSync(path.join(OUT, name + ".svg"),
      svgPlot("Sajben diffuser, " + name + " shock (p_exit/p0 = " + cs.exitPressureRatio + "): wall pressure", comp, ex, walls));
    log(name + ": " + JSON.stringify(r));
  }
  if (o.write) fs.writeFileSync(path.join(OUT, "report.json"), JSON.stringify(report, null, 2) + "\n");
  return report;
}

if (require.main === module) {
  try { main(); } catch (e) { console.error("sajben.js: " + e.message); process.exit(1); }
} else {
  module.exports = { main };
}
