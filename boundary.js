/* ============================================================
   boundary.js  —  boundary-layer displacement correction for
   quasi-1D duct flow (math layer, no DOM)
   ------------------------------------------------------------
   Viscous–inviscid interaction: the inviscid core sees an effective
   area reduced by the displacement thickness δ* of the wall boundary
   layers, and δ* is computed from the core's edge conditions:

       A_eff(x) = span · (h(x) − nWalls · δ*(x))          ("channel")
       A_eff(x) = π (r(x) − δ*(x))²,  r = √(A/π)          ("axisymmetric")

   Boundary layer: turbulent, adiabatic wall, integral method
     momentum integral (compressible von Kármán):
       dθ/dx = Cf/2 − (θ/uₑ)(2 + H − Mₑ²) duₑ/dx
     entrainment (Head 1958, Green's compressible form):
       d(ρₑ uₑ θ H₁)/dx = ρₑ uₑ F(H₁)
       H₁ = 3.3 + 0.8234 (Hᵢ − 1.1)^−1.287   (Hᵢ ≤ 1.6)
          = 3.3 + 1.5501 (Hᵢ − 0.6778)^−3.064 (Hᵢ > 1.6)
       F  = 0.0306 (H₁ − 3)^−0.6169
     compressible shape factor (adiabatic wall, recovery factor r):
       H + 1 = (Hᵢ + 1)(1 + r (γ−1)/2 Mₑ²)
     skin friction: Ludwieg–Tillmann, Cf = 0.246·10^(−0.678 Hᵢ) Re_θ^−0.268,
       evaluated at Eckert's reference temperature
       T* /Tₑ = 0.5 + 0.039 Mₑ² + 0.5 T_aw/Tₑ,  Cf = Cf_LT(Re_θ*) · ρ* /ρₑ
     viscosity: Sutherland's law for air.
   Separation: when Hᵢ reaches HI_SEP the layer is flagged separated and
   Hᵢ is held there (Cf → 0) until entrainment brings it back — the usual
   crude continuation of an integral method, reported, not hidden.

   Coupling: exact steady quasi-1D (nozzle.js quasi1DSolution) on A_eff with
   under-relaxed δ*; with a shock in the duct the shock position is solved
   for by bisection (see viscousQuasi1D), because direct iteration locks
   into a limit cycle there.

   What it does NOT do: wall friction does not remove momentum from the
   core (displacement effect only), sidewalls are ignored in "channel"
   mode, the shock / boundary-layer interaction is not resolved (the edge
   velocity jumps at the shock and the integral equations take the jump
   as given), and the inlet boundary-layer thickness θ₀ is an input.
   References: White, Viscous Fluid Flow (3rd ed.) §6-7; Cebeci & Bradshaw,
   Momentum Transfer in Boundary Layers; Green, Aero. Res. Council R&M 3791.

   Public API: integrateBoundaryLayer(edge, opts),
               viscousQuasi1D(nozzleParams, blOpts) → { solution, bl, Aeff, deltaStar, xShock, reliable, … },
               viscousQuasi1DGen (the same, as a generator for time-sliced browser runs)
   ============================================================ */

if (typeof quasi1DSolution === "undefined" && typeof require === "function") {
  var _blNozzle       = require("./nozzle");
  var quasi1DSolution = _blNozzle.quasi1DSolution;
  var initNozzleState = _blNozzle.initNozzleState;
  var nozzleGeometry  = _blNozzle.nozzleGeometry;
}

const BL_RECOVERY = 0.89;          // turbulent recovery factor r ≈ Pr^(1/3)
const BL_HI_ZPG = 1.4;             // incompressible shape factor, zero pressure gradient
const BL_HI_SEP = 2.4;             // Head's separation criterion on Hᵢ
const BL_THETA_MIN = 1e-7;         // [m] floor so Re_θ stays finite at a sharp leading edge
const BL_DLNU_STEP = 0.004;        // max |Δ ln uₑ| per integration sub-step
const SUTH_MU0 = 1.716e-5, SUTH_T0 = 273.15, SUTH_S = 110.4;   // Sutherland's law, air
const VI_RELAX = 0.15;               // under-relaxation of δ* between coupling iterations
const VI_MAX_ITERS = 150;
const VI_TOL = 1e-5;               // max |Δδ*| / h_min to stop iterating
const VI_SAMPLES = 300;            // stations for the effective-area table
const VI_BISECT_ITERS = 40;        // bisection steps on the shock position
const VI_PEXIT_GAP = 0.02;         // |p_exit − p_b|/p_b above this after bisection: no matching state
                                   // (the discrete coupling leaves jumps of ~1%; separation onset ~5%+)
const SHOCK_JUMP_DX = 1e-9;       // width of the inserted shock jump, fraction of the duct length
const VI_SHOCK_TOL = 1e-5;         // shock-position tolerance, fraction of the duct length

function sutherlandAir(T) { return SUTH_MU0 * Math.pow(T / SUTH_T0, 1.5) * (SUTH_T0 + SUTH_S) / (T + SUTH_S); }

// Head's shape-factor relations.
function headH1(Hi) {
  return Hi <= 1.6 ? 3.3 + 0.8234 * Math.pow(Hi - 1.1, -1.287) : 3.3 + 1.5501 * Math.pow(Hi - 0.6778, -3.064);
}
function headF(H1) { return 0.0306 * Math.pow(Math.max(H1 - 3, 1e-6), -0.6169); }
// H₁ decreases monotonically with Hᵢ: invert by bisection on (1.1, BL_HI_SEP].
function headHi(H1) {
  if (H1 >= headH1(1.1 + 1e-6)) return 1.1 + 1e-6;
  if (H1 <= headH1(BL_HI_SEP)) return BL_HI_SEP;
  let lo = 1.1 + 1e-6, hi = BL_HI_SEP;
  for (let it = 0; it < 60; it++) { const m = 0.5 * (lo + hi); if (headH1(m) > H1) lo = m; else hi = m; }
  return 0.5 * (lo + hi);
}

// Integrate the boundary layer along edge = { x, ue, Me, rhoe, Te, gamma } (arrays).
// opts: { theta0 [m] (inlet momentum thickness), Hi0 }
function integrateBoundaryLayer(edge, opts) {
  const o = Object.assign({ theta0: 0, Hi0: BL_HI_ZPG }, opts || {});
  const n = edge.x.length;
  const out = { x: edge.x.slice(), theta: new Array(n), H: new Array(n), Hi: new Array(n), deltaStar: new Array(n),
                Cf: new Array(n), separated: new Array(n), separation: [] };
  const at = (i, w) => {                 // edge state linearly interpolated inside interval i → i+1
    const L = (a) => a[i] + w * (a[i + 1] - a[i]);
    return { ue: Math.exp(Math.log(edge.ue[i]) + w * (Math.log(edge.ue[i + 1]) - Math.log(edge.ue[i]))),
             Me: L(edge.Me), rhoe: L(edge.rhoe), Te: L(edge.Te), g: L(edge.gamma) };
  };
  const props = (e, theta, Hi) => {
    const m = 1 + BL_RECOVERY * 0.5 * (e.g - 1) * e.Me * e.Me;          // T_aw / Tₑ
    const H = (Hi + 1) * m - 1;
    const Tstar = e.Te * (0.5 + 0.039 * e.Me * e.Me + 0.5 * m);          // Eckert reference temperature
    const rhoStar = e.rhoe * e.Te / Tstar;
    const Re = Math.max(rhoStar * e.ue * theta / sutherlandAir(Tstar), 1);
    const Cf = Hi >= BL_HI_SEP ? 0 : 0.246 * Math.pow(10, -0.678 * Hi) * Math.pow(Re, -0.268) * rhoStar / e.rhoe;
    return { H, Cf };
  };
  // state: theta and Y = ρₑ uₑ θ H₁
  let theta = Math.max(o.theta0, BL_THETA_MIN), Hi = o.Hi0;
  const e0 = at(0, 0);
  let Y = e0.rhoe * e0.ue * theta * headH1(Hi);
  const record = (i, e, sep) => {
    const p = props(e, theta, Hi);
    out.theta[i] = theta; out.Hi[i] = Hi; out.H[i] = p.H; out.Cf[i] = p.Cf; out.deltaStar[i] = p.H * theta; out.separated[i] = sep;
  };
  let sep = false;
  record(0, e0, sep);
  for (let i = 0; i < n - 1; i++) {
    const dx = edge.x[i + 1] - edge.x[i];
    const dlnu = Math.log(edge.ue[i + 1] / edge.ue[i]);
    const nSub = Math.max(1, Math.ceil(Math.abs(dlnu) / BL_DLNU_STEP));
    const h = 1 / nSub;
    for (let k = 0; k < nSub; k++) {
      const rhs = (w, th, y) => {
        const e = at(i, w);
        const hi = headHi(y / (e.rhoe * e.ue * th));
        const p = props(e, th, hi);
        const dlnudx = dlnu / dx;                         // uₑ varies exponentially inside the interval
        return { dth: p.Cf / 2 - th * (2 + p.H - e.Me * e.Me) * dlnudx, dy: e.rhoe * e.ue * headF(headH1(hi)) };
      };
      const w0 = k * h, w1 = (k + 1) * h, ddx = h * dx;
      const a = rhs(w0, theta, Y);
      const thP = Math.max(BL_THETA_MIN, theta + ddx * a.dth), yP = Y + ddx * a.dy;
      const b = rhs(w1, thP, yP);
      theta = Math.max(BL_THETA_MIN, theta + 0.5 * ddx * (a.dth + b.dth));
      Y = Y + 0.5 * ddx * (a.dy + b.dy);
      // separation cap: hold Hᵢ at BL_HI_SEP (H₁ at its floor)
      const e = at(i, w1), H1 = Y / (e.rhoe * e.ue * theta), H1sep = headH1(BL_HI_SEP);
      if (H1 <= H1sep) { Y = e.rhoe * e.ue * theta * H1sep; Hi = BL_HI_SEP; } else Hi = headHi(H1);
    }
    const nowSep = Hi >= BL_HI_SEP - 1e-9;
    if (nowSep && !sep) out.separation.push({ x: edge.x[i + 1], event: "separation" });
    if (!nowSep && sep) out.separation.push({ x: edge.x[i + 1], event: "reattachment" });
    sep = nowSep;
    record(i + 1, at(i, 1), sep);
  }
  return out;
}

// ---------- Viscous–inviscid coupling for a quasi-1D duct ----------
// With a shock in the duct, iterating δ* ↔ inviscid flow directly does not
// converge: the shock position is very sensitive to the downstream effective
// area, and the δ* jump at the shock moves with it (a limit cycle). So the
// shock position is the unknown instead:
//   for a trial x_s, iterate δ* with the shock *held* at x_s (stable), and
//   read off the exit pressure that would hold it there; then bisect on x_s
//   until that exit pressure equals p_b (it falls as x_s moves aft). Where the
//   layer starts to separate at the shock, p_exit(x_s) jumps; if p_b lies in
//   that jump no steady state of this model matches it, and the result says so
//   (pressureGap) rather than inventing one. Once the layer separates, p_exit(x_s)
//   is no longer even monotone and the answer depends on the grid. That is the
//   limit of an integral method in direct mode: it cannot follow a separated
//   shock–boundary-layer interaction, so such results come back reliable: false.
// Without a shock in the duct the plain iteration is used.
//
// nozzleParams: as for initNozzleState (a gas.js gas and SI units).
// blOpts: { geometry: "channel" | "axisymmetric", span [m] (channel), nWalls (channel, default 2),
//           theta0 [m], relax, maxIters, tol, N (stations) }
// viscousQuasi1DGen is the same computation as a generator that yields
// { step, xShock, change } after every inviscid + boundary-layer pass, so the
// browser can run it in time slices (LabUI.runSliced); its return value is the result.
function* viscousQuasi1DGen(nozzleParams, blOpts) {
  const o = Object.assign({ geometry: "channel", span: 1, nWalls: 2, theta0: 0, relax: VI_RELAX,
                            maxIters: VI_MAX_ITERS, tol: VI_TOL, N: VI_SAMPLES });
  for (const k in blOpts || {}) if (blOpts[k] != null) o[k] = blOpts[k];   // unset options keep their defaults
  const geo = nozzleGeometry(nozzleParams), L = geo.L, M = o.N;
  const xs = [], A0 = [];
  for (let k = 0; k <= M; k++) { const x = (k / M) * L; xs.push(x); A0.push(geo.A(x)); }
  // wall-normal size: channel height or pipe radius
  const size = A0.map(a => (o.geometry === "axisymmetric" ? Math.sqrt(a / Math.PI) : a / o.span));
  const sizeMin = Math.min(...size);
  const areaOf = (i, ds) => (o.geometry === "axisymmetric" ? Math.PI * Math.pow(size[i] - ds, 2)
                                                            : o.span * (size[i] - o.nWalls * ds));
  const base = Object.assign({}, nozzleParams, { N: M });
  const pb = nozzleParams.pb;

  // One inviscid solve on δ*, then the boundary layer along it.
  const pass = (ds, mode) => {
    const Aeff = ds.map((d, i) => areaOf(i, d));
    if (Aeff.some(a => !(a > 0))) throw new Error("Displacement thickness closes the duct: the boundary layer is too thick for this geometry.");
    const state = initNozzleState(Object.assign({}, base, { type: "table", areaTable: { x: xs, A: Aeff } }));
    const sol = quasi1DSolution(state, mode);
    const R = state.thermo.R;
    const edge = { x: sol.x.slice(), ue: sol.u.map(u => Math.max(u, 1e-3)), Me: sol.M.slice(), rhoe: sol.rho.slice(), Te: sol.T.slice(),
                   gamma: sol.a.map((a, i) => a * a / (R * sol.T[i])) };
    // Put the shock jump exactly at x_s, each side's state extrapolated linearly
    // from its own two nearest cells, so δ* moves continuously with x_s instead
    // of snapping as the shock crosses a cell centre.
    const k = sol.regime === "shock" ? sol.x.findIndex(x => x > sol.xShock) : -1;
    if (k >= 2 && k + 1 < sol.x.length) {
      const xsh = sol.xShock, X = edge.x;
      for (const key of ["ue", "Me", "rhoe", "Te", "gamma"]) {
        const a = edge[key];
        const ext = (i, j) => a[i] + (xsh - X[i]) / (X[j] - X[i]) * (a[j] - a[i]);
        a.splice(k, 0, ext(k - 1, k - 2), ext(k, k + 1));
      }
      X.splice(k, 0, xsh, xsh + SHOCK_JUMP_DX * L);
    }
    const bl = integrateBoundaryLayer(edge, { theta0: o.theta0 });
    const dsNew = xs.map(x => {
      const X = bl.x; if (x <= X[0]) return bl.deltaStar[0]; if (x >= X[X.length - 1]) return bl.deltaStar[X.length - 1];
      let a = 0, b = X.length - 1; while (b - a > 1) { const m = (a + b) >> 1; if (X[m] <= x) a = m; else b = m; }
      return bl.deltaStar[a] + (x - X[a]) / (X[b] - X[a]) * (bl.deltaStar[b] - bl.deltaStar[a]);
    });
    let change = 0;
    for (let i = 0; i <= M; i++) change = Math.max(change, Math.abs(dsNew[i] - ds[i]) / sizeMin);
    return { state, sol, bl, dsNew, change, Aeff };
  };
  // Fixed-mode iteration (mode: "auto", or { shockAt }) from a starting δ*.
  const history = [];
  const iterate = function* (ds0, mode) {
    let ds = ds0.slice(), r = null, it = 0;
    for (; it < o.maxIters; it++) {
      r = pass(ds, mode);
      yield { step: history.length, xShock: mode.shockAt != null ? mode.shockAt : null, change: r.change };
      if (r.change < o.tol && it > 0) break;
      ds = ds.map((d, i) => d + o.relax * (r.dsNew[i] - d));
    }
    return Object.assign(r, { ds, iters: it + 1, converged: r.change < o.tol });
  };

  const zero = new Array(M + 1).fill(0);
  const inviscid = pass(zero, "auto");
  let res = null, xShock = null, pressureGap = false, pExitMismatch = 0;
  if (inviscid.sol.regime !== "shock") {
    res = (yield* iterate(zero, "auto"));                       // no shock in the duct: direct iteration is stable
    history.push({ step: "direct", iters: res.iters, change: res.change });
  } else {
    // Bracket the shock between the throat of the effective duct and the exit.
    const exitFor = function* (xsTry, ds0) {
      const r = yield* iterate(ds0, { shockAt: xsTry });
      history.push({ step: "shock at", xShock: xsTry, pExit: r.sol.pExit, iters: r.iters, change: r.change });
      return r;
    };
    let lo = inviscid.sol.xThroat + 1e-3 * L, hi = L * (1 - 1e-3);
    let rLo = yield* exitFor(lo, zero);
    lo = Math.max(lo, rLo.sol.xThroat + 1e-3 * L);
    if (rLo.sol.pExit < pb) {
      // even a shock at the throat leaves the exit pressure below p_b: the flow unchokes
      res = (yield* iterate(rLo.ds, "auto"));
      history.push({ step: "unchoked", iters: res.iters, change: res.change });
    } else {
      const rHi = yield* exitFor(hi, rLo.ds);
      if (rHi.sol.pExit > pb) {
        res = (yield* iterate(rHi.ds, "auto"));                 // shock blown out of the duct
        history.push({ step: "shock leaves duct", iters: res.iters, change: res.change });
      } else {
        let warm = rLo.ds;
        for (let k = 0; k < VI_BISECT_ITERS && hi - lo > VI_SHOCK_TOL * L; k++) {
          const mid = 0.5 * (lo + hi), r = yield* exitFor(mid, warm);
          warm = r.ds;
          if (r.sol.pExit > pb) lo = mid; else hi = mid;
        }
        xShock = 0.5 * (lo + hi);
        res = (yield* iterate(warm, { shockAt: xShock }));
        // p_exit(x_s) is discontinuous where the boundary layer starts to separate at
        // the shock; if p_b falls in that gap, no attached/separated state matches it.
        pExitMismatch = (res.sol.pExit - pb) / pb;
        pressureGap = Math.abs(pExitMismatch) > VI_PEXIT_GAP;
      }
    }
  }
  return { state: res.state, solution: res.sol, bl: res.bl, x: xs, A: A0, Aeff: res.Aeff, deltaStar: res.ds,
           xShock: res.sol.regime === "shock" ? res.sol.xShock : null, xShockInviscid: inviscid.sol.xShock,
           converged: res.converged, pressureGap, pExitMismatch,
           separated: res.bl.separation.length > 0,
           // attached and matched: the only results this direct-mode method can stand behind
           reliable: res.converged && !pressureGap && res.bl.separation.length === 0, history, options: o };
}

function viscousQuasi1D(nozzleParams, blOpts) {
  const gen = viscousQuasi1DGen(nozzleParams, blOpts);
  for (;;) { const r = gen.next(); if (r.done) return r.value; }
}

if (typeof module !== "undefined" && module.exports) {
  module.exports = { sutherlandAir, headH1, headHi, headF, integrateBoundaryLayer, viscousQuasi1D, viscousQuasi1DGen,
                     BL_HI_SEP };
}
