/* ============================================================
   nozzle.js  —  Quasi-1D Euler nozzle solver
   ------------------------------------------------------------
   Equations:  dU/dt + dF/dx = H
               H = -(dA/dx)/A * [rho*u, rho*u^2, u*(E+p)]
   Gas:  legacy non-dimensional ideal gas (gamma, R = 1, A_throat = 1), or
         any gas.js model in SI (params.gas): calorically or thermally perfect.
   A(x): parabolic profiles ("convergent" | "divergent" | "CD") scaled by
         Athroat, or a tabulated profile (type "table", params.areaTable
         { x: [...], A: [...] }) interpolated with monotone cubic PCHIP.
   Reference solution: exact steady quasi-1D flow for every regime —
         unchoked, choked with a normal shock in the divergent part, or
         shock-free supersonic — for any gas (quasi1DSolution).
   Depends on solver.js (primToCons, FLUXES, hllc, computeRhs, stableDt,
   extract, idealEos, gasEos) and gas.js (thermo helpers).
   ============================================================ */

// ---------- Dependencies on solver.js / gas.js ----------
// Browser: solver.js and gas.js are loaded first and their functions are page
// globals, so `typeof primToCons` is "function" and this block is skipped (a
// bare `var` redeclaration of a global *function* leaves it untouched).
// Node: module scope does not see another module's functions, so pull them in
// explicitly. Only function declarations may be aliased with `var` here:
// solver.js declares FLUXES with `const`, and a `var FLUXES` anywhere in this
// file would be a redeclaration error in the browser — see _fluxTable().
if (typeof primToCons === "undefined" && typeof require === "function") {
  var _solver    = require("./solver");
  var primToCons = _solver.primToCons;
  var hllc       = _solver.hllc;
  var computeRhs = _solver.computeRhs;
  var stableDt   = _solver.stableDt;
  var extract    = _solver.extract;
  var idealEos   = _solver.idealEos;
  var gasEos     = _solver.gasEos;
  var _gasN      = require("./gas");
  var gasSound   = _gasN.gasSound;
  var gasGamma   = _gasN.gasGamma;
  var gasTfromH  = _gasN.gasTfromH;
  var gasIsentropeT = _gasN.gasIsentropeT;
  var gasStagnation = _gasN.gasStagnation;
  var gasModel   = _gasN.gasModel;
}
function _fluxTable() {
  return typeof FLUXES !== "undefined" ? FLUXES : _solver.FLUXES;
}

const BISECT_ITERS = 100;      // bisection cap for the steady-flow root finds
const INIT_U_FRACTION = 0.01;  // initial guess: u = 1% of the stagnation sound speed
const AREA_SAMPLES = 2000;     // samples used to locate the throat of an A(x)

// ---------- Thermo for the stagnation relations ----------
// The legacy path is a calorically perfect gas with R = 1, cp = gamma/(gamma-1);
// using the same h / phi interface lets one code path serve both.
function _legacyThermo(gamma) {
  const cp = gamma / (gamma - 1);
  return { name: "ideal", R: 1, cp: () => cp, h: T => cp * T, phi: T => cp * Math.log(T) };
}

// ---------- Area profiles ----------
// Analytic profiles: A_throat = 1 (legacy) and AR = A_max / A_throat.
function nozzleArea(x, L, AR, type) {
  const s = x / L;
  if (type === "convergent") {
    return 1 + (AR - 1) * (1 - s) * (1 - s);
  }
  if (type === "divergent") {
    return 1 + (AR - 1) * s * s;
  }
  // CD: symmetric, throat at x = L/2
  const t = 2 * s - 1;  // -1 at inlet, 0 at throat, +1 at exit
  return 1 + (AR - 1) * t * t;
}

function nozzleAreaDeriv(x, L, AR, type) {
  const s = x / L;
  if (type === "convergent") {
    return -2 * (AR - 1) * (1 - s) / L;
  }
  if (type === "divergent") {
    return  2 * (AR - 1) * s / L;
  }
  // CD: dA/dx = 4*(AR-1)*(2s-1)/L
  const t = 2 * s - 1;
  return 4 * (AR - 1) * t / L;
}

// Tabulated A(x): monotone piecewise-cubic Hermite (Fritsch–Carlson 1980), so
// dA/dx is continuous and the curve never overshoots the data. x is shifted so
// the duct starts at 0. Returns { L, A(x), dA(x) }.
function makeAreaTable(xs, As) {
  if (!xs || !As || xs.length !== As.length || xs.length < 2) throw new Error("Area table needs matching x and A arrays (≥ 2 points).");
  const n = xs.length, x0 = xs[0];
  const X = xs.map(v => v - x0), A = As.slice();
  for (let i = 1; i < n; i++) if (!(X[i] > X[i - 1])) throw new Error("Area table x must be strictly increasing.");
  for (let i = 0; i < n; i++) if (!(A[i] > 0)) throw new Error("Area table A must be positive.");
  const d = [], m = new Array(n);
  for (let i = 0; i < n - 1; i++) d.push((A[i + 1] - A[i]) / (X[i + 1] - X[i]));
  m[0] = d[0]; m[n - 1] = d[n - 2];
  for (let i = 1; i < n - 1; i++) {
    if (d[i - 1] * d[i] <= 0) { m[i] = 0; continue; }
    const w1 = 2 * (X[i + 1] - X[i]) + (X[i] - X[i - 1]), w2 = (X[i + 1] - X[i]) + 2 * (X[i] - X[i - 1]);
    m[i] = (w1 + w2) / (w1 / d[i - 1] + w2 / d[i]);      // weighted harmonic mean
  }
  const seg = x => {
    let lo = 0, hi = n - 2;
    while (lo < hi) { const mid = (lo + hi + 1) >> 1; if (X[mid] <= x) lo = mid; else hi = mid - 1; }
    return lo;
  };
  const clampX = x => Math.max(0, Math.min(X[n - 1], x));
  return {
    L: X[n - 1],
    A: x => {
      x = clampX(x); const i = seg(x), h = X[i + 1] - X[i], t = (x - X[i]) / h;
      const h00 = (1 + 2 * t) * (1 - t) * (1 - t), h10 = t * (1 - t) * (1 - t), h01 = t * t * (3 - 2 * t), h11 = t * t * (t - 1);
      return h00 * A[i] + h10 * h * m[i] + h01 * A[i + 1] + h11 * h * m[i + 1];
    },
    dA: x => {
      x = clampX(x); const i = seg(x), h = X[i + 1] - X[i], t = (x - X[i]) / h;
      const d00 = 6 * t * t - 6 * t, d10 = 3 * t * t - 4 * t + 1, d01 = -6 * t * t + 6 * t, d11 = 3 * t * t - 2 * t;
      return (d00 * A[i] + d01 * A[i + 1]) / h + d10 * m[i] + d11 * m[i + 1];
    }
  };
}

// Area function for a nozzle geometry { type, L, AR, Athroat, areaTable }.
function nozzleGeometry(g) {
  if (g.type === "table") {
    const tab = makeAreaTable(g.areaTable.x, g.areaTable.A);
    return { L: tab.L, A: tab.A, dA: tab.dA };
  }
  const s = g.Athroat || 1;
  return { L: g.L, A: x => s * nozzleArea(x, g.L, g.AR, g.type), dA: x => s * nozzleAreaDeriv(x, g.L, g.AR, g.type) };
}
// Area at x for a nozzle state or geometry (used by the UI to draw the duct).
function nozzleAreaAt(g, x) { return (g.geom || nozzleGeometry(g)).A(x); }

// Build area + dA/dx arrays for all N+4 cells (including ghost cells).
function buildAreaArrays(N, dx, L, AR, type, geom) {
  const G = geom || nozzleGeometry({ type, L, AR });
  const tot = N + 4;
  const A = new Array(tot), dA = new Array(tot);
  for (let i = 0; i < tot; i++) {
    const xc = (i - 1.5) * dx;                      // ghost cell x positions
    const xc_clamped = Math.max(0, Math.min(L, xc));
    A[i]  = G.A(xc_clamped);
    dA[i] = G.dA(xc_clamped);
  }
  return { A, dA };
}

// ---------- Quasi-1D RHS with geometric source ----------
// Reuses computeRhs from solver.js for the flux divergence part,
// then adds the area source term on interior cells.
function computeNozzleRhs(rho, mom, ene, N, dx, eos, fluxFn, A, dA) {
  const E = typeof eos === "number" ? idealEos(eos) : eos;
  const [r1, r2, r3] = computeRhs(rho, mom, ene, N, dx, E, fluxFn);
  for (let i = 0; i < N; i++) {
    const k = i + 2;
    const q = E.prim(rho[k], mom[k], ene[k]);
    const u = q.u;
    const src = -dA[k] / A[k];
    r1[i] += src * mom[k];            // -(dA/dx)/A * rho*u
    r2[i] += src * rho[k] * u * u;    // -(dA/dx)/A * rho*u^2  (NOT +p)
    r3[i] += src * u * (ene[k] + q.p); // -(dA/dx)/A * u*(E+p)
  }
  return [r1, r2, r3];
}

// Sonic speed of the flow expanded isentropically from (T0): u = a where
// 2 (h0 - h(T)) = gamma(T) R T. Inflow from a reservoir at rest can't exceed it.
function _sonicFromStagnation(th, T0) {
  const h0 = th.h(T0);
  let lo = 1e-3 * T0, hi = T0;                 // f(lo) > 0, f(T0) < 0
  const f = T => 2 * (h0 - th.h(T)) - gasGamma(th, T) * th.R * T;
  for (let it = 0; it < BISECT_ITERS; it++) {
    const mid = 0.5 * (lo + hi);
    if (f(mid) > 0) lo = mid; else hi = mid;
  }
  const T = 0.5 * (lo + hi);
  return { T, u: Math.sqrt(2 * (h0 - th.h(T))) };
}

// ---------- Inlet BC: isentropic reservoir ----------
// Subsonic inflow from a stagnation reservoir (Anderson, CFD: The Basics, Ch. 7):
// one characteristic leaves the domain, so the velocity is extrapolated from
// the first interior cell; the two incoming ones fix the stagnation state, so
// T and p follow from (T0, p0) along the steady isentrope
//     h(T) = h(T0) - u^2 / 2,     p = p0 exp((phi(T) - phi(T0)) / R)
// (for a calorically perfect gas: T = T0 - u^2/(2 cp), p = p0 (T/T0)^(gamma/(gamma-1))).
// This holds total enthalpy and total pressure at exactly (T0, p0). (A simple-
// wave Riemann invariant u + 2a/(gamma-1) = 2a0/(gamma-1) does not: it is an
// unsteady-expansion relation and bleeds ~gamma*M*p0 of total pressure.)
function applyNozzleInletBC(rho, mom, ene, nz) {
  const { eos, thermo: th, p0, T0, uMax } = nz;
  // No backflow into the reservoir, and at most sonic
  const u_b = Math.min(Math.max(mom[2] / rho[2], 0), uMax);
  const T_b = gasTfromH(th, th.h(T0) - 0.5 * u_b * u_b, T0);
  const p_b = p0 * Math.exp((th.phi(T_b) - th.phi(T0)) / th.R);
  const [rg, mg, eg] = eos.cons(p_b / (th.R * T_b), u_b, p_b);
  rho[0] = rg; mom[0] = mg; ene[0] = eg;
  rho[1] = rg; mom[1] = mg; ene[1] = eg;
}

// ---------- Exit BC: back pressure or supersonic extrapolation ----------
function applyNozzleExitBC(rho, mom, ene, nz) {
  const { eos, pb, N } = nz;
  const k = N + 1;  // last interior cell
  const r = rho[k], m = mom[k], e = ene[k];
  const q = eos.prim(r, m, e);
  const Mexit = q.u / q.a;

  if (Mexit >= 1.0) {
    // Supersonic exit: all characteristics leave, extrapolate everything
    rho[N + 2] = r; mom[N + 2] = m; ene[N + 2] = e;
    rho[N + 3] = r; mom[N + 3] = m; ene[N + 3] = e;
  } else {
    // Subsonic exit: specified back pressure, zero-gradient rho and u
    const E_b = eos.cons(r, q.u, pb)[2];
    rho[N + 2] = r; mom[N + 2] = m; ene[N + 2] = E_b;
    rho[N + 3] = r; mom[N + 3] = m; ene[N + 3] = E_b;
  }
}

// ---------- Initialize nozzle state ----------
// params: { type, AR, L, Athroat?, areaTable?, p0, T0, pb, N, cfl, flux,
//           gamma (legacy) | gas (gas.js model, SI) }
function initNozzleState(p) {
  const { N, p0, pb, flux, cfl } = p;
  const geom = nozzleGeometry(p);
  const L = geom.L;
  const dx = L / N;
  const tot = N + 4;
  const T0 = Math.max(p.T0, 1e-6);
  const gas = typeof p.gas === "string" ? gasModel(p.gas, p.thermo) : p.gas;   // name or model
  const eos = gas ? gasEos(gas) : idealEos(p.gamma);
  const thermo = gas || _legacyThermo(p.gamma);
  const sonic = _sonicFromStagnation(thermo, T0);
  const a0 = gasSound(thermo, T0);

  const rho = new Array(tot), mom = new Array(tot), ene = new Array(tot);
  // Linear pressure gradient from p0 at inlet to pb at exit, isentropic temperature.
  for (let i = 0; i < N; i++) {
    const frac = (i + 0.5) / N;
    const pi   = p0 + (pb - p0) * frac;
    const Ti   = gasIsentropeT(thermo, T0, p0, pi);
    const [r, m, e] = eos.cons(pi / (thermo.R * Ti), INIT_U_FRACTION * a0, pi);
    rho[i + 2] = r; mom[i + 2] = m; ene[i + 2] = e;
  }
  const { A, dA } = buildAreaArrays(N, dx, L, p.AR, p.type, geom);
  const state = {
    rho, mom, ene, t: 0, step: 0,
    dx, N, gamma: p.gamma, eos, thermo, gas: gas || null, L, A, dA, geom,
    p0, T0, pb, AR: p.AR, type: p.type, Athroat: p.Athroat || 1, areaTable: p.areaTable || null,
    uMax: sonic.u,
    flux: flux || "HLLC", cfl: cfl || 0.4
  };
  applyNozzleInletBC(rho, mom, ene, state);
  applyNozzleExitBC(rho, mom, ene, state);
  return state;
}

// ---------- Step forward to tTarget (SSP-RK2 + adaptive CFL) ----------
function nozzleStepUntil(state, tTarget) {
  const { N, dx, cfl, A, dA, eos } = state;
  const fluxFn = _fluxTable()[state.flux] || hllc;
  let { rho, mom, ene, t, step } = state;
  let budget = 80000;

  while (t < tTarget && budget-- > 0) {
    let dt = stableDt({ rho, mom, ene, N, dx, cfl, eos });
    if (!(dt < Infinity)) break;
    if (t + dt > tTarget) dt = tTarget - t;

    // Stage 1
    const L1 = computeNozzleRhs(rho, mom, ene, N, dx, eos, fluxFn, A, dA);
    const rho1 = rho.slice(), mom1 = mom.slice(), ene1 = ene.slice();
    for (let i = 0; i < N; i++) {
      rho1[i + 2] = rho[i + 2] + dt * L1[0][i];
      mom1[i + 2] = mom[i + 2] + dt * L1[1][i];
      ene1[i + 2] = ene[i + 2] + dt * L1[2][i];
    }
    applyNozzleInletBC(rho1, mom1, ene1, state);
    applyNozzleExitBC(rho1, mom1, ene1, state);

    // Stage 2 (SSP)
    const L2 = computeNozzleRhs(rho1, mom1, ene1, N, dx, eos, fluxFn, A, dA);
    for (let i = 0; i < N; i++) {
      rho[i + 2] = 0.5 * rho[i + 2] + 0.5 * (rho1[i + 2] + dt * L2[0][i]);
      mom[i + 2] = 0.5 * mom[i + 2] + 0.5 * (mom1[i + 2] + dt * L2[1][i]);
      ene[i + 2] = 0.5 * ene[i + 2] + 0.5 * (ene1[i + 2] + dt * L2[2][i]);
    }
    applyNozzleInletBC(rho, mom, ene, state);
    applyNozzleExitBC(rho, mom, ene, state);

    t += dt; step++;
  }
  return Object.assign({}, state, { rho, mom, ene, t, step });
}

// ---------- Extract primitive + derived variables ----------
// Extends solver.js extract() with total pressure pt, area A and mass flow
// mdot = rho u A (kg/s when A is in m^2).
function extractNozzle(state) {
  const base = extract(state);
  const { A, thermo } = state;
  const N = base.x.length;
  const pt = new Array(N), Ar = new Array(N), mdot = new Array(N);
  for (let i = 0; i < N; i++) {
    pt[i] = gasStagnation(thermo, base.T[i], base.p[i], base.u[i]).p0;
    Ar[i] = A[i + 2];
    mdot[i] = base.rho[i] * base.u[i] * Ar[i];
  }
  base.pt = pt;
  base.A  = Ar;
  base.mdot = mdot;
  return base;
}

// ---------- Analytical: A/A* relation and inverse (calorically perfect) ----------
function areaMachRelation(M, gamma) {
  const gm1 = gamma - 1, gp1 = gamma + 1;
  const base = (2 / gp1) * (1 + 0.5 * gm1 * M * M);
  return (1 / M) * Math.pow(base, gp1 / (2 * gm1));
}

// Invert A/A* to find M. Uses bisection for robustness.
function machFromAreaRatio(AR, gamma, supersonic) {
  if (AR <= 1 + 1e-10) return 1.0;
  let lo, hi;
  if (supersonic) {
    lo = 1.0; hi = 50.0;
    while (areaMachRelation(hi, gamma) < AR) hi *= 2;
  } else {
    lo = 1e-9; hi = 1.0;
  }
  const f = M => areaMachRelation(M, gamma) - AR;
  for (let i = 0; i < 80; i++) {
    const mid = 0.5 * (lo + hi);
    const fm = f(mid);
    if (supersonic ? fm < 0 : fm > 0) lo = mid; else hi = mid;
    if (hi - lo < 1e-12) break;
  }
  return 0.5 * (lo + hi);
}

// ---------- Steady isentropic flow from a stagnation state, any gas ----------
// Along the isentrope from (T0, p0): u = sqrt(2 (h0 - h)), p = p0 exp((phi - phi0)/R),
// mass flux G = rho u. G peaks at the sonic point (G*); A/A* = G*/G.
function _isentrope(th, T0, p0) {
  const h0 = th.h(T0), phi0 = th.phi(T0);
  const at = T => {
    const u = Math.sqrt(Math.max(0, 2 * (h0 - th.h(T))));
    const p = p0 * Math.exp((th.phi(T) - phi0) / th.R);
    const rho = p / (th.R * T), a = gasSound(th, T);
    return { T, u, p, rho, a, M: u / a, G: rho * u };
  };
  const star = at(_sonicFromStagnation(th, T0).T);
  const Tmin = 1e-3 * T0;
  // Temperature on the chosen branch where the mass flux equals G (<= G*).
  const withG = (G, supersonic) => {
    if (G >= star.G) return star;
    let lo = supersonic ? Tmin : star.T, hi = supersonic ? star.T : T0;
    for (let it = 0; it < BISECT_ITERS; it++) {
      const mid = 0.5 * (lo + hi), g = at(mid).G;
      // subsonic branch: G falls as T rises; supersonic branch: G rises with T
      if ((g > G) === !supersonic) lo = mid; else hi = mid;
    }
    return at(0.5 * (lo + hi));
  };
  const withP = p => at(gasIsentropeT(th, T0, p0, p));
  return { T0, p0, at, star, withG, withP };
}

// Normal shock in a steady stream (state 1 on an isentrope with total T0):
// same h0, mass flux m = rho1 u1; solve the momentum balance for the subsonic
// downstream temperature T2 in (T*, T0). (Anderson, Modern Compressible Flow §3.6.)
function _normalShock(th, T0, s1) {
  const h0 = th.h(T0), m = s1.rho * s1.u, I1 = s1.p + m * s1.u;
  const sonicT = _sonicFromStagnation(th, T0).T;
  const state2 = T => {
    const u = Math.sqrt(Math.max(1e-30, 2 * (h0 - th.h(T))));
    const rho = m / u, p = rho * th.R * T;
    return { T, u, rho, p, I: p + m * u };
  };
  let lo = sonicT, hi = T0 * (1 - 1e-12);
  for (let it = 0; it < BISECT_ITERS; it++) {
    const mid = 0.5 * (lo + hi);
    if (state2(mid).I < I1) lo = mid; else hi = mid;     // impulse rises toward T0
  }
  const s2 = state2(0.5 * (lo + hi));
  s2.p0 = s2.p * Math.exp((th.phi(T0) - th.phi(s2.T)) / th.R);
  return s2;
}

// Exact steady quasi-1D solution on the cell centres of a nozzle state.
// mode "auto" picks the regime from the back pressure; "subsonic" /
// "supersonic" force a choked isentropic branch downstream of the throat.
function quasi1DSolution(state, mode) {
  const { N, dx, geom, p0, T0, pb, thermo: th } = state;
  const L = geom.L;
  // Throat: minimum of A(x) on a fine sampling.
  let xt = 0, At = Infinity;
  for (let k = 0; k <= AREA_SAMPLES; k++) {
    const x = (k / AREA_SAMPLES) * L, a = geom.A(x);
    if (a < At) { At = a; xt = x; }
  }
  const Ae = geom.A(L);
  const S0 = _isentrope(th, T0, p0);
  const mdotChoked = S0.star.G * At;

  const pSubExit = S0.withG(mdotChoked / Ae, false).p;   // choked, subsonic exit
  const pSupExit = S0.withG(mdotChoked / Ae, true).p;    // choked, supersonic exit
  const pShockAtExit = Ae > At * (1 + 1e-9) ? _normalShock(th, T0, S0.withG(mdotChoked / Ae, true)).p : pSubExit;

  let regime, xShock = null, mdot, S2 = null;
  if (mode === "subsonic" || mode === "supersonic") {
    regime = mode === "subsonic" ? "choked-subsonic" : "choked-supersonic";
    mdot = mdotChoked;
  } else if (pb >= pSubExit) {
    regime = "subsonic";
    mdot = S0.withP(Math.min(pb, p0)).G * Ae;              // exit sets the flow rate
  } else if (pb <= pShockAtExit) {
    regime = "supersonic";
    mdot = mdotChoked;
  } else {
    // Shock inside: exit pressure behind a shock at xs falls as xs moves aft.
    regime = "shock";
    mdot = mdotChoked;
    const exitP = xs => {
      const sh = _normalShock(th, T0, S0.withG(mdot / geom.A(xs), true));
      return _isentrope(th, T0, sh.p0).withG(mdot / Ae, false).p;
    };
    let lo = xt, hi = L;
    for (let it = 0; it < 60; it++) {
      const mid = 0.5 * (lo + hi);
      if (exitP(mid) > pb) lo = mid; else hi = mid;
    }
    xShock = 0.5 * (lo + hi);
    S2 = _isentrope(th, T0, _normalShock(th, T0, S0.withG(mdot / geom.A(xShock), true)).p0);
  }

  const out = { x: [], M: [], p: [], rho: [], T: [], a: [], u: [], pt: [], mdot: [],
                regime, xShock, xThroat: xt, Athroat: At, mdotChoked, pSubExit, pSupExit, pShockAtExit };
  for (let i = 0; i < N; i++) {
    const x = (i + 0.5) * dx, A = geom.A(x), G = mdot / A;
    const downstream = x > xt;
    let s, pt = p0;
    if (regime === "subsonic" || regime === "choked-subsonic") s = S0.withG(G, false);
    else if (regime === "shock" && x > xShock) { s = S2.withG(G, false); pt = S2.p0; }
    else s = S0.withG(G, downstream);
    out.x.push(x); out.M.push(s.M); out.p.push(s.p); out.rho.push(s.rho); out.T.push(s.T);
    out.a.push(s.a); out.u.push(s.u); out.pt.push(pt); out.mdot.push(mdot);
  }
  return out;
}

// Backward-compatible name for the analytic overlay. branch: "auto" (exact
// quasi-1D regime from p_b, the default), "subsonic" or "supersonic".
function isentropicSolution(state, branch) {
  return quasi1DSolution(state, branch || "auto");
}

// ---------- Presets ----------
const NOZZLE_PRESETS = {
  // Non-dimensional (gamma = 1.4, R = 1, A_throat = 1). For AR = 2:
  // Choked, subsonic branch: M_exit ≈ 0.306, p_exit/p0 ≈ 0.937  →  pb above this never chokes
  // Design point:            M_exit ≈ 2.197, p_exit/p0 ≈ 0.094
  // Normal shock at exit:    p2/p0 ≈ 0.513  →  pb in (0.513, 0.937) puts a shock inside the
  //                          divergent section (pb = 0.93 gives a weak one just past the throat)
  // tEnd is long enough for each regime to settle to steady state from the initial guess.
  "Subsonic CD":       { type: "CD",          AR: 2.0, p0: 1, T0: 1, pb: 0.97,  N: 300, L: 1, gamma: 1.4, cfl: 0.4,  tEnd: 60, branch: "auto" },
  "Choked CD":         { type: "CD",          AR: 2.0, p0: 1, T0: 1, pb: 0.93,  N: 300, L: 1, gamma: 1.4, cfl: 0.4,  tEnd: 30, branch: "auto" },
  "Normal Shock":      { type: "CD",          AR: 2.0, p0: 1, T0: 1, pb: 0.75,  N: 300, L: 1, gamma: 1.4, cfl: 0.35, tEnd: 15, branch: "auto" },
  "Full Supersonic CD":{ type: "CD",          AR: 2.0, p0: 1, T0: 1, pb: 0.09,  N: 300, L: 1, gamma: 1.4, cfl: 0.35, tEnd: 20, branch: "auto" },
  "Convergent":        { type: "convergent",  AR: 3.0, p0: 1, T0: 1, pb: 0.60,  N: 300, L: 1, gamma: 1.4, cfl: 0.4,  tEnd: 5,  branch: "auto" },
  "Divergent":         { type: "divergent",   AR: 2.5, p0: 1, T0: 1, pb: 0.10,  N: 300, L: 1, gamma: 1.4, cfl: 0.35, tEnd: 10, branch: "auto" },
  // SI cases (gas.js models). L in m, p in Pa, T in K, Athroat in m^2, tEnd in s.
  "Air tunnel M2 (SI)":  { type: "CD", AR: 1.6875, Athroat: 1e-3, L: 0.5, gas: "air", thermo: "thermally-perfect",
                           p0: 3e5, T0: 300, pb: 2.4e4, N: 300, cfl: 0.4, tEnd: 0.006, branch: "auto" },
  "Hot exhaust CH₄–air (SI)": { type: "CD", AR: 4, Athroat: 1e-3, L: 0.3, gas: "CH4-air products", thermo: "thermally-perfect",
                           p0: 20e5, T0: 2200, pb: 1e5, N: 300, cfl: 0.4, tEnd: 0.002, branch: "auto" }
};

if (typeof module !== "undefined" && module.exports) {
  module.exports = {
    nozzleArea, nozzleAreaDeriv, makeAreaTable, nozzleGeometry, nozzleAreaAt, buildAreaArrays,
    computeNozzleRhs,
    applyNozzleInletBC, applyNozzleExitBC,
    initNozzleState, nozzleStepUntil, extractNozzle,
    areaMachRelation, machFromAreaRatio, quasi1DSolution, isentropicSolution,
    NOZZLE_PRESETS
  };
}
