/* ============================================================
   nozzle.js  —  Quasi-1D Euler nozzle solver
   ------------------------------------------------------------
   Equations:  dU/dt + dF/dx = H
               H = -(dA/dx)/A * [rho*u, rho*u^2, u*(E+p)]
   A(x): parabolic profiles, A_throat = 1 by convention.
   Depends on solver.js: primToCons, FLUXES, hllc, computeRhs, stableDt, extract
   ============================================================ */

// ---------- Dependencies on solver.js ----------
// Browser: solver.js is loaded first and its functions are page globals, so
// `typeof primToCons` is "function" and this block is skipped (a bare `var`
// redeclaration of a global *function* leaves it untouched).
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
}
function _fluxTable() {
  return typeof FLUXES !== "undefined" ? FLUXES : _solver.FLUXES;
}

// ---------- Area profiles ----------
// AR = A_max / A_throat  (A_throat = 1 by convention)

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

// Build area + dA/dx arrays for all N+4 cells (including ghost cells).
function buildAreaArrays(N, dx, L, AR, type) {
  const tot = N + 4;
  const A = new Array(tot), dA = new Array(tot);
  for (let i = 0; i < tot; i++) {
    const xc = (i - 1.5) * dx;                      // ghost cell x positions
    const xc_clamped = Math.max(0, Math.min(L, xc));
    A[i]  = nozzleArea(xc_clamped, L, AR, type);
    dA[i] = nozzleAreaDeriv(xc_clamped, L, AR, type);
  }
  return { A, dA };
}

// ---------- Quasi-1D RHS with geometric source ----------
// Reuses computeRhs from solver.js for the flux divergence part,
// then adds the area source term on interior cells.
function computeNozzleRhs(rho, mom, ene, N, dx, gamma, fluxFn, A, dA) {
  const [r1, r2, r3] = computeRhs(rho, mom, ene, N, dx, gamma, fluxFn);
  for (let i = 0; i < N; i++) {
    const k = i + 2;
    const u = mom[k] / rho[k];
    const p = (gamma - 1) * (ene[k] - 0.5 * rho[k] * u * u);
    const src = -dA[k] / A[k];
    r1[i] += src * mom[k];          // -(dA/dx)/A * rho*u
    r2[i] += src * rho[k] * u * u; // -(dA/dx)/A * rho*u^2  (NOT +p)
    r3[i] += src * u * (ene[k] + p); // -(dA/dx)/A * u*(E+p)
  }
  return [r1, r2, r3];
}

// ---------- Inlet BC: isentropic reservoir ----------
// Subsonic inflow from a stagnation reservoir (Anderson, CFD: The Basics, Ch. 7):
// one characteristic leaves the domain, so the velocity is extrapolated from
// the first interior cell; the two incoming ones fix the stagnation state, so
// T and p follow from (T0, p0) along the steady isentrope
//     T = T0 - u^2 / (2 cp),   cp = gamma / (gamma - 1)   (R = 1)
//     p = p0 (T / T0)^(gamma / (gamma - 1))
// This holds total enthalpy and total pressure at exactly (T0, p0). (A simple-
// wave Riemann invariant u + 2a/(gamma-1) = 2a0/(gamma-1) does not: it is an
// unsteady-expansion relation and bleeds ~gamma*M*p0 of total pressure.)
function applyNozzleInletBC(rho, mom, ene, gamma, p0, T0) {
  const gm1 = gamma - 1;
  const cp  = gamma / gm1;
  // Inflow from a reservoir at rest is at most sonic: u^2 <= a*^2 = 2 gamma T0 / (gamma + 1)
  const uMax = Math.sqrt(2 * gamma * T0 / (gamma + 1));
  // No backflow into the reservoir
  const u_b = Math.min(Math.max(mom[2] / rho[2], 0), uMax);

  const T_b   = T0 - u_b * u_b / (2 * cp);
  const p_b   = p0 * Math.pow(T_b / T0, gamma / gm1);
  const rho_b = p_b / T_b;
  const [rg, mg, eg] = primToCons(rho_b, u_b, p_b, gamma);
  rho[0] = rg; mom[0] = mg; ene[0] = eg;
  rho[1] = rg; mom[1] = mg; ene[1] = eg;
}

// ---------- Exit BC: back pressure or supersonic extrapolation ----------
function applyNozzleExitBC(rho, mom, ene, gamma, pb, N) {
  const k = N + 1;  // last interior cell
  const r = rho[k], m = mom[k], e = ene[k];
  const u = m / r;
  const p = (gamma - 1) * (e - 0.5 * r * u * u);
  const a = Math.sqrt(gamma * Math.max(p, 1e-12) / Math.max(r, 1e-12));
  const Mexit = u / a;

  if (Mexit >= 1.0) {
    // Supersonic exit: all characteristics leave, extrapolate everything
    rho[N + 2] = r; mom[N + 2] = m; ene[N + 2] = e;
    rho[N + 3] = r; mom[N + 3] = m; ene[N + 3] = e;
  } else {
    // Subsonic exit: specified back pressure, zero-gradient rho and u
    const E_b = pb / (gamma - 1) + 0.5 * r * u * u;
    rho[N + 2] = r; mom[N + 2] = m; ene[N + 2] = E_b;
    rho[N + 3] = r; mom[N + 3] = m; ene[N + 3] = E_b;
  }
}

// ---------- Initialize nozzle state ----------
function initNozzleState(p) {
  const { N, L, gamma, p0, T0, pb, AR, type, flux, cfl } = p;
  const dx = L / N;
  const tot = N + 4;
  const T0_safe = Math.max(T0, 1e-6);
  const a0 = Math.sqrt(gamma * T0_safe);

  const rho = new Array(tot), mom = new Array(tot), ene = new Array(tot);
  // Linear pressure gradient from p0 at inlet to pb at exit, isentropic temperature.
  for (let i = 0; i < N; i++) {
    const frac = (i + 0.5) / N;
    const pi   = p0 + (pb - p0) * frac;
    const Ti   = T0_safe * Math.pow(Math.max(pi / p0, 1e-9), (gamma - 1) / gamma);
    const ri   = pi / Math.max(Ti, 1e-12);
    const [r, m, e] = primToCons(ri, 0.01 * a0, pi, gamma);
    rho[i + 2] = r; mom[i + 2] = m; ene[i + 2] = e;
  }
  const { A, dA } = buildAreaArrays(N, dx, L, AR, type);
  applyNozzleInletBC(rho, mom, ene, gamma, p0, T0_safe);
  applyNozzleExitBC(rho, mom, ene, gamma, pb, N);
  return {
    rho, mom, ene, t: 0, step: 0,
    dx, N, gamma, L, A, dA,
    p0, T0: T0_safe, pb, AR, type,
    flux: flux || "HLLC", cfl: cfl || 0.4
  };
}

// ---------- Step forward to tTarget (SSP-RK2 + adaptive CFL) ----------
function nozzleStepUntil(state, tTarget) {
  const { N, dx, gamma, cfl, A, dA, p0, T0, pb } = state;
  const fluxFn = _fluxTable()[state.flux] || hllc;
  let { rho, mom, ene, t, step } = state;
  let budget = 80000;

  while (t < tTarget && budget-- > 0) {
    let dt = stableDt({ rho, mom, ene, N, dx, gamma, cfl });
    if (!(dt < Infinity)) break;
    if (t + dt > tTarget) dt = tTarget - t;

    // Stage 1
    const L1 = computeNozzleRhs(rho, mom, ene, N, dx, gamma, fluxFn, A, dA);
    const rho1 = rho.slice(), mom1 = mom.slice(), ene1 = ene.slice();
    for (let i = 0; i < N; i++) {
      rho1[i + 2] = rho[i + 2] + dt * L1[0][i];
      mom1[i + 2] = mom[i + 2] + dt * L1[1][i];
      ene1[i + 2] = ene[i + 2] + dt * L1[2][i];
    }
    applyNozzleInletBC(rho1, mom1, ene1, gamma, p0, T0);
    applyNozzleExitBC(rho1, mom1, ene1, gamma, pb, N);

    // Stage 2 (SSP)
    const L2 = computeNozzleRhs(rho1, mom1, ene1, N, dx, gamma, fluxFn, A, dA);
    for (let i = 0; i < N; i++) {
      rho[i + 2] = 0.5 * rho[i + 2] + 0.5 * (rho1[i + 2] + dt * L2[0][i]);
      mom[i + 2] = 0.5 * mom[i + 2] + 0.5 * (mom1[i + 2] + dt * L2[1][i]);
      ene[i + 2] = 0.5 * ene[i + 2] + 0.5 * (ene1[i + 2] + dt * L2[2][i]);
    }
    applyNozzleInletBC(rho, mom, ene, gamma, p0, T0);
    applyNozzleExitBC(rho, mom, ene, gamma, pb, N);

    t += dt; step++;
  }
  return Object.assign({}, state, { rho, mom, ene, t, step });
}

// ---------- Extract primitive + derived variables ----------
// Extends solver.js extract() with total pressure pt and area A.
function extractNozzle(state) {
  const base = extract(state);
  const { gamma, A } = state;
  const N = base.x.length;
  const pt = new Array(N), Ar = new Array(N);
  for (let i = 0; i < N; i++) {
    const Mi = base.M[i];
    const f = 1 + 0.5 * (gamma - 1) * Mi * Mi;
    pt[i] = base.p[i] * Math.pow(Math.max(f, 1e-12), gamma / (gamma - 1));
    Ar[i] = A[i + 2];
  }
  base.pt = pt;
  base.A  = Ar;
  return base;
}

// ---------- Analytical: A/A* relation and inverse ----------
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

// ---------- Isentropic analytical solution ----------
// Returns {x, M, p, rho, T, a, pt} arrays for overlay.
// branch: "subsonic" or "supersonic" (only relevant in divergent/CD sections)
function isentropicSolution(state, branch) {
  const { N, dx, L, AR, type, gamma, p0, T0 } = state;
  const gm1 = gamma - 1;
  const rho0 = p0 / T0;
  const x = new Array(N), M = new Array(N);
  const p = new Array(N), rho = new Array(N), T = new Array(N);
  const a = new Array(N), pt = new Array(N);

  for (let i = 0; i < N; i++) {
    const xc = (i + 0.5) * dx;
    const Ax = nozzleArea(xc, L, AR, type);
    // For CD nozzle: subsonic in convergent section, chosen branch in divergent section.
    let sup = (branch === "supersonic");
    if (type === "CD") {
      // Force subsonic in the convergent half regardless of branch choice
      if (xc < 0.5 * L) sup = false;
    }
    const Mi = machFromAreaRatio(Ax, gamma, sup);
    const f  = 1 + 0.5 * gm1 * Mi * Mi;
    x[i]   = xc;
    M[i]   = Mi;
    T[i]   = T0 / f;
    p[i]   = p0 * Math.pow(T[i] / T0, gamma / gm1);
    rho[i] = rho0 * Math.pow(T[i] / T0, 1 / gm1);
    a[i]   = Math.sqrt(gamma * T[i]);
    pt[i]  = p0; // isentropic: total pressure is constant
  }
  return { x, M, p, rho, T, a, pt };
}

// ---------- Presets ----------
const NOZZLE_PRESETS = {
  // pb values are chosen to produce the labeled flow regime for AR=2, gamma=1.4.
  // Choked, subsonic branch: M_exit ≈ 0.306, p_exit/p0 ≈ 0.937  →  pb above this never chokes
  // Design point:            M_exit ≈ 2.197, p_exit/p0 ≈ 0.094
  // Normal shock at exit:    p2/p0 ≈ 0.513  →  pb in (0.513, 0.937) puts a shock inside the
  //                          divergent section (pb = 0.93 gives a weak one just past the throat)
  // tEnd is long enough for each regime to settle to steady state from the initial guess.
  "Subsonic CD":       { type: "CD",          AR: 2.0, p0: 1, T0: 1, pb: 0.97,  N: 300, L: 1, gamma: 1.4, cfl: 0.4,  tEnd: 60, branch: "subsonic"    },
  "Choked CD":         { type: "CD",          AR: 2.0, p0: 1, T0: 1, pb: 0.93,  N: 300, L: 1, gamma: 1.4, cfl: 0.4,  tEnd: 30, branch: "subsonic"    },
  "Normal Shock":      { type: "CD",          AR: 2.0, p0: 1, T0: 1, pb: 0.75,  N: 300, L: 1, gamma: 1.4, cfl: 0.35, tEnd: 15, branch: "subsonic"    },
  "Full Supersonic CD":{ type: "CD",          AR: 2.0, p0: 1, T0: 1, pb: 0.09,  N: 300, L: 1, gamma: 1.4, cfl: 0.35, tEnd: 20, branch: "supersonic"  },
  "Convergent":        { type: "convergent",  AR: 3.0, p0: 1, T0: 1, pb: 0.60,  N: 300, L: 1, gamma: 1.4, cfl: 0.4,  tEnd: 5,  branch: "subsonic"    },
  "Divergent":         { type: "divergent",   AR: 2.5, p0: 1, T0: 1, pb: 0.10,  N: 300, L: 1, gamma: 1.4, cfl: 0.35, tEnd: 10, branch: "supersonic"  }
};

if (typeof module !== "undefined" && module.exports) {
  module.exports = {
    nozzleArea, nozzleAreaDeriv, buildAreaArrays,
    computeNozzleRhs,
    applyNozzleInletBC, applyNozzleExitBC,
    initNozzleState, nozzleStepUntil, extractNozzle,
    areaMachRelation, machFromAreaRatio, isentropicSolution,
    NOZZLE_PRESETS
  };
}
