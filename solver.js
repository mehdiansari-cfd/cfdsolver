/* ============================================================
   solver.js
   1D Euler shock tube solver — pure math layer (no DOM, no UI)
   ------------------------------------------------------------
   Methods : HLLC / HLL / Roe / Rusanov flux + MUSCL (minmod)
             + SSP-RK2 + adaptive CFL
   Equations: dU/dt + dF(U)/dx = 0,  U = [rho, rho*u, E (, rho*Y)]
             E = rho e(T) + (1/2) rho u^2,   p = rho R T

   Equation of state. Every function that used to take `gamma` takes an
   "eos" instead, which can be:
     - a number gamma        legacy non-dimensional ideal gas, R = 1
                             (bit-identical to the original solver)
     - gasEos(gas)           one gas from gas.js, SI units, calorically
                             or thermally perfect (gamma = gamma(T))
     - binaryEos(gasA, gasB) two gases (driver / driven) separated by a
                             contact; a mass fraction Y of gasA rides along
                             as a passive scalar, U gets a 4th row rho*Y

   Public API (attached to globals once this script is loaded):
     - idealEos(gamma), gasEos(gas), binaryEos(gasA, gasB)
     - primToCons(rho, u, p, gamma)
     - fluxFromCons(rho, mom, ene, gamma)
     - hllc(UL, UR, eos) / hll / roe / rusanov      UL = [rho, mom, E (, rhoY)]
     - FLUXES                                     name -> flux function
     - minmod(a, b)
     - applyBcs(rho, mom, ene, bcType, rY?)       'transmissive' | 'reflective' | 'periodic'
     - computeRhs(rho, mom, ene, N, dx, eos, fluxFn, rY?)
     - stableDt(state)                             CFL-limited time step
     - initState(params)                          params.flux ('HLLC' | 'HLL' | 'Roe' | 'Rusanov')
                                                  params.gas, or params.gasLeft + params.gasRight
     - stepUntil(state, tTarget)
     - extract(state)                              -> { x, rho, u, p, e, M, T, s, a (, gamma, Y) }
     - exactRiemann(left, right, gamma, x0, t, xArr) -> { x, rho, u, p } | null
     - riemannGasSolve(gasL, left, gasR, right)    real-gas exact Riemann problem
     - exactRiemannGas(gasL, left, gasR, right, x0, t, xArr)
     - localVelocity(state, x)
     - PRESETS                                     classic shock tube setups
   ============================================================ */

// ---------- Gas library (gas.js) ----------
// Browser: gas.js loads first and these are page globals, so this block is
// skipped. Node: pull them in. Only gas.js *function declarations* may be
// aliased with var (its constants are `const` and would collide).
if (typeof gasTfromE === "undefined" && typeof require === "function") {
  var _gaslib       = require("./gas");
  var gasE          = _gaslib.gasE;
  var gasGamma      = _gaslib.gasGamma;
  var gasSound      = _gaslib.gasSound;
  var gasTfromE     = _gaslib.gasTfromE;
  var gasIsentropeT = _gaslib.gasIsentropeT;
  var mixGases      = _gaslib.mixGases;
  var gasTfromH     = _gaslib.gasTfromH;
}

const P_FLOOR = 1e-12, RHO_FLOOR = 1e-12;   // floors for sound speed (CLAUDE.md §8.3)
const EOS_TREF = 298.15;                    // gas.js reference temperature [K]
const EOS_PREF = 101325;                    // entropy reference pressure [Pa]
const EOS_T_GUESS_MIN = 20;                 // lower clamp for the Newton start [K]

// ---------- Equation-of-state adapters ----------
// prim(rho, mom, ene, rY) -> { u, p, a, T, gamma, R }
// cons(rho, u, p, Y)      -> [rho, rho u, E (, rho Y)]
const _idealCache = {};
function idealEos(gamma) {
  if (_idealCache[gamma]) return _idealCache[gamma];
  const eos = {
    kind: "ideal", gamma, binary: false, R: 1,
    prim: function (rho, mom, ene) {
      const u = mom / rho;
      const p = (gamma - 1) * (ene - 0.5 * rho * u * u);
      const a = Math.sqrt(gamma * Math.max(p, P_FLOOR) / Math.max(rho, RHO_FLOOR));
      return { u, p, a, T: p / Math.max(rho, RHO_FLOOR), gamma, R: 1 };
    },
    cons: function (rho, u, p) { return [rho, rho * u, p / (gamma - 1) + 0.5 * rho * u * u]; },
    gasAt: function () { return null; }
  };
  _idealCache[gamma] = eos;
  return eos;
}

// Temperature from specific internal energy, Newton started from a
// calorically perfect estimate so it converges in a few iterations.
function _tFromE(gas, e) {
  const cvRef = gas.cp(EOS_TREF) - gas.R;
  const guess = Math.max(EOS_T_GUESS_MIN, EOS_TREF + (e - gasE(gas, EOS_TREF)) / cvRef);
  return gasTfromE(gas, e, guess);
}
function _gasPrim(gas, rho, mom, ene) {
  const u = mom / rho;
  const T = _tFromE(gas, ene / rho - 0.5 * u * u);
  const cp = gas.cp(T), gamma = cp / (cp - gas.R);     // a^2 = gamma R T
  return { u, p: rho * gas.R * T, a: Math.sqrt(gamma * gas.R * T), T, gamma, R: gas.R };
}

function gasEos(gas) {
  return {
    kind: "gas", gas, binary: false, R: gas.R,
    prim: function (rho, mom, ene) { return _gasPrim(gas, rho, mom, ene); },
    cons: function (rho, u, p) {
      const T = p / (rho * gas.R);
      return [rho, rho * u, rho * (gasE(gas, T) + 0.5 * u * u)];
    },
    gasAt: function () { return gas; }
  };
}

// Two gases; Y = mass fraction of gasA (the left / driver gas).
function binaryEos(gasA, gasB) {
  const at = Y => mixGases(gasA, gasB, Math.min(1, Math.max(0, Y)));
  return {
    kind: "binary", gasA, gasB, binary: true,
    prim: function (rho, mom, ene, rY) { return _gasPrim(at(rY / rho), rho, mom, ene); },
    cons: function (rho, u, p, Y) {
      const g = at(Y), T = p / (rho * g.R);
      return [rho, rho * u, rho * (gasE(g, T) + 0.5 * u * u), rho * Y];
    },
    gasAt: at
  };
}

function _eos(g) { return typeof g === "number" ? idealEos(g) : g; }

// ---------- Variable conversions (legacy ideal-gas forms) ----------
function primToCons(rho, u, p, gamma) {
  return [rho, rho * u, p / (gamma - 1) + 0.5 * rho * u * u];
}

function fluxFromCons(rho, mom, ene, gamma) {
  const u = mom / rho;
  const p = (gamma - 1) * (ene - 0.5 * rho * u * u);
  return [mom, mom * u + p, u * (ene + p)];
}

// Physical flux F(U) given U and its primitive state q; a passive scalar
// rho*Y (4th row) is carried with the mass flux.
function _physFlux(U, q) {
  const F = [U[1], U[1] * q.u + q.p, q.u * (U[2] + q.p)];
  if (U.length > 3) F.push(U[1] * U[3] / U[0]);
  return F;
}

// ---------- HLLC approximate Riemann solver ----------
function hllc(L, R, gamma) {
  const eos = _eos(gamma);
  const qL = eos.prim(L[0], L[1], L[2], L[3]), qR = eos.prim(R[0], R[1], R[2], R[3]);
  const rhoL = L[0], rhoR = R[0];
  const uL = qL.u, uR = qR.u, pL = qL.p, pR = qR.p, aL = qL.a, aR = qR.a;

  // Davis wave-speed estimates
  const SL = Math.min(uL - aL, uR - aR);
  const SR = Math.max(uL + aL, uR + aR);

  // Contact wave speed (Toro Eq. 10.37)
  const sStar = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) /
                (rhoL * (SL - uL) - rhoR * (SR - uR));

  const FL = _physFlux(L, qL);
  const FR = _physFlux(R, qR);

  if (SL >= 0) return FL;
  if (SR <= 0) return FR;

  // Star state (Toro Eq. 10.39); the scalar Y is constant across S_K.
  const U = sStar >= 0 ? L : R, q = sStar >= 0 ? qL : qR, F = sStar >= 0 ? FL : FR;
  const S = sStar >= 0 ? SL : SR;
  const f = U[0] * (S - q.u) / (S - sStar);
  const rS = f, mS = f * sStar;
  const eS = f * (U[2] / U[0] + (sStar - q.u) * (sStar + q.p / (U[0] * (S - q.u))));
  const out = [F[0] + S * (rS - U[0]), F[1] + S * (mS - U[1]), F[2] + S * (eS - U[2])];
  if (U.length > 3) out.push(F[3] + S * (f * U[3] / U[0] - U[3]));
  return out;
}

// ---------- HLL approximate Riemann solver ----------
// No contact wave; just averaged star state between SL and SR.
function hll(L, R, gamma) {
  const eos = _eos(gamma);
  const qL = eos.prim(L[0], L[1], L[2], L[3]), qR = eos.prim(R[0], R[1], R[2], R[3]);
  const SL = Math.min(qL.u - qL.a, qR.u - qR.a);
  const SR = Math.max(qL.u + qL.a, qR.u + qR.a);
  const FL = _physFlux(L, qL);
  const FR = _physFlux(R, qR);
  if (SL >= 0) return FL;
  if (SR <= 0) return FR;
  const denom = SR - SL;
  return FL.map((_, k) => (SR * FL[k] - SL * FR[k] + SL * SR * (R[k] - L[k])) / denom);
}

// ---------- Rusanov (local Lax-Friedrichs) flux ----------
function rusanov(L, R, gamma) {
  const eos = _eos(gamma);
  const qL = eos.prim(L[0], L[1], L[2], L[3]), qR = eos.prim(R[0], R[1], R[2], R[3]);
  const s  = Math.max(Math.abs(qL.u) + qL.a, Math.abs(qR.u) + qR.a);
  const FL = _physFlux(L, qL);
  const FR = _physFlux(R, qR);
  return FL.map((_, k) => 0.5 * (FL[k] + FR[k]) - 0.5 * s * (R[k] - L[k]));
}

// ---------- Roe-Pike flux with Harten entropy fix ----------
// For a gas-model EOS (SI, enthalpy with a reference offset) two things change:
//  - the Roe-averaged sound speed is the Roe-weighted mean of aL^2 and aR^2,
//    an approximate linearisation (Glaister 1988 gives the exact one);
//  - the entropy wave carries dE/drho|_p = u^2/2 + chi with chi = h - cp T,
//    and, for two gases, the species wave carries de/dY|_{rho,p}
//    = (hA - hB) - cp T (RA - RB)/R. Both vanish for the legacy ideal gas.
// HLLC remains the recommended flux for real gases.
function roe(L, R, gamma) {
  const eos = _eos(gamma);
  const qL = eos.prim(L[0], L[1], L[2], L[3]), qR = eos.prim(R[0], R[1], R[2], R[3]);
  const rhoL = L[0], rhoR = R[0], eneL = L[2], eneR = R[2];
  const uL = qL.u, uR = qR.u, pL = qL.p, pR = qR.p;
  const HL = (eneL + pL) / rhoL;
  const HR = (eneR + pR) / rhoR;

  // Roe averages
  const sL = Math.sqrt(Math.max(rhoL, RHO_FLOOR));
  const sR = Math.sqrt(Math.max(rhoR, RHO_FLOOR));
  const sum = sL + sR;
  const uH = (sL * uL + sR * uR) / sum;
  const HH = (sL * HL + sR * HR) / sum;
  const rhoH = sL * sR;
  const a2 = eos.kind === "ideal" ? (eos.gamma - 1) * (HH - 0.5 * uH * uH)
                                  : (sL * qL.a * qL.a + sR * qR.a * qR.a) / sum;
  const aH = Math.sqrt(Math.max(a2, P_FLOOR));

  const drho = rhoR - rhoL;
  const du   = uR - uL;
  const dp   = pR - pL;

  // Wave strengths
  const alpha1 = (dp - rhoH * aH * du) / (2 * aH * aH);   // left acoustic
  const alpha2 =  drho - dp / (aH * aH);                  // entropy/contact
  const alpha3 = (dp + rhoH * aH * du) / (2 * aH * aH);   // right acoustic

  // Eigenvalues (with Harten fix on the genuinely nonlinear fields)
  let lam1 = uH - aH, lam2 = uH, lam3 = uH + aH;
  const eps = 0.1 * aH;
  if (Math.abs(lam1) < eps) lam1 = (lam1 * lam1 + eps * eps) / (2 * eps);
  if (Math.abs(lam3) < eps) lam3 = (lam3 * lam3 + eps * eps) / (2 * eps);
  const al1 = Math.abs(lam1), al2 = Math.abs(lam2), al3 = Math.abs(lam3);

  // chi = h - cp T: zero for the legacy ideal gas (h = cp T)
  let chiH = 0;
  if (eos.kind !== "ideal") {
    const chi = (q, H) => H - 0.5 * q.u * q.u - q.gamma * q.R / (q.gamma - 1) * q.T;
    chiH = (sL * chi(qL, HL) + sR * chi(qR, HR)) / sum;
  }

  // Right eigenvectors
  const r1 = [1, uH - aH, HH - uH * aH];
  const r2 = [1, uH,      0.5 * uH * uH + chiH];
  const r3 = [1, uH + aH, HH + uH * aH];

  const FL = _physFlux(L, qL);
  const FR = _physFlux(R, qR);

  const out = [0, 0, 0];
  for (let k = 0; k < 3; k++) {
    out[k] = 0.5 * (FL[k] + FR[k])
           - 0.5 * (al1 * alpha1 * r1[k] + al2 * alpha2 * r2[k] + al3 * alpha3 * r3[k]);
  }
  if (L.length > 3) {
    // Passive scalar: every wave carries the Roe-averaged Y, plus a scalar
    // wave of strength rhoH * dY moving with uH.
    const YL = L[3] / rhoL, YR = R[3] / rhoR, YH = (sL * YL + sR * YR) / sum;
    const alpha4 = rhoH * (YR - YL);
    // energy carried by the species wave: de/dY at fixed rho, p (Roe-averaged state)
    const A = eos.gasA, B = eos.gasB, TH = (sL * qL.T + sR * qR.T) / sum, g = eos.gasAt(YH);
    const deY = (A.h(TH) - B.h(TH)) - g.cp(TH) * TH * (A.R - B.R) / g.R;
    out[2] -= 0.5 * al2 * alpha4 * deY;
    out.push(0.5 * (FL[3] + FR[3])
             - 0.5 * (YH * (al1 * alpha1 + al2 * alpha2 + al3 * alpha3) + al2 * alpha4));
  }
  return out;
}

// ---------- Flux dispatch table ----------
const FLUXES = { HLLC: hllc, HLL: hll, Roe: roe, Rusanov: rusanov };

// ---------- Slope limiter ----------
function minmod(a, b) {
  if (a * b <= 0) return 0;
  return Math.sign(a) * Math.min(Math.abs(a), Math.abs(b));
}

// ---------- Boundary conditions ----------
// rY (optional) is the passive scalar rho*Y; it is copied like rho.
function applyBcs(rho, mom, ene, bc, rY) {
  const n = rho.length;
  const arrs = rY ? [rho, mom, ene, rY] : [rho, mom, ene];
  if (bc === "periodic") {
    // ghost cells copy the interior cells one period away
    for (const a of arrs) { a[0] = a[n - 4]; a[1] = a[n - 3]; a[n - 2] = a[2]; a[n - 1] = a[3]; }
  } else if (bc === "reflective") {
    rho[1] = rho[2]; mom[1] = -mom[2]; ene[1] = ene[2];
    rho[0] = rho[3]; mom[0] = -mom[3]; ene[0] = ene[3];
    rho[n - 2] = rho[n - 3]; mom[n - 2] = -mom[n - 3]; ene[n - 2] = ene[n - 3];
    rho[n - 1] = rho[n - 4]; mom[n - 1] = -mom[n - 4]; ene[n - 1] = ene[n - 4];
    if (rY) { rY[1] = rY[2]; rY[0] = rY[3]; rY[n - 2] = rY[n - 3]; rY[n - 1] = rY[n - 4]; }
  } else {
    // transmissive (zero-gradient)
    rho[1] = rho[2]; mom[1] = mom[2]; ene[1] = ene[2];
    rho[0] = rho[2]; mom[0] = mom[2]; ene[0] = ene[2];
    rho[n - 2] = rho[n - 3]; mom[n - 2] = mom[n - 3]; ene[n - 2] = ene[n - 3];
    rho[n - 1] = rho[n - 3]; mom[n - 1] = mom[n - 3]; ene[n - 1] = ene[n - 3];
    if (rY) { rY[1] = rY[2]; rY[0] = rY[2]; rY[n - 2] = rY[n - 3]; rY[n - 1] = rY[n - 3]; }
  }
}

// ---------- Spatial operator: dU/dt = -dF/dx (interior cells only) ----------
// fluxFn defaults to HLLC for backward compatibility. With rY (two-gas runs)
// the scalar is reconstructed and fluxed like the other rows and a 4th RHS
// row is returned.
function computeRhs(rho, mom, ene, N, dx, gamma, fluxFn, rY) {
  const flux = fluxFn || hllc;
  const nVar = rY ? 4 : 3;
  const Q = rY ? [rho, mom, ene, rY] : [rho, mom, ene];
  // Slopes in cells 1..N+2 (N+2 of them)
  const sl = [];
  for (let v = 0; v < nVar; v++) {
    const q = Q[v], s = new Array(N + 2);
    for (let i = 0; i < N + 2; i++) {
      const k = i + 1;
      s[i] = minmod(q[k] - q[k - 1], q[k + 1] - q[k]);
    }
    sl.push(s);
  }
  // Fluxes at N+1 interfaces (MUSCL-reconstructed left/right states)
  const F = [];
  for (let v = 0; v < nVar; v++) F.push(new Array(N + 1));
  for (let i = 0; i <= N; i++) {
    const kL = i + 1, kR = i + 2;
    const UL = new Array(nVar), UR = new Array(nVar);
    for (let v = 0; v < nVar; v++) {
      UL[v] = Q[v][kL] + 0.5 * sl[v][i];
      UR[v] = Q[v][kR] - 0.5 * sl[v][i + 1];
    }
    const f = flux(UL, UR, gamma);
    for (let v = 0; v < nVar; v++) F[v][i] = f[v];
  }
  // Divergence in interior cells
  const out = [];
  for (let v = 0; v < nVar; v++) {
    const r = new Array(N), Fv = F[v];
    for (let i = 0; i < N; i++) r[i] = -(Fv[i + 1] - Fv[i]) / dx;
    out.push(r);
  }
  return out;
}

// ---------- Initial state from a parameter object ----------
// Legacy:  { gamma, leftState: {rho, u, p}, ... }                (non-dim, R = 1)
// One gas: { gas, leftState: {p, T, u} or {rho, u, p}, ... }      (SI)
// Two gas: { gasLeft, gasRight, ... } — left of the diaphragm is gasLeft.
function _eosFromParams(p) {
  if (p.gasLeft && p.gasRight) return binaryEos(p.gasLeft, p.gasRight);
  if (p.gas) return gasEos(p.gas);
  return idealEos(p.gamma);
}
function _densityOf(s, R) { return s.rho != null ? s.rho : s.p / (R * s.T); }

function initState(p) {
  const { N, L, gamma, diaphragm, leftState, rightState, bcType } = p;
  const eos = _eosFromParams(p);
  const dx = L / N;
  const tot = N + 4; // 2 ghost cells per side
  const rho = new Array(tot).fill(0), mom = new Array(tot).fill(0), ene = new Array(tot).fill(0);
  const rY = eos.binary ? new Array(tot).fill(0) : null;
  const Rl = eos.kind === "ideal" ? 1 : eos.binary ? p.gasLeft.R : p.gas.R;
  const Rr = eos.kind === "ideal" ? 1 : eos.binary ? p.gasRight.R : p.gas.R;
  for (let i = 0; i < N; i++) {
    const x = (i + 0.5) * dx;
    const left = x < diaphragm;
    const s = left ? leftState : rightState;
    const U = eos.cons(_densityOf(s, left ? Rl : Rr), s.u, s.p, left ? 1 : 0);
    rho[i + 2] = U[0]; mom[i + 2] = U[1]; ene[i + 2] = U[2];
    if (rY) rY[i + 2] = U[3];
  }
  applyBcs(rho, mom, ene, bcType, rY);
  const fluxName = p.flux || "HLLC";
  const snap = s => Object.assign({}, s);
  return {
    rho, mom, ene, rY, t: 0, step: 0, dx, N, gamma, eos, bcType,
    cfl: p.cfl, L, flux: fluxName,
    // Snapshot the initial states + diaphragm so we can compute the
    // exact Riemann solution later without depending on UI params.
    leftState: snap(leftState), rightState: snap(rightState),
    diaphragm
  };
}

// ---------- Stable time step from the CFL condition ----------
// dt = CFL * dx / max_i(|u_i| + a_i) over the interior cells. Recomputed every
// step, since wave speeds change as shocks form (CLAUDE.md §8.6). Works for any
// state with { rho, mom, ene, N, dx, gamma, cfl } (+ eos, rY) — the nozzle reuses it.
function stableDt(state) {
  const { rho, mom, ene, N, dx, gamma, cfl, eos, rY } = state;
  let maxSpeed = 0;
  if (!eos || eos.kind === "ideal") {
    const g = eos ? eos.gamma : gamma;
    for (let i = 2; i < N + 2; i++) {
      const u = mom[i] / rho[i];
      const p = (g - 1) * (ene[i] - 0.5 * rho[i] * u * u);
      const a = Math.sqrt(g * Math.max(p, 1e-12) / Math.max(rho[i], 1e-12));
      const s = Math.abs(u) + a;
      if (s > maxSpeed) maxSpeed = s;
    }
  } else {
    for (let i = 2; i < N + 2; i++) {
      const q = eos.prim(rho[i], mom[i], ene[i], rY ? rY[i] : undefined);
      const s = Math.abs(q.u) + q.a;
      if (s > maxSpeed) maxSpeed = s;
    }
  }
  return cfl * dx / maxSpeed;
}

// ---------- Step forward in time to tTarget (SSP-RK2 + adaptive CFL) ----------
function stepUntil(state, tTarget) {
  const { N, dx, gamma, bcType, cfl, rY } = state;
  const eos = state.eos || gamma;
  const fluxFn = FLUXES[state.flux] || hllc;
  let { rho, mom, ene, t, step } = state;
  const Q = rY ? [rho, mom, ene, rY] : [rho, mom, ene];
  let budget = 5000; // safety cap per call
  while (t < tTarget && budget-- > 0) {
    let dt = stableDt({ rho, mom, ene, N, dx, gamma, cfl, eos: state.eos, rY });
    if (t + dt > tTarget) dt = tTarget - t;

    // Stage 1
    const L1 = computeRhs(rho, mom, ene, N, dx, eos, fluxFn, rY);
    const Q1 = Q.map(q => q.slice());
    for (let v = 0; v < Q.length; v++) {
      for (let i = 0; i < N; i++) Q1[v][i + 2] = Q[v][i + 2] + dt * L1[v][i];
    }
    applyBcs(Q1[0], Q1[1], Q1[2], bcType, Q1[3]);

    // Stage 2 (SSP)
    const L2 = computeRhs(Q1[0], Q1[1], Q1[2], N, dx, eos, fluxFn, Q1[3]);
    for (let v = 0; v < Q.length; v++) {
      for (let i = 0; i < N; i++) {
        Q[v][i + 2] = 0.5 * Q[v][i + 2] + 0.5 * (Q1[v][i + 2] + dt * L2[v][i]);
      }
    }
    applyBcs(rho, mom, ene, bcType, rY);

    t += dt; step++;
  }
  return Object.assign({}, state, { rho, mom, ene, rY, t, step });
}

// ---------- Extract primitive + derived variables on interior cells ----------
// Legacy (non-dim): rho, u, p, e (internal energy), M (Mach), T (= p/rho),
//   s (= ln(p/rho^gamma)), a (sound speed).
// Gas (SI): T [K], e = e(T) - e(298.15 K) [J/kg], s = phi(T) - R ln(p/101325)
//   [J/(kg K)], plus the local gamma and, for two-gas runs, Y.
function extract(state) {
  const { rho, mom, ene, N, dx, gamma, eos, rY } = state;
  const x = new Array(N), R = new Array(N), U = new Array(N), P = new Array(N), E = new Array(N);
  const M = new Array(N), T = new Array(N), S = new Array(N), A = new Array(N);
  if (!eos || eos.kind === "ideal") {
    const g = eos ? eos.gamma : gamma;
    for (let i = 0; i < N; i++) {
      x[i] = (i + 0.5) * dx;
      R[i] = rho[i + 2];
      U[i] = mom[i + 2] / R[i];
      P[i] = (g - 1) * (ene[i + 2] - 0.5 * R[i] * U[i] * U[i]);
      E[i] = P[i] / ((g - 1) * R[i]);
      const aa = Math.sqrt(g * Math.max(P[i], 1e-12) / Math.max(R[i], 1e-12));
      A[i] = aa;
      M[i] = U[i] / aa;
      T[i] = P[i] / Math.max(R[i], 1e-12); // RT (R=1 implicit, non-dim)
      S[i] = Math.log(Math.max(P[i], 1e-12) / Math.pow(Math.max(R[i], 1e-12), g));
    }
    return { x, rho: R, u: U, p: P, e: E, M, T, s: S, a: A };
  }
  const G = new Array(N), Y = rY ? new Array(N) : null;
  for (let i = 0; i < N; i++) {
    const k = i + 2;
    const q = eos.prim(rho[k], mom[k], ene[k], rY ? rY[k] : undefined);
    const gas = eos.gasAt(rY ? rY[k] / rho[k] : 1);
    x[i] = (i + 0.5) * dx;
    R[i] = rho[k]; U[i] = q.u; P[i] = q.p; T[i] = q.T; A[i] = q.a; G[i] = q.gamma;
    M[i] = q.u / q.a;
    E[i] = gasE(gas, q.T) - gasE(gas, EOS_TREF);
    S[i] = gas.phi(q.T) - gas.R * Math.log(Math.max(q.p, P_FLOOR) / EOS_PREF);
    if (Y) Y[i] = Math.min(1, Math.max(0, rY[k] / rho[k]));
  }
  const out = { x, rho: R, u: U, p: P, e: E, M, T, s: S, a: A, gamma: G };
  if (Y) out.Y = Y;
  return out;
}

// ---------- Exact Riemann solver for ideal gas ----------
// Returns the analytical solution sampled at the cell centers xArr at time t.
// Returns null if a vacuum forms (the iterative procedure isn't valid there).
function exactRiemann(left, right, gamma, x0, t, xArr) {
  const rhoL = left.rho,  uL = left.u,  pL = left.p;
  const rhoR = right.rho, uR = right.u, pR = right.p;
  const aL = Math.sqrt(gamma * pL / rhoL);
  const aR = Math.sqrt(gamma * pR / rhoR);
  const g  = gamma, gm1 = g - 1, gp1 = g + 1;

  // Trivial t=0: return the step function exactly.
  if (!(t > 0)) {
    const N = xArr.length;
    const rho = new Array(N), u = new Array(N), p = new Array(N);
    for (let i = 0; i < N; i++) {
      if (xArr[i] < x0) { rho[i] = rhoL; u[i] = uL; p[i] = pL; }
      else              { rho[i] = rhoR; u[i] = uR; p[i] = pR; }
    }
    return { x: xArr, rho, u, p };
  }

  // Vacuum check (Toro Eq. 4.40)
  if ((2 * aL / gm1) + (2 * aR / gm1) <= (uR - uL)) return null;

  // f and f' for each side (shock vs rarefaction branches).
  const f  = (p, rho_s, p_s, a_s) => {
    if (p > p_s) {
      const A = 2 / (gp1 * rho_s), B = (gm1 / gp1) * p_s;
      return (p - p_s) * Math.sqrt(A / (p + B));
    }
    return (2 * a_s / gm1) * (Math.pow(p / p_s, gm1 / (2 * g)) - 1);
  };
  const fp = (p, rho_s, p_s, a_s) => {
    if (p > p_s) {
      const A = 2 / (gp1 * rho_s), B = (gm1 / gp1) * p_s;
      const sq = Math.sqrt(A / (p + B));
      return sq * (1 - (p - p_s) / (2 * (p + B)));
    }
    return (1 / (rho_s * a_s)) * Math.pow(p / p_s, -gp1 / (2 * g));
  };

  // Newton iteration on p*  (Toro Eq. 4.5, two-rarefaction initial guess)
  let p = 0.5 * (pL + pR) - 0.125 * (uR - uL) * (rhoL + rhoR) * (aL + aR);
  p = Math.max(p, 1e-8 * Math.min(pL, pR));
  for (let iter = 0; iter < 60; iter++) {
    const F  = f(p, rhoL, pL, aL) + f(p, rhoR, pR, aR) + (uR - uL);
    const Fp = fp(p, rhoL, pL, aL) + fp(p, rhoR, pR, aR);
    let pNew = p - F / Fp;
    if (pNew < 0) pNew = 0.5 * p;
    if (Math.abs(2 * (pNew - p) / (pNew + p)) < 1e-9) { p = pNew; break; }
    p = pNew;
  }
  const pStar = p;
  const uStar = 0.5 * (uL + uR) + 0.5 * (f(pStar, rhoR, pR, aR) - f(pStar, rhoL, pL, aL));

  // Sample the self-similar solution at each x_i.
  const N = xArr.length;
  const rhoOut = new Array(N), uOut = new Array(N), pOut = new Array(N);
  for (let i = 0; i < N; i++) {
    const sVar = (xArr[i] - x0) / t;
    let R, U, P;
    if (sVar < uStar) {
      // ---- LEFT WAVE ----
      if (pStar > pL) {
        // Left shock
        const SL = uL - aL * Math.sqrt(gp1 / (2 * g) * pStar / pL + gm1 / (2 * g));
        if (sVar < SL) { R = rhoL; U = uL; P = pL; }
        else {
          const r = (pStar / pL + gm1 / gp1) / ((gm1 / gp1) * (pStar / pL) + 1);
          R = rhoL * r; U = uStar; P = pStar;
        }
      } else {
        // Left rarefaction fan
        const aSL = aL * Math.pow(pStar / pL, gm1 / (2 * g));
        const SHL = uL - aL;
        const STL = uStar - aSL;
        if (sVar < SHL) { R = rhoL; U = uL; P = pL; }
        else if (sVar > STL) {
          R = rhoL * Math.pow(pStar / pL, 1 / g); U = uStar; P = pStar;
        } else {
          const af = (2 / gp1) * (aL + (gm1 / 2) * (uL - sVar));
          U = (2 / gp1) * (aL + (gm1 / 2) * uL + sVar);
          R = rhoL * Math.pow(af / aL, 2 / gm1);
          P = pL   * Math.pow(af / aL, 2 * g / gm1);
        }
      }
    } else {
      // ---- RIGHT WAVE ----
      if (pStar > pR) {
        // Right shock
        const SR = uR + aR * Math.sqrt(gp1 / (2 * g) * pStar / pR + gm1 / (2 * g));
        if (sVar > SR) { R = rhoR; U = uR; P = pR; }
        else {
          const r = (pStar / pR + gm1 / gp1) / ((gm1 / gp1) * (pStar / pR) + 1);
          R = rhoR * r; U = uStar; P = pStar;
        }
      } else {
        // Right rarefaction fan
        const aSR = aR * Math.pow(pStar / pR, gm1 / (2 * g));
        const SHR = uR + aR;
        const STR = uStar + aSR;
        if (sVar > SHR) { R = rhoR; U = uR; P = pR; }
        else if (sVar < STR) {
          R = rhoR * Math.pow(pStar / pR, 1 / g); U = uStar; P = pStar;
        } else {
          const af = (2 / gp1) * (aR - (gm1 / 2) * (uR - sVar));
          U = (2 / gp1) * (-aR + (gm1 / 2) * uR + sVar);
          R = rhoR * Math.pow(af / aR, 2 / gm1);
          P = pR   * Math.pow(af / aR, 2 * g / gm1);
        }
      }
    }
    rhoOut[i] = R; uOut[i] = U; pOut[i] = P;
  }
  return { x: xArr, rho: rhoOut, u: uOut, p: pOut };
}

// ---------- Exact Riemann solver for real (calorically / thermally perfect) gases ----------
// Same wave structure as Toro Ch. 4, with the ideal-gas closed forms replaced
// by their general ideal-gas (p = rho R T, cp = cp(T)) counterparts:
//   shock:        Rankine–Hugoniot  h2 - h1 = (p2 - p1)(v1 + v2)/2  solved for T2,
//                 du = sqrt((p2 - p1)(v1 - v2))
//   rarefaction:  isentrope phi(T) - R ln p = const,
//                 du = ∫ dp / (rho a)  (Riemann invariant), integrated in ln p
// Left and right may be different gases (driver / driven).
// States: { p, u, T } or { p, u, rho }.
const RAREF_NODES = 64;      // Simpson intervals for the rarefaction integral (even)
const FAN_TABLE   = 400;     // pressure samples tabulating each rarefaction fan
const PSTAR_ITERS = 200;     // bisection cap on p* (log-space)
const PSTAR_RTOL  = 1e-12;
const RH_ITERS    = 60;

function _side(gas, s) {
  const T = s.T != null ? s.T : s.p / (s.rho * gas.R);
  const rho = s.rho != null ? s.rho : s.p / (gas.R * T);
  return { gas, p: s.p, u: s.u, T, rho, a: gasSound(gas, T) };
}

// State behind a shock of pressure p2 running into state K.
function _shock(K, p2) {
  const g = K.gas, v1 = 1 / K.rho, dp = p2 - K.p, h1 = g.h(K.T);
  // calorically perfect Rankine–Hugoniot as the starting guess (Anderson Eq. 3.58)
  const gm = gasGamma(g, K.T), c = (gm + 1) / (gm - 1), r = p2 / K.p;
  let T = K.T * r * (c + r) / (1 + c * r);
  for (let it = 0; it < RH_ITERS; it++) {
    const f  = g.h(T) - h1 - 0.5 * dp * (v1 + g.R * T / p2);
    const fp = g.cp(T) - 0.5 * dp * g.R / p2;          // > 0 since cp > R
    const Tn = T - f / fp;
    if (Math.abs(Tn - T) < 1e-12 * T) { T = Tn; break; }
    T = Tn > 0 ? Tn : 0.5 * T;
  }
  const rho2 = p2 / (g.R * T);
  const m = Math.sqrt(dp / (v1 - 1 / rho2));           // mass flux through the shock
  return { T, rho: rho2, du: dp / m, m };
}

// Rarefaction from K down to p: magnitude of the velocity change.
function _rarefDu(K, p) {
  const g = K.gas, n = RAREF_NODES, a = Math.log(p), b = Math.log(K.p), h = (b - a) / n;
  let sum = 0;
  for (let k = 0; k <= n; k++) {
    const pk = Math.exp(a + k * h);
    const T = gasIsentropeT(g, K.T, K.p, pk);
    const w = k === 0 || k === n ? 1 : (k % 2 ? 4 : 2);
    sum += w * g.R * T / gasSound(g, T);                // dp/(rho a) = (R T / a) d ln p
  }
  return sum * h / 3;
}

// Velocity function f_K(p) (Toro Eq. 4.6 generalised): positive for a shock.
function _fK(K, p) { return p > K.p ? _shock(K, p).du : -_rarefDu(K, p); }

// Tabulate a rarefaction fan from pK down to p*: u(p), T, rho, a along the isentrope.
function _fanTable(K, pStar, sign) {
  const g = K.gas, n = FAN_TABLE, a = Math.log(K.p), b = Math.log(pStar), h = (b - a) / n;
  const tab = [];
  let u = K.u, prev = null;
  for (let k = 0; k <= n; k++) {
    const p = Math.exp(a + k * h), T = gasIsentropeT(g, K.T, K.p, p), c = gasSound(g, T);
    const w = g.R * T / c;                              // integrand in ln p
    if (prev) u += sign * 0.5 * (prev.w + w) * h;       // du = -/+ dp/(rho a): left fan speeds up as p drops
    prev = { p, T, rho: p / (g.R * T), a: c, u, w };
    tab.push(prev);
  }
  return tab;
}

// Solve the Riemann problem; returns the star region and the wave pattern.
function riemannGasSolve(gasL, left, gasR, right) {
  const L = _side(gasL, left), R = _side(gasR, right);
  const F = p => _fK(L, p) + _fK(R, p) + (R.u - L.u);

  // Bracket p* (F is monotone increasing); F(0+) > 0 means a vacuum forms.
  let lo = 1e-9 * Math.min(L.p, R.p), hi = Math.max(L.p, R.p);
  if (F(lo) > 0) return null;
  while (F(hi) < 0) hi *= 2;
  for (let it = 0; it < PSTAR_ITERS && hi / lo - 1 > PSTAR_RTOL; it++) {
    const mid = Math.sqrt(lo * hi);
    if (F(mid) < 0) lo = mid; else hi = mid;
  }
  const pStar = Math.sqrt(lo * hi);
  const fL = _fK(L, pStar), fR = _fK(R, pStar);
  const uStar = 0.5 * (L.u + R.u) + 0.5 * (fR - fL);

  const wave = (K, sign) => {        // sign = -1 left wave, +1 right wave
    if (pStar > K.p) {
      const s = _shock(K, pStar);
      return { kind: "shock", T: s.T, rho: s.rho, speed: K.u + sign * s.m / K.rho, mach: s.m / (K.rho * K.a) };
    }
    const tab = _fanTable(K, pStar, sign);
    const last = tab[tab.length - 1];
    return { kind: "rarefaction", T: last.T, rho: last.rho, table: tab,
             head: K.u + sign * K.a, tail: uStar + sign * last.a };
  };
  return { pStar, uStar, L, R, left: wave(L, -1), right: wave(R, +1) };
}

// Sample the self-similar solution at xi = (x - x0)/t.
function _sampleGas(sol, xi) {
  const { L, R, pStar, uStar } = sol;
  const onLeft = xi < uStar;
  const K = onLeft ? L : R, W = onLeft ? sol.left : sol.right, sign = onLeft ? -1 : 1;
  const out = (gas, rho, u, p, T) => ({ gas, rho, u, p, T, left: onLeft });
  if (W.kind === "shock") {
    const outside = sign < 0 ? xi < W.speed : xi > W.speed;
    return outside ? out(K.gas, K.rho, K.u, K.p, K.T) : out(K.gas, W.rho, uStar, pStar, W.T);
  }
  const outsideHead = sign < 0 ? xi < W.head : xi > W.head;
  if (outsideHead) return out(K.gas, K.rho, K.u, K.p, K.T);
  const insideTail = sign < 0 ? xi > W.tail : xi < W.tail;
  if (insideTail) return out(K.gas, W.rho, uStar, pStar, W.T);
  // inside the fan: find where the characteristic speed u -/+ a equals xi
  const tab = W.table, lam = e => e.u + sign * e.a;
  let a = 0, b = tab.length - 1;
  const up = lam(tab[b]) > lam(tab[a]);
  while (b - a > 1) {
    const m = (a + b) >> 1;
    if ((lam(tab[m]) < xi) === up) a = m; else b = m;
  }
  const A = tab[a], B = tab[b], w = Math.max(0, Math.min(1, (xi - lam(A)) / (lam(B) - lam(A) || 1)));
  const T = A.T + w * (B.T - A.T), p = A.p + w * (B.p - A.p);
  return out(K.gas, p / (K.gas.R * T), A.u + w * (B.u - A.u), p, T);
}

// Exact solution sampled at xArr, time t, with the derived fields the UI plots.
// Returns null if a vacuum forms.
function exactRiemannGas(gasL, left, gasR, right, x0, t, xArr, solved) {
  const sol = solved || riemannGasSolve(gasL, left, gasR, right);
  if (!sol) return null;
  const n = xArr.length;
  const o = { x: xArr, rho: new Array(n), u: new Array(n), p: new Array(n), T: new Array(n),
              a: new Array(n), M: new Array(n), e: new Array(n), s: new Array(n), gamma: new Array(n), solution: sol };
  for (let i = 0; i < n; i++) {
    const xi = t > 0 ? (xArr[i] - x0) / t : (xArr[i] < x0 ? -Infinity : Infinity);
    const q = _sampleGas(sol, xi), g = q.gas;
    o.rho[i] = q.rho; o.u[i] = q.u; o.p[i] = q.p; o.T[i] = q.T;
    o.a[i] = gasSound(g, q.T); o.M[i] = q.u / o.a[i]; o.gamma[i] = gasGamma(g, q.T);
    o.e[i] = gasE(g, q.T) - gasE(g, EOS_TREF);
    o.s[i] = g.phi(q.T) - g.R * Math.log(q.p / EOS_PREF);
  }
  return o;
}

// ---------- Local fluid velocity at position x (for particle advection) ----------
function localVelocity(state, x) {
  const { N, dx, rho, mom, L } = state;
  if (x < 0 || x > L) return 0;
  const i = Math.min(N - 1, Math.max(0, Math.floor(x / dx)));
  return mom[i + 2] / rho[i + 2];
}

// ---------- Classic shock tube setups ----------
const PRESETS = {
  "Sod":          { L: { rho: 1.000, u:  0.0,   p: 1.000  }, R: { rho: 0.125, u:  0.0, p: 0.100 }, tEnd: 0.200, diaphragm: 0.5 },
  "Lax":          { L: { rho: 0.445, u:  0.698, p: 3.528  }, R: { rho: 0.500, u:  0.0, p: 0.571 }, tEnd: 0.160, diaphragm: 0.5 },
  "123 problem":  { L: { rho: 1.000, u: -2.0,   p: 0.400  }, R: { rho: 1.000, u:  2.0, p: 0.400 }, tEnd: 0.150, diaphragm: 0.5 },
  "Strong shock": { L: { rho: 1.000, u:  0.0,   p: 1000.0 }, R: { rho: 1.000, u:  0.0, p: 0.010 }, tEnd: 0.012, diaphragm: 0.5 },
  "Blast wave":   { L: { rho: 1.000, u:  0.0,   p: 0.010  }, R: { rho: 1.000, u:  0.0, p: 100.0 }, tEnd: 0.035, diaphragm: 0.5 }
};

// ---------- Node.js export (lets you unit-test from the command line) ----------
if (typeof module !== "undefined" && module.exports) {
  module.exports = {
    idealEos, gasEos, binaryEos,
    primToCons, fluxFromCons,
    hllc, hll, roe, rusanov, FLUXES,
    minmod, applyBcs,
    computeRhs, stableDt, initState, stepUntil, extract, localVelocity,
    exactRiemann, riemannGasSolve, exactRiemannGas,
    PRESETS
  };
}
