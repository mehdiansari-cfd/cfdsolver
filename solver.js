/* ============================================================
   solver.js
   1D Euler shock tube solver — pure math layer (no DOM, no UI)
   ------------------------------------------------------------
   Methods : HLLC / HLL / Roe / Rusanov flux + MUSCL (minmod)
             + SSP-RK2 + adaptive CFL
   Equations: dU/dt + dF(U)/dx = 0,  U = [rho, rho*u, E]
             E = p/(gamma - 1) + (1/2) rho u^2

   Public API (attached to globals once this script is loaded):
     - primToCons(rho, u, p, gamma)
     - fluxFromCons(rho, mom, ene, gamma)
     - hllc(UL, UR, gamma)
     - hll(UL, UR, gamma)
     - roe(UL, UR, gamma)
     - rusanov(UL, UR, gamma)
     - FLUXES                                     name -> flux function
     - minmod(a, b)
     - applyBcs(rho, mom, ene, bcType)            'transmissive' | 'reflective'
     - computeRhs(rho, mom, ene, N, dx, gamma, fluxFn)
     - initState(params)                          params.flux ('HLLC' | 'HLL' | 'Roe' | 'Rusanov')
     - stepUntil(state, tTarget)
     - extract(state)                              -> { x, rho, u, p, e, M, T, s, a }
     - exactRiemann(left, right, gamma, x0, t, xArr) -> { x, rho, u, p } | null
     - localVelocity(state, x)
     - PRESETS                                     classic shock tube setups
   ============================================================ */

// ---------- Variable conversions ----------
function primToCons(rho, u, p, gamma) {
  return [rho, rho * u, p / (gamma - 1) + 0.5 * rho * u * u];
}

function fluxFromCons(rho, mom, ene, gamma) {
  const u = mom / rho;
  const p = (gamma - 1) * (ene - 0.5 * rho * u * u);
  return [mom, mom * u + p, u * (ene + p)];
}

// ---------- HLLC approximate Riemann solver ----------
function hllc(L, R, gamma) {
  const [rhoL, momL, eneL] = L;
  const [rhoR, momR, eneR] = R;
  const uL = momL / rhoL, uR = momR / rhoR;
  const pL = (gamma - 1) * (eneL - 0.5 * rhoL * uL * uL);
  const pR = (gamma - 1) * (eneR - 0.5 * rhoR * uR * uR);
  const aL = Math.sqrt(gamma * Math.max(pL, 1e-12) / Math.max(rhoL, 1e-12));
  const aR = Math.sqrt(gamma * Math.max(pR, 1e-12) / Math.max(rhoR, 1e-12));

  // Davis wave-speed estimates
  const SL = Math.min(uL - aL, uR - aR);
  const SR = Math.max(uL + aL, uR + aR);

  // Contact wave speed (Toro Eq. 10.37)
  const sStar = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) /
                (rhoL * (SL - uL) - rhoR * (SR - uR));

  const FL = fluxFromCons(rhoL, momL, eneL, gamma);
  const FR = fluxFromCons(rhoR, momR, eneR, gamma);

  if (SL >= 0) return FL;
  if (SR <= 0) return FR;

  if (sStar >= 0) {
    const f = rhoL * (SL - uL) / (SL - sStar);
    const rS = f, mS = f * sStar;
    const eS = f * (eneL / rhoL + (sStar - uL) * (sStar + pL / (rhoL * (SL - uL))));
    return [FL[0] + SL * (rS - rhoL), FL[1] + SL * (mS - momL), FL[2] + SL * (eS - eneL)];
  } else {
    const f = rhoR * (SR - uR) / (SR - sStar);
    const rS = f, mS = f * sStar;
    const eS = f * (eneR / rhoR + (sStar - uR) * (sStar + pR / (rhoR * (SR - uR))));
    return [FR[0] + SR * (rS - rhoR), FR[1] + SR * (mS - momR), FR[2] + SR * (eS - eneR)];
  }
}

// ---------- HLL approximate Riemann solver ----------
// No contact wave; just averaged star state between SL and SR.
function hll(L, R, gamma) {
  const [rhoL, momL, eneL] = L;
  const [rhoR, momR, eneR] = R;
  const uL = momL / rhoL, uR = momR / rhoR;
  const pL = (gamma - 1) * (eneL - 0.5 * rhoL * uL * uL);
  const pR = (gamma - 1) * (eneR - 0.5 * rhoR * uR * uR);
  const aL = Math.sqrt(gamma * Math.max(pL, 1e-12) / Math.max(rhoL, 1e-12));
  const aR = Math.sqrt(gamma * Math.max(pR, 1e-12) / Math.max(rhoR, 1e-12));
  const SL = Math.min(uL - aL, uR - aR);
  const SR = Math.max(uL + aL, uR + aR);
  const FL = fluxFromCons(rhoL, momL, eneL, gamma);
  const FR = fluxFromCons(rhoR, momR, eneR, gamma);
  if (SL >= 0) return FL;
  if (SR <= 0) return FR;
  const denom = SR - SL;
  return [
    (SR * FL[0] - SL * FR[0] + SL * SR * (rhoR - rhoL)) / denom,
    (SR * FL[1] - SL * FR[1] + SL * SR * (momR - momL)) / denom,
    (SR * FL[2] - SL * FR[2] + SL * SR * (eneR - eneL)) / denom
  ];
}

// ---------- Rusanov (local Lax-Friedrichs) flux ----------
function rusanov(L, R, gamma) {
  const [rhoL, momL, eneL] = L;
  const [rhoR, momR, eneR] = R;
  const uL = momL / rhoL, uR = momR / rhoR;
  const pL = (gamma - 1) * (eneL - 0.5 * rhoL * uL * uL);
  const pR = (gamma - 1) * (eneR - 0.5 * rhoR * uR * uR);
  const aL = Math.sqrt(gamma * Math.max(pL, 1e-12) / Math.max(rhoL, 1e-12));
  const aR = Math.sqrt(gamma * Math.max(pR, 1e-12) / Math.max(rhoR, 1e-12));
  const s  = Math.max(Math.abs(uL) + aL, Math.abs(uR) + aR);
  const FL = fluxFromCons(rhoL, momL, eneL, gamma);
  const FR = fluxFromCons(rhoR, momR, eneR, gamma);
  return [
    0.5 * (FL[0] + FR[0]) - 0.5 * s * (rhoR - rhoL),
    0.5 * (FL[1] + FR[1]) - 0.5 * s * (momR - momL),
    0.5 * (FL[2] + FR[2]) - 0.5 * s * (eneR - eneL)
  ];
}

// ---------- Roe-Pike flux with Harten entropy fix ----------
function roe(L, R, gamma) {
  const [rhoL, momL, eneL] = L;
  const [rhoR, momR, eneR] = R;
  const uL = momL / rhoL, uR = momR / rhoR;
  const pL = (gamma - 1) * (eneL - 0.5 * rhoL * uL * uL);
  const pR = (gamma - 1) * (eneR - 0.5 * rhoR * uR * uR);
  const HL = (eneL + pL) / rhoL;
  const HR = (eneR + pR) / rhoR;

  // Roe averages
  const sL = Math.sqrt(Math.max(rhoL, 1e-12));
  const sR = Math.sqrt(Math.max(rhoR, 1e-12));
  const sum = sL + sR;
  const uH = (sL * uL + sR * uR) / sum;
  const HH = (sL * HL + sR * HR) / sum;
  const rhoH = sL * sR;
  const a2 = (gamma - 1) * (HH - 0.5 * uH * uH);
  const aH = Math.sqrt(Math.max(a2, 1e-12));

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

  // Right eigenvectors
  const r1 = [1, uH - aH, HH - uH * aH];
  const r2 = [1, uH,      0.5 * uH * uH];
  const r3 = [1, uH + aH, HH + uH * aH];

  const FL = fluxFromCons(rhoL, momL, eneL, gamma);
  const FR = fluxFromCons(rhoR, momR, eneR, gamma);

  const out = [0, 0, 0];
  for (let k = 0; k < 3; k++) {
    out[k] = 0.5 * (FL[k] + FR[k])
           - 0.5 * (al1 * alpha1 * r1[k] + al2 * alpha2 * r2[k] + al3 * alpha3 * r3[k]);
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
function applyBcs(rho, mom, ene, bc) {
  const n = rho.length;
  if (bc === "reflective") {
    rho[1] = rho[2]; mom[1] = -mom[2]; ene[1] = ene[2];
    rho[0] = rho[3]; mom[0] = -mom[3]; ene[0] = ene[3];
    rho[n - 2] = rho[n - 3]; mom[n - 2] = -mom[n - 3]; ene[n - 2] = ene[n - 3];
    rho[n - 1] = rho[n - 4]; mom[n - 1] = -mom[n - 4]; ene[n - 1] = ene[n - 4];
  } else {
    // transmissive (zero-gradient)
    rho[1] = rho[2]; mom[1] = mom[2]; ene[1] = ene[2];
    rho[0] = rho[2]; mom[0] = mom[2]; ene[0] = ene[2];
    rho[n - 2] = rho[n - 3]; mom[n - 2] = mom[n - 3]; ene[n - 2] = ene[n - 3];
    rho[n - 1] = rho[n - 3]; mom[n - 1] = mom[n - 3]; ene[n - 1] = ene[n - 3];
  }
}

// ---------- Spatial operator: dU/dt = -dF/dx (interior cells only) ----------
// fluxFn defaults to HLLC for backward compatibility.
function computeRhs(rho, mom, ene, N, dx, gamma, fluxFn) {
  const flux = fluxFn || hllc;
  // Slopes in cells 1..N+2 (N+2 of them)
  const sR = new Array(N + 2), sM = new Array(N + 2), sE = new Array(N + 2);
  for (let i = 0; i < N + 2; i++) {
    const k = i + 1;
    sR[i] = minmod(rho[k] - rho[k - 1], rho[k + 1] - rho[k]);
    sM[i] = minmod(mom[k] - mom[k - 1], mom[k + 1] - mom[k]);
    sE[i] = minmod(ene[k] - ene[k - 1], ene[k + 1] - ene[k]);
  }
  // Fluxes at N+1 interfaces (MUSCL-reconstructed left/right states)
  const F1 = new Array(N + 1), F2 = new Array(N + 1), F3 = new Array(N + 1);
  for (let i = 0; i <= N; i++) {
    const kL = i + 1, kR = i + 2;
    const UL = [rho[kL] + 0.5 * sR[i],     mom[kL] + 0.5 * sM[i],     ene[kL] + 0.5 * sE[i]];
    const UR = [rho[kR] - 0.5 * sR[i + 1], mom[kR] - 0.5 * sM[i + 1], ene[kR] - 0.5 * sE[i + 1]];
    const F = flux(UL, UR, gamma);
    F1[i] = F[0]; F2[i] = F[1]; F3[i] = F[2];
  }
  // Divergence in interior cells
  const r1 = new Array(N), r2 = new Array(N), r3 = new Array(N);
  for (let i = 0; i < N; i++) {
    r1[i] = -(F1[i + 1] - F1[i]) / dx;
    r2[i] = -(F2[i + 1] - F2[i]) / dx;
    r3[i] = -(F3[i + 1] - F3[i]) / dx;
  }
  return [r1, r2, r3];
}

// ---------- Initial state from a parameter object ----------
function initState(p) {
  const { N, L, gamma, diaphragm, leftState, rightState, bcType } = p;
  const dx = L / N;
  const tot = N + 4; // 2 ghost cells per side
  const rho = new Array(tot).fill(0), mom = new Array(tot).fill(0), ene = new Array(tot).fill(0);
  for (let i = 0; i < N; i++) {
    const x = (i + 0.5) * dx;
    const s = x < diaphragm ? leftState : rightState;
    const [r, m, e] = primToCons(s.rho, s.u, s.p, gamma);
    rho[i + 2] = r; mom[i + 2] = m; ene[i + 2] = e;
  }
  applyBcs(rho, mom, ene, bcType);
  const fluxName = p.flux || "HLLC";
  return {
    rho, mom, ene, t: 0, step: 0, dx, N, gamma, bcType,
    cfl: p.cfl, L, flux: fluxName,
    // Snapshot the initial states + diaphragm so we can compute the
    // exact Riemann solution later without depending on UI params.
    leftState: { rho: leftState.rho, u: leftState.u, p: leftState.p },
    rightState: { rho: rightState.rho, u: rightState.u, p: rightState.p },
    diaphragm
  };
}

// ---------- Step forward in time to tTarget (SSP-RK2 + adaptive CFL) ----------
function stepUntil(state, tTarget) {
  const { N, dx, gamma, bcType, cfl } = state;
  const fluxFn = FLUXES[state.flux] || hllc;
  let { rho, mom, ene, t, step } = state;
  let budget = 5000; // safety cap per call
  while (t < tTarget && budget-- > 0) {
    let maxSpeed = 0;
    for (let i = 2; i < N + 2; i++) {
      const u = mom[i] / rho[i];
      const p = (gamma - 1) * (ene[i] - 0.5 * rho[i] * u * u);
      const a = Math.sqrt(gamma * Math.max(p, 1e-12) / Math.max(rho[i], 1e-12));
      const s = Math.abs(u) + a;
      if (s > maxSpeed) maxSpeed = s;
    }
    let dt = cfl * dx / maxSpeed;
    if (t + dt > tTarget) dt = tTarget - t;

    // Stage 1
    const L1 = computeRhs(rho, mom, ene, N, dx, gamma, fluxFn);
    const rho1 = rho.slice(), mom1 = mom.slice(), ene1 = ene.slice();
    for (let i = 0; i < N; i++) {
      rho1[i + 2] = rho[i + 2] + dt * L1[0][i];
      mom1[i + 2] = mom[i + 2] + dt * L1[1][i];
      ene1[i + 2] = ene[i + 2] + dt * L1[2][i];
    }
    applyBcs(rho1, mom1, ene1, bcType);

    // Stage 2 (SSP)
    const L2 = computeRhs(rho1, mom1, ene1, N, dx, gamma, fluxFn);
    for (let i = 0; i < N; i++) {
      rho[i + 2] = 0.5 * rho[i + 2] + 0.5 * (rho1[i + 2] + dt * L2[0][i]);
      mom[i + 2] = 0.5 * mom[i + 2] + 0.5 * (mom1[i + 2] + dt * L2[1][i]);
      ene[i + 2] = 0.5 * ene[i + 2] + 0.5 * (ene1[i + 2] + dt * L2[2][i]);
    }
    applyBcs(rho, mom, ene, bcType);

    t += dt; step++;
  }
  return Object.assign({}, state, { rho, mom, ene, t, step });
}

// ---------- Extract primitive + derived variables on interior cells ----------
// Returns rho, u, p, e (internal energy), M (Mach), T (= p/rho), s (= ln(p/rho^gamma)), a (sound speed).
function extract(state) {
  const { rho, mom, ene, N, dx, gamma } = state;
  const x = new Array(N), R = new Array(N), U = new Array(N), P = new Array(N), E = new Array(N);
  const M = new Array(N), T = new Array(N), S = new Array(N), A = new Array(N);
  for (let i = 0; i < N; i++) {
    x[i] = (i + 0.5) * dx;
    R[i] = rho[i + 2];
    U[i] = mom[i + 2] / R[i];
    P[i] = (gamma - 1) * (ene[i + 2] - 0.5 * R[i] * U[i] * U[i]);
    E[i] = P[i] / ((gamma - 1) * R[i]);
    const aa = Math.sqrt(gamma * Math.max(P[i], 1e-12) / Math.max(R[i], 1e-12));
    A[i] = aa;
    M[i] = U[i] / aa;
    T[i] = P[i] / Math.max(R[i], 1e-12); // RT (R=1 implicit, non-dim)
    S[i] = Math.log(Math.max(P[i], 1e-12) / Math.pow(Math.max(R[i], 1e-12), gamma));
  }
  return { x, rho: R, u: U, p: P, e: E, M, T, s: S, a: A };
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
    primToCons, fluxFromCons,
    hllc, hll, roe, rusanov, FLUXES,
    minmod, applyBcs,
    computeRhs, initState, stepUntil, extract, localVelocity,
    exactRiemann,
    PRESETS
  };
}
