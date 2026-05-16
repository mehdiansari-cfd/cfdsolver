/* ============================================================
   solver.js
   1D Euler shock tube solver — pure math layer (no DOM, no UI)
   ------------------------------------------------------------
   Methods : HLLC flux + MUSCL (minmod) + SSP-RK2 + adaptive CFL
   Equations: dU/dt + dF(U)/dx = 0,  U = [rho, rho*u, E]
             E = p/(gamma - 1) + (1/2) rho u^2

   Public API (attached to globals once this script is loaded):
     - primToCons(rho, u, p, gamma)
     - fluxFromCons(rho, mom, ene, gamma)
     - hllc(UL, UR, gamma)
     - minmod(a, b)
     - applyBcs(rho, mom, ene, bcType)            'transmissive' | 'reflective'
     - computeRhs(rho, mom, ene, N, dx, gamma)
     - initState(params)
     - stepUntil(state, tTarget)
     - extract(state)                              -> { x, rho, u, p, e }
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
function computeRhs(rho, mom, ene, N, dx, gamma) {
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
    const F = hllc(UL, UR, gamma);
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
  return { rho, mom, ene, t: 0, step: 0, dx, N, gamma, bcType, cfl: p.cfl, L };
}

// ---------- Step forward in time to tTarget (SSP-RK2 + adaptive CFL) ----------
function stepUntil(state, tTarget) {
  const { N, dx, gamma, bcType, cfl } = state;
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
    const L1 = computeRhs(rho, mom, ene, N, dx, gamma);
    const rho1 = rho.slice(), mom1 = mom.slice(), ene1 = ene.slice();
    for (let i = 0; i < N; i++) {
      rho1[i + 2] = rho[i + 2] + dt * L1[0][i];
      mom1[i + 2] = mom[i + 2] + dt * L1[1][i];
      ene1[i + 2] = ene[i + 2] + dt * L1[2][i];
    }
    applyBcs(rho1, mom1, ene1, bcType);

    // Stage 2 (SSP)
    const L2 = computeRhs(rho1, mom1, ene1, N, dx, gamma);
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

// ---------- Extract primitive variables on interior cells ----------
function extract(state) {
  const { rho, mom, ene, N, dx, gamma } = state;
  const x = new Array(N), R = new Array(N), U = new Array(N), P = new Array(N), E = new Array(N);
  for (let i = 0; i < N; i++) {
    x[i] = (i + 0.5) * dx;
    R[i] = rho[i + 2];
    U[i] = mom[i + 2] / R[i];
    P[i] = (gamma - 1) * (ene[i + 2] - 0.5 * R[i] * U[i] * U[i]);
    E[i] = P[i] / ((gamma - 1) * R[i]);
  }
  return { x, rho: R, u: U, p: P, e: E };
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
    primToCons, fluxFromCons, hllc, minmod, applyBcs,
    computeRhs, initState, stepUntil, extract, localVelocity,
    PRESETS
  };
}
