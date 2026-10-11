/* ============================================================
   solver2d.js
   2D compressible Euler solver — pure math layer (no DOM, no UI)
   ------------------------------------------------------------
   Methods : HLLC / HLL / Roe / Rusanov flux (with the tangential
             momentum carried as a passive scalar across the wave
             pattern) + MUSCL (minmod) + SSP-RK2 + adaptive CFL.
   Equations:
       dU/dt + dF(U)/dx + dG(U)/dy = 0
       U = [rho, rho*u, rho*v, E]
       F = [rho*u, rho*u^2 + p, rho*u*v,     u*(E+p)]
       G = [rho*v, rho*u*v,     rho*v^2 + p, v*(E+p)]
       E = p/(gamma - 1) + (1/2) rho (u^2 + v^2)

   Storage: flat Float64Array of length Nxt*Nyt where
            Nxt = Nx + 4, Nyt = Ny + 4 (2 ghost cells per side).
            idx(i,j) = j*Nxt + i.   Interior: i in [2, Nx+1], j in [2, Ny+1].

   Public API (attached to globals once this script is loaded):
     - hllc2d(L, R, gamma, axis)        axis: 0 = x, 1 = y
     - hll2d  (L, R, gamma, axis)
     - roe2d  (L, R, gamma, axis)
     - rusanov2d(L, R, gamma, axis)
     - FLUXES2D                          name -> flux function
     - applyBcs2D(state)                 transmissive on all four walls
     - initState2D(params)
     - stableDt2D(state)                 CFL-limited time step
     - stepUntil2D(state, tTarget)
     - extract2D(state)                  -> { Nx, Ny, xs, ys, rho, u, v, p, M, e }
     - slice2D(field, Nx, Ny, axis, idx) -> 1D array along x or y
     - PRESETS2D
   ============================================================ */

// ---------- Equation-of-state helpers ----------
function _pressure2D(rho, u, v, E, gamma) {
  return (gamma - 1) * (E - 0.5 * rho * (u * u + v * v));
}
function _sound2D(p, rho, gamma) {
  return Math.sqrt(gamma * Math.max(p, 1e-12) / Math.max(rho, 1e-12));
}

// ---------- Physical flux in x (axis=0) or y (axis=1) ----------
// state[0..3] = rho, rho*u, rho*v, E.
function _flux2D(s, gamma, axis) {
  const rho = s[0], mu = s[1], mv = s[2], E = s[3];
  const u = mu / rho, v = mv / rho;
  const p = _pressure2D(rho, u, v, E, gamma);
  if (axis === 0) {
    return [mu, mu * u + p, mu * v, u * (E + p)];
  }
  return [mv, mv * u, mv * v + p, v * (E + p)];
}

// ---------- Rusanov / local Lax-Friedrichs ----------
function rusanov2d(L, R, gamma, axis) {
  const rhoL = L[0], uL = L[1] / rhoL, vL = L[2] / rhoL;
  const rhoR = R[0], uR = R[1] / rhoR, vR = R[2] / rhoR;
  const pL = _pressure2D(rhoL, uL, vL, L[3], gamma);
  const pR = _pressure2D(rhoR, uR, vR, R[3], gamma);
  const aL = _sound2D(pL, rhoL, gamma);
  const aR = _sound2D(pR, rhoR, gamma);
  const nL = axis === 0 ? uL : vL;
  const nR = axis === 0 ? uR : vR;
  const s  = Math.max(Math.abs(nL) + aL, Math.abs(nR) + aR);
  const FL = _flux2D(L, gamma, axis);
  const FR = _flux2D(R, gamma, axis);
  return [
    0.5 * (FL[0] + FR[0]) - 0.5 * s * (R[0] - L[0]),
    0.5 * (FL[1] + FR[1]) - 0.5 * s * (R[1] - L[1]),
    0.5 * (FL[2] + FR[2]) - 0.5 * s * (R[2] - L[2]),
    0.5 * (FL[3] + FR[3]) - 0.5 * s * (R[3] - L[3])
  ];
}

// ---------- HLL ----------
function hll2d(L, R, gamma, axis) {
  const rhoL = L[0], uL = L[1] / rhoL, vL = L[2] / rhoL;
  const rhoR = R[0], uR = R[1] / rhoR, vR = R[2] / rhoR;
  const pL = _pressure2D(rhoL, uL, vL, L[3], gamma);
  const pR = _pressure2D(rhoR, uR, vR, R[3], gamma);
  const aL = _sound2D(pL, rhoL, gamma);
  const aR = _sound2D(pR, rhoR, gamma);
  const nL = axis === 0 ? uL : vL;
  const nR = axis === 0 ? uR : vR;
  const SL = Math.min(nL - aL, nR - aR);
  const SR = Math.max(nL + aL, nR + aR);
  const FL = _flux2D(L, gamma, axis);
  const FR = _flux2D(R, gamma, axis);
  if (SL >= 0) return FL;
  if (SR <= 0) return FR;
  const denom = SR - SL;
  return [
    (SR * FL[0] - SL * FR[0] + SL * SR * (R[0] - L[0])) / denom,
    (SR * FL[1] - SL * FR[1] + SL * SR * (R[1] - L[1])) / denom,
    (SR * FL[2] - SL * FR[2] + SL * SR * (R[2] - L[2])) / denom,
    (SR * FL[3] - SL * FR[3] + SL * SR * (R[3] - L[3])) / denom
  ];
}

// ---------- HLLC ----------
// 2D extension of Toro Ch. 10: tangential momentum is advected with the
// contact wave, so the only change vs. 1D is one extra row in the star
// state and the use of |velocity|^2 in the pressure.
function hllc2d(L, R, gamma, axis) {
  const rhoL = L[0], uL = L[1] / rhoL, vL = L[2] / rhoL;
  const rhoR = R[0], uR = R[1] / rhoR, vR = R[2] / rhoR;
  const pL = _pressure2D(rhoL, uL, vL, L[3], gamma);
  const pR = _pressure2D(rhoR, uR, vR, R[3], gamma);
  const aL = _sound2D(pL, rhoL, gamma);
  const aR = _sound2D(pR, rhoR, gamma);
  const nL = axis === 0 ? uL : vL;
  const nR = axis === 0 ? uR : vR;
  const tL = axis === 0 ? vL : uL;       // tangential
  const tR = axis === 0 ? vR : uR;
  const SL = Math.min(nL - aL, nR - aR);
  const SR = Math.max(nL + aL, nR + aR);
  const sStar = (pR - pL + rhoL * nL * (SL - nL) - rhoR * nR * (SR - nR)) /
                (rhoL * (SL - nL) - rhoR * (SR - nR));
  const FL = _flux2D(L, gamma, axis);
  const FR = _flux2D(R, gamma, axis);
  if (SL >= 0) return FL;
  if (SR <= 0) return FR;

  // Build the star state for whichever side the contact is on.
  const star = (rho, n, t, p, S, E) => {
    const f  = rho * (S - n) / (S - sStar);
    const Es = f * (E / rho + (sStar - n) * (sStar + p / (rho * (S - n))));
    if (axis === 0) return [f, f * sStar, f * t, Es];   // normal = u, tangential = v
    return [f, f * t, f * sStar, Es];                    // normal = v, tangential = u
  };

  if (sStar >= 0) {
    const US = star(rhoL, nL, tL, pL, SL, L[3]);
    return [
      FL[0] + SL * (US[0] - L[0]),
      FL[1] + SL * (US[1] - L[1]),
      FL[2] + SL * (US[2] - L[2]),
      FL[3] + SL * (US[3] - L[3])
    ];
  }
  const US = star(rhoR, nR, tR, pR, SR, R[3]);
  return [
    FR[0] + SR * (US[0] - R[0]),
    FR[1] + SR * (US[1] - R[1]),
    FR[2] + SR * (US[2] - R[2]),
    FR[3] + SR * (US[3] - R[3])
  ];
}

// ---------- Roe with Harten entropy fix ----------
// Uses the rotational invariance of the Euler equations: work in
// (normal, tangential) coordinates, then assemble the flux in (x,y).
function roe2d(L, R, gamma, axis) {
  const rhoL = L[0], uL = L[1] / rhoL, vL = L[2] / rhoL, EL = L[3];
  const rhoR = R[0], uR = R[1] / rhoR, vR = R[2] / rhoR, ER = R[3];
  const pL = _pressure2D(rhoL, uL, vL, EL, gamma);
  const pR = _pressure2D(rhoR, uR, vR, ER, gamma);
  const HL = (EL + pL) / rhoL;
  const HR = (ER + pR) / rhoR;

  // Roe averages
  const sL = Math.sqrt(Math.max(rhoL, 1e-12));
  const sR = Math.sqrt(Math.max(rhoR, 1e-12));
  const sum = sL + sR;
  const uH = (sL * uL + sR * uR) / sum;
  const vH = (sL * vL + sR * vR) / sum;
  const HH = (sL * HL + sR * HR) / sum;
  const qH2 = uH * uH + vH * vH;
  const a2  = (gamma - 1) * (HH - 0.5 * qH2);
  const aH  = Math.sqrt(Math.max(a2, 1e-12));
  const rhoH = sL * sR;

  // Pick normal / tangential
  const nH = axis === 0 ? uH : vH;
  const tH = axis === 0 ? vH : uH;
  const dnL = axis === 0 ? uL : vL;
  const dnR = axis === 0 ? uR : vR;
  const dtL = axis === 0 ? vL : uL;
  const dtR = axis === 0 ? vR : uR;

  const drho = rhoR - rhoL;
  const dn   = dnR - dnL;
  const dt   = dtR - dtL;
  const dp   = pR - pL;

  // Wave strengths
  const alpha1 = (dp - rhoH * aH * dn) / (2 * aH * aH);
  const alpha2 = drho - dp / (aH * aH);
  const alpha3 = (dp + rhoH * aH * dn) / (2 * aH * aH);
  const alpha4 = rhoH;                   // shear wave strength factor

  // Eigenvalues with Harten fix on the genuinely nonlinear fields
  let lam1 = nH - aH, lam2 = nH, lam3 = nH + aH, lam4 = nH;
  const eps = 0.1 * aH;
  if (Math.abs(lam1) < eps) lam1 = (lam1 * lam1 + eps * eps) / (2 * eps);
  if (Math.abs(lam3) < eps) lam3 = (lam3 * lam3 + eps * eps) / (2 * eps);
  const al1 = Math.abs(lam1), al2 = Math.abs(lam2), al3 = Math.abs(lam3), al4 = Math.abs(lam4);

  // Right eigenvectors in (rho, rho*n, rho*t, E) coordinates
  const r1 = [1, nH - aH, tH, HH - nH * aH];
  const r2 = [1, nH,      tH, 0.5 * qH2];
  const r3 = [1, nH + aH, tH, HH + nH * aH];
  const r4 = [0, 0,       1,  tH];               // shear wave

  const FL = _flux2D(L, gamma, axis);
  const FR = _flux2D(R, gamma, axis);

  // Combine, then unrotate the (normal, tangential) momentum entries back to (x, y).
  const diss = [
    al1 * alpha1 * r1[0] + al2 * alpha2 * r2[0] + al3 * alpha3 * r3[0] + al4 * alpha4 * dt * r4[0],
    al1 * alpha1 * r1[1] + al2 * alpha2 * r2[1] + al3 * alpha3 * r3[1] + al4 * alpha4 * dt * r4[1],
    al1 * alpha1 * r1[2] + al2 * alpha2 * r2[2] + al3 * alpha3 * r3[2] + al4 * alpha4 * dt * r4[2],
    al1 * alpha1 * r1[3] + al2 * alpha2 * r2[3] + al3 * alpha3 * r3[3] + al4 * alpha4 * dt * r4[3]
  ];

  // diss[1] is the dissipation on the *normal* momentum row, diss[2] on the
  // *tangential* row. Place them into the x- and y-momentum slots correctly.
  let dissMx, dissMy;
  if (axis === 0) { dissMx = diss[1]; dissMy = diss[2]; }
  else            { dissMx = diss[2]; dissMy = diss[1]; }

  return [
    0.5 * (FL[0] + FR[0]) - 0.5 * diss[0],
    0.5 * (FL[1] + FR[1]) - 0.5 * dissMx,
    0.5 * (FL[2] + FR[2]) - 0.5 * dissMy,
    0.5 * (FL[3] + FR[3]) - 0.5 * diss[3]
  ];
}

// ---------- Flux dispatch ----------
const FLUXES2D = { HLLC: hllc2d, HLL: hll2d, Roe: roe2d, Rusanov: rusanov2d };

// ---------- minmod slope limiter (same as 1D) ----------
function _minmod2D(a, b) {
  if (a * b <= 0) return 0;
  return Math.sign(a) * Math.min(Math.abs(a), Math.abs(b));
}

// ---------- Transmissive BC on all four walls (2 ghost cells per side) ----------
function applyBcs2D(state) {
  const { rho, mu, mv, E, Nx, Ny } = state;
  const Nxt = Nx + 4;
  // Left and right columns
  for (let j = 0; j < Ny + 4; j++) {
    const rowL = j * Nxt;
    const i2 = rowL + 2, iN = rowL + Nx + 1;
    rho[rowL    ] = rho[i2]; mu[rowL    ] = mu[i2]; mv[rowL    ] = mv[i2]; E[rowL    ] = E[i2];
    rho[rowL + 1] = rho[i2]; mu[rowL + 1] = mu[i2]; mv[rowL + 1] = mv[i2]; E[rowL + 1] = E[i2];
    rho[rowL + Nx + 2] = rho[iN]; mu[rowL + Nx + 2] = mu[iN]; mv[rowL + Nx + 2] = mv[iN]; E[rowL + Nx + 2] = E[iN];
    rho[rowL + Nx + 3] = rho[iN]; mu[rowL + Nx + 3] = mu[iN]; mv[rowL + Nx + 3] = mv[iN]; E[rowL + Nx + 3] = E[iN];
  }
  // Bottom and top rows
  for (let i = 0; i < Nxt; i++) {
    const idxB = 2 * Nxt + i;            // first interior row
    const idxT = (Ny + 1) * Nxt + i;     // last interior row
    rho[              i] = rho[idxB]; mu[              i] = mu[idxB]; mv[              i] = mv[idxB]; E[              i] = E[idxB];
    rho[    Nxt    + i] = rho[idxB]; mu[    Nxt    + i] = mu[idxB]; mv[    Nxt    + i] = mv[idxB]; E[    Nxt    + i] = E[idxB];
    rho[(Ny + 2) * Nxt + i] = rho[idxT]; mu[(Ny + 2) * Nxt + i] = mu[idxT]; mv[(Ny + 2) * Nxt + i] = mv[idxT]; E[(Ny + 2) * Nxt + i] = E[idxT];
    rho[(Ny + 3) * Nxt + i] = rho[idxT]; mu[(Ny + 3) * Nxt + i] = mu[idxT]; mv[(Ny + 3) * Nxt + i] = mv[idxT]; E[(Ny + 3) * Nxt + i] = E[idxT];
  }
}

// ---------- Compute RHS = -dF/dx - dG/dy on interior cells ----------
// Slopes and flux arrays are allocated in state.work and reused each call.
function computeRhs2D(state, fluxFn, rho, mu, mv, E) {
  const { Nx, Ny, dx, dy, gamma, work } = state;
  const Nxt = Nx + 4, Nyt = Ny + 4;
  const { sxR, sxMu, sxMv, sxE, syR, syMu, syMv, syE,
          Fx0, Fx1, Fx2, Fx3, Gy0, Gy1, Gy2, Gy3,
          r1, r2, r3, r4 } = work;

  // ---- x-slopes (over j in 0..Nyt-1, i in 1..Nxt-2) ----
  for (let j = 0; j < Nyt; j++) {
    const row = j * Nxt;
    for (let i = 1; i < Nxt - 1; i++) {
      const k = row + i;
      sxR [k] = _minmod2D(rho[k] - rho[k - 1], rho[k + 1] - rho[k]);
      sxMu[k] = _minmod2D(mu [k] - mu [k - 1], mu [k + 1] - mu [k]);
      sxMv[k] = _minmod2D(mv [k] - mv [k - 1], mv [k + 1] - mv [k]);
      sxE [k] = _minmod2D(E  [k] - E  [k - 1], E  [k + 1] - E  [k]);
    }
  }
  // ---- y-slopes (over j in 1..Nyt-2, i in 0..Nxt-1) ----
  for (let j = 1; j < Nyt - 1; j++) {
    const rowN = j * Nxt, rowP = (j - 1) * Nxt, rowF = (j + 1) * Nxt;
    for (let i = 0; i < Nxt; i++) {
      const k = rowN + i;
      syR [k] = _minmod2D(rho[k] - rho[rowP + i], rho[rowF + i] - rho[k]);
      syMu[k] = _minmod2D(mu [k] - mu [rowP + i], mu [rowF + i] - mu [k]);
      syMv[k] = _minmod2D(mv [k] - mv [rowP + i], mv [rowF + i] - mv [k]);
      syE [k] = _minmod2D(E  [k] - E  [rowP + i], E  [rowF + i] - E  [k]);
    }
  }

  // ---- x-faces: between interior cells i, i+1 for i in [1, Nx+2], j in [2, Ny+1] ----
  // We need fluxes at x-faces bordering each interior cell, so faces i=1..Nx+1 (Nx+1 of them).
  const UL = [0, 0, 0, 0], UR = [0, 0, 0, 0];
  for (let j = 2; j < Ny + 2; j++) {
    const row = j * Nxt;
    for (let i = 1; i <= Nx + 1; i++) {
      const kL = row + i, kR = row + i + 1;
      UL[0] = rho[kL] + 0.5 * sxR [kL]; UR[0] = rho[kR] - 0.5 * sxR [kR];
      UL[1] = mu [kL] + 0.5 * sxMu[kL]; UR[1] = mu [kR] - 0.5 * sxMu[kR];
      UL[2] = mv [kL] + 0.5 * sxMv[kL]; UR[2] = mv [kR] - 0.5 * sxMv[kR];
      UL[3] = E  [kL] + 0.5 * sxE [kL]; UR[3] = E  [kR] - 0.5 * sxE [kR];
      const F = fluxFn(UL, UR, gamma, 0);
      const fi = row + i;
      Fx0[fi] = F[0]; Fx1[fi] = F[1]; Fx2[fi] = F[2]; Fx3[fi] = F[3];
    }
  }

  // ---- y-faces: between cells j, j+1 for j in [1, Ny+1], i in [2, Nx+1] ----
  for (let j = 1; j <= Ny + 1; j++) {
    const rowB = j * Nxt, rowT = (j + 1) * Nxt;
    for (let i = 2; i < Nx + 2; i++) {
      const kB = rowB + i, kT = rowT + i;
      UL[0] = rho[kB] + 0.5 * syR [kB]; UR[0] = rho[kT] - 0.5 * syR [kT];
      UL[1] = mu [kB] + 0.5 * syMu[kB]; UR[1] = mu [kT] - 0.5 * syMu[kT];
      UL[2] = mv [kB] + 0.5 * syMv[kB]; UR[2] = mv [kT] - 0.5 * syMv[kT];
      UL[3] = E  [kB] + 0.5 * syE [kB]; UR[3] = E  [kT] - 0.5 * syE [kT];
      const G = fluxFn(UL, UR, gamma, 1);
      const fi = rowB + i;
      Gy0[fi] = G[0]; Gy1[fi] = G[1]; Gy2[fi] = G[2]; Gy3[fi] = G[3];
    }
  }

  // ---- Divergence on interior cells (i in [2, Nx+1], j in [2, Ny+1]) ----
  // Flux indexing convention (matches loops above):
  //   Fx[row + i]      = x-flux at the face between cells (i,j) and (i+1,j)
  //                      → east face of cell i, west face of cell (i+1)
  //   Gy[j*Nxt + i]    = y-flux at the face between cells (i,j) and (i,j+1)
  //                      → north face of cell j, south face of cell (j+1)
  const invDx = 1 / dx, invDy = 1 / dy;
  for (let j = 2; j < Ny + 2; j++) {
    const row  = j * Nxt;
    const rowS = (j - 1) * Nxt;    // south face of cell (i, j) — flux between (i, j-1) and (i, j)
    const rowN = j * Nxt;          // north face of cell (i, j) — flux between (i, j) and (i, j+1)
    for (let i = 2; i < Nx + 2; i++) {
      const k  = row + i;
      const fW = row + i - 1;      // west face of cell i — flux between (i-1, j) and (i, j)
      const fE = row + i;          // east face of cell i — flux between (i, j) and (i+1, j)
      const gS = rowS + i, gN = rowN + i;
      r1[k] = -(Fx0[fE] - Fx0[fW]) * invDx - (Gy0[gN] - Gy0[gS]) * invDy;
      r2[k] = -(Fx1[fE] - Fx1[fW]) * invDx - (Gy1[gN] - Gy1[gS]) * invDy;
      r3[k] = -(Fx2[fE] - Fx2[fW]) * invDx - (Gy2[gN] - Gy2[gS]) * invDy;
      r4[k] = -(Fx3[fE] - Fx3[fW]) * invDx - (Gy3[gN] - Gy3[gS]) * invDy;
    }
  }
}

// ---------- Initial state ----------
// params: { Nx, Ny, Lx, Ly, gamma, cfl, flux, initFn(x,y) -> {rho,u,v,p} }
function initState2D(params) {
  const { Nx, Ny, Lx, Ly, gamma, initFn } = params;
  const dx = Lx / Nx, dy = Ly / Ny;
  const Nxt = Nx + 4, Nyt = Ny + 4, n = Nxt * Nyt;

  const rho = new Float64Array(n);
  const mu  = new Float64Array(n);
  const mv  = new Float64Array(n);
  const E   = new Float64Array(n);

  for (let j = 0; j < Ny; j++) {
    const yc = (j + 0.5) * dy;
    for (let i = 0; i < Nx; i++) {
      const xc = (i + 0.5) * dx;
      const s = initFn(xc, yc);
      const k = (j + 2) * Nxt + (i + 2);
      rho[k] = s.rho;
      mu [k] = s.rho * s.u;
      mv [k] = s.rho * s.v;
      E  [k] = s.p / (gamma - 1) + 0.5 * s.rho * (s.u * s.u + s.v * s.v);
    }
  }

  // Work arrays — allocate once, reuse forever.
  const work = {
    sxR:  new Float64Array(n), sxMu: new Float64Array(n), sxMv: new Float64Array(n), sxE:  new Float64Array(n),
    syR:  new Float64Array(n), syMu: new Float64Array(n), syMv: new Float64Array(n), syE:  new Float64Array(n),
    Fx0:  new Float64Array(n), Fx1:  new Float64Array(n), Fx2:  new Float64Array(n), Fx3:  new Float64Array(n),
    Gy0:  new Float64Array(n), Gy1:  new Float64Array(n), Gy2:  new Float64Array(n), Gy3:  new Float64Array(n),
    r1:   new Float64Array(n), r2:   new Float64Array(n), r3:   new Float64Array(n), r4:   new Float64Array(n),
    rho1: new Float64Array(n), mu1:  new Float64Array(n), mv1:  new Float64Array(n), E1:   new Float64Array(n)
  };

  const state = {
    rho, mu, mv, E,
    Nx, Ny, dx, dy, Lx, Ly, gamma,
    cfl:  params.cfl  || 0.4,
    flux: params.flux || "HLLC",
    t: 0, step: 0,
    work
  };
  applyBcs2D(state);
  return state;
}

// ---------- Stable time step from the CFL condition ----------
// dt = CFL * min(dx, dy) / max(|u| + |v| + a) over the interior, recomputed
// every step (CLAUDE.md §4, §8.6).
function stableDt2D(state) {
  const { rho, mu, mv, E, Nx, Ny, dx, dy, gamma, cfl } = state;
  const Nxt = Nx + 4;
  let maxSpeed = 0;
  for (let j = 2; j < Ny + 2; j++) {
    const row = j * Nxt;
    for (let i = 2; i < Nx + 2; i++) {
      const k = row + i;
      const r = rho[k];
      const u = mu[k] / r, v = mv[k] / r;
      const p = _pressure2D(r, u, v, E[k], gamma);
      const a = _sound2D(p, r, gamma);
      const s = Math.abs(u) + Math.abs(v) + a;
      if (s > maxSpeed) maxSpeed = s;
    }
  }
  return cfl * Math.min(dx, dy) / maxSpeed;
}

// ---------- Time-stepping (SSP-RK2 + adaptive CFL) ----------
function stepUntil2D(state, tTarget) {
  const { Nx, Ny } = state;
  const Nxt = Nx + 4;
  const fluxFn = FLUXES2D[state.flux] || hllc2d;
  const { rho, mu, mv, E, work } = state;
  const { rho1, mu1, mv1, E1, r1, r2, r3, r4 } = work;
  let t = state.t, step = state.step;
  let budget = 50000;

  while (t < tTarget && budget-- > 0) {
    let dt = stableDt2D(state);
    if (!(dt < Infinity)) break;
    if (t + dt > tTarget) dt = tTarget - t;

    // Stage 1: U1 = U + dt * L(U)
    computeRhs2D(state, fluxFn, rho, mu, mv, E);
    rho1.set(rho); mu1.set(mu); mv1.set(mv); E1.set(E);
    for (let j = 2; j < Ny + 2; j++) {
      const row = j * Nxt;
      for (let i = 2; i < Nx + 2; i++) {
        const k = row + i;
        rho1[k] = rho[k] + dt * r1[k];
        mu1 [k] = mu [k] + dt * r2[k];
        mv1 [k] = mv [k] + dt * r3[k];
        E1  [k] = E  [k] + dt * r4[k];
      }
    }
    applyBcs2D({ rho: rho1, mu: mu1, mv: mv1, E: E1, Nx, Ny });

    // Stage 2: U^{n+1} = 0.5*U + 0.5*(U1 + dt*L(U1))
    computeRhs2D(state, fluxFn, rho1, mu1, mv1, E1);
    for (let j = 2; j < Ny + 2; j++) {
      const row = j * Nxt;
      for (let i = 2; i < Nx + 2; i++) {
        const k = row + i;
        rho[k] = 0.5 * rho[k] + 0.5 * (rho1[k] + dt * r1[k]);
        mu [k] = 0.5 * mu [k] + 0.5 * (mu1 [k] + dt * r2[k]);
        mv [k] = 0.5 * mv [k] + 0.5 * (mv1 [k] + dt * r3[k]);
        E  [k] = 0.5 * E  [k] + 0.5 * (E1  [k] + dt * r4[k]);
      }
    }
    applyBcs2D(state);

    t += dt; step++;
  }

  state.t = t; state.step = step;
  return state;
}

// ---------- Extract primitive + derived fields on interior cells ----------
// Returns a snapshot of flat arrays of length Nx*Ny (no ghost cells).
function extract2D(state) {
  const { rho, mu, mv, E, Nx, Ny, dx, dy, gamma } = state;
  const Nxt = Nx + 4;
  const n = Nx * Ny;
  const R = new Float64Array(n), U = new Float64Array(n), V = new Float64Array(n);
  const P = new Float64Array(n), eF = new Float64Array(n), M = new Float64Array(n);
  const xs = new Float64Array(Nx), ys = new Float64Array(Ny);
  for (let i = 0; i < Nx; i++) xs[i] = (i + 0.5) * dx;
  for (let j = 0; j < Ny; j++) ys[j] = (j + 0.5) * dy;
  for (let j = 0; j < Ny; j++) {
    const rowI = j * Nx;
    const rowG = (j + 2) * Nxt;
    for (let i = 0; i < Nx; i++) {
      const k = rowG + i + 2;
      const ki = rowI + i;
      const r = rho[k];
      const u = mu[k] / r, v = mv[k] / r;
      const p = _pressure2D(r, u, v, E[k], gamma);
      const a = _sound2D(p, r, gamma);
      R[ki] = r; U[ki] = u; V[ki] = v;
      P[ki] = p;
      eF[ki] = p / ((gamma - 1) * r);
      M[ki] = Math.sqrt(u * u + v * v) / Math.max(a, 1e-12);
    }
  }
  return { Nx, Ny, xs, ys, rho: R, u: U, v: V, p: P, e: eF, M };
}

// ---------- Slice a 2D field (flat Nx*Ny array) along x or y ----------
// axis = "x" returns the row at j = jIdx (length Nx)
// axis = "y" returns the column at i = iIdx (length Ny)
function slice2D(field, Nx, Ny, axis, idx) {
  if (axis === "x") {
    const j = Math.min(Ny - 1, Math.max(0, idx));
    const out = new Float64Array(Nx);
    for (let i = 0; i < Nx; i++) out[i] = field[j * Nx + i];
    return out;
  } else {
    const i = Math.min(Nx - 1, Math.max(0, idx));
    const out = new Float64Array(Ny);
    for (let j = 0; j < Ny; j++) out[j] = field[j * Nx + i];
    return out;
  }
}

// ---------- 2D presets ----------
const PRESETS2D = {
  "Sod-radial": {
    Lx: 1, Ly: 1, N: 200, gamma: 1.4, cfl: 0.4, tEnd: 0.25,
    init: function (x, y) {
      const dx = x - 0.5, dy = y - 0.5;
      if (dx * dx + dy * dy < 0.16) return { rho: 1.0,   u: 0, v: 0, p: 1.0 };
      return                                 { rho: 0.125, u: 0, v: 0, p: 0.1 };
    }
  },
  "2D Riemann (config 3)": {
    // Lax-Liu / Schulz-Rinne 1993, configuration 3: four shocks meet.
    Lx: 1, Ly: 1, N: 200, gamma: 1.4, cfl: 0.35, tEnd: 0.30,
    init: function (x, y) {
      if (x >= 0.5 && y >= 0.5) return { rho: 1.5,    u: 0,     v: 0,     p: 1.5   };
      if (x <  0.5 && y >= 0.5) return { rho: 0.5323, u: 1.206, v: 0,     p: 0.3   };
      if (x <  0.5 && y <  0.5) return { rho: 0.138,  u: 1.206, v: 1.206, p: 0.029 };
      return                            { rho: 0.5323, u: 0,     v: 1.206, p: 0.3   };
    }
  },
  "Sedov blast": {
    Lx: 2, Ly: 2, N: 200, gamma: 1.4, cfl: 0.35, tEnd: 1.0,
    init: function (x, y) {
      // Energy spike in the center: a small disk gets very high pressure.
      const dx = x - 1.0, dy = y - 1.0;
      const r2 = dx * dx + dy * dy;
      if (r2 < 0.0025) return { rho: 1.0, u: 0, v: 0, p: 4e2 };
      return                   { rho: 1.0, u: 0, v: 0, p: 1e-4 };
    }
  },
  "Kelvin-Helmholtz": {
    Lx: 1, Ly: 1, N: 256, gamma: 1.4, cfl: 0.4, tEnd: 2.0,
    init: function (x, y) {
      // Two-layer shear with a small v perturbation. Transmissive BCs degrade
      // the canonical periodic setup but the rollup still appears clearly.
      const dy = (y > 0.25 && y < 0.75) ? 1 : 0;
      const u  = dy ?  0.5 : -0.5;
      const r  = dy ?  2.0 :  1.0;
      const v  = 0.01 * Math.sin(4 * Math.PI * x) *
                 (Math.exp(-((y - 0.25) * (y - 0.25)) / 0.002) +
                  Math.exp(-((y - 0.75) * (y - 0.75)) / 0.002));
      return { rho: r, u, v, p: 2.5 };
    }
  }
};

// ---------- Node CommonJS export ----------
if (typeof module !== "undefined" && module.exports) {
  module.exports = {
    hllc2d, hll2d, roe2d, rusanov2d, FLUXES2D,
    applyBcs2D, computeRhs2D,
    initState2D, stableDt2D, stepUntil2D, extract2D, slice2D,
    PRESETS2D
  };
}
