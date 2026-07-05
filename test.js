/* ============================================================
   test.js  —  zero-dependency verification suite
   ------------------------------------------------------------
   Run with:  node test.js
   Uses only Node's built-in `assert`. No test framework, no
   dev-dependencies — consistent with the buildless project rule.

   Checks:
     1. Accuracy    — numerical Sod matches the exact Riemann
                      solution (L1 density error below tolerance).
     2. Positivity  — density and pressure stay > 0 on the
                      strong-shock (1D) and Sedov blast (2D)
                      problems, which stress positivity hardest.
     3. Conservation — total mass and energy are conserved to
                      near machine precision while the waves stay
                      inside a transmissive domain.
   ============================================================ */

const assert = require("assert");
const S  = require("./solver");
require("./solver");            // nozzle reuses solver's globals if needed
const S2 = require("./solver2d");

let passed = 0;
function check(name, fn) {
  try {
    fn();
    console.log("  ok   " + name);
    passed++;
  } catch (e) {
    console.error("  FAIL " + name);
    console.error("       " + e.message);
    process.exitCode = 1;
  }
}

// ---- helpers ---------------------------------------------------------------

// Integrate a conserved variable over the interior cells (indices 2 .. N+1).
function integrate(arr, N, dx) {
  let sum = 0;
  for (let i = 2; i < N + 2; i++) sum += arr[i];
  return sum * dx;
}

// L1 error between two equal-length arrays, normalized by cell count.
function l1(a, b) {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += Math.abs(a[i] - b[i]);
  return s / a.length;
}

// ---- 1. Accuracy: Sod vs exact Riemann -------------------------------------

console.log("Accuracy (numerical vs exact Riemann):");

check("Sod density L1 error < 5e-3 at N=800", () => {
  const params = {
    N: 800, L: 1.0, gamma: 1.4, cfl: 0.4,
    diaphragm: 0.5, bcType: "transmissive", flux: "HLLC",
    leftState:  { rho: 1.0,   u: 0.0, p: 1.0 },
    rightState: { rho: 0.125, u: 0.0, p: 0.1 }
  };
  let state = S.initState(params);
  state = S.stepUntil(state, 0.2);
  const sol   = S.extract(state);
  const exact = S.exactRiemann(params.leftState, params.rightState,
                               params.gamma, params.diaphragm, 0.2, sol.x);
  assert(exact, "exact Riemann returned null (unexpected vacuum)");
  const err = l1(sol.rho, exact.rho);
  assert(err < 5e-3, "L1 density error " + err.toExponential(3) + " exceeds 5e-3");
});

check("Sod error decreases as the grid is refined", () => {
  const base = {
    L: 1.0, gamma: 1.4, cfl: 0.4, diaphragm: 0.5,
    bcType: "transmissive", flux: "HLLC",
    leftState:  { rho: 1.0,   u: 0.0, p: 1.0 },
    rightState: { rho: 0.125, u: 0.0, p: 0.1 }
  };
  const errAt = N => {
    let s = S.initState(Object.assign({ N }, base));
    s = S.stepUntil(s, 0.2);
    const sol = S.extract(s);
    const ex  = S.exactRiemann(base.leftState, base.rightState, base.gamma, 0.5, 0.2, sol.x);
    return l1(sol.rho, ex.rho);
  };
  const e100 = errAt(100), e400 = errAt(400);
  assert(e400 < e100, "error did not drop from N=100 (" + e100.toExponential(2) +
                      ") to N=400 (" + e400.toExponential(2) + ")");
});

// ---- 2. Positivity ---------------------------------------------------------

console.log("Positivity (density and pressure stay > 0):");

check("1D strong shock keeps rho, p > 0", () => {
  const p = S.PRESETS["Strong shock"];
  const params = {
    N: 400, L: 1.0, gamma: 1.4, cfl: 0.4,
    diaphragm: p.diaphragm, bcType: "transmissive", flux: "HLLC",
    leftState: p.L, rightState: p.R
  };
  let state = S.initState(params);
  state = S.stepUntil(state, p.tEnd);
  const sol = S.extract(state);
  for (let i = 0; i < sol.rho.length; i++) {
    assert(sol.rho[i] > 0 && Number.isFinite(sol.rho[i]), "rho non-positive/NaN at i=" + i);
    assert(sol.p[i]   > 0 && Number.isFinite(sol.p[i]),   "p non-positive/NaN at i=" + i);
  }
});

check("2D Sedov blast keeps rho, p > 0", () => {
  const preset = S2.PRESETS2D["Sedov blast"];
  const params = {
    Nx: 120, Ny: 120, Lx: preset.Lx, Ly: preset.Ly,
    gamma: preset.gamma, cfl: preset.cfl, flux: "HLLC",
    initFn: preset.init
  };
  let state = S2.initState2D(params);
  state = S2.stepUntil2D(state, 0.05);   // short march; the blast stresses positivity immediately
  const sol = S2.extract2D(state);
  for (let k = 0; k < sol.rho.length; k++) {
    assert(sol.rho[k] > 0 && Number.isFinite(sol.rho[k]), "rho non-positive/NaN at k=" + k);
    assert(sol.p[k]   > 0 && Number.isFinite(sol.p[k]),   "p non-positive/NaN at k=" + k);
  }
});

// ---- 3. Conservation -------------------------------------------------------
// A finite-volume scheme is conservative by construction: the only way a
// conserved quantity changes is through the domain boundary fluxes. For the Sod
// problem the boundary gas is at rest (u = 0), so the mass flux (rho*u) and
// energy flux (u(E+p)) through both ends are zero as long as the waves have not
// reached the boundaries. Total mass and energy are therefore conserved to
// machine precision. (Momentum is NOT conserved here — the boundary pressure
// exerts a net impulse — so we deliberately do not assert it.)

console.log("Conservation (mass and energy, transmissive, waves interior):");

check("Sod conserves mass and energy to ~1e-10", () => {
  const params = {
    N: 400, L: 1.0, gamma: 1.4, cfl: 0.4,
    diaphragm: 0.5, bcType: "transmissive", flux: "HLLC",
    leftState:  { rho: 1.0,   u: 0.0, p: 1.0 },
    rightState: { rho: 0.125, u: 0.0, p: 0.1 }
  };
  let state = S.initState(params);
  const m0 = integrate(state.rho, state.N, state.dx);
  const e0 = integrate(state.ene, state.N, state.dx);
  state = S.stepUntil(state, 0.2);       // by t=0.2 the waves are still well inside [0,1]
  const m1 = integrate(state.rho, state.N, state.dx);
  const e1 = integrate(state.ene, state.N, state.dx);
  const dm = Math.abs(m1 - m0) / Math.abs(m0);
  const de = Math.abs(e1 - e0) / Math.abs(e0);
  assert(dm < 1e-10, "mass drift " + dm.toExponential(3) + " exceeds 1e-10");
  assert(de < 1e-10, "energy drift " + de.toExponential(3) + " exceeds 1e-10");
});

// ---- summary ---------------------------------------------------------------

console.log("");
if (process.exitCode) {
  console.error("Some checks failed.");
} else {
  console.log("All " + passed + " checks passed.");
}
