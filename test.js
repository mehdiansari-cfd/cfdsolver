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
     4. Nozzle      — the quasi-1D nozzle loads standalone in Node,
                      the reservoir inlet holds total pressure p0,
                      the choked mass flow matches theory, and the
                      design-point flow matches the isentropic
                      solution.
   ============================================================ */

const assert = require("assert");
const S  = require("./solver");
const G  = require("./gas");
const Z  = require("./nozzle");   // pulls in what it needs from ./solver itself
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

// ---- 4. Nozzle -------------------------------------------------------------
// Steady quasi-1D theory (Anderson, Modern Compressible Flow Ch. 5), R = 1:
// a choked nozzle passes mdot = p0 A* sqrt(gamma/T0) (2/(gamma+1))^((gamma+1)/(2(gamma-1))),
// and wherever the flow is shock-free the total pressure stays at p0.

console.log("Nozzle (reservoir inlet, choking, isentropic design point):");

function runNozzle(name) {
  const q = Z.NOZZLE_PRESETS[name];
  let s = Z.initNozzleState(q);
  s = Z.nozzleStepUntil(s, q.tEnd);
  return { q, s, sol: Z.extractNozzle(s) };
}

check("nozzle.js loads in Node without solver globals", () => {
  for (const fn of ["initNozzleState", "nozzleStepUntil", "extractNozzle", "isentropicSolution"]) {
    assert.strictEqual(typeof Z[fn], "function", fn + " not exported");
  }
});

check("Choked CD: inlet total pressure = p0, throat sonic, mdot = choked theory", () => {
  const { q, sol } = runNozzle("Choked CD");
  const g = q.gamma, N = sol.x.length;
  const ptIn = sol.pt[0];
  assert(Math.abs(ptIn - q.p0) / q.p0 < 5e-3, "inlet total pressure " + ptIn.toFixed(4) + " != p0");
  const Mthroat = sol.M[N >> 1];
  assert(Math.abs(Mthroat - 1) < 0.03, "throat Mach " + Mthroat.toFixed(3) + " not sonic");
  const mdotTheory = q.p0 * Math.sqrt(g / q.T0) * Math.pow(2 / (g + 1), (g + 1) / (2 * (g - 1)));
  for (const i of [0, N >> 1, N - 1]) {
    const mdot = sol.rho[i] * sol.u[i] * sol.A[i];
    assert(Math.abs(mdot - mdotTheory) / mdotTheory < 0.01,
           "mdot " + mdot.toFixed(4) + " at i=" + i + " vs theory " + mdotTheory.toFixed(4));
  }
});

check("Full Supersonic CD matches the isentropic solution", () => {
  const { q, s, sol } = runNozzle("Full Supersonic CD");
  const iso = Z.isentropicSolution(s, q.branch);
  const N = sol.x.length;
  for (let i = 0; i < N; i++) {
    assert(Math.abs(sol.M[i] - iso.M[i]) < 0.03, "Mach off by " + Math.abs(sol.M[i] - iso.M[i]).toFixed(3) + " at i=" + i);
    assert(Math.abs(sol.p[i] - iso.p[i]) / iso.p[i] < 0.03, "pressure off by " +
           (100 * Math.abs(sol.p[i] - iso.p[i]) / iso.p[i]).toFixed(1) + "% at i=" + i);
  }
});

// ---- 5. Gas library ---------------------------------------------------------
// Reference values: JANAF / NIST cp at 298.15 K [J/(mol K)], textbook air gamma.

console.log("Gas library (NASA-7 thermo, SI):");

check("cp(298.15 K) matches JANAF for N2, O2, Ar, H2, CO2 to 0.1%", () => {
  const ref = { N2: 29.124, O2: 29.376, Ar: 20.786, H2: 28.836, CO2: 37.135 };
  for (const k of Object.keys(ref)) {
    const g = G.GASES[k], cpMol = g.cp(298.15) * g.W;
    assert(Math.abs(cpMol / ref[k] - 1) < 1e-3, k + " cp = " + cpMol.toFixed(3) + " vs " + ref[k]);
  }
});

check("air: R = 287.05, gamma = 1.400 at 298 K and 1.336 at 1000 K", () => {
  const air = G.GASES.air;
  assert(Math.abs(air.R - 287.05) < 0.05, "R = " + air.R);
  assert(Math.abs(G.gasGamma(air, 298.15) - 1.400) < 1e-3, "gamma(298) = " + G.gasGamma(air, 298.15));
  assert(Math.abs(G.gasGamma(air, 1000) - 1.336) < 2e-3, "gamma(1000) = " + G.gasGamma(air, 1000));
});

check("thermo inverses round-trip: T(e(T)), isentrope, stagnation", () => {
  const air = G.GASES.air;
  for (const T of [120, 300, 1500, 5000, 8000]) {
    assert(Math.abs(G.gasTfromE(air, G.gasE(air, T)) / T - 1) < 1e-9, "T(e) at " + T);
  }
  const st = G.gasStagnation(air, 250, 5e4, 400), back = G.gasStaticFromStagnation(air, st.T0, st.p0, 400);
  assert(Math.abs(back.T - 250) < 1e-6 && Math.abs(back.p / 5e4 - 1) < 1e-9, "stagnation round trip");
});

// ---- 6. Real-gas exact Riemann and numerics -----------------------------------

console.log("Real gases (exact Riemann, thermally perfect, two-gas):");

// Calorically perfect gas with R = 1, cp = 3.5 (gamma = 1.4): the general solver
// must reproduce the closed-form ideal-gas one.
const CPG_TEST = { name: "test", R: 1, cp: () => 3.5, h: T => 3.5 * (T - 298.15), phi: T => 3.5 * Math.log(T / 298.15) };

check("general exact Riemann = ideal-gas exact Riemann on all five presets", () => {
  const x = []; for (let i = 0; i < 400; i++) x.push((i + 0.5) / 400);
  for (const name of Object.keys(S.PRESETS)) {
    const pr = S.PRESETS[name];
    const a = S.exactRiemann(pr.L, pr.R, 1.4, 0.5, pr.tEnd, x);
    const b = S.exactRiemannGas(CPG_TEST, pr.L, CPG_TEST, pr.R, 0.5, pr.tEnd, x);
    for (let i = 0; i < x.length; i++) {
      const d = Math.max(Math.abs(a.rho[i] / b.rho[i] - 1), Math.abs(a.p[i] / b.p[i] - 1));
      assert(d < 5e-5, name + ": rel. difference " + d.toExponential(2) + " at i=" + i);
    }
  }
});

check("calorically perfect air reproduces the normal-shock T2 formula", () => {
  const air = G.gasModel("air", "calorically-perfect"), He = G.gasModel("He");
  const sol = S.riemannGasSolve(He, { p: 20e5, T: 300, u: 0 }, air, { p: 1e4, T: 300, u: 0 });
  const g = G.gasGamma(air, 300), M = sol.right.mach;
  const T2 = 300 * (2 * g * M * M - (g - 1)) * ((g - 1) * M * M + 2) / ((g + 1) * (g + 1) * M * M);
  assert(Math.abs(sol.right.T / T2 - 1) < 1e-6, "T2 " + sol.right.T + " vs " + T2);
});

function runRealGas(P, tEnd) {
  let s = S.initState(P);
  const sum = a => { let t = 0; for (let i = 2; i < P.N + 2; i++) t += a[i]; return t; };
  const m0 = sum(s.rho), e0 = sum(s.ene);
  while (s.t < tEnd) s = S.stepUntil(s, tEnd);
  const sol = S.extract(s);
  const ex = S.exactRiemannGas(P.gasLeft || P.gas, P.leftState, P.gasRight || P.gas, P.rightState, P.diaphragm, tEnd, sol.x);
  let l1 = 0, rmax = 0;
  for (let i = 0; i < P.N; i++) { l1 += Math.abs(sol.rho[i] - ex.rho[i]); rmax = Math.max(rmax, ex.rho[i]); }
  return { err: l1 / P.N / rmax, dm: Math.abs(sum(s.rho) / m0 - 1), de: Math.abs(sum(s.ene) / e0 - 1), sol };
}

check("thermally perfect air shock tube (10 : 1 bar) matches exact, conserves mass & energy", () => {
  const r = runRealGas({ N: 400, L: 1, cfl: 0.4, diaphragm: 0.5, bcType: "transmissive", flux: "HLLC",
                         gas: G.GASES.air, leftState: { p: 10e5, T: 300, u: 0 }, rightState: { p: 1e5, T: 300, u: 0 } }, 6e-4);
  assert(r.err < 5e-3, "L1 rho error " + r.err.toExponential(2));
  assert(r.dm < 1e-10 && r.de < 1e-10, "drift mass " + r.dm.toExponential(1) + " energy " + r.de.toExponential(1));
});

check("two-gas He -> air (Ms ~ 4.1) with every flux matches exact, Y stays in [0, 1]", () => {
  for (const flux of Object.keys(S.FLUXES)) {
    const r = runRealGas({ N: 400, L: 1, cfl: 0.4, diaphragm: 0.5, bcType: "transmissive", flux,
                           gasLeft: G.GASES.He, gasRight: G.GASES.air,
                           leftState: { p: 20e5, T: 300, u: 0 }, rightState: { p: 1e4, T: 300, u: 0 } }, 2e-4);
    assert(r.err < 6e-3, flux + ": L1 rho error " + r.err.toExponential(2));
    assert(r.dm < 1e-10 && r.de < 1e-10, flux + ": conservation drift");
    for (const y of r.sol.Y) assert(y >= 0 && y <= 1, flux + ": Y out of range " + y);
  }
});

// ---- 7. Nozzle: exact quasi-1D in every regime, any gas ----------------------

console.log("Nozzle vs exact quasi-1D theory (all regimes, real gases):");

function nozzleVsExact(name) {
  const q = Z.NOZZLE_PRESETS[name];
  let st = Z.initNozzleState(q);
  while (st.t < q.tEnd) st = Z.nozzleStepUntil(st, q.tEnd);
  const e = Z.extractNozzle(st), ex = Z.quasi1DSolution(st, "auto");
  let dM = 0, dp = 0, xsNum = null, jump = 0;
  for (let i = 0; i < e.x.length; i++) {
    if (i > 0 && e.p[i] - e.p[i - 1] > jump) { jump = e.p[i] - e.p[i - 1]; xsNum = e.x[i]; }
    if (ex.xShock != null && Math.abs(e.x[i] - ex.xShock) < 0.04 * st.L) continue;   // skip the shock itself
    dM = Math.max(dM, Math.abs(e.M[i] - ex.M[i]));
    dp = Math.max(dp, Math.abs(e.p[i] - ex.p[i]) / st.p0);
  }
  return { ex, dM, dp, xsNum, mdotRatio: e.mdot[e.x.length >> 1] / ex.mdot[0], dx: st.dx };
}

check("Normal Shock preset: shock within 2 cells of the exact quasi-1D location", () => {
  const r = nozzleVsExact("Normal Shock");
  assert.strictEqual(r.ex.regime, "shock");
  assert(Math.abs(r.xsNum - r.ex.xShock) < 2 * r.dx, "shock at " + r.xsNum + " vs exact " + r.ex.xShock.toFixed(4));
  assert(r.dM < 0.01 && r.dp < 0.01, "away from the shock: |dM| " + r.dM.toFixed(4) + ", |dp|/p0 " + r.dp.toFixed(4));
});

check("SI air M2 tunnel (thermally perfect): matches exact quasi-1D, mdot to 0.2%", () => {
  const r = nozzleVsExact("Air tunnel M2 (SI)");
  assert.strictEqual(r.ex.regime, "supersonic");
  assert(r.dM < 0.01 && r.dp < 0.005, "|dM| " + r.dM.toFixed(4) + ", |dp|/p0 " + r.dp.toFixed(4));
  assert(Math.abs(r.mdotRatio - 1) < 2e-3, "mdot ratio " + r.mdotRatio);
});

check("tabulated area profile (PCHIP) reproduces the analytic CD nozzle", () => {
  const xs = [], As = [];
  for (let k = 0; k <= 40; k++) { const x = k / 40; xs.push(x); As.push(Z.nozzleArea(x, 1, 2, "CD")); }
  const tab = Z.makeAreaTable(xs, As);
  for (let k = 0; k <= 400; k++) {
    const x = k / 400;
    assert(Math.abs(tab.A(x) - Z.nozzleArea(x, 1, 2, "CD")) < 2e-3, "A mismatch at x=" + x);
  }
});

// ---- summary ---------------------------------------------------------------

console.log("");
if (process.exitCode) {
  console.error("Some checks failed.");
} else {
  console.log("All " + passed + " checks passed.");
}
