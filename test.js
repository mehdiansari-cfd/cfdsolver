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
const Z  = require("./nozzle");
const FA = require("./facility");   // pulls in what it needs from ./solver itself
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

// ---- 8. Shock-tube facility ------------------------------------------------
// Calorically perfect gases must reproduce the classical closed forms
// (Anderson, Modern Compressible Flow Ch. 7; Gaydon & Hurle 1963); the
// thermally perfect analysis is cross-checked against the two-gas solver.

console.log("Shock-tube facility (classical theory and two-gas simulation):");

check("calorically perfect He -> air: shock-tube equation, p5/p2 and Mr relations hold", () => {
  const He = G.gasModel("He", "calorically-perfect"), air = G.gasModel("air", "calorically-perfect");
  const r = FA.facilityAnalysis({ driver: He, driven: air, p4: 20e5, T4: 300, p1: 5e3, T1: 300, L4: 2, L1: 5 });
  const g1 = G.gasGamma(air, 300), g4 = G.gasGamma(He, 300), a1 = G.gasSound(air, 300), a4 = G.gasSound(He, 300);
  const P = r.p2p1, Ms = r.Ms, Mr = r.reflected.Mr;
  const p4p1 = P * Math.pow(1 - (g4 - 1) * (a1 / a4) * (P - 1) / Math.sqrt(2 * g1 * (2 * g1 + (g1 + 1) * (P - 1))), -2 * g4 / (g4 - 1));
  assert(Math.abs(p4p1 / r.p4p1 - 1) < 1e-6, "p4/p1 " + p4p1 + " vs " + r.p4p1);
  const p5p2 = (((3 * g1 - 1) / (g1 - 1)) * P - 1) / (P + (g1 + 1) / (g1 - 1));
  assert(Math.abs(r.states[5].p / r.states[2].p / p5p2 - 1) < 1e-6, "p5/p2");
  const rhs = Ms / (Ms * Ms - 1) * Math.sqrt(1 + 2 * (g1 - 1) / ((g1 + 1) * (g1 + 1)) * (Ms * Ms - 1) * (g1 + 1 / (Ms * Ms)));
  assert(Math.abs(Mr / (Mr * Mr - 1) / rhs - 1) < 1e-6, "Mr relation");
});

check("facilityDesignP4 inverts the analysis", () => {
  const base = { driver: "He", driven: "air", T4: 300, p1: 5e3, T1: 300, L4: 2, L1: 5 };
  const p4 = FA.facilityDesignP4(base, 3.5);
  assert(Math.abs(FA.facilityAnalysis(Object.assign({}, base, { p4 })).Ms - 3.5) < 1e-6);
});

check("two-gas simulation reproduces t_reflection, p5, T5 and the under-tailored wave", () => {
  const cfg = { driver: "He", driven: "N2", T4: 300, p1: 5e3, T1: 300, L4: 3, L1: 6 };
  cfg.p4 = FA.facilityDesignP4(cfg, 4.2);
  const a = FA.facilityAnalysis(cfg), P = FA.facilitySimParams(cfg, 600);
  const tW = a.times.shockAtWall, tA = a.interaction.arriveT;
  assert.strictEqual(a.interaction.tailoring, "under-tailored");
  let s = S.initState(P);
  const wall = () => { const k = s.N + 1; return s.eos.prim(s.rho[k], s.mom[k], s.ene[k], s.rY[k]); };
  const hist = [];
  const tStop = tA + 0.3 * (a.expansion.arriveT - tA);
  while (s.t < tStop) { s = S.stepUntil(s, Math.min(tStop, s.t + 2e-5)); const q = wall(); hist.push({ t: s.t, p: q.p, T: q.T }); }
  const p5 = a.states[5].p, T5 = a.states[5].T;
  const tRef = hist.find(h => h.p > 0.5 * (a.states[2].p + p5)).t;
  assert(Math.abs(tRef / tW - 1) < 0.01, "reflection at " + tRef + " vs " + tW);
  const mid = hist.filter(h => h.t > tW + 0.3 * (tA - tW) && h.t < tW + 0.7 * (tA - tW));
  const pm = mid.reduce((acc, h) => acc + h.p, 0) / mid.length, Tm = mid.reduce((acc, h) => acc + h.T, 0) / mid.length;
  assert(Math.abs(pm / p5 - 1) < 0.01, "p5 plateau " + pm + " vs " + p5);
  assert(Math.abs(Tm / T5 - 1) < 0.01, "T5 plateau " + Tm + " vs " + T5);
  // the reflected shock from the contact raises the end-wall pressure afterwards
  assert(hist[hist.length - 1].p > 1.2 * p5, "no pressure rise after the interaction wave");
});

// ---- 9. Case files and the command-line runner --------------------------------

console.log("Case files, CSV and the CLI (node run.js):");

const CASE = require("./case");
const fs = require("fs"), os = require("os"), path = require("path");
const { execFileSync } = require("child_process");

check("case validation reports bad input clearly", () => {
  assert.throws(() => CASE.normalizeCase({ kind: "rocket" }), /kind/);
  assert.throws(() => CASE.normalizeCase({ kind: "shocktube", gas: "air", L: 1, diaphragm: 0.5, tEnd: 1e-3,
                                            left: { p: 1e5, u: 0 }, right: { p: 1e5, T: 300, u: 0 } }), /left: 'T'/);
  assert.throws(() => CASE.resolveGas("unobtainium"), /Unknown gas/);
});

check("run.js: sod.json → CSV + summary; --set and --sweep work", () => {
  const out = fs.mkdtempSync(path.join(os.tmpdir(), "cfd-run-"));
  const summary = JSON.parse(execFileSync(process.execPath, ["run.js", "cases/sod.json", "--out", out, "--quiet", "--set", "numerics.N=200"]));
  assert(summary.L1rhoRelative < 1e-2, "L1 " + summary.L1rhoRelative);
  const csv = fs.readFileSync(path.join(out, "sod.csv"), "utf8").trim().split("\n").filter(l => !l.startsWith("#"));
  assert.strictEqual(csv.length, 201, "200 rows + header");
  assert(/^x,rho,u,p,T,M,a,e,s$/.test(csv[0]), "header " + csv[0]);
  const sweep = execFileSync(process.execPath, ["run.js", "cases/sod.json", "--out", out, "--quiet",
                                                "--set", "numerics.N=100", "--sweep", "numerics.flux=\"HLLC\",\"Rusanov\""]).toString();
  const rows = sweep.trim().split("\n");
  assert.strictEqual(rows.length, 3, sweep);
  assert(fs.existsSync(path.join(out, "sod.sweep.csv")));
});

// ---- 10. Verification tooling: order of accuracy and GCI ---------------------

console.log("Verification (order of accuracy, grid convergence index):");

const V = require("./verify");

check("GCI recovers p = 2 and the exact limit from f = 1 + 0.3 h^2 + 0.05 h^3 (r = 2 and 1.5)", () => {
  const f = h => 1 + 0.3 * h * h + 0.05 * h * h * h;
  for (const r of [2, 1.5]) {
    const g = V.gci3(f(0.01), f(0.01 * r), f(0.01 * r * r), 0.01, 0.01 * r, 0.01 * r * r);
    assert(Math.abs(g.p - 2) < 0.01, "p = " + g.p);
    assert(Math.abs(g.fExt - 1) < 1e-6, "fExt = " + g.fExt);
    assert.strictEqual(g.convergence, "monotone");
  }
});

check("smooth periodic wave: MUSCL + SSP-RK2 reaches 2nd order in L1 (pairwise > 1.8)", () => {
  const o = V.orderStudy("smooth-advection", [50, 100, 200, 400]);
  const last = o.order.L1.pairwise[o.order.L1.pairwise.length - 1];
  assert(last > 1.8, "finest pairwise L1 order " + last.toFixed(3));
});

check("GCI bands contain the exact answer: shock position (nozzle) and rho_2 (real-gas tube)", () => {
  for (const name of ["nozzle-shock-position", "shocktube-plateau"]) {
    const r = V.gciStudy(name, [100, 200, 400]);
    assert(r.referenceInsideBand, name + ": exact " + r.reference + " outside " + r.band.join(" … "));
    assert(r.gci.convergence === "monotone", name + ": " + r.gci.convergence);
  }
});

// ---- 11. Boundary layer: integral turbulent method --------------------------
// Flat plate at low Mach: momentum thickness against the 1/7-power law
// θ = 0.036 x Re_x^-0.2, and Hᵢ settling near its zero-pressure-gradient 1.4.

const BLM = require("./boundary");
console.log("Boundary layer (boundary.js):");

check("flat plate: θ(x) within 15% of the 1/7-power law, Hi → 1.4", () => {
  const n = 400, Lp = 2, ue = 30, Te = 300, rhoe = 1.2, gam = 1.4, Me = ue / Math.sqrt(gam * 287 * Te);
  const x = Array.from({ length: n + 1 }, (_, i) => (i / n) * Lp);
  const r = BLM.integrateBoundaryLayer({ x, ue: x.map(() => ue), Me: x.map(() => Me), rhoe: x.map(() => rhoe),
                                         Te: x.map(() => Te), gamma: x.map(() => gam) });
  const nu = BLM.sutherlandAir(Te) / rhoe;
  for (const i of [n / 4, n / 2, n]) {
    const law = 0.036 * x[i] * Math.pow(ue * x[i] / nu, -0.2);
    assert(Math.abs(r.theta[i] / law - 1) < 0.15, "θ/θ_1/7 at x = " + x[i] + ": " + (r.theta[i] / law).toFixed(3));
  }
  assert(Math.abs(r.Hi[n] - 1.4) < 0.1, "Hi at the end " + r.Hi[n].toFixed(3));
  assert(r.separation.length === 0, "no separation on a flat plate");
});

check("coupled CD nozzle, no shock: converges, δ* > 0, choked mass flow drops", () => {
  const gas = G.gasModel("air"), p0 = 5e5, T0 = 300;
  const prm = { type: "CD", gas, p0, T0, pb: 0.05 * p0, L: 0.3, AR: 2.5, Athroat: 1e-3 };
  const v = BLM.viscousQuasi1D(prm, { geometry: "channel", span: 0.05, N: 150 });
  assert(v.converged, "coupling converged");
  assert(v.deltaStar.every(d => d >= 0) && v.deltaStar[v.deltaStar.length - 1] > 0, "δ* positive");
  const inv = Z.quasi1DSolution(Z.initNozzleState(Object.assign({}, prm, { N: 150 })), "auto");
  assert(v.solution.mdotChoked < inv.mdotChoked, "boundary layer blocks the throat");
});

// ---- 12. Experimental validation: Sajben transonic diffuser -------------------
// NPARC wall-pressure data (validation/sajben). Quasi-1D inviscid flow must
// match the isentropic acceleration upstream of the shock and the recovered
// pressure downstream. It puts the shock too far aft (no boundary-layer
// displacement or separation); the displacement correction must pull the
// weak shock back to within VISC_SHOCK_TOL of the data. The strong shock
// separates the layer, which a direct integral method cannot follow: there it
// must say so (reliable: false) rather than report its shock position as an answer.

console.log("Validation (Sajben diffuser, NPARC data):");

check("Sajben: upstream & recovery match data; inviscid shock aft; displacement correction moves it to the data", () => {
  const SJ = require("./validation/sajben.js");
  const r = SJ.main({ N: 200, write: false, quiet: true, theta0: [0], blN: 150 });
  if (r === null) { console.log("       (data files absent: skipped)"); return; }
  const cellH = 0.13;                                            // two cells at N = 200, in throat heights
  const VISC_SHOCK_TOL = 0.4;                                     // throat heights
  for (const [name, c] of Object.entries(r)) {
    assert(Math.abs(c.shock_xH_computed - c.shock_xH_exactQuasi1D) < cellH, name + ": numerical shock vs exact quasi-1D");
    for (const [w, v] of Object.entries(c.walls)) {
      assert(v.upstream.rmsDeltaP < 0.03, name + "/" + w + ": upstream RMS Δp/p0 " + v.upstream.rmsDeltaP.toFixed(4));
      assert(v.downstream.rmsDeltaP < 0.05, name + "/" + w + ": recovered-flow RMS Δp/p0 " + v.downstream.rmsDeltaP.toFixed(4));
      assert(v.shock_xH_error > 0, name + "/" + w + ": inviscid shock should sit aft of the measured one");
    }
  }
  const weak = r.weak.viscous[0];
  assert(weak.reliable, "weak: attached and matched");
  assert(weak.shock_xH < r.weak.shock_xH_exactQuasi1D, "weak: displacement moves the shock upstream");
  for (const [w, v] of Object.entries(weak.walls))
    assert(Math.abs(v.shock_xH_error) < VISC_SHOCK_TOL, "weak/" + w + ": viscous shock error " + v.shock_xH_error.toFixed(3) + " H");
  assert(!r.strong.viscous[0].reliable && r.strong.viscous[0].separated, "strong: separation flagged, result marked unreliable");
});

// ---- summary ---------------------------------------------------------------

console.log("");
if (process.exitCode) {
  console.error("Some checks failed.");
} else {
  console.log("All " + passed + " checks passed.");
}
