/* ============================================================
   verify.js  —  code verification and solution verification
   ------------------------------------------------------------
   Order of accuracy (code verification)
     orderStudy(name, Ns)        errors vs an exact solution on a grid
                                 sequence + observed order (pairwise and
                                 least-squares)
       "smooth-advection"   rho = 1 + 0.2 sin(2 pi x), u = 1, p = 1,
                            periodic, one period: exact = initial state.
                            Smooth, so it measures the scheme's design
                            order (2 for MUSCL + SSP-RK2, a little lower
                            in L-inf where minmod clips the extrema).
       "sod"                Sod vs the exact Riemann solution. Shocks and
                            contacts cap the L1 order near 1 for *every*
                            scheme, so this is a correctness check, not an
                            order test.

   Grid convergence index (solution verification)
     gci3(f1, f2, f3, h1, h2, h3)
       Celik, Ghia, Roache, Freitas, Coleman & Raad, "Procedure for
       Estimation and Reporting of Uncertainty Due to Discretization in
       CFD Applications", J. Fluids Eng. 130 (2008) 078001. Grid 1 is the
       finest. Returns the observed order p, the Richardson-extrapolated
       value, relative errors and GCI_fine (safety factor 1.25).
     gciStudy(name, Ns)          runs a scalar functional on three grids
       "nozzle-mdot-unchoked"   mass flow, unchoked CD nozzle (non-dim)
       "nozzle-shock-position"  normal-shock location, CD nozzle (non-dim)
       "shocktube-plateau"      contact-plateau density, real-gas air tube
   Every run is a generator that yields progress (0..1) between solver
   chunks, so a browser page can drive it in time slices without freezing
   (labui.js runSliced); orderStudy / gciStudy drive them synchronously.
   Math layer: no DOM. Node: require("./verify").
   ============================================================ */

if (typeof initState === "undefined" && typeof require === "function") {
  var _vS = require("./solver");
  var initState = _vS.initState, stepUntil = _vS.stepUntil, extract = _vS.extract, exactRiemann = _vS.exactRiemann;
  var exactRiemannGas = _vS.exactRiemannGas;
  var _vZ = require("./nozzle");
  var initNozzleState = _vZ.initNozzleState, nozzleStepUntil = _vZ.nozzleStepUntil;
  var extractNozzle = _vZ.extractNozzle, quasi1DSolution = _vZ.quasi1DSolution;
  var gasModel = require("./gas").gasModel;
}

const GCI_SAFETY = 1.25;        // Celik et al. 2008: Fs = 1.25 for three-grid studies
const CHUNK_T = 0.02;           // fraction of a run between progress yields
const STEADY_MAX_TRANSITS = 400;
const STEADY_SUBCHUNKS = 40;    // solver calls per flow-through time between yields
const SHOCK_TOL = 1e-6;         // steady tolerance for the shock-position study
const GCI_P_ITERS = 100;
const SMOOTH_AMPLITUDE = 0.2;

// ---------- Norms ----------
function errorNorms(num, ex, dx) {
  let l1 = 0, l2 = 0, li = 0;
  for (let i = 0; i < num.length; i++) {
    const d = Math.abs(num[i] - ex[i]);
    l1 += d; l2 += d * d; li = Math.max(li, d);
  }
  const n = num.length;
  return { L1: l1 / n, L2: Math.sqrt(l2 / n), Linf: li };
}

// Least-squares slope of log(err) against log(h), and the pairwise orders.
function observedOrder(hs, errs) {
  const X = hs.map(Math.log), Y = errs.map(Math.log), n = X.length;
  const mx = X.reduce((a, b) => a + b) / n, my = Y.reduce((a, b) => a + b) / n;
  let sxy = 0, sxx = 0;
  for (let i = 0; i < n; i++) { sxy += (X[i] - mx) * (Y[i] - my); sxx += (X[i] - mx) * (X[i] - mx); }
  const pairwise = [];
  for (let i = 1; i < n; i++) pairwise.push(Math.log(errs[i - 1] / errs[i]) / Math.log(hs[i - 1] / hs[i]));
  return { slope: sxy / sxx, pairwise };
}

// ---------- Three-grid GCI (Celik et al. 2008, Eqs. 1–10) ----------
function gci3(f1, f2, f3, h1, h2, h3) {
  const r21 = h2 / h1, r32 = h3 / h2;
  const e21 = f2 - f1, e32 = f3 - f2;
  if (e21 === 0 || e32 === 0) return { p: NaN, fExt: f1, ea21: 0, eext21: 0, gciFine: 0, convergence: "exact" };
  const R = e21 / e32;                                       // convergence ratio
  const s = Math.sign(e32 / e21);
  // fixed-point iteration for the apparent order (Eqs. 3a–3c)
  let p = Math.abs(Math.log(Math.abs(e32 / e21))) / Math.log(r21);
  for (let it = 0; it < GCI_P_ITERS; it++) {
    const q = Math.log((Math.pow(r21, p) - s) / (Math.pow(r32, p) - s));
    const pn = Math.abs(Math.log(Math.abs(e32 / e21)) + q) / Math.log(r21);
    if (!Number.isFinite(pn)) break;
    if (Math.abs(pn - p) < 1e-12) { p = pn; break; }
    p = pn;
  }
  const rp = Math.pow(r21, p);
  const fExt = (rp * f1 - f2) / (rp - 1);                    // Richardson extrapolation (Eq. 4)
  const ea21 = Math.abs((f1 - f2) / f1);
  const eext21 = Math.abs((fExt - f1) / fExt);
  const gciFine = GCI_SAFETY * ea21 / (rp - 1);
  const ea32 = Math.abs((f2 - f3) / f2), gciCoarse = GCI_SAFETY * ea32 / (Math.pow(r32, p) - 1);
  return {
    p, fExt, ea21, eext21, gciFine, gciCoarse,
    // ≈ 1 when the grids are in the asymptotic range (Roache 1998)
    asymptoticRatio: gciCoarse / (Math.pow(r21, p) * gciFine),
    convergence: R > 0 && R < 1 ? "monotone" : R < 0 && Math.abs(R) < 1 ? "oscillatory" : "divergent",
    r21, r32
  };
}

// ---------- Order-of-accuracy problems ----------
const ORDER_STUDIES = {
  "smooth-advection": {
    label: "Smooth density wave, periodic (design order)",
    expected: 2,
    run: function* (N, flux) {
      const P = { N, L: 1, gamma: 1.4, cfl: 0.4, diaphragm: 2, bcType: "periodic", flux: flux || "HLLC",
                  leftState: { rho: 1, u: 1, p: 1 }, rightState: { rho: 1, u: 1, p: 1 } };
      let s = initState(P);
      const x0 = i => (i + 0.5) / N;
      for (let i = 0; i < N; i++) {                           // overwrite with the smooth profile
        const rho = 1 + SMOOTH_AMPLITUDE * Math.sin(2 * Math.PI * x0(i));
        s.rho[i + 2] = rho; s.mom[i + 2] = rho; s.ene[i + 2] = 1 / 0.4 + 0.5 * rho;
      }
      while (s.t < 1) { s = stepUntil(s, Math.min(1, s.t + CHUNK_T)); yield s.t; }   // one period
      const sol = extract(s);
      const exact = sol.x.map(x => 1 + SMOOTH_AMPLITUDE * Math.sin(2 * Math.PI * x));
      return { x: sol.x, num: sol.rho, exact, steps: s.step };
    }
  },
  "sod": {
    label: "Sod shock tube vs exact Riemann (discontinuous: order ≤ 1)",
    expected: 1,
    run: function* (N, flux) {
      const L = { rho: 1, u: 0, p: 1 }, R = { rho: 0.125, u: 0, p: 0.1 };
      let s = initState({ N, L: 1, gamma: 1.4, cfl: 0.4, diaphragm: 0.5, bcType: "transmissive",
                          flux: flux || "HLLC", leftState: L, rightState: R });
      while (s.t < 0.2) { s = stepUntil(s, Math.min(0.2, s.t + 0.2 * CHUNK_T)); yield s.t / 0.2; }
      const sol = extract(s);
      return { x: sol.x, num: sol.rho, exact: exactRiemann(L, R, 1.4, 0.5, 0.2, sol.x).rho, steps: s.step };
    }
  }
};

// Run a generator to completion synchronously; returns its return value.
function runToEnd(gen) {
  for (;;) { const r = gen.next(); if (r.done) return r.value; }
}

// Generator: yields { grid, of, progress } and returns the study result.
function* orderStudyGen(name, Ns, flux) {
  const st = ORDER_STUDIES[name];
  if (!st) throw new Error("Unknown order study '" + name + "'");
  const rows = [];
  for (let g = 0; g < Ns.length; g++) {
    const N = Ns[g], it = st.run(N, flux);
    let r = it.next();
    while (!r.done) { yield { grid: g, of: Ns.length, N, progress: r.value }; r = it.next(); }
    const e = errorNorms(r.value.num, r.value.exact);
    rows.push({ N, h: 1 / N, L1: e.L1, L2: e.L2, Linf: e.Linf, steps: r.value.steps });
  }
  const hs = rows.map(r => r.h);
  const ord = k => rows.length > 1 ? observedOrder(hs, rows.map(r => r[k])) : null;
  return { name, label: st.label, expected: st.expected, rows, order: { L1: ord("L1"), L2: ord("L2"), Linf: ord("Linf") } };
}
function orderStudy(name, Ns, flux) { return runToEnd(orderStudyGen(name, Ns, flux)); }

// ---------- GCI problems (scalar functionals) ----------
// March a nozzle to steady state: stop when the mean |Δρ|/ρ over one
// flow-through time falls below tol. Yields progress as -log10(residual)/-log10(tol).
function* _nozzleRunSteady(q, N, tol) {
  let s = initNozzleState(Object.assign({}, q, { N }));
  const transit = s.L / Math.sqrt(1.4 * s.T0);
  let prev = s.rho.slice(), t = transit, progress = 0;
  for (let k = 0; k < STEADY_MAX_TRANSITS; k++) {
    for (let j = 1; j <= STEADY_SUBCHUNKS; j++) {             // yield often: keeps a page responsive
      s = nozzleStepUntil(s, t - transit + j * transit / STEADY_SUBCHUNKS);
      yield progress;
    }
    let d = 0, m = 0;
    for (let i = 2; i < N + 2; i++) { d += Math.abs(s.rho[i] - prev[i]); m += s.rho[i]; }
    progress = Math.max(0, Math.min(0.99, Math.log10(m / d) / Math.log10(1 / tol)));
    yield progress;
    if (d / m < tol) break;
    prev = s.rho.slice(); t += transit;
  }
  return s;
}

const GCI_STUDIES = {
  "nozzle-mdot-unchoked": {
    label: "Mass flow, unchoked CD nozzle (AR 2, pb/p0 0.97)",
    quantity: "mdot", unit: "",
    reference: () => {
      const s = initNozzleState({ type: "CD", AR: 2, p0: 1, T0: 1, pb: 0.97, N: 50, L: 1, gamma: 1.4 });
      return quasi1DSolution(s, "auto").mdot[0];
    },
    run: function* (N) {
      const s = yield* _nozzleRunSteady({ type: "CD", AR: 2, p0: 1, T0: 1, pb: 0.97, L: 1, gamma: 1.4, cfl: 0.4 }, N, 1e-6);
      const e = extractNozzle(s);
      return e.mdot.reduce((a, b) => a + b) / e.mdot.length;   // duct average
    }
  },
  "nozzle-shock-position": {
    label: "Normal-shock position, CD nozzle (AR 2, pb/p0 0.75)",
    quantity: "x_shock", unit: "L",
    reference: () => {
      const s = initNozzleState({ type: "CD", AR: 2, p0: 1, T0: 1, pb: 0.75, N: 50, L: 1, gamma: 1.4 });
      return quasi1DSolution(s, "auto").xShock;
    },
    run: function* (N) {
      const s = yield* _nozzleRunSteady({ type: "CD", AR: 2, p0: 1, T0: 1, pb: 0.75, L: 1, gamma: 1.4, cfl: 0.35 }, N, SHOCK_TOL);
      const e = extractNozzle(s);
      // shock = where p crosses the mean of the pre- and post-shock plateaus,
      // located to sub-cell accuracy by linear interpolation
      let k = 1, best = 0;
      for (let i = 1; i < N; i++) if (e.p[i] - e.p[i - 1] > best) { best = e.p[i] - e.p[i - 1]; k = i; }
      const lo = e.p[Math.max(0, k - 4)], hi = e.p[Math.min(N - 1, k + 4)], mid = 0.5 * (lo + hi);
      for (let i = Math.max(1, k - 4); i <= Math.min(N - 1, k + 4); i++) {
        if (e.p[i - 1] < mid && e.p[i] >= mid) return e.x[i - 1] + (mid - e.p[i - 1]) / (e.p[i] - e.p[i - 1]) * (e.x[i] - e.x[i - 1]);
      }
      return e.x[k];
    }
  },
  "shocktube-plateau": {
    label: "Density between contact and shock, thermally perfect air 10:1 bar",
    quantity: "rho_2", unit: "kg/m³",
    reference: () => {
      const air = gasModel("air");
      return exactRiemannGas(air, { p: 1e6, T: 300, u: 0 }, air, { p: 1e5, T: 300, u: 0 }, 0.5, 6e-4, [0.75]).rho[0];
    },
    run: function* (N) {
      const air = gasModel("air");
      let s = initState({ N, L: 1, cfl: 0.4, diaphragm: 0.5, bcType: "transmissive", flux: "HLLC", gas: air,
                          leftState: { p: 1e6, T: 300, u: 0 }, rightState: { p: 1e5, T: 300, u: 0 } });
      while (s.t < 6e-4) { s = stepUntil(s, Math.min(6e-4, s.t + 6e-4 * CHUNK_T)); yield s.t / 6e-4; }
      const e = extract(s);
      // average over the middle of the plateau (x in [0.72, 0.78] at t = 6e-4)
      let sum = 0, n = 0;
      for (let i = 0; i < N; i++) if (e.x[i] > 0.72 && e.x[i] < 0.78) { sum += e.rho[i]; n++; }
      return sum / n;
    }
  }
};

// Generator over three (or more) grids; the GCI uses the three finest.
function* gciStudyGen(name, Ns) {
  const st = GCI_STUDIES[name];
  if (!st) throw new Error("Unknown GCI study '" + name + "'");
  const grids = Ns.slice().sort((a, b) => a - b);
  if (grids.length < 3) throw new Error("A GCI study needs at least three grids.");
  const rows = [];
  for (let g = 0; g < grids.length; g++) {
    const N = grids[g], it = st.run(N);
    let r = it.next();
    while (!r.done) { yield { grid: g, of: grids.length, N, progress: r.value }; r = it.next(); }
    rows.push({ N, h: 1 / N, f: r.value });
  }
  const [g3, g2, g1] = rows.slice(-3);                         // coarse, medium, fine
  const g = gci3(g1.f, g2.f, g3.f, g1.h, g2.h, g3.h);
  const ref = st.reference ? st.reference() : null;
  return { name, label: st.label, quantity: st.quantity, unit: st.unit, rows, gci: g,
           fine: g1.f, band: [g1.f - g.gciFine * Math.abs(g1.f), g1.f + g.gciFine * Math.abs(g1.f)], reference: ref,
           referenceInsideBand: ref == null ? null : Math.abs(ref - g1.f) <= g.gciFine * Math.abs(g1.f) };
}
function gciStudy(name, Ns) { return runToEnd(gciStudyGen(name, Ns)); }

if (typeof module !== "undefined" && module.exports) {
  module.exports = { errorNorms, observedOrder, gci3, runToEnd,
                     ORDER_STUDIES, orderStudyGen, orderStudy, GCI_STUDIES, gciStudyGen, gciStudy };
}
