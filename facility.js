/* ============================================================
   facility.js  —  shock-tube facility analysis (math layer, no DOM)
   ------------------------------------------------------------
   Ideal shock-tube theory for real (calorically or thermally
   perfect) gases, built from the exact Riemann solver in solver.js:

       driver (4)          | driven (1)                 end wall
     x = -L4 ............. 0 ........................... x = L1

   1. Diaphragm burst  = Riemann problem (4 | 1)
        → incident shock (Ms, W), shocked driven gas (2),
          contact, expanded driver gas (3), expansion fan.
   2. Shock reflection = Riemann problem (2 | mirror of 2)
        → reflected shock, stagnant test gas (5) at the end wall.
   3. Reflected shock meets the contact = Riemann problem (3 | 5)
        → the wave sent back to the end wall: shock (under-tailored),
          rarefaction (over-tailored) or none (tailored).
   4. Expansion head reflected off the driver end wall, integrated
        as a C+ characteristic through the exact fan (an estimate:
        later wave interactions are not followed).
   Test time = first disturbance (3 or 4) reaching the end wall
               minus the reflected-shock start time.

   Assumptions: inviscid, ideal (instant) diaphragm burst, no boundary
   layer growth (real test times are shorter), frozen composition.
   For the full unsteady flow, including every wave interaction, run the
   two-gas numerical solver (see facilitySimParams) and read its x–t diagram.

   References: Gaydon & Hurle, "The Shock Tube in High-Temperature Chemical
   Physics" (1963); Anderson, "Modern Compressible Flow" Ch. 7.

   Public API: facilityAnalysis(cfg), facilityDesignP4(cfg, MsTarget),
               facilitySimParams(cfg, N), FACILITY_PRESETS
   cfg = { driver, driven: gas.js models, p4, T4, p1, T1 [Pa, K],
           L4, L1 [m] }
   ============================================================ */

if (typeof riemannGasSolve === "undefined" && typeof require === "function") {
  var _solverF        = require("./solver");
  var riemannGasSolve = _solverF.riemannGasSolve;
  var _gasF           = require("./gas");
  var gasSound        = _gasF.gasSound;
  var gasGamma        = _gasF.gasGamma;
  var gasModel        = _gasF.gasModel;
}

const TAILORED_TOL = 0.02;     // |p*/p5 - 1| below this counts as tailored
const CHAR_STEPS   = 4000;     // RK2 steps for the reflected-expansion characteristic
const DESIGN_ITERS = 100;

function _gasOf(g, thermo) { return typeof g === "string" ? gasModel(g, thermo) : g; }

// Full thermodynamic description of a uniform state.
function _describe(gas, p, T, u) {
  const rho = p / (gas.R * T), a = gasSound(gas, T);
  return { p, T, rho, u, a, M: u / a, gamma: gasGamma(gas, T), h: gas.h(T), gas: gas.name };
}

function facilityAnalysis(cfgIn) {
  const cfg = Object.assign({}, cfgIn);
  const driver = _gasOf(cfg.driver, cfg.thermo), driven = _gasOf(cfg.driven, cfg.thermo);
  const { p4, T4, p1, T1, L4, L1 } = cfg;
  if (!(p4 > p1)) throw new Error("The driver pressure p4 must exceed the driven pressure p1.");

  // 1. Diaphragm burst
  const inc = riemannGasSolve(driver, { p: p4, T: T4, u: 0 }, driven, { p: p1, T: T1, u: 0 });
  if (!inc) throw new Error("The incident Riemann problem produces a vacuum.");
  const W = inc.right.speed, Ms = inc.right.mach, u2 = inc.uStar, p2 = inc.pStar;
  const s1 = _describe(driven, p1, T1, 0), s4 = _describe(driver, p4, T4, 0);
  const s2 = _describe(driven, p2, inc.right.T, u2), s3 = _describe(driver, p2, inc.left.T, u2);
  const tWall = L1 / W;                                     // incident shock reaches the end wall

  // 2. Reflection at the end wall: state 2 against its mirror image → u5 = 0
  const ref = riemannGasSolve(driven, { p: p2, T: s2.T, u: u2 }, driven, { p: p2, T: s2.T, u: -u2 });
  const Wr = ref.left.speed;                                // lab frame, < 0
  const s5 = _describe(driven, ref.pStar, ref.left.T, 0);
  const Mr = (u2 - Wr) / s2.a;                              // reflected-shock Mach into gas 2

  // 3. Reflected shock meets the contact x_c = u2 t
  const tInt = (L1 - Wr * tWall) / (u2 - Wr), xInt = u2 * tInt;
  const intr = riemannGasSolve(driver, { p: p2, T: s3.T, u: u2 }, driven, { p: s5.p, T: s5.T, u: 0 });
  const ratio = intr.pStar / s5.p;
  const tailoring = Math.abs(ratio - 1) < TAILORED_TOL ? "tailored"
                  : ratio > 1 ? "under-tailored" : "over-tailored";
  // the wave sent toward the end wall runs through gas 5 (at rest)
  const sendSpeed = intr.right.kind === "shock" ? intr.right.speed : intr.right.head;
  const tIntArrive = tInt + (L1 - xInt) / sendSpeed;
  const s5b = _describe(driven, intr.pStar, intr.right.T, 0);   // gas 5 just behind the interaction wave

  // 4. Expansion head: reaches x = -L4 at L4/a4, reflects, and runs back as a
  // C+ characteristic dx/dt = u + a through the undisturbed exact solution.
  const tHeadWall = L4 / s4.a;
  const ua = (x, t) => {                                    // u + a of the incident-wave pattern
    const xi = x / t;
    if (xi <= inc.left.head) return s4.a;                   // region 4 (at rest)
    if (xi < inc.left.tail) {                               // inside the fan
      const tab = inc.left.table;
      let lo = 0, hi = tab.length - 1;
      while (hi - lo > 1) { const m = (lo + hi) >> 1; if (tab[m].u - tab[m].a < xi) lo = m; else hi = m; }
      const A = tab[lo], B = tab[hi], w = (xi - (A.u - A.a)) / ((B.u - B.a) - (A.u - A.a) || 1);
      return A.u + A.a + w * ((B.u + B.a) - (A.u + A.a));
    }
    if (x < u2 * Math.min(t, tInt)) return u2 + s3.a;        // region 3 (driver gas)
    const xr = t > tWall ? L1 + Wr * (t - tWall) : Infinity;
    return x < xr ? u2 + s2.a : s5.a;                       // region 2, then 5
  };
  const path = [{ x: -L4, t: tHeadWall }];
  let x = -L4, t = tHeadWall;
  const dt = (2 * (L4 + L1) / s4.a) / CHAR_STEPS;
  for (let k = 0; k < CHAR_STEPS && x < L1; k++) {
    const k1 = ua(x, t), k2 = ua(x + dt * k1, t + dt);
    x += 0.5 * dt * (k1 + k2); t += dt;
    if (k % 20 === 0) path.push({ x, t });
  }
  const tExpArrive = x >= L1 ? t - (x - L1) / ua(x, t) : Infinity;
  path.push({ x: Math.min(x, L1), t: x >= L1 ? tExpArrive : t });

  // Test time at the end wall
  const intSignificant = tailoring !== "tailored";
  const limit = intSignificant && tIntArrive < tExpArrive ? "reflected wave from contact" : "reflected expansion";
  const tEnd = intSignificant ? Math.min(tIntArrive, tExpArrive) : tExpArrive;
  const testTime = tEnd - tWall;

  // Wave diagram (x [m], t [s]) for plotting
  const lines = [
    { name: "incident shock",   kind: "shock",   pts: [{ x: 0, t: 0 }, { x: L1, t: tWall }] },
    { name: "contact surface",  kind: "contact", pts: [{ x: 0, t: 0 }, { x: xInt, t: tInt }] },
    { name: "expansion head",   kind: "fan",     pts: [{ x: 0, t: 0 }, { x: -L4, t: tHeadWall }] },
    { name: "expansion tail",   kind: "fan",     pts: [{ x: 0, t: 0 }, { x: Math.max(-L4, inc.left.tail * tHeadWall), t: tHeadWall }] },
    { name: "reflected shock",  kind: "shock",   pts: [{ x: L1, t: tWall }, { x: xInt, t: tInt }] },
    { name: "interaction wave", kind: intr.right.kind === "shock" ? "shock" : "fan",
      pts: [{ x: xInt, t: tInt }, { x: L1, t: tIntArrive }] },
    { name: "reflected expansion head (est.)", kind: "fan", pts: path }
  ];

  return {
    states: { 1: s1, 2: s2, 3: s3, 4: s4, 5: s5 },
    Ms, W, u2, p2p1: p2 / p1, p4p1: p4 / p1,
    reflected: { speed: Wr, Mr, p5p1: s5.p / p1, T5T1: s5.T / T1 },
    interaction: { t: tInt, x: xInt, pStarOverP5: ratio, tailoring, wave: intr.right.kind,
                   arriveT: tIntArrive, endWallAfter: s5b },
    expansion: { headAtDriverWall: tHeadWall, arriveT: tExpArrive },
    times: { shockAtWall: tWall, testEnd: tEnd, testTime, limitedBy: limit },
    lines
  };
}

// Driver pressure p4 giving a target incident shock Mach number (bisection on
// log p4; Ms rises monotonically with p4/p1).
function facilityDesignP4(cfgIn, MsTarget) {
  const cfg = Object.assign({}, cfgIn);
  const driver = _gasOf(cfg.driver, cfg.thermo), driven = _gasOf(cfg.driven, cfg.thermo);
  const Ms = p4 => riemannGasSolve(driver, { p: p4, T: cfg.T4, u: 0 }, driven, { p: cfg.p1, T: cfg.T1, u: 0 }).right.mach;
  let lo = cfg.p1 * (1 + 1e-6), hi = cfg.p1 * 10;
  while (Ms(hi) < MsTarget) {
    hi *= 10;
    if (hi > cfg.p1 * 1e12) return null;                    // unreachable with this driver
  }
  for (let it = 0; it < DESIGN_ITERS; it++) {
    const mid = Math.sqrt(lo * hi);
    if (Ms(mid) < MsTarget) lo = mid; else hi = mid;
  }
  return Math.sqrt(lo * hi);
}

// Parameters for the two-gas numerical run of the same facility: one domain
// [0, L4 + L1] with reflective end walls, diaphragm at L4 (shift x by -L4 to
// match the analysis coordinates).
function facilitySimParams(cfgIn, N, opts) {
  const cfg = Object.assign({}, cfgIn), o = opts || {};
  return {
    N: N || 800, L: cfg.L4 + cfg.L1, diaphragm: cfg.L4, bcType: "reflective",
    cfl: o.cfl || 0.4, flux: o.flux || "HLLC",
    gasLeft: _gasOf(cfg.driver, cfg.thermo), gasRight: _gasOf(cfg.driven, cfg.thermo),
    leftState: { p: cfg.p4, T: cfg.T4, u: 0 }, rightState: { p: cfg.p1, T: cfg.T1, u: 0 }
  };
}

// Presets: p4 chosen with facilityDesignP4 for the labelled Ms (T4 = T1 = 300 K),
// thermally perfect gases. T5 is kept below ~2500 K, where frozen,
// non-dissociating composition is still a fair model.
const FACILITY_PRESETS = {
  "He → N₂, tailored (Ms 3.34)": { driver: "He",  driven: "N2",  p4: 376e3, T4: 300, p1: 5e3, T1: 300, L4: 3, L1: 6 },
  "He → Ar, kinetics (Ms 2.8)":  { driver: "He",  driven: "Ar",  p4: 140e3, T4: 300, p1: 5e3, T1: 300, L4: 3, L1: 6 },
  "H₂ → N₂ (Ms 4.0)":            { driver: "H2",  driven: "N2",  p4: 351e3, T4: 300, p1: 5e3, T1: 300, L4: 3, L1: 6 },
  "air → air, weak (Ms 1.4)":    { driver: "air", driven: "air", p4: 50e3,  T4: 300, p1: 1e4, T1: 300, L4: 2, L1: 4 }
};

if (typeof module !== "undefined" && module.exports) {
  module.exports = { facilityAnalysis, facilityDesignP4, facilitySimParams, FACILITY_PRESETS };
}
