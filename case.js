/* ============================================================
   case.js  —  case files, batch runs and CSV (driver layer)
   ------------------------------------------------------------
   No DOM and no numerics: turns a JSON case into solver params,
   runs it to completion, and formats results as CSV. Shared by the
   browser pages (save / load / export) and the command-line runner
   (run.js), so a case saved from the UI runs unchanged in a script.

   Case format (SI units unless the gas is { "gamma": ... }, which
   selects the legacy non-dimensional ideal gas with R = 1):

   { "kind": "shocktube",
     "name": "air 10:1",
     "gas": "air" | { "name": "air", "model": "thermally-perfect" | "calorically-perfect" }
            | { "gamma": 1.4 },
     ("gasLeft" + "gasRight" instead of "gas": two gases, e.g. a He driver into air)
     "left":  { "p": 1e6, "T": 300, "u": 0 },   (non-dim: { "rho", "u", "p" })
     "right": { "p": 1e5, "T": 300, "u": 0 },
     "L": 1, "diaphragm": 0.5, "bc": "transmissive" | "reflective" | "periodic",
     "numerics": { "N": 400, "cfl": 0.4, "flux": "HLLC" },
     "tEnd": 6e-4,
     "output": { "times": [2e-4, 4e-4], "probes": [0.75] } }

   { "kind": "nozzle", "gas": ..., "geometry": { "type": "CD", "AR": 2, "L": 0.5, "Athroat": 1e-3 }
                                  | { "type": "table", "x": [...], "A": [...] },
     "p0": 3e5, "T0": 300, "pb": 2.4e4, "numerics": {...}, "tEnd": 0.006,
     "steady": { "tol": 1e-6, "tMax": 0.5 } }      optional: run until steady instead
       (stops when the mean |Δρ|/ρ over one flow-through time L/a0 drops below tol)

   Nozzle cases may add a displacement correction (boundary.js), in SI:
   "boundaryLayer": { "geometry": "channel", "span": 0.1, "nWalls": 2, "theta0": 0 }
   ("axisymmetric" takes A = πr²; theta0 = inlet momentum thickness, m).
   { "kind": "facility", "driver": "He", "driven": "N2", "thermo": "thermally-perfect",
     "p4": 376e3, "T4": 300, "p1": 5e3, "T1": 300, "L4": 3, "L1": 6,
     "simulate": true, "numerics": { "N": 600 } }
   ============================================================ */

if (typeof initState === "undefined" && typeof require === "function") {
  var _caseSolver = require("./solver");
  var initState = _caseSolver.initState, stepUntil = _caseSolver.stepUntil, extract = _caseSolver.extract;
  var exactRiemannGas = _caseSolver.exactRiemannGas, exactRiemann = _caseSolver.exactRiemann;
  var _caseNozzle = require("./nozzle");
  var initNozzleState = _caseNozzle.initNozzleState, nozzleStepUntil = _caseNozzle.nozzleStepUntil;
  var extractNozzle = _caseNozzle.extractNozzle, quasi1DSolution = _caseNozzle.quasi1DSolution;
  var _caseFac = require("./facility");
  var facilityAnalysis = _caseFac.facilityAnalysis, facilitySimParams = _caseFac.facilitySimParams;
  var gasModel = require("./gas").gasModel;
}
if (typeof viscousQuasi1D === "undefined" && typeof require === "function") {
  var viscousQuasi1D = require("./boundary").viscousQuasi1D;
}

const CASE_KINDS = ["shocktube", "nozzle", "facility"];
const CASE_DEFAULT_NUMERICS = { N: 400, cfl: 0.4, flux: "HLLC" };
const CSV_DIGITS = 9;

// Gas spec → { gas } (SI model) or { gamma } (legacy non-dimensional).
function resolveGas(spec, thermo) {
  if (spec == null) throw new Error("Missing gas.");
  if (typeof spec === "string") return { gas: gasModel(spec, thermo) };
  if (typeof spec.gamma === "number") return { gamma: spec.gamma };
  if (typeof spec.name === "string") return { gas: gasModel(spec.name, spec.model || thermo) };
  throw new Error("Gas must be a name (\"air\"), { name, model } or { gamma }.");
}

function _need(obj, keys, where) {
  for (const k of keys) {
    if (typeof obj[k] !== "number" || !Number.isFinite(obj[k])) throw new Error(where + ": '" + k + "' must be a number.");
  }
}

// Validate a case and fill defaults. Returns a normalised deep copy.
function normalizeCase(input) {
  const c = JSON.parse(JSON.stringify(input));
  if (!CASE_KINDS.includes(c.kind)) throw new Error("Case 'kind' must be one of " + CASE_KINDS.join(", ") + ".");
  c.name = c.name || c.kind;
  c.numerics = Object.assign({}, CASE_DEFAULT_NUMERICS, c.numerics || {});
  _need(c.numerics, ["N", "cfl"], "numerics");
  if (!(c.numerics.cfl > 0 && c.numerics.cfl <= 0.9)) throw new Error("numerics.cfl must be in (0, 0.9].");
  c.output = c.output || {};
  if (c.kind === "shocktube") {
    _need(c, ["L", "diaphragm", "tEnd"], "shocktube");
    c.bc = c.bc || "transmissive";
    if (c.gasLeft || c.gasRight) {
      if (!c.gasLeft || !c.gasRight) throw new Error("shocktube: give both gasLeft and gasRight (or just gas).");
      if (resolveGas(c.gasLeft, c.thermo).gamma != null || resolveGas(c.gasRight, c.thermo).gamma != null)
        throw new Error("shocktube: gasLeft / gasRight must be real gases, not { gamma }.");
    }
    const legacy = !c.gasLeft && resolveGas(c.gas, c.thermo).gamma != null;
    _need(c.left, legacy ? ["rho", "u", "p"] : ["p", "T", "u"], "left");
    _need(c.right, legacy ? ["rho", "u", "p"] : ["p", "T", "u"], "right");
  } else if (c.kind === "nozzle") {
    if (c.steady) { _need(c.steady, ["tol", "tMax"], "steady"); c.tEnd = c.steady.tMax; }
    _need(c, ["p0", "T0", "pb", "tEnd"], "nozzle");
    c.geometry = c.geometry || {};
    if (c.geometry.type === "table") {
      if (!Array.isArray(c.geometry.x) || !Array.isArray(c.geometry.A)) throw new Error("geometry table needs arrays x and A.");
    } else {
      _need(c.geometry, ["AR", "L"], "geometry");
      c.geometry.type = c.geometry.type || "CD";
    }
    if (c.boundaryLayer) {
      if (resolveGas(c.gas, c.thermo).gamma != null) throw new Error("boundaryLayer needs a real gas in SI units, not { gamma }.");
      const bl = c.boundaryLayer = Object.assign({ geometry: "channel", nWalls: 2, theta0: 0 }, c.boundaryLayer);
      if (!["channel", "axisymmetric"].includes(bl.geometry)) throw new Error("boundaryLayer.geometry must be \"channel\" or \"axisymmetric\".");
      if (bl.geometry === "channel") _need(bl, ["span"], "boundaryLayer (channel: A = span × height)");
      _need(bl, ["theta0"], "boundaryLayer");
    }
  } else {
    _need(c, ["p4", "T4", "p1", "T1", "L4", "L1"], "facility");
    if (c.simulate && !(c.tEnd > 0)) c.tEnd = null;          // filled from the analysis
  }
  return c;
}

// Solver parameters for a normalised case.
function caseToParams(c) {
  const num = c.numerics;
  if (c.kind === "shocktube") {
    const g = c.gasLeft ? { gasLeft: resolveGas(c.gasLeft, c.thermo).gas, gasRight: resolveGas(c.gasRight, c.thermo).gas }
                        : resolveGas(c.gas, c.thermo);
    return Object.assign({ N: num.N, L: c.L, cfl: num.cfl, flux: num.flux, diaphragm: c.diaphragm,
                           bcType: c.bc, leftState: c.left, rightState: c.right }, g);
  }
  if (c.kind === "nozzle") {
    const g = resolveGas(c.gas, c.thermo), geo = c.geometry;
    const p = Object.assign({ N: num.N, cfl: num.cfl, flux: num.flux, p0: c.p0, T0: c.T0, pb: c.pb,
                              type: geo.type, AR: geo.AR, L: geo.L, Athroat: geo.Athroat }, g);
    if (geo.type === "table") { p.areaTable = { x: geo.x, A: geo.A }; p.L = geo.x[geo.x.length - 1] - geo.x[0]; }
    return p;
  }
  const fc = { driver: c.driver, driven: c.driven, thermo: c.thermo, p4: c.p4, T4: c.T4, p1: c.p1, T1: c.T1, L4: c.L4, L1: c.L1 };
  return { facility: fc, sim: facilitySimParams(fc, num.N, { cfl: num.cfl, flux: num.flux }) };
}

// Columns of a field result in CSV order, with units.
function fieldColumns(sol, si) {
  const u = (s, n) => (si ? s + " [" + n + "]" : s);
  const cols = [["x", u("x", "m")], ["rho", u("rho", "kg/m3")], ["u", u("u", "m/s")], ["p", u("p", "Pa")],
                ["T", u("T", "K")], ["M", "M"], ["a", u("a", "m/s")], ["e", u("e", "J/kg")], ["s", u("s", "J/(kg K)")]];
  if (sol.gamma) cols.push(["gamma", "gamma"]);
  if (sol.Y) cols.push(["Y", "Y_driver"]);
  if (sol.A) cols.push(["A", u("A", "m2")], ["pt", u("p0", "Pa")], ["mdot", u("mdot", "kg/s")]);
  return cols.filter(([k]) => sol[k]);
}

function fieldsToCSV(sol, opts) {
  const o = opts || {}, cols = fieldColumns(sol, o.si !== false);
  const fmt = v => (Number.isFinite(v) ? Number(v.toPrecision(CSV_DIGITS)) : "");
  const lines = [];
  if (o.header) for (const h of o.header) lines.push("# " + h);
  lines.push(cols.map(c => c[1]).join(","));
  for (let i = 0; i < sol.x.length; i++) lines.push(cols.map(([k]) => fmt(sol[k][i])).join(","));
  return lines.join("\n") + "\n";
}

// Generic table (array of plain objects with the same keys) → CSV.
function tableToCSV(rows) {
  if (!rows.length) return "";
  const keys = Object.keys(rows[0]);
  const cell = v => (typeof v === "number" ? Number(v.toPrecision(CSV_DIGITS)) : JSON.stringify(v == null ? "" : v));
  return keys.join(",") + "\n" + rows.map(r => keys.map(k => cell(r[k])).join(",")).join("\n") + "\n";
}

// Run a case to completion (synchronous; for Node and small browser jobs).
// onProgress(t, tEnd) is called between chunks.
function runCase(input, onProgress) {
  const c = normalizeCase(input);
  const P = caseToParams(c);
  const out = { case: c, snapshots: [], probes: [], summary: {} };
  const times = (c.output.times || []).slice().sort((a, b) => a - b);
  const probes = c.output.probes || [];

  let s, tEnd, ext, step;
  if (c.kind === "shocktube") {
    s = initState(P); tEnd = c.tEnd; ext = extract; step = stepUntil;
  } else if (c.kind === "nozzle") {
    s = initNozzleState(P); tEnd = c.tEnd; ext = extractNozzle; step = nozzleStepUntil;
  } else {
    const a = facilityAnalysis(P.facility);
    out.analysis = a;
    out.summary = {
      Ms: a.Ms, W: a.W, u2: a.u2, p2: a.states[2].p, T2: a.states[2].T, p3: a.states[3].p, T3: a.states[3].T,
      p5: a.states[5].p, T5: a.states[5].T, rho5: a.states[5].rho, Mr: a.reflected.Mr, reflectedSpeed: a.reflected.speed,
      tailoring: a.interaction.tailoring, pStarOverP5: a.interaction.pStarOverP5,
      shockAtWall: a.times.shockAtWall, testTime: a.times.testTime, testLimitedBy: a.times.limitedBy
    };
    if (!c.simulate) return out;
    s = initState(P.sim); ext = extract; step = stepUntil;
    tEnd = c.tEnd || a.times.testEnd + 0.25 * a.times.testTime;
    c.tEnd = tEnd;
  }

  const probeAt = (sol, x) => {
    const i = Math.max(0, Math.min(sol.x.length - 1, Math.floor(x / s.dx)));
    return { x: sol.x[i], p: sol.p[i], T: sol.T[i], u: sol.u[i], rho: sol.rho[i] };
  };
  const PROBE_CHUNK = probes.length ? tEnd / 2000 : Infinity;   // probe sampling interval
  // Steady-state stop (nozzle): compare rho over one flow-through time.
  const steady = c.kind === "nozzle" && c.steady ? c.steady : null;
  const transit = steady ? s.L / Math.sqrt(s.thermo.cp(s.T0) / (s.thermo.cp(s.T0) - s.thermo.R) * s.thermo.R * s.T0) : 0;
  let tCheck = transit, rhoPrev = steady ? s.rho.slice() : null;
  out.converged = steady ? false : null;
  let ti = 0;
  while (s.t < tEnd) {
    if (steady && s.t >= tCheck) {
      let d = 0, m = 0;
      for (let i = 2; i < s.N + 2; i++) { d += Math.abs(s.rho[i] - rhoPrev[i]); m += s.rho[i]; }
      out.residual = d / m;
      if (out.residual < steady.tol) { out.converged = true; break; }
      rhoPrev = s.rho.slice(); tCheck = s.t + transit;
    }
    const next = Math.min(tEnd, ti < times.length ? times[ti] : Infinity, s.t + PROBE_CHUNK, steady ? tCheck : Infinity);
    s = step(s, next);
    if (probes.length) {
      const sol = ext(s);
      for (let k = 0; k < probes.length; k++) out.probes.push(Object.assign({ t: s.t, probe: k }, probeAt(sol, probes[k])));
    }
    while (ti < times.length && s.t >= times[ti] * (1 - 1e-12)) { out.snapshots.push({ t: s.t, fields: ext(s) }); ti++; }
    if (onProgress) onProgress(s.t, tEnd);
  }
  out.final = ext(s);
  out.steps = s.step;
  out.t = s.t;

  // Reference solution and scalar outputs
  if (c.kind === "shocktube") {
    const gL = P.gasLeft || P.gas, gR = P.gasRight || P.gas;
    if (c.bc === "transmissive") {
      out.exact = gL
        ? exactRiemannGas(gL, c.left, gR, c.right, c.diaphragm, s.t, out.final.x)
        : exactRiemann(c.left, c.right, P.gamma, c.diaphragm, s.t, out.final.x);
      if (out.exact) {
        let l1 = 0, rmax = 0;
        for (let i = 0; i < out.final.x.length; i++) { l1 += Math.abs(out.final.rho[i] - out.exact.rho[i]); rmax = Math.max(rmax, out.exact.rho[i]); }
        out.summary.L1rhoRelative = l1 / out.final.x.length / rmax;
        if (out.exact.solution) { out.summary.pStar = out.exact.solution.pStar; out.summary.uStar = out.exact.solution.uStar; }
      }
    }
  } else if (c.kind === "nozzle") {
    const ex = quasi1DSolution(s, "auto"), f = out.final, n = f.x.length;
    out.exact = ex;
    Object.assign(out.summary, {
      regime: ex.regime, mdot: f.mdot[n >> 1], mdotExact: ex.mdot[0], exitMach: f.M[n - 1], exitMachExact: ex.M[n - 1],
      exitPressure: f.p[n - 1], xShockExact: ex.xShock
    });
    if (c.boundaryLayer) {
      // Displacement-corrected steady solution (boundary.js), next to the inviscid one.
      const v = viscousQuasi1D(P, c.boundaryLayer);
      const at = (X, Y, x) => { let k = 1; while (k < X.length - 1 && X[k] < x) k++; return Y[k - 1] + (x - X[k - 1]) / (X[k] - X[k - 1]) * (Y[k] - Y[k - 1]); };
      out.viscous = {
        x: v.solution.x, p: v.solution.p, M: v.solution.M, T: v.solution.T, u: v.solution.u, rho: v.solution.rho,
        Aeff: v.solution.x.map(x => at(v.x, v.Aeff, x)), deltaStar: v.solution.x.map(x => at(v.x, v.deltaStar, x)),
        theta: v.solution.x.map(x => at(v.bl.x, v.bl.theta, x)), H: v.solution.x.map(x => at(v.bl.x, v.bl.H, x)),
        Cf: v.solution.x.map(x => at(v.bl.x, v.bl.Cf, x))
      };
      Object.assign(out.summary, {
        viscousRegime: v.solution.regime, xShockViscous: v.xShock, mdotViscous: v.solution.mdot[0],
        exitMachViscous: v.solution.M[v.solution.M.length - 1], viscousReliable: v.reliable,
        viscousConverged: v.converged, separation: v.bl.separation
      });
    }
  }
  out.summary.steps = s.step; out.summary.t = s.t;
  if (steady) { out.summary.converged = out.converged; out.summary.residual = out.residual; }
  return out;
}

// Set a dotted path in an object, creating intermediate objects ("numerics.N").
function setPath(obj, path, value) {
  const keys = path.split(".");
  let o = obj;
  for (let i = 0; i < keys.length - 1; i++) { if (typeof o[keys[i]] !== "object" || o[keys[i]] == null) o[keys[i]] = {}; o = o[keys[i]]; }
  o[keys[keys.length - 1]] = value;
  return obj;
}
// "300" → 300, "true" → true, "air" → "air", "[1,2]" → [1, 2]
function parseValue(v) {
  try { return JSON.parse(v); } catch (e) { return v; }
}

if (typeof module !== "undefined" && module.exports) {
  module.exports = { CASE_KINDS, resolveGas, normalizeCase, caseToParams, runCase,
                     fieldColumns, fieldsToCSV, tableToCSV, setPath, parseValue };
}
