#!/usr/bin/env node
/* ============================================================
   validation/sajben.js — quasi-1D comparison with the Sajben diffuser
   (see validation/sajben/README.md). Node built-ins only.
   Skips cleanly (exit 0) when the experimental files are absent.
   ============================================================ */
"use strict";
const fs = require("fs");
const path = require("path");
const Z = require("../nozzle");
const C = require("../case");

const HERE = path.join(__dirname, "sajben");
const OUT = path.join(__dirname, "..", "results", "sajben");

// Rows of a whitespace/comma separated file whose cells are all numbers.
function numericRows(file) {
  const rows = [];
  for (const line of fs.readFileSync(file, "utf8").split(/\r?\n/)) {
    const cells = line.trim().split(/[\s,;]+/).filter(Boolean);
    if (!cells.length) continue;
    const nums = cells.map(Number);
    if (nums.every(Number.isFinite)) rows.push(nums);
  }
  return rows;
}

function interp(xs, ys, x) {
  if (x <= xs[0]) return ys[0];
  if (x >= xs[xs.length - 1]) return ys[ys.length - 1];
  let lo = 0, hi = xs.length - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
  return ys[lo] + (x - xs[lo]) / (xs[hi] - xs[lo]) * (ys[hi] - ys[lo]);
}

function wall(file, cols) {
  const rows = numericRows(file).map(r => [r[cols.x], r[cols.y]]).sort((a, b) => a[0] - b[0]);
  if (rows.length < 2) throw new Error(file + ": fewer than two numeric (x, y) rows — check geometry.columns");
  return { x: rows.map(r => r[0]), y: rows.map(r => r[1]) };
}

function main() {
  const cfg = JSON.parse(fs.readFileSync(path.join(HERE, "config.json"), "utf8"));
  const need = [cfg.geometry.upper, cfg.geometry.lower].concat(Object.values(cfg.cases).map(c => c.file));
  const missing = need.filter(f => !fs.existsSync(path.join(HERE, f)));
  if (missing.length) {
    console.log("SKIPPED: Sajben data not present. Missing:\n  " + missing.map(f => path.join("validation/sajben", f)).join("\n  "));
    console.log("See validation/sajben/README.md for where to get them.");
    return;
  }
  const unset = ["p0", "T0"].filter(k => cfg[k] == null)
    .concat(Object.entries(cfg.cases).filter(([, c]) => c.exitPressureRatio == null).map(([k]) => "cases." + k + ".exitPressureRatio"));
  if (unset.length) throw new Error("Fill these in validation/sajben/config.json from the NPARC case page: " + unset.join(", "));

  // Geometry: h(x) = y_upper - y_lower on the union of the two x grids.
  const up = wall(path.join(HERE, cfg.geometry.upper), cfg.geometry.columns);
  const lo = wall(path.join(HERE, cfg.geometry.lower), cfg.geometry.columns);
  const xs = Array.from(new Set(up.x.concat(lo.x))).sort((a, b) => a - b)
    .filter(x => x >= Math.max(up.x[0], lo.x[0]) && x <= Math.min(up.x[up.x.length - 1], lo.x[lo.x.length - 1]));
  const h = xs.map(x => interp(up.x, up.y, x) - interp(lo.x, lo.y, x));
  const hThroat = Math.min(...h), xThroat = xs[h.indexOf(hThroat)];
  console.log("geometry: " + xs.length + " stations, x in [" + xs[0] + ", " + xs[xs.length - 1] + "], throat h = " + hThroat + " at x = " + xThroat);

  const gas = C.resolveGas(cfg.gas).gas;
  fs.mkdirSync(OUT, { recursive: true });
  const report = {};
  for (const [name, cs] of Object.entries(cfg.cases)) {
    const N = cfg.numerics.N;
    let s = Z.initNozzleState({ type: "table", areaTable: { x: xs, A: h }, gas, p0: cfg.p0, T0: cfg.T0,
                                pb: cs.exitPressureRatio * cfg.p0, N, cfl: cfg.numerics.cfl, flux: cfg.numerics.flux });
    const transit = s.L / Math.sqrt(1.4 * gas.R * cfg.T0);
    let prev = s.rho.slice(), t = transit, res = 1;
    for (let k = 0; k < cfg.numerics.maxTransits && res > cfg.numerics.steadyTol; k++) {
      s = Z.nozzleStepUntil(s, t);
      let d = 0, m = 0;
      for (let i = 2; i < N + 2; i++) { d += Math.abs(s.rho[i] - prev[i]); m += s.rho[i]; }
      res = d / m; prev = s.rho.slice(); t += transit;
    }
    const sol = Z.extractNozzle(s), x0 = xs[0];
    const xPhys = sol.x.map(x => x + x0);
    const pRatio = sol.p.map(p => p / cfg.p0);
    const toPhys = xd => (cs.xScale === "throat-heights" ? xThroat + xd * hThroat : xd);

    const data = numericRows(path.join(HERE, cs.file)).map(r => ({ x: toPhys(r[cs.columns.x]), p: cs.pIsRatio ? r[cs.columns.p] : r[cs.columns.p] / cfg.p0 }))
      .filter(d => d.x >= xPhys[0] && d.x <= xPhys[xPhys.length - 1]);
    let sq = 0;
    for (const d of data) { const pn = interp(xPhys, pRatio, d.x); sq += (pn - d.p) * (pn - d.p); }
    const shockAt = (X, P) => { let k = 1, b = -Infinity; for (let i = 1; i < P.length; i++) if (P[i] - P[i - 1] > b) { b = P[i] - P[i - 1]; k = i; } return 0.5 * (X[k] + X[k - 1]); };
    report[name] = {
      residual: res, dataPoints: data.length, rmsDeltaPoverP0: Math.sqrt(sq / Math.max(1, data.length)),
      shockX_computed: shockAt(xPhys, pRatio),
      shockX_measured: data.length > 2 ? shockAt(data.map(d => d.x), data.map(d => d.p)) : null
    };
    const rows = xPhys.map((x, i) => ({ x, x_over_h_throat: (x - xThroat) / hThroat, p_over_p0_computed: pRatio[i], M_computed: sol.M[i] }));
    fs.writeFileSync(path.join(OUT, name + ".computed.csv"), C.tableToCSV(rows));
    fs.writeFileSync(path.join(OUT, name + ".measured.csv"), C.tableToCSV(data.map(d => ({ x: d.x, x_over_h_throat: (d.x - xThroat) / hThroat, p_over_p0_measured: d.p }))));
    console.log(name + ": " + JSON.stringify(report[name]));
  }
  fs.writeFileSync(path.join(OUT, "report.json"), JSON.stringify(report, null, 2) + "\n");
}

try { main(); } catch (e) { console.error("sajben.js: " + e.message); process.exit(1); }
