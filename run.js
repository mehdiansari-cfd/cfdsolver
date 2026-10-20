#!/usr/bin/env node
/* ============================================================
   run.js  —  headless command-line runner (Node built-ins only)
   ------------------------------------------------------------
   node run.js case.json [options]

     --out DIR            write results here (default: ./results)
     --set key=value      override a case field; dotted paths allowed
                          (--set pb=2e4 --set numerics.N=800); repeatable
     --sweep key=v1,v2,…  run once per value; writes one field CSV per
                          run and a <name>.sweep.csv of scalar outputs
     --area FILE.csv      nozzle: read the area table from a CSV file
                          (two columns x, A; header lines are skipped)
     --quiet              no progress on stderr

   Outputs (per run):
     <name>.csv            final fields (SI units in the header)
     <name>.exact.csv      reference solution, when one exists
     <name>.viscous.csv    displacement-corrected solution (nozzle cases with boundaryLayer)
     <name>.t=<t>.csv      snapshots requested in case.output.times
     <name>.probes.csv     time histories at case.output.probes
     <name>.summary.json   scalar results (mdot, Ms, p5, test time, …)
   ============================================================ */
"use strict";
const fs = require("fs");
const path = require("path");
const C = require("./case");

function usage(msg) {
  if (msg) process.stderr.write("run.js: " + msg + "\n\n");
  process.stderr.write("usage: node run.js case.json [--out DIR] [--set key=value]... [--sweep key=v1,v2,...] [--area FILE.csv] [--quiet]\n");
  process.exit(msg ? 2 : 0);
}

function parseArgs(argv) {
  const a = { sets: [], out: "results", quiet: false };
  for (let i = 0; i < argv.length; i++) {
    const k = argv[i];
    const val = () => { if (i + 1 >= argv.length) usage("missing value after " + k); return argv[++i]; };
    if (k === "--out") a.out = val();
    else if (k === "--set") a.sets.push(val());
    else if (k === "--sweep") a.sweep = val();
    else if (k === "--area") a.area = val();
    else if (k === "--quiet") a.quiet = true;
    else if (k === "-h" || k === "--help") usage();
    else if (k.startsWith("--")) usage("unknown option " + k);
    else if (!a.file) a.file = k;
    else usage("more than one case file given");
  }
  if (!a.file) usage("no case file given");
  return a;
}

function splitKV(s) {
  const i = s.indexOf("=");
  if (i < 1) usage("expected key=value, got '" + s + "'");
  return [s.slice(0, i), s.slice(i + 1)];
}

// Two-column numeric CSV → { x, A }; non-numeric lines (headers, comments) skipped.
function readAreaCSV(file) {
  const x = [], A = [];
  for (const line of fs.readFileSync(file, "utf8").split(/\r?\n/)) {
    const cells = line.split(/[,;\s]+/).filter(Boolean).map(Number);
    if (cells.length >= 2 && cells.every(Number.isFinite)) { x.push(cells[0]); A.push(cells[1]); }
  }
  if (x.length < 2) throw new Error(file + ": no numeric x, A rows found");
  return { x, A };
}

const safe = s => String(s).replace(/[^A-Za-z0-9._=+-]+/g, "_");

function runOne(caseObj, outDir, label, quiet) {
  const t0 = Date.now();
  let last = 0;
  const res = C.runCase(caseObj, (t, tEnd) => {
    if (quiet) return;
    const now = Date.now();
    if (now - last > 500) { last = now; process.stderr.write("\r  " + label + "  t = " + t.toExponential(3) + " / " + tEnd.toExponential(3) + "   "); }
  });
  if (!quiet) process.stderr.write("\r  " + label + "  done in " + ((Date.now() - t0) / 1000).toFixed(1) + " s" + " ".repeat(30) + "\n");

  const base = path.join(outDir, safe(label));
  const header = ["case: " + res.case.name + " (" + res.case.kind + ")", "t = " + (res.t != null ? res.t : "analysis only")];
  if (res.final) fs.writeFileSync(base + ".csv", C.fieldsToCSV(res.final, { header, si: !res.case.gas || res.case.gas.gamma == null }));
  if (res.exact && res.exact.x) {
    const ex = Object.assign({}, res.exact);
    for (const k of Object.keys(ex)) if (!Array.isArray(ex[k])) delete ex[k];
    fs.writeFileSync(base + ".exact.csv", C.fieldsToCSV(ex, { header: header.concat(["reference solution"]) }));
  }
  if (res.viscous) fs.writeFileSync(base + ".viscous.csv", C.tableToCSV(res.viscous.x.map((x, i) => {
    const row = {};
    for (const [k, unit] of [["x", "m"], ["p", "Pa"], ["M", ""], ["T", "K"], ["u", "m/s"], ["rho", "kg/m3"], ["Aeff", "m2"],
                             ["deltaStar", "m"], ["theta", "m"], ["H", ""], ["Cf", ""]]) row[unit ? k + " [" + unit + "]" : k] = res.viscous[k][i];
    return row;
  })));
  for (const snap of res.snapshots) fs.writeFileSync(base + ".t=" + snap.t.toExponential(4) + ".csv", C.fieldsToCSV(snap.fields, { header }));
  if (res.probes.length) fs.writeFileSync(base + ".probes.csv", C.tableToCSV(res.probes));
  fs.writeFileSync(base + ".summary.json", JSON.stringify({ case: res.case, summary: res.summary }, null, 2) + "\n");
  return res;
}

function main() {
  const a = parseArgs(process.argv.slice(2));
  let caseObj;
  try { caseObj = JSON.parse(fs.readFileSync(a.file, "utf8")); }
  catch (e) { usage("cannot read " + a.file + ": " + e.message); }
  for (const kv of a.sets) { const [k, v] = splitKV(kv); C.setPath(caseObj, k, C.parseValue(v)); }
  if (a.area) caseObj.geometry = Object.assign({ type: "table" }, readAreaCSV(a.area));
  fs.mkdirSync(a.out, { recursive: true });

  const name = caseObj.name || path.basename(a.file, ".json");
  try {
    if (!a.sweep) {
      const res = runOne(caseObj, a.out, name, a.quiet);
      process.stdout.write(JSON.stringify(res.summary, null, 2) + "\n");
      return;
    }
    const [key, list] = splitKV(a.sweep);
    const rows = [];
    for (const raw of list.split(",")) {
      const v = C.parseValue(raw);
      const c = C.setPath(JSON.parse(JSON.stringify(caseObj)), key, v);
      const res = runOne(c, a.out, name + "_" + key + "=" + raw, a.quiet);
      const row = { [key]: v };
      for (const [k, x] of Object.entries(res.summary)) if (typeof x !== "object") row[k] = x;
      rows.push(row);
    }
    const file = path.join(a.out, safe(name) + ".sweep.csv");
    fs.writeFileSync(file, C.tableToCSV(rows));
    process.stdout.write(C.tableToCSV(rows));
  } catch (e) {
    process.stderr.write("\nrun.js: " + e.message + "\n");
    process.exit(1);
  }
}

main();
