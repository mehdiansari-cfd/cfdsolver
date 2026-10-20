# Sajben transonic diffuser — experimental validation harness

Quasi-1D comparison against the Sajben diffuser wall-pressure measurements
(Bogar, Sajben & Kroutil; distributed by the NPARC Alliance Validation Archive,
<https://www.grc.nasa.gov/www/wind/valid/transdif/transdif.html>).

**No experimental data is committed here yet.** The archive could not be
reached from the environment this harness was written in, so the files must be
added by hand. Nothing in this folder is a stand-in for the real data.

## 1. Add the files

Download `sajben.tar` (or the individual files) from the archive and put these
in `validation/sajben/data/`:

| File | What it is |
|------|------------|
| `upper.pts`, `lower.pts` | upper / lower wall coordinates (geometry) |
| `sajwpres.data` | wall static pressures, weak-shock case |
| `sajspres.data` | wall static pressures, strong-shock case |

## 2. Fill in `config.json`

Copy the flow conditions from the archive's case description:
inlet total pressure and temperature, and the exit static-to-inlet-total
pressure ratio of each case. Check the column layout of the `.pts` and `.data`
files and set `geometry.columns` / `cases.*.columns` and the `xScale` (whether
x is in the same units as the geometry or in throat heights) to match.
The loader skips any line that is not all numbers, so headers are fine.

## 3. Run

```sh
node validation/sajben.js
```

It builds A(x) ∝ h(x) = y_upper − y_lower (unit depth), marches the quasi-1D
nozzle to steady state with thermally perfect air, and writes
`results/sajben/<case>.csv` (computed and measured p/p0 vs x) plus the RMS
difference and the shock locations. If the data files are missing it says so
and exits without failing.

## What to expect

Quasi-1D inviscid flow has no boundary layer and no shock-induced separation.
The weak-shock case is the fair comparison (shock position and post-shock
recovery); the strong-shock case separates in the experiment, so a quasi-1D
code will put its shock too far downstream. Report both, with that caveat.
