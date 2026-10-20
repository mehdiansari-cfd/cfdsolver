# Sajben transonic diffuser — experimental validation

Quasi-1D comparison against the Sajben diffuser wall-pressure measurements
(Bogar, Sajben & Kroutil; distributed by the NPARC Alliance Validation Archive,
<https://www.grc.nasa.gov/www/wind/valid/transdif/transdif.html>).

## Data (in `data/`, unmodified NPARC files)

| File | What it is |
|------|------------|
| `upper.pts`, `lower.pts` | upper / lower wall coordinates, inches (flat lower wall) |
| `sajwpres.data` | wall static pressure p/p0 vs x/H*, bottom and top walls, weak-shock case |
| `sajspres.data` | the same, strong-shock case |

Conditions, from the data-file headers (copied into `config.json`): T0 = 300 K,
p0/p_atm = 1.195 (p0 = 121.08 kPa), p_exit/p0 = 0.82 (weak) and 0.72 (strong).
x is measured from the throat; H* = 1.732 in = 44.0 mm; A_exit/A* = 1.50.

## Run

```sh
node validation/sajben.js
```

A(x) ∝ h(x) = y_upper − y_lower (per metre of span), thermally perfect air, the
quasi-1D nozzle marched to steady state (mean |Δρ|/ρ per flow-through < 1e-6).
Writes to `results/sajben/`: `<case>.computed.csv`, `<case>.measured.csv`,
`<case>.svg`, `<case>.case.json` (runs with `run.js` or loads in the Nozzle tab),
`<case>.measured_duct_x.csv` (measurements in that case's x, for the Nozzle tab's
*Load p/p₀ CSV*) and `report.json`. Copies of the case files and plots are
committed in `cases/` and `docs/validation/`.

## Results (N = 400)

| Case | RMS Δ(p/p0) ahead of the shock | Shock x/H*: computed vs measured | RMS Δ(p/p0) in the recovered flow |
|------|------|------|------|
| weak (0.82)   | 0.023 (bottom) · 0.022 (top) | 2.53 vs 1.49 (bottom) · 1.35 (top) | 0.012 · 0.007 |
| strong (0.72) | 0.022 · 0.017 | 3.89 vs 2.35 · 2.12 | 0.033 · 0.035 |

The numerical shock sits within 0.02 H* of the exact quasi-1D solution, so the
remaining differences are the model's, not the discretisation's:

- **Upstream of the shock** the isentropic acceleration through the throat
  matches the measurements to about 2% of p0 on both walls.
- **The shock is placed 1.0–1.2 H* too far aft (weak) and 1.5–1.8 H* (strong)**,
  and the computed jump is sharper and stronger than measured. Inviscid
  quasi-1D flow has no boundary layer: in the experiment the displacement
  thickness, and in the strong case shock-induced separation, reduce the
  effective area ratio and pull the shock upstream.
- **Downstream** the pressure recovers to the measured exit value (to about 1%
  weak, 3–4% strong, where separated flow recovers more slowly).

This is the expected behaviour of an inviscid quasi-1D model and a useful
calibration of what it can and cannot predict.

## With the boundary-layer displacement correction

`config.json` → `boundaryLayer` couples the exact quasi-1D solution with a
turbulent integral boundary layer on the top and bottom walls (`boundary.js`,
started at the duct inlet). The inlet momentum thickness θ0 is not part of the
data, so it is swept (0, 0.2, 0.4 mm) rather than fitted. The teal curves in the
plots show the result.

| Case | Shock x/H*: measured | inviscid | + boundary layer (θ0 = 0 / 0.2 / 0.4 mm) | RMS Δ(p/p0), all points |
|------|------|------|------|------|
| weak (0.82)   | 1.49 · 1.35 | 2.53 | 1.54 / 1.54 / 1.54 | 0.114 · 0.109 → 0.022 · 0.024 |
| strong (0.72) | 2.35 · 2.12 | 3.89 | 3.32 / 3.32 / 3.83, **separated** | 0.085 · 0.087 → 0.071–0.087 |

- **Weak shock.** The displacement thickness (about 7% of the exit height per
  wall) moves the shock 1.0 H* upstream, to within 0.06 H* (bottom) and
  0.19 H* (top) of the measured jump. It also cuts the overall RMS error
  five-fold and does not depend on θ0. The layer stays attached.
- **Strong shock.** The layer separates at the shock. A direct-mode integral
  method cannot follow separated flow: the exit pressure is no longer a
  monotone function of the shock position, and the answer moves with θ0 and
  the grid. `report.json` marks these runs `reliable: false` and the plot
  labels them "sep.". Placing the strong shock needs an inverse /
  interaction method or a 2D viscous solver.

`test.js` checks the upstream and recovery agreement and that the inviscid
shock is aft of the measured one. With the correction it also checks that the
weak shock lies within 0.4 H* of the data, and that the strong case is flagged
as separated.
