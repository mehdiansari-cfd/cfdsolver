# 1D Euler Shock Tube Laboratory

A self-contained, browser-based playground for solving the 1D compressible Euler equations with several of the standard finite-volume flux schemes, side-by-side with the analytical Riemann solution. No build step, no dependencies to install — open `index.html` in a modern browser and the lab is running.

The math lives in `solver.js` as a pure module with no DOM coupling, so you can also `require` it from Node for headless experiments, convergence studies, or test scripts.

---

## What it solves

The 1D Euler equations for an ideal calorically-perfect gas:

```
dU/dt + dF(U)/dx = 0,    U = [ρ, ρu, E],    E = p/(γ-1) + ½ρu²
```

with a Riemann (shock-tube) initial condition — two constant states separated by a diaphragm at `x = x₀`. At `t = 0` the diaphragm bursts and the resulting wave pattern (some combination of shocks, rarefactions, and a contact discontinuity) is evolved on a finite-volume grid.

---

## Try it

```text
git clone https://github.com/<you>/cfdsolver.git
cd cfdsolver
open index.html        # or just double-click it
```

That's it. The page loads React, ReactDOM, and Babel-standalone from a CDN, then `solver.js` for the math, then renders the lab into `#root`. There's no toolchain, no `npm install`, no bundler. Everything is a single page.

---

## Features

**Four approximate Riemann solvers**, switchable from a dropdown without restarting:

- **HLLC** — three-wave estimate with a contact wave; sharp on contact discontinuities.
- **HLL** — two-wave estimate; simple and robust, smears contacts.
- **Roe** — linearised solver on the Roe-averaged state, with a Harten entropy fix on the genuinely nonlinear fields.
- **Rusanov** (local Lax-Friedrichs) — maximally dissipative one-wave estimate; the safety net.

**Second-order in space and time.** MUSCL reconstruction with a minmod slope limiter gives second-order spatial accuracy where the solution is smooth and drops to first order at extrema (TVD). Time integration is SSP-RK2 (Shu-Osher) with an adaptive CFL-limited time step.

**Exact Riemann reference.** A Newton iteration on the star pressure plus full self-similar wave-structure sampling gives the analytical solution at the current time. Toggle "exact overlay" to draw it as a dashed grey curve on every plot — solver error becomes visible at a glance.

**Two boundary conditions.** Transmissive (zero-gradient) or reflective (mirror with sign flip on momentum). The exact overlay auto-disables for reflective BCs, since the analytical solution assumes an infinite domain.

**Four tube views.** The horizontal "inside the tube" strip can render density, schlieren-style `|∂ρ/∂x|`, pressure, or Mach number as a grey-scale field.

**Eight scalar plots.** Top row: ρ, u, p, e. Bottom row (toggleable): Mach M = u/a, temperature T = p/ρ (non-dim, R=1), entropy s = ln(p/ρ^γ), and sound speed a = √(γp/ρ).

**Live animation.** Hit "animate flow" to march the solution forward in 90 frames while a few hundred Lagrangian tracer particles drift with the local velocity field, drawing short streaks across the tube.

**Five classic presets.** Sod, Lax, the 1-2-3 problem, the strong shock, and the Woodward-Colella blast wave — one click each.

---

## Numerical methods at a glance

| Stage              | Scheme                                              |
|--------------------|-----------------------------------------------------|
| Spatial discretisation | Finite-volume, cell-centred, 2 ghost cells per side |
| Reconstruction     | MUSCL with minmod slope limiter (piecewise linear)  |
| Riemann solver     | HLLC / HLL / Roe / Rusanov (selectable)             |
| Time integration   | SSP-RK2 (two-stage, strong-stability-preserving)    |
| Time step          | Δt = CFL · Δx / max(|u| + a), recomputed each step  |
| Boundary conditions| Transmissive or reflective                          |

---

## Canonical test problems

| Preset         | Left (ρ, u, p)        | Right (ρ, u, p)       | t_end | Notes |
|----------------|-----------------------|-----------------------|-------|-------|
| Sod            | (1.000, 0.0, 1.0)     | (0.125, 0.0, 0.1)     | 0.20  | Mild left rarefaction + contact + right shock. |
| Lax            | (0.445, 0.698, 3.528) | (0.500, 0.0, 0.571)   | 0.16  | Asymmetric, stronger contact. |
| 123 problem    | (1.000, −2.0, 0.4)    | (1.000, 2.0, 0.4)     | 0.15  | Two strong rarefactions, near vacuum at the centre. |
| Strong shock   | (1.000, 0.0, 1000.0)  | (1.000, 0.0, 0.01)    | 0.012 | Stress test for the flux limiter. |
| Blast wave     | (1.000, 0.0, 0.01)    | (1.000, 0.0, 100.0)   | 0.035 | Mirrored stress test. |

γ = 1.4 (air), diaphragm at x₀ = 0.5, domain L = 1.0 by default.

---

## File layout

```
cfdsolver/
├── index.html      UI (React via Babel-standalone), canvas + SVG rendering
├── solver.js       Pure math layer: fluxes, Riemann, time integration
└── README.md       This file
```

`solver.js` knows nothing about the DOM. `index.html` knows nothing about HLLC. They communicate through plain JS globals.

---

## Programmatic API

`solver.js` exports everything as both browser globals and a CommonJS module:

```js
const S = require("./solver");

const params = {
  N: 400, L: 1.0, gamma: 1.4, cfl: 0.4,
  diaphragm: 0.5, bcType: "transmissive", flux: "HLLC",
  leftState:  { rho: 1.0,   u: 0.0, p: 1.0 },
  rightState: { rho: 0.125, u: 0.0, p: 0.1 }
};

let state = S.initState(params);
state = S.stepUntil(state, 0.2);

const sol   = S.extract(state);                                  // { x, rho, u, p, e, M, T, s, a }
const exact = S.exactRiemann(params.leftState, params.rightState,
                             params.gamma, params.diaphragm, 0.2, sol.x);
```

The public surface:

- `primToCons(ρ, u, p, γ)` / `fluxFromCons(ρ, m, E, γ)` — variable conversions.
- `hllc`, `hll`, `roe`, `rusanov` — flux functions. All have the same signature `(UL, UR, γ) → F`.
- `FLUXES` — `{ HLLC, HLL, Roe, Rusanov }` name-to-function table.
- `minmod(a, b)` — the slope limiter.
- `applyBcs(ρ, m, E, "transmissive" | "reflective")`.
- `computeRhs(ρ, m, E, N, Δx, γ, fluxFn?)` — MUSCL + flux divergence on interior cells.
- `initState(params)` / `stepUntil(state, tTarget)` — drive the simulation forward.
- `extract(state)` — primitive + derived fields on interior cells.
- `exactRiemann(left, right, γ, x₀, t, xArr)` — analytical solution sampled at `xArr`; returns `null` if vacuum forms.
- `localVelocity(state, x)` — used by the particle tracer.
- `PRESETS` — the five canonical problems above.

---

## UI quick reference

| Control            | What it does                                                                  |
|--------------------|-------------------------------------------------------------------------------|
| Preset buttons     | Load a canonical Riemann problem.                                             |
| Left / Right state | ρ, u, p on each side of the diaphragm.                                        |
| γ, CFL, N, t_end, L, x₀ | Numerics and domain.                                                     |
| Flux               | Pick which approximate Riemann solver feeds the finite-volume update.         |
| BCs                | Transmissive (zero-gradient) or reflective (wall).                            |
| Tube               | What scalar to render in the horizontal strip: ρ, schlieren, p, or M.         |
| Exact overlay      | Dashed grey reference curve on every plot. Skipped for reflective BCs.        |
| Derived fields     | Show or hide the second row of plots (M, T, s, a).                            |
| Run to t_end       | One-shot simulation to the final time.                                        |
| Animate flow       | Frame-by-frame march with Lagrangian tracer streaks.                          |

---

## References

The choice of methods follows Toro's textbook closely:

- E. F. Toro, *Riemann Solvers and Numerical Methods for Fluid Dynamics*, 3rd ed., Springer, 2009.
  Chapters 4 (exact Riemann), 10 (HLL/HLLC), 11 (Roe), 13 (MUSCL).
- A. Harten, "High resolution schemes for hyperbolic conservation laws," *J. Comp. Phys.*, 49, 1983 — the entropy fix on the Roe solver.
- C.-W. Shu and S. Osher, "Efficient implementation of essentially non-oscillatory shock-capturing schemes," *J. Comp. Phys.*, 77, 1988 — SSP-RK time integration.
- G. A. Sod, "A survey of several finite difference methods for systems of nonlinear hyperbolic conservation laws," *J. Comp. Phys.*, 27, 1978 — the Sod test problem.

---

## Limitations and possible extensions

The lab is deliberately small and pedagogical. A non-exhaustive list of things it does *not* do (and could):

- Higher-order reconstruction (PPM, WENO5, MP5).
- Source terms (gravity, friction, geometric source for cylindrical / spherical 1D).
- Real-gas equations of state, multi-species, or reactive flow.
- Multiple spatial dimensions.
- Convergence studies / order-of-accuracy plots against the exact Riemann.
- Probes with time-history recording and a scrub bar.
- Custom-preset save/load and CSV/JSON export of the solution.

PRs and forks welcome.
