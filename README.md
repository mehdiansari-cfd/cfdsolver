# 1D Euler Laboratory

A self-contained, browser-based playground for the 1D compressible Euler equations, covering two canonical problems on the same finite-volume engine:

- **Shock tube** — the classic Riemann problem with shocks, contacts, and rarefactions.
- **Quasi-1D nozzle** — convergent, divergent, or convergent-divergent (CD) duct driven by a stagnation reservoir at the inlet and a static back-pressure at the exit; reproduces subsonic, choked, normal-shock-in-divergent, and fully supersonic flow regimes.

Both modules share four selectable Riemann fluxes, MUSCL reconstruction, SSP-RK2 time integration, and a matching analytical reference overlay (exact Riemann for the tube, isentropic A/A* theory for the nozzle).

No build step, no dependencies to install — open `index.html` in a modern browser and the lab is running.

The math lives in `solver.js` and `nozzle.js` as pure modules with no DOM coupling, so you can also `require` them from Node for headless experiments, convergence studies, or test scripts.

---

## What it solves

The 1D Euler equations for an ideal calorically-perfect gas, in two flavours.

**Shock tube** (uniform area):

```
dU/dt + dF(U)/dx = 0,    U = [ρ, ρu, E],    E = p/(γ-1) + ½ρu²
```

with a Riemann initial condition — two constant states separated by a diaphragm at `x = x₀`. At `t = 0` the diaphragm bursts and the resulting wave pattern (some combination of shocks, rarefactions, and a contact discontinuity) is evolved on a finite-volume grid.

**Quasi-1D nozzle** (variable area `A(x)`):

```
dU/dt + dF(U)/dx = H,    H = -(1/A) (dA/dx) · [ρu, ρu², u(E+p)]
```

driven by a stagnation reservoir `(p₀, T₀)` at the inlet (closed with the outgoing Riemann invariant) and a static back pressure `p_b` at the exit (switching automatically between subsonic and supersonic behaviour as the local exit Mach crosses unity). Time-marched to steady state.

---

## Try it

```text
git clone https://github.com/<you>/cfdsolver.git
cd cfdsolver
open index.html        # or just double-click it
```

That's it. The page loads React, ReactDOM, and Babel-standalone from a CDN, then `solver.js` for the math, then renders the lab into `#root`. There's no toolchain, no `npm install`, no bundler. Everything is a single page.

---

## Shock-tube module

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
| Boundary conditions| Shock tube: transmissive or reflective. Nozzle: reservoir `(p₀, T₀)` inlet + static back-pressure exit (with auto sub/supersonic switch). |

---

## Shock-tube test problems

| Preset         | Left (ρ, u, p)        | Right (ρ, u, p)       | t_end | Notes |
|----------------|-----------------------|-----------------------|-------|-------|
| Sod            | (1.000, 0.0, 1.0)     | (0.125, 0.0, 0.1)     | 0.20  | Mild left rarefaction + contact + right shock. |
| Lax            | (0.445, 0.698, 3.528) | (0.500, 0.0, 0.571)   | 0.16  | Asymmetric, stronger contact. |
| 123 problem    | (1.000, −2.0, 0.4)    | (1.000, 2.0, 0.4)     | 0.15  | Two strong rarefactions, near vacuum at the centre. |
| Strong shock   | (1.000, 0.0, 1000.0)  | (1.000, 0.0, 0.01)    | 0.012 | Stress test for the flux limiter. |
| Blast wave     | (1.000, 0.0, 0.01)    | (1.000, 0.0, 100.0)   | 0.035 | Mirrored stress test. |

γ = 1.4 (air), diaphragm at x₀ = 0.5, domain L = 1.0 by default.

---

## Nozzle module

Quasi-1D Euler in a duct with cross-section `A(x)`. The flux divergence is reused unchanged from `solver.js` (so all four Riemann solvers, MUSCL, and SSP-RK2 carry over); on top of that, `nozzle.js` adds the geometric source term, the reservoir / back-pressure boundary closures, and the isentropic analytical reference.

**Area profiles** — three families, all with `A_throat = 1` by convention and an area ratio `AR = A_max / A_throat`:

- **Convergent** — quadratic taper from `A(0) = AR` down to `A(L) = 1`.
- **Divergent** — quadratic flare from `A(0) = 1` out to `A(L) = AR`.
- **CD** (convergent–divergent) — symmetric parabola with the throat at `x = L/2`.

**Boundary conditions** — physically motivated rather than zero-gradient:

- **Inlet:** isentropic reservoir at `(p₀, T₀)`. The incoming Riemann invariant is fixed from the reservoir; the outgoing invariant is extrapolated from the first interior cell. This naturally handles subsonic inflow and chokes cleanly when the throat reaches `M = 1`.
- **Exit:** static back pressure `p_b` when the exit is subsonic; full extrapolation when the exit is supersonic. The switch is driven by the local exit Mach, so a CD nozzle can transition through choked → shock-in-divergent → fully supersonic without changing the BC code.

**Isentropic analytical overlay.** Inverts the `A/A*` area-Mach relation by bisection on each cell. In a CD nozzle the convergent half is always subsonic; the divergent half follows the user-selected `branch` (subsonic for a Venturi, supersonic for a fully-expanded design point). Useful as a sanity reference, but note that the analytical solution is *not* valid when a normal shock sits in the divergent section.

### Nozzle test problems

| Preset             | Type       | AR  | p_b  | branch     | Notes |
|--------------------|------------|-----|------|------------|-------|
| Subsonic CD        | CD         | 2.0 | 0.97 | subsonic   | Venturi regime, isentropic acceleration then deceleration. |
| Choked CD          | CD         | 2.0 | 0.93 | subsonic   | Throat at M = 1, subsonic recompression downstream. |
| Normal Shock       | CD         | 2.0 | 0.75 | subsonic   | Shock stands inside the divergent section. |
| Full Supersonic CD | CD         | 2.0 | 0.09 | supersonic | Design-point expansion to M ≈ 2.2 at the exit. |
| Convergent         | convergent | 3.0 | 0.60 | subsonic   | Pure acceleration up to (at most) sonic at the exit. |
| Divergent          | divergent  | 2.5 | 0.10 | supersonic | Supersonic expansion from a sonic inlet plenum. |

Stagnation conditions are `p₀ = 1`, `T₀ = 1` (non-dimensional, R = 1), γ = 1.4, `N = 300`, `L = 1` for every preset.

---

## File layout

```
cfdsolver/
├── index.html      UI (React via Babel-standalone), canvas + SVG rendering
├── solver.js       Shock-tube math: fluxes, MUSCL, SSP-RK2, exact Riemann
├── nozzle.js       Quasi-1D nozzle math: area source, reservoir + back-pressure BCs,
│                   isentropic A/A* overlay. Builds on solver.js's flux machinery.
└── README.md       This file
```

`solver.js` and `nozzle.js` know nothing about the DOM. `index.html` knows nothing about HLLC. They communicate through plain JS globals. `nozzle.js` depends on `solver.js` (it reuses `primToCons`, `FLUXES`, `computeRhs`, and `extract`), so `solver.js` must be loaded first.

---

## Programmatic API

Both `solver.js` and `nozzle.js` export everything as browser globals and as CommonJS modules.

### Shock tube

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

Public surface:

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
- `PRESETS` — the five canonical shock-tube problems above.

### Nozzle

```js
require("./solver");                 // must load first — nozzle reuses its globals
const Z = require("./nozzle");

const params = {
  type: "CD", AR: 2.0,               // convergent-divergent, area ratio 2
  p0: 1.0, T0: 1.0, pb: 0.09,        // reservoir + back pressure (design point)
  N: 300, L: 1.0, gamma: 1.4, cfl: 0.35,
  flux: "HLLC"
};

let state = Z.initNozzleState(params);
state = Z.nozzleStepUntil(state, 20);            // march to steady state

const sol   = Z.extractNozzle(state);            // adds pt (total pressure) and A(x)
const ideal = Z.isentropicSolution(state, "supersonic");
```

Public surface:

- `nozzleArea(x, L, AR, type)` / `nozzleAreaDeriv(...)` — the three area profiles `"convergent" | "divergent" | "CD"` and their derivatives.
- `buildAreaArrays(N, dx, L, AR, type)` — precomputed `A` and `dA/dx` on the ghosted grid.
- `computeNozzleRhs(ρ, m, E, N, Δx, γ, fluxFn, A, dA)` — flux divergence plus the geometric source.
- `applyNozzleInletBC(ρ, m, E, γ, p₀, T₀)` — reservoir closure via outgoing Riemann invariant.
- `applyNozzleExitBC(ρ, m, E, γ, p_b, N)` — back-pressure / supersonic-extrapolation switch.
- `initNozzleState(params)` / `nozzleStepUntil(state, tTarget)` — drive the nozzle simulation.
- `extractNozzle(state)` — primitives + derived fields, *plus* `pt` (total pressure) and `A(x)`.
- `areaMachRelation(M, γ)` and `machFromAreaRatio(AR, γ, supersonic)` — the A/A* relation and its bisection inverse.
- `isentropicSolution(state, "subsonic" | "supersonic")` — analytical isentropic profile for overlay.
- `NOZZLE_PRESETS` — the six regimes tabulated above.

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

The choice of methods follows Toro's and Anderson's textbooks closely:

- E. F. Toro, *Riemann Solvers and Numerical Methods for Fluid Dynamics*, 3rd ed., Springer, 2009. Chapters 4 (exact Riemann), 10 (HLL/HLLC), 11 (Roe), 13 (MUSCL).
- J. D. Anderson, *Modern Compressible Flow*, 3rd ed., McGraw-Hill, 2003. Chapters on quasi-1D nozzle flow, choking, and shock-in-nozzle behaviour.
- A. Harten, "High resolution schemes for hyperbolic conservation laws," *J. Comp. Phys.*, 49, 1983 — the entropy fix on the Roe solver.
- C.-W. Shu and S. Osher, "Efficient implementation of essentially non-oscillatory shock-capturing schemes," *J. Comp. Phys.*, 77, 1988 — SSP-RK time integration.
- G. A. Sod, "A survey of several finite difference methods for systems of nonlinear hyperbolic conservation laws," *J. Comp. Phys.*, 27, 1978 — the Sod test problem.

---

## Limitations and possible extensions

The lab is deliberately small and pedagogical. A non-exhaustive list of things it does *not* do (and could):

- Higher-order reconstruction (PPM, WENO5, MP5).
- Other source terms (gravity, friction, cylindrical / spherical 1D — the quasi-1D area source is already in `nozzle.js`).
- Real-gas equations of state, multi-species, or reactive flow.
- Multiple spatial dimensions.
- Convergence studies / order-of-accuracy plots against the exact Riemann or isentropic reference.
- Probes with time-history recording and a scrub bar.
- Custom-preset save/load and CSV/JSON export of the solution.
- Unsteady nozzle scenarios (back-pressure ramps, starting-shock transients).

PRs and forks welcome.
