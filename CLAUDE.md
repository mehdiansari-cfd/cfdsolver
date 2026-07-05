# CLAUDE.md

Project guidance for a compressible Euler CFD laboratory supporting **1D shock tube, quasi-1D nozzle, and 2D Cartesian** finite-volume simulations. This file orients Claude when reading, modifying, or extending the codebase.

The numerical framework follows the structure laid out in *Write Your First CFD Solver* by Tom-Robin Teschner (CFD.University, 2024), realized here as a browser-based JavaScript port. When in doubt about a method's intent, that book — together with Toro's *Riemann Solvers* and Anderson's *Modern Compressible Flow* — is the canonical reference.

---

## 0. The hard rule: this project is buildless and dependency-free

**Do not add a build step, bundler, package manager install, transpile step, or framework migration.** This is the single most important constraint.

- The pages (`index.html`, `index2d.html`) load React, ReactDOM, and Babel-standalone from a CDN, then the solver `.js` files, then render. Open the HTML in a browser and it runs. There is no `npm install`, no `node_modules` in the runtime path, no CMake, no webpack/vite.
- The solver math (`solver.js`, `nozzle.js`, `solver2d.js`) is plain, ES5-compatible JavaScript with **no imports and no third-party libraries**. Each file exposes browser globals *and* a CommonJS `module.exports` so it can be `require`d from Node for headless runs and tests.
- The test suite (`test.js`) and CI use **only Node's built-in modules** (`assert`). No test framework, no dev-dependencies.
- If a dependency ever seems genuinely warranted, **flag it and make the case first** — the default answer is "keep it vanilla." See `CONTRIBUTING.md`.

Everything below assumes this constraint.

---

## 1. What this project is

A second-order, explicit, finite-volume solver for the compressible Euler equations on Cartesian grids, with an interactive browser UI.

**Equations solved:**
```
∂U/∂t + ∇·F(U) = H
```

- **1D conserved state:** `U = [ρ, ρu, E]ᵀ`, flux `F = [ρu, ρu² + p, u(E+p)]ᵀ`
- **2D conserved state:** `U = [ρ, ρu, ρv, E]ᵀ`, fluxes `F` (x-direction) and `G` (y-direction)
- **Closure:** ideal gas, `p = (γ−1)(E − ½ρ‖u‖²)`, γ = 1.4 unless overridden
- **Source term `H`:** zero for the shock tube and the 2D solver. The **quasi-1D nozzle** (`nozzle.js`) adds a geometric area source `H = -(1/A)(dA/dx)·[ρu, ρu², u(E+p)]`. That is the only source term in the codebase.
- **Inviscid only.** No viscous terms, no turbulence model, no reactions.

The solver is meant to be educationally clear first, optimized second. Favor readability over micro-optimization unless a measurable hot path is identified (the 2D inner loop is the one place where reusing preallocated buffers matters — see §3).

---

## 2. Solver framework

Every operation belongs to one of three phases. **Preserve this separation** when editing or adding features.

### Pre-processing (runs once)
1. Read parameters (from the UI or a params object)
2. Allocate memory (ghosted arrays)
3. Build the mesh / area profile
4. Initialize the solution (`initState` / `initNozzleState` / `initState2D`)

### Solving (time loop — `stepUntil` / `nozzleStepUntil` / `stepUntil2D`)
1. Prepare solution update (compute stable Δt from CFL)
2. Solve equations
   - Reconstruct states at faces (MUSCL + minmod)
   - Compute fluxes; solve the Riemann problem at each face
   - Integrate in time (SSP-RK2), adding the nozzle source term where applicable
3. Update boundary conditions (ghost cells)
4. Per-step post-processing (extract fields for the UI, advance tracer particles)
5. Check stop condition

### Post-processing (runs once)
1. Extract final primitive + derived fields (`extract` / `extractNozzle` / `extract2D`)
2. Render (canvas heatmap / SVG plots). No file output — see §9.

If a feature doesn't fit cleanly into one of these slots, that's a signal to re-think where it belongs.

---

## 3. Code organization

### Actual file layout
```
cfdsolver/
├── index.html      1D UI: shock tube + nozzle (React via Babel-standalone, canvas + SVG)
├── index2d.html    2D UI: heatmap + centerline line-outs
├── solver.js       1D shock-tube math: fluxes, MUSCL, SSP-RK2, exact Riemann
├── nozzle.js       Quasi-1D nozzle math: area source, reservoir + back-pressure BCs,
│                   isentropic A/A* overlay. Reuses solver.js's flux machinery.
├── solver2d.js     2D Cartesian Euler: 2D-aware HLLC/HLL/Roe/Rusanov, MUSCL, SSP-RK2
├── test.js         Zero-dependency Node verification suite
└── docs/           Screenshot/GIF capture guide + suggested repo metadata
```

**Separation of concerns (do not break this):** the `.js` math files know nothing about the DOM; the HTML files know nothing about HLLC. They communicate through plain JS globals (and CommonJS exports for Node). `nozzle.js` reuses `solver.js` globals (`primToCons`, `FLUXES`, `computeRhs`, `extract`), so in the browser `solver.js` must load **before** `nozzle.js`. `solver2d.js` is self-contained.

Keep the time-stepping loop visible and linear, top-to-bottom like prose. Helper functions are fine; deep class hierarchies are not. This is the book's deliberately flat style, kept for the same pedagogical reason.

### Named constants, not raw literals
There is no C++-style `enum` here. Scheme, flux, and BC selections are **string constants** dispatched through lookup tables:

```js
const FLUXES = { HLLC: hllc, HLL: hll, Roe: roe, Rusanov: rusanov };   // solver.js
// BCs: "transmissive" | "reflective" (1D); nozzle: reservoir inlet + back-pressure exit
// nozzle types: "convergent" | "divergent" | "CD"
```

Compare against these named strings, never a raw index. When adding a scheme, register it in the relevant `FLUXES` / `FLUXES2D` table rather than branching on integers.

### Data structures
- **1D** (`solver.js`, `nozzle.js`): separate parallel arrays `rho`, `mom`, `ene`, each length `N + 4` (**two ghost cells per side**; interior is indices `2 … N+1`). This is simpler in JS than an array-of-triples and keeps the inner loop flat.
- **2D** (`solver2d.js`): flat `Float64Array`s of length `(Nx+4)·(Ny+4)`, indexed `j*(Nx+4) + i`. All slope/flux/RHS buffers are **allocated once in `initState2D` and reused** every step — this keeps the JS engine off the heap during the time loop and makes 200² interactive. Do not allocate inside the 2D time loop.

---

## 4. Numerical methods

### Spatial reconstruction (cell-centered → face-reconstructed)
- **Piecewise constant (1st-order Godunov):** `ϕ_face = ϕ_cell`. Robust, cheap, smears discontinuities. It is the implicit fallback at boundary-adjacent cells (the MUSCL stencil is clipped there via the ghost cells).
- **MUSCL (2nd-order, the shipped default):** piecewise-linear reconstruction with a **minmod** slope limiter. `minmod(a, b)` returns 0 on sign change, else the smaller-magnitude argument. Slopes are limited component-wise on the conserved variables.

**Hard rule:** MUSCL without a limiter is unstable for discontinuous flows (Godunov's theorem — linear higher-order schemes are not monotone). Any higher-order reconstruction must be paired with a TVD limiter. Only **minmod** is currently shipped; superbee, van Albada, MC, and van Leer are valid extensions.

### Riemann solver (consolidating left/right face states into a single flux)
Four are shipped and selectable at runtime from the UI dropdown (`FLUXES` / `FLUXES2D`); the **default is HLLC**:
- **HLLC** (HLL + Contact wave) — sharp on contact discontinuities; the default and what most production compressible solvers use. Toro Ch. 10.
- **HLL** — two-wave estimate; robust, smears contacts.
- **Roe** — linearized on the Roe-averaged state, with a Harten entropy fix on the genuinely nonlinear fields. Toro Ch. 11.
- **Rusanov / local Lax-Friedrichs** — one-wave, maximally dissipative; the debugging safety net. When a new feature breaks, regressing to Rusanov clarifies whether the issue is in the Riemann solver or elsewhere.

`solver2d.js` extends these to 2D by carrying the tangential momentum as a passive scalar across the wave pattern (rotational invariance of the Euler equations makes this exact for axis-aligned faces).

### Time integration
- **Default: SSP-RK2** (two-stage, strong-stability-preserving; Shu-Osher). SSP is required because MUSCL is in the spatial reconstruction — plain RK can violate TVD.
- **CFL limit:** ~1.0 in theory for explicit + linear advection; ~0.5 in practice for nonlinear Euler with MUSCL+limiter. Use **CFL ≤ 0.4** as the default for new cases. Divergence sets in around CFL ≈ 0.58 for Sod with MUSCL.

### Adaptive timestep
Recompute every step: `Δt = CFL · Δx / maxᵢ(|uᵢ| + aᵢ)` in 1D, `Δt = CFL · min(Δx, Δy) / max(|u| + |v| + a)` in 2D. The max wavespeed must be recomputed each step — local conditions change as shocks form. (See pitfall §8.6.)

---

## 5. Boundary conditions

Three primitive types, everything else is built from these:
- **Periodic:** wrap-around.
- **Dirichlet:** prescribe value at boundary.
- **Neumann:** prescribe gradient (most commonly zero, giving `ϕ_B = ϕ_neighbor`).

Composite conditions in this codebase:

| BC type             | Where | Normal vel. | Tangential vel. | ρ, p, E |
|---------------------|-------|-------------|-----------------|---------|
| Transmissive (open) | 1D shock tube, all 2D walls | Neumann | Neumann | Neumann |
| Reflective / inviscid wall | 1D shock tube option | Dirichlet (0) | Neumann | Neumann |
| Reservoir inlet     | nozzle | from Riemann invariant | — | isentropic from `(p₀, T₀)` |
| Back-pressure / supersonic outlet | nozzle | Neumann | — | `p_b` (subsonic) or full extrapolation (supersonic) |

Reflective walls: copy the conserved state from the interior cell and **negate the normal momentum component**; tangential momentum is copied as-is. With ghost cells, mirror across the boundary index. For 2D, south/north treat y-momentum as normal, east/west treat x-momentum as normal.

The nozzle inlet closes the subsonic boundary by fixing the incoming Riemann invariant from the reservoir and extrapolating the outgoing one from the first interior cell; the exit switches between imposed back pressure and full extrapolation based on the local exit Mach. See `nozzle.js` (`applyNozzleInletBC` / `applyNozzleExitBC`).

---

## 6. From 1D to 2D — what actually changes

The 2D solver is structurally the same as 1D. Three concrete differences:

### 6.1 State vector grows
1D `U = [ρ, ρu, E]` → 2D `U = [ρ, ρu, ρv, E]`. `E = p/(γ−1) + ½ρ(u² + v²)`.

### 6.2 Two flux vectors
```
F = [ρu, ρu²+p, ρuv, u(E+p)]      G = [ρv, ρuv, ρv²+p, v(E+p)]
∂U/∂t + (F̃_{i+½,j} − F̃_{i−½,j})/Δx + (G̃_{i,j+½} − G̃_{i,j−½})/Δy = 0
```
Apply the 1D Riemann solver **dimension by dimension**: solve the x-Riemann problem with `u` as the normal velocity, then the y-Riemann problem with `v`. The flux formulas don't change; just feed them the correct normal/tangential velocities.

### 6.3 An extra loop everywhere
Every 1D `for i` gains a `for j`. Reconstruction, flux assembly, time update, BC application, output — all gain a y-loop, and allocations gain a dimension. See pitfall §8.7 (index swapping).

---

## 7. Test cases

Every change to numerics or boundary handling must reproduce these before merging. A subset is automated in `test.js` (§ Validation in the README); the rest are available as UI presets.

### 1D verification (shock tube)

| Case | Initial state (left / right) | t_end | What to check |
|------|------------------------------|-------|---------------|
| **Sod** | (1, 0, 1) / (0.125, 0, 0.1) | 0.20 | Rarefaction + contact + shock. Compare against the exact Riemann solver for ρ, u, p, e. |
| **Lax** | (0.445, 0.698, 3.528) / (0.5, 0, 0.571) | 0.16 | Robustness with non-zero initial velocity. |
| **123 problem** | (1, −2, 0.4) / (1, 2, 0.4) | 0.15 | Near-vacuum at center — positivity-failing schemes crash. |
| **Strong shock** | (1, 0, 1000) / (1, 0, 0.01) | 0.012 | Very strong shock; CFL robustness and positivity. |
| **Blast wave** | (1, 0, 0.01) / (1, 0, 100) | 0.035 | Mirrored strong shock. |

The exact Sod solution is implemented (`exactRiemann`, Toro Ch. 4) and used both as the UI's dashed overlay and as the accuracy check in `test.js`.

### Nozzle verification
The six `NOZZLE_PRESETS` (Subsonic CD, Choked CD, Normal Shock, Full Supersonic CD, Convergent, Divergent) cover the regimes; compare against the `isentropicSolution` overlay where no shock is present. The total-pressure plot dropping across a shock is the key qualitative check.

### 2D verification

| Case | Setup | What to check |
|------|-------|---------------|
| **Sod-radial** | circular diaphragm at r = 0.4 | Radial symmetry preserved; anisotropy is a numerical artifact. |
| **2D Riemann (config 3)** | Lax–Liu four-quadrant | Shock-vortex interaction. |
| **Sedov blast** | point energy deposition, ambient p = 10⁻⁴ | Positivity at extreme ratios. |
| **Kelvin–Helmholtz** | shear layer with small perturbation | Sensitive to numerical dissipation; good scheme comparator. |

### Convergence study
Run Sod at N = 100, 200, 400, 800 and verify the L1 density error drops at the expected order (~1 for constant + Euler-time, ~2 for MUSCL + RK2). **Order recovery is the strongest evidence the solver is correct.** `test.js` checks the error is below a bound at N=800 and decreases under refinement.

---

## 8. Common pitfalls

If you suspect a bug, check these first.

1. **MUSCL without a limiter diverges** on any discontinuity. Godunov's theorem, not a bug.
2. **CFL > 0.5 with explicit time + nonlinear flux** is unsafe. Divergence around CFL ≈ 0.58 for Sod. Upgrade to a higher-order SSP scheme if you need bigger steps.
3. **Division-by-zero in ratios/limiters** on flat regions. Keep the `ε` floors (`Math.max(p, 1e-12)` etc. appear throughout the wave-speed and sound-speed calculations).
4. **Negative pressure/density** in star regions of HLLC for very strong shocks. Detection + fallback to HLL/Rusanov is the standard mitigation. The strong-shock and Sedov cases in `test.js` guard this.
5. **Wrong ghost-cell sign for reflective walls.** Tangential momentum is **copied**; only the **normal** component negates. Getting this wrong tunnels flow through walls or spawns spurious vorticity.
6. **Forgetting to recompute Δt each step.** A Δt fixed from the initial condition violates CFL once shocks form.
7. **In 2D, swapping x and y indices** in nested loops (especially flux assembly and the flat `j*(Nx+4)+i` indexing). Output looks plausible but is wrong. Test 2D against a known 1D solution by collapsing one direction first.
8. **Allocating inside the 2D time loop.** The buffers in `initState2D` are reused deliberately; new allocations each step tank performance.
9. **Loading `nozzle.js` before `solver.js`.** It reuses solver globals and will throw `ReferenceError` otherwise. In the browser, script order matters; in Node, `require("./solver")` first.
10. **Assuming `extract`/`extract2D` return ghosted arrays.** They return interior-only fields for display; the raw ghosted state lives on the state object (`rho`/`mom`/`ene`).

---

## 9. Output and post-processing

There is **no file output** — this is an interactive browser tool, not a batch solver.
- **1D (`index.html`):** a canvas "inside the tube" strip (density / schlieren / pressure / Mach) plus SVG line plots (ρ, u, p, e and derived M, T, s, a), with the exact/isentropic reference drawn as a dashed overlay. Live animation advects Lagrangian tracer particles.
- **2D (`index2d.html`):** a canvas heatmap (density / schlieren / pressure / Mach) rendered via `ImageData`, plus SVG line-outs cut through the x- and y-midlines.
- For headless/batch work, `require` the solver from Node and read the arrays from `extract*` directly (see the README's Programmatic API and `test.js`).

If you add export (CSV/JSON/PNG), it belongs in the UI layer, not the math files.

---

## 10. Style and conventions

- **No build, no dependencies** (see §0). Solver files stay ES5-compatible plain JS with browser-global + CommonJS dual export.
- **Variable naming:** physical quantities use their symbols (`rho`, `u`, `v`, `p`, `E`, `gamma`, `cfl`); loop indices `i, j, k`; wave speeds `SL, SR, sStar`. Readers of Toro/Teschner should be able to scan the code.
- **Comments** point at the *equation reference* (`// Toro Eq. 10.37: contact speed`) rather than restating the code.
- **Magic numbers:** every numerical constant gets a name (`gamma`, `cfl`, `eps = 1e-8/1e-12`). No bare floats mid-expression.
- **UI:** React via Babel-standalone, function components with hooks. Keep the visual style (Fraunces + JetBrains Mono, `#efe9d8` background, `#8b3a2f` accent). The pages register a **classic** JSX runtime preset (`react-classic`) so Babel-standalone emits `React.createElement` — do not switch to the automatic runtime, which emits an ES-module `import` a non-module page can't resolve.

---

## 11. Extension roadmap

Already shipped: HLL, HLLC, Roe, Rusanov; SSP-RK2; the quasi-1D nozzle with area source; the 2D Cartesian solver. Remaining ideas, roughly by effort:

- [ ] Additional limiters (superbee, van Albada, MC, van Leer) — add to the reconstruction, keep them TVD.
- [ ] SSP-RK3 time integration; verify temporal order on a smooth problem.
- [ ] WENO5 reconstruction. Needs a lower-order fallback for the two boundary-adjacent cells. Significantly more involved.
- [ ] Reflective / periodic / inflow BCs for the 2D module (currently transmissive only).
- [ ] More nozzle source terms (cylindrical/spherical 1D reuse the same geometric-source machinery).
- [ ] Convergence-study / order-of-accuracy tooling beyond the single `test.js` check.
- [ ] Solution export (CSV/JSON/PNG) from the UI.
- [ ] Adaptive mesh refinement (AMR). Major rewrite; only after the above are solid.
- [ ] Viscous terms (Navier–Stokes). Separate diffusive flux; more restrictive CFL.

When implementing any of these, add the corresponding test/preset before merging, and keep it buildless.

---

## 12. When in doubt

The books are the reference for the framework, equations, and conventions; the constraints in §0 are the reference for *how* to implement. If a request conflicts with either, surface the conflict rather than silently going either way. Examples:

- "Add a 4th-order scheme without a limiter." → flag; Godunov's theorem.
- "Use CFL = 1.5 to speed things up." → flag; explicit schemes have a hard stability limit.
- "Add a bundler / npm build / TypeScript compile / a charting library." → flag; this project is buildless and dependency-free by design (§0).
- "Implement a finite-element solver in this file." → flag; this codebase is FVM by design.

The goal is a clear, correct, extensible lab — not the most features, not the cleverest code. Optimize for the reader.
