# CLAUDE.md

Project guidance for a compressible Euler CFD laboratory supporting **1D shock tube, quasi-1D nozzle, shock-tube facility design, and 2D Cartesian** finite-volume simulations, with real gases in SI units, case files, a CLI runner and verification tooling. This file orients Claude when reading, modifying, or extending the codebase.

The numerical framework follows the structure laid out in *Write Your First CFD Solver* by Tom-Robin Teschner (CFD.University, 2024), realized here as a browser-based JavaScript port. When in doubt about a method's intent, that book — together with Toro's *Riemann Solvers* and Anderson's *Modern Compressible Flow* — is the canonical reference.

---

## 0. The hard rule: this project is buildless and dependency-free

**Do not add a build step, bundler, package manager install, transpile step, or framework migration.** This is the single most important constraint.

- The start page (`index.html`) is static HTML that loads only `lab.css`, `gas.js` and `solver.js` (its hero runs a live Sod tube). The lab pages (`lab.html`, `index2d.html`) load React, ReactDOM, and Babel-standalone from a CDN (exact versions pinned), then `labui.js` and the solver `.js` files, then render. Open the HTML in a browser and it runs. There is no `npm install`, no `node_modules` in the runtime path, no CMake, no webpack/vite.
- The math (`gas.js`, `solver.js`, `nozzle.js`, `facility.js`, `verify.js`, `solver2d.js`) and the driver layer (`case.js`) are plain JavaScript with **no third-party libraries** (only each other). Each file exposes browser globals *and* a CommonJS `module.exports` so it can be `require`d from Node for headless runs and tests.
- The test suite (`test.js`), the CLI (`run.js`) and `validation/` use **only Node's built-in modules**. No test framework, no dev-dependencies.
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
- **Closure (1D):** ideal gas `p = ρRT` with an equation-of-state adapter (§3): the legacy non-dimensional γ-law (`p = (γ−1)(E − ½ρu²)`, R = 1), or a `gas.js` model in SI — calorically perfect (cₚ fixed) or thermally perfect (cₚ(T) from NASA-7 polynomials, so γ = γ(T)) — or two gases with a transported mass fraction Y (row `ρY`). T is recovered from e by Newton iteration.
- **Closure (2D):** non-dimensional γ-law only.
- **Source term `H`:** zero for the shock tube and the 2D solver. The **quasi-1D nozzle** (`nozzle.js`) adds a geometric area source `H = -(1/A)(dA/dx)·[ρu, ρu², u(E+p)]`. That is the only source term in the codebase.
- **Inviscid solvers.** No viscous terms in any time-marching solver, no reactions. The one viscous model is `boundary.js`: a turbulent integral boundary layer whose displacement thickness is coupled to the *exact steady* quasi-1D solution (A_eff = A − δ*), never to `nozzleStepUntil`. Gas composition is **frozen** (no dissociation): fine to ~2500 K, increasingly wrong above.

The solver is meant to be educationally clear first, optimized second. Favor readability over micro-optimization unless a measurable hot path is identified (the 2D inner loop is the one place where reusing preallocated buffers matters — see §3).

---

## 2. Solver framework

Every operation belongs to one of three phases. **Preserve this separation** when editing or adding features.

### Pre-processing (runs once)
1. Read parameters (from the UI, a params object, or a JSON case via `case.js`); resolve gas names to `gas.js` models
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
2. Render (canvas fields and plots); optionally write CSV / JSON — see §9.

If a feature doesn't fit cleanly into one of these slots, that's a signal to re-think where it belongs.

---

## 3. Code organization

### Actual file layout
```
cfdsolver/
├── index.html      Start page: static HTML, live Sod hero (gas.js + solver.js), links into the lab; redirects old index.html#tab links
├── lab.html        1D UI: shock tube, nozzle, facility, verify workspaces (React via Babel-standalone, canvas)
├── index2d.html    2D UI: heatmap + line-outs along movable cuts
├── labui.js        Shared UI layer: Player (time-budgeted compute + playback), canvas
│                   plots, controls. Plain JS + React.createElement, no math, no JSX
├── lab.css         Shared stylesheet: light + dark colour tokens on :root (canvases read them via LabUI.THEME)
├── gas.js          Gas models in SI: NASA-7 thermo, mixtures, T(e)/T(h), isentropes, stagnation
├── solver.js       1D Euler math: EOS adapters, fluxes, MUSCL, SSP-RK2, exact Riemann (ideal + real gas)
├── nozzle.js       Quasi-1D nozzle math: area source, PCHIP area tables, reservoir + back-pressure
│                   BCs, exact quasi-1D solution (incl. a forced shock position). Reuses solver.js's flux machinery.
├── boundary.js     Turbulent integral boundary layer (Head/Green, Ludwieg–Tillmann) + displacement
│                   coupling with quasi1DSolution; generator form for time-sliced browser runs
├── facility.js     Shock-tube facility analysis (incident/reflected shock, tailoring, test time)
├── case.js         Driver layer: JSON cases → params → run → CSV. Shared by the UI and run.js
├── run.js          Headless CLI: node run.js case.json [--set] [--sweep] [--area]
├── verify.js       Order-of-accuracy studies, Celik et al. GCI; runs as generators
├── solver2d.js     2D Cartesian Euler: 2D-aware HLLC/HLL/Roe/Rusanov, MUSCL, SSP-RK2
├── cases/          Example case files
├── validation/     Sajben diffuser comparison with the unmodified NPARC data (never edit or invent data)
├── test.js         Zero-dependency Node verification suite
└── docs/           Screenshot/GIF capture guide + suggested repo metadata
```

**Separation of concerns (do not break this):** the `.js` math files know nothing about the DOM; the HTML files know nothing about HLLC. They communicate through plain JS globals (and CommonJS exports for Node). Browser load order: `labui.js`, then **`gas.js → solver.js → nozzle.js → boundary.js → facility.js → case.js → verify.js`**; each file uses the globals of the ones before it. In Node every file `require`s what it needs through a guarded `if (typeof X === "undefined" && typeof require === "function") { var … }` block at its top (§8.9). `solver2d.js` is self-contained.

**Equation-of-state adapters (`solver.js`).** Every 1D flux / RHS function takes an `eos` argument that is either a **number γ** (legacy non-dimensional path — keep it bit-identical; `test.js` and the presets depend on it) or an adapter from `idealEos(γ)`, `gasEos(gas)` or `binaryEos(gasA, gasB)` exposing `prim(ρ, ρu, E, ρY) → {u, p, a, T, gamma, R}` and `cons(ρ, u, p, Y)`. Never compute `p` from `E` with a γ formula outside the adapters. `labui.js` is UI-only: the pages hand its `Player` plain callbacks (`init`, `advance` = `stepUntil*`, `stableDt`, `record`, `check`), so it never touches the equations. Shared UI code lives there as plain JS (not JSX) because Babel-standalone can't fetch external JSX over `file://`.

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
| Periodic            | 1D shock tube (`"periodic"`) | wrap | wrap | wrap |
| Reflective / inviscid wall | 1D shock tube option | Dirichlet (0) | Neumann | Neumann |
| Reservoir inlet     | nozzle | extrapolated from interior | — | isentropic from `(p₀, T₀)` |
| Back-pressure / supersonic outlet | nozzle | Neumann | — | `p_b` (subsonic) or full extrapolation (supersonic) |

Reflective walls: copy the conserved state from the interior cell and **negate the normal momentum component**; tangential momentum is copied as-is. With ghost cells, mirror across the boundary index. For 2D, south/north treat y-momentum as normal, east/west treat x-momentum as normal.

The nozzle inlet closes the subsonic boundary by extrapolating the velocity from the first interior cell and fixing T and p from the stagnation state (`T = T₀ − u²/2c_p`, isentropic p), which holds total pressure at exactly p₀. Do not swap this for the simple-wave invariant `u + 2a/(γ−1) = 2a₀/(γ−1)`: that is an unsteady-expansion relation and bleeds ~γ·M·p₀ of total pressure, so the nozzle never chokes; the exit switches between imposed back pressure and full extrapolation based on the local exit Mach. See `nozzle.js` (`applyNozzleInletBC` / `applyNozzleExitBC`).

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

### Real gases and two gases
`exactRiemannGas` must reproduce `exactRiemann` for a calorically perfect R = 1, cₚ = 3.5 gas (all five presets, 5×10⁻⁵), and the calorically perfect limit must reproduce the closed-form normal-shock relations. Numerics: thermally perfect air and two-gas He → air vs exact, with every flux, conserving mass and energy.

### Nozzle verification
The `NOZZLE_PRESETS` (six non-dimensional regimes plus two SI cases) are compared against `quasi1DSolution` — the exact steady quasi-1D solution for any gas, including the normal-shock position. The total-pressure plot dropping across a shock is the key qualitative check.

### Facility verification
In the calorically perfect limit `facilityAnalysis` must reproduce the shock-tube equation, p₅/p₂ and the Mr–Mₛ relation (Gaydon & Hurle; Anderson Ch. 7) to 10⁻⁶; the two-gas simulation must match the reflection time, p₅ and T₅ to < 1%.

### 2D verification

| Case | Setup | What to check |
|------|-------|---------------|
| **Sod-radial** | circular diaphragm at r = 0.4 | Radial symmetry preserved; anisotropy is a numerical artifact. |
| **2D Riemann (config 3)** | Lax–Liu four-quadrant | Shock-vortex interaction. |
| **Sedov blast** | point energy deposition, ambient p = 10⁻⁴ | Positivity at extreme ratios. |
| **Kelvin–Helmholtz** | shear layer with small perturbation | Sensitive to numerical dissipation; good scheme comparator. |

### Convergence study and GCI (`verify.js`, Verify tab)
**Design order must be measured on a smooth problem**: the periodic density wave (`"smooth-advection"`) reaches L1 order ≈ 1.9 (2 by design; minmod clips extrema, so L∞ ≈ 1.3). Sod converges at L1 ≈ 0.85 for *any* scheme because of its discontinuities — it is a correctness check, not an order test. Solution verification uses the three-grid GCI of Celik et al. (2008); `test.js` checks that the GCI bands contain the exact answer for the nozzle shock position and the real-gas plateau density.

### Experimental validation
`validation/sajben.js` compares the quasi-1D nozzle with the Sajben diffuser wall pressures (unmodified NPARC files in `validation/sajben/data/`). Expected: RMS Δ(p/p₀) ≈ 0.02 ahead of the shock and ≤ 0.035 in the recovered flow; the inviscid shock lies 1.0–1.8 H* **aft** of the measured one (no boundary layer / separation). That is a model limit, not a bug, and must not be "tuned away". With the boundary-layer displacement correction (`config.json` → `boundaryLayer`) the weak shock moves to 1.54 H* (measured 1.35–1.49) for every inlet θ₀ swept. The strong case separates and is reported `reliable: false`. θ₀ is not in the data: sweep it, never fit it. `test.js` checks these. Never edit the data files or commit invented or "representative" experimental numbers.

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
9. **Loading `nozzle.js` before `solver.js` in the browser.** It reuses solver globals and will throw `ReferenceError` otherwise, so script order matters. In Node, `nozzle.js` `require`s `./solver` itself (the `var` block at its top). Only alias solver.js *function declarations* there, and only with `var`: `const` would collide with the browser globals, and so would `var FLUXES`, because solver.js declares `FLUXES` with `const` (nozzle.js reads it through `_fluxTable()` instead). Load both pages in a browser after touching this block.
10. **Assuming `extract`/`extract2D` return ghosted arrays.** They return interior-only fields for display; the raw ghosted state lives on the state object (`rho`/`mom`/`ene`).
11. **Assuming ρe = p/(γ−1) with a gas model.** SI energies carry an enthalpy reference offset (e(298 K) = −RT). Roe's entropy-wave eigenvector needs dE/dρ|ₚ = u²/2 + (h − cₚT), and its species wave de/dY|ρ,ₚ; dropping them made Roe blow up. Any new flux must work through `eos.prim`.
12. **`const` name collisions across browser scripts.** All page scripts share one global lexical scope: a top-level `const` reused in two files (or `var X` where another file has `const X`) kills the page with a SyntaxError. Give module constants distinctive names; alias only *function declarations* with `var`. Load both pages in a browser after adding a file.
13. **Mutating arrays that the UI memoises.** `PlotCard` memoises on array identity; append by building a new array (`concat`), never `push` into one already passed to a plot.
14. **Two-gas interface blip.** Conservative multi-gas schemes produce a few-percent spurious pressure oscillation where the driver/driven interface is smeared (Abgrall 1996). It is a known artefact, documented in the UI; don't "fix" it by tuning, and don't mistake it for a tailoring wave.
15. **Frozen chemistry.** Above ~2500 K real gases dissociate; the frozen gas models over-predict T. The UI warns; presets stay below it.
16. **Viscous–inviscid coupling with a shock.** Iterating δ* ↔ the inviscid flow directly locks into a limit cycle, because the shock position is hypersensitive to the downstream area. `viscousQuasi1D` therefore bisects on the shock position (`quasi1DSolution(state, { shockAt })`) and inserts the shock as an exact jump in the boundary-layer edge data, so δ* moves continuously with x_s. Once the layer separates, p_exit(x_s) is not monotone and the result is grid- and θ₀-dependent. Such runs come back `reliable: false` and must be shown as such, never quoted as predictions.

---

## 9. Output and post-processing

File output is **CSV and JSON only**, and lives outside the math files: the UI writes case files and CSV through `case.js` + `labui.js` download helpers; `run.js` writes field CSVs, reference-solution CSVs, probe histories, sweep tables and summary JSON.
- **1D (`lab.html`):** a canvas "inside the tube" strip (density / schlieren / pressure / Mach) plus canvas line plots (ρ, u, p, e and derived M, T, s, a), with the exact/isentropic reference drawn as a dashed overlay. Playback advects Lagrangian tracer particles. The nozzle adds ṁ(x) and a steady-state convergence (residual) plot.
- **2D (`index2d.html`):** a canvas heatmap (density / schlieren / pressure / Mach) rendered via `ImageData`, plus line-outs along two cuts the user places by clicking.
- **Runs are recorded and replayable.** `LabUI.Player` computes in ~10 ms slices per animation frame (never blocking the page), records a frame each time the solution passes the next multiple of `t_end / maxFrames`, and plays them back at a steady rate with a scrubbable timeline. It never clips a step to land on a frame time — that would add a sliver step per frame.
- For headless/batch work, `require` the solver from Node and read the arrays from `extract*` directly (see the README's Programmatic API and `test.js`).

New export formats belong in `case.js` / the UI layer, not the math files.

---

## 10. Style and conventions

- **No build, no dependencies** (see §0). Math files stay plain script JavaScript (ES2015+ syntax is fine: `const`, arrows, generators — but no `import`/`export`, no transpile) with browser-global + CommonJS dual export.
- **Variable naming:** physical quantities use their symbols (`rho`, `u`, `v`, `p`, `E`, `gamma`, `cfl`); loop indices `i, j, k`; wave speeds `SL, SR, sStar`. Readers of Toro/Teschner should be able to scan the code.
- **Comments** point at the *equation reference* (`// Toro Eq. 10.37: contact speed`) rather than restating the code.
- **Magic numbers:** every numerical constant gets a name (`gamma`, `cfl`, `eps = 1e-8/1e-12`). No bare floats mid-expression.
- **UI:** React via Babel-standalone, function components with hooks. Keep the workspace layout (top bar with pill tabs, control rail, transport bar, stage, plot grid) and the clean style: system sans-serif UI text with tabular figures (mono only for code and axis ticks), sentence-case labels, cards with hairline borders, generous spacing. Two themes, light and dark (system default, toggled in the top bar, remembered in `localStorage["eulerlab-theme"]`): every colour is a token on `:root` in `lab.css`, defined for both themes (rust accent, amber `--accent-2` for live/cursor, `--series` / `--ref` / `--ref2` / `--data` for plot ink). **Canvases never hard-code colours**: they draw with `LabUI.THEME`, which `readTheme()` fills from those tokens, and `useCanvas` redraws on a theme switch. Fields keep the warm cream → rust colormap in both themes. Numerical solution = solid accent line, analytical reference = dashed grey, in every plot. Draw plots and fields on `<canvas>` (not SVG) so per-frame updates stay cheap. The pages register a **classic** JSX runtime preset (`react-classic`) so Babel-standalone emits `React.createElement` — do not switch to the automatic runtime, which emits an ES-module `import` a non-module page can't resolve.

---

## 11. Extension roadmap

Already shipped: HLL, HLLC, Roe, Rusanov; SSP-RK2; the quasi-1D nozzle with area source; the 2D Cartesian solver. Remaining ideas, roughly by effort:

- [ ] Additional limiters (superbee, van Albada, MC, van Leer) — add to the reconstruction, keep them TVD.
- [ ] SSP-RK3 time integration; verify temporal order on a smooth problem.
- [ ] WENO5 reconstruction. Needs a lower-order fallback for the two boundary-adjacent cells. Significantly more involved.
- [ ] Reflective / periodic / inflow BCs for the 2D module (currently transmissive only).
- [ ] More nozzle source terms (cylindrical/spherical 1D reuse the same geometric-source machinery).
- [x] Real gases (calorically / thermally perfect), SI units, two-gas shock tubes.
- [x] Case files, CSV export, headless CLI with sweeps.
- [x] Order-of-accuracy studies and GCI.
- [x] Shock-tube facility analysis.
- [x] Sajben experimental comparison (NPARC data).
- [x] Boundary-layer displacement correction for quasi-1D (integral turbulent layer; weak Sajben shock to within 0.2 H*).
- [ ] Separated shock / boundary-layer interaction (inverse or quasi-simultaneous integral coupling); wall friction in the core.
- [ ] Equilibrium / finite-rate chemistry (frozen composition only today).
- [ ] Fanno (friction) and Rayleigh (heat addition) source terms.
- [ ] Quasi-conservative two-gas formulation (removes the interface blip).
- [ ] Nozzle design outputs (thrust, Isp, C_F), back-pressure ramps, MOC contours.
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
