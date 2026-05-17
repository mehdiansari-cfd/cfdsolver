# CLAUDE.md

Project guidance for a compressible Euler CFD solver supporting **1D and 2D** finite-volume simulations. This file orients Claude when reading, modifying, or extending the codebase.

The numerical framework follows the structure laid out in *Write Your First CFD Solver* by Tom-Robin Teschner (CFD.University, 2024). When in doubt about a method's intent, that book is the canonical reference.

---

## 1. What this project is

A second-order, explicit, finite-volume solver for the compressible Euler equations on Cartesian grids.

**Equations solved** (conservative form, no source terms):
```
∂U/∂t + ∇·F(U) = 0
```

- **1D conserved state:** `U = [ρ, ρu, E]ᵀ`, flux `F = [ρu, ρu² + p, u(E+p)]ᵀ`
- **2D conserved state:** `U = [ρ, ρu, ρv, E]ᵀ`, fluxes `F` (x-direction) and `G` (y-direction)
- **Closure:** ideal gas, `p = (γ−1)(E − ½ρ‖u‖²)`, γ = 1.4 unless overridden
- **Inviscid only.** No viscous terms, no source terms, no turbulence model.

The solver is meant to be educationally clear first, optimized second. Favor readability over micro-optimization unless a measurable hot path is identified.

---

## 2. Solver framework

Every operation in this codebase belongs to one of three phases. **Preserve this separation** when editing or adding features.

### Pre-processing (runs once)
1. Read parameters
2. Allocate memory
3. Create / read mesh
4. Initialize solution

### Solving (time loop)
1. Prepare solution update (store previous state, compute stable Δt from CFL)
2. Solve equations
   - Reconstruct states at faces
   - Compute fluxes; solve Riemann problem at each face
   - Integrate solution in time
3. Update boundary conditions
4. Custom per-step post-processing (CSV write, residuals, etc.)
5. Check stop condition

### Post-processing (runs once)
1. Write final solution
2. Deallocate (only relevant in languages with manual memory)

If a feature doesn't fit cleanly into one of these slots, that's a signal to re-think where it belongs.

---

## 3. Code organization

### File layout
Keep one source file per logical unit. Suggested layout:

```
solver/
├── euler1d.{cpp,py,js}      # 1D driver / main
├── euler2d.{cpp,py,js}      # 2D driver / main
├── solver.{h,py,js}         # core math, shared between 1D and 2D where possible
├── CMakeLists.txt           # build (C++ only)
└── tests/                   # verification cases
```

For **C++**, follow the book's deliberately flat style: a single `main()` function that reads top-to-bottom like prose. Avoid premature abstraction (no classes, no templates beyond what STL provides) **unless the abstraction visibly reduces duplication between 1D and 2D**. The book justifies this choice on pedagogical grounds; we keep it for the same reason.

For **Python/JS** ports, the same principle applies: keep the time-stepping loop visible and linear. Helper functions are fine; deep class hierarchies are not.

### Enums and named constants
Use enums (or string constants in Python) to make face/direction/scheme references self-documenting. The book uses these in C++:

```cpp
enum SCHEME  { CONSTANT = 0, MUSCL };
enum LIMITER { NONE = 0, MINMOD, VANLEER };
enum FACE    { WEST = 0, EAST };
```

For 2D, extend with:

```cpp
enum FACE      { WEST = 0, EAST, SOUTH, NORTH };
enum DIRECTION { X = 0, Y };
```

Never use raw integers (`if (face == 0)`) or raw strings (`if (face == "west")`) for these comparisons. The first is opaque; the second silently fails on typos.

### Data structures
- 1D primary array shape: `U[i][k]` where `i` = cell index, `k ∈ {0,1,2}` = conserved variable.
- 2D primary array shape: `U[i][j][k]` where `k ∈ {0,1,2,3}`.
- Face-reconstructed states: `Ufaces[i][face][k]` (1D) or `Ufaces[i][j][face][k]` (2D).
- Face fluxes: `Ffaces` shaped the same as `Ufaces`.

Prefer `std::array` for the fixed inner dimension (the variable index) and `std::vector` for the cell dimension(s). Stack allocation beats heap for inner-loop access.

---

## 4. Numerical methods

### Spatial reconstruction (cell-centered → face-reconstructed)
Two options live in the codebase, selectable via the `SCHEME` enum:

1. **Piecewise constant (1st-order Godunov):** `ϕ_face = ϕ_cell`. Robust and cheap; smears discontinuities aggressively. Use as the boundary-adjacent fallback for any higher-order scheme, since MUSCL/WENO need wider stencils than are available next to a boundary.

2. **MUSCL (2nd-order):** Linear reconstruction with a slope limiter.
   ```
   ϕᴸ_{i+½} = ϕᵢ + (Ψ(rᵢ₊½)/2)·(ϕᵢ − ϕᵢ₋₁)
   ϕᴿ_{i−½} = ϕᵢ − (Ψ(rᵢ₋½)/2)·(ϕᵢ₊₁ − ϕᵢ)
   ```
   The smoothness indicator `r` includes an `ε = 1e−8` floor in the denominator to prevent NaN on flat regions:
   ```
   rᵢ₊½ = (ϕᵢ₊₁ − ϕᵢ) / (ϕᵢ − ϕᵢ₋₁ + ε)
   ```

**Hard rule:** MUSCL without a limiter is unstable for discontinuous flows (Godunov's theorem — linear higher-order schemes are not monotone). If you implement MUSCL, you must pair it with a TVD limiter. The currently shipped limiters:
- `MINMOD`: `Ψ(r) = max(0, min(1, r))` — robust, slightly diffusive at smooth extrema.
- `VANLEER`: `Ψ(r) = (r + |r|) / (1 + |r|)` — smoother, sharper resolution.

Other limiters (superbee, van Albada, MC) are valid extensions and should be added to the `LIMITER` enum.

### Riemann solver (consolidating Lᴿ and Rᴿ face states into a single flux)
The default is the **Rusanov / local Lax-Friedrichs** flux because it's the simplest robust choice:
```
F̃_{i+½} = ½(Fᴸ + Fᴿ) − S·(Uᴿ − Uᴸ),    S = max(|uᴸ| + aᴸ, |uᴿ| + aᴿ)
```
This is diffusive — particularly at contact discontinuities, where it smears badly. Acceptable upgrades, in order of recommendation:
- **HLL** (Harten–Lax–van Leer): better shock resolution, still smears contacts.
- **HLLC** (HLL + Contact wave): correctly captures contact discontinuities. This is what most production compressible solvers use. See Toro Ch. 10.

If swapping in HLL/HLLC, keep Rusanov available as a fallback for debugging — when a new feature breaks, regressing to Rusanov often clarifies whether the issue is in the Riemann solver or elsewhere.

### Time integration
- **Default:** 1st-order forward Euler. `Uⁿ⁺¹ = Uⁿ − (Δt/Δx)·(F̃ᴱ − F̃ᵂ)` (1D); add y-flux divergence for 2D.
- **CFL limit:** 1.0 in theory for explicit Euler + linear advection; ~0.5 in practice for the nonlinear Euler equations with MUSCL+limiter. Use **CFL ≤ 0.4** as the default for new test cases. The book shows divergence sets in around CFL ≈ 0.58 for the Sod problem with MUSCL+vanLeer.
- **Upgrade path:** SSP-RK2 or SSP-RK3 (Shu-Osher). Strong-stability-preserving is required if MUSCL or WENO is in the spatial reconstruction; plain RK schemes can violate TVD.

### Adaptive timestep
Every step, compute `Δt = CFL · Δx / max_i(|u_i| + a_i)` in 1D, and `Δt = CFL · min(Δx, Δy) / max_{i,j}(|u| + |v| + a)` in 2D. The max wavespeed must be recomputed each step — local conditions change.

---

## 5. Boundary conditions

Three primitive types, everything else is built from these:
- **Periodic:** wrap-around, `ϕ_{B_L−1} = ϕ_{B_R−1}` etc.
- **Dirichlet:** prescribe value at boundary.
- **Neumann:** prescribe gradient (most commonly zero, giving `ϕ_B = ϕ_neighbor`).

Composite conditions used in this solver:

| BC type             | Velocity (normal) | Velocity (tangential) | ρ, p, E       |
|---------------------|-------------------|-----------------------|---------------|
| Transmissive (open) | Neumann (zero)    | Neumann (zero)        | Neumann       |
| Reflective / inviscid wall | Dirichlet (0) | Neumann (zero)    | Neumann       |
| Symmetry            | Dirichlet (0)     | Neumann (zero)        | Neumann       |
| Viscous wall (no-slip) | Dirichlet (0) | Dirichlet (0)         | Neumann       |
| Supersonic inlet    | Dirichlet         | Dirichlet             | Dirichlet (all) |
| Subsonic inlet      | Dirichlet         | Dirichlet             | Dirichlet, one Neumann |
| Supersonic outlet   | Neumann           | Neumann               | Neumann (all) |
| Subsonic outlet     | Neumann           | Neumann               | one Dirichlet (usually p) |

Implement reflective walls by copying the conserved state from the interior cell and **negating the normal momentum component**. Tangential momentum is copied as-is. With ghost cells, mirror across the boundary index.

For 2D, the south and north boundaries get the same treatment as east and west, with the y-momentum being the "normal" component for south/north and x-momentum for east/west.

---

## 6. From 1D to 2D — what actually changes

The 2D solver is structurally the same as 1D. Three concrete differences:

### 6.1 State vector grows
1D `U = [ρ, ρu, E]` → 2D `U = [ρ, ρu, ρv, E]`. The energy expression becomes:
```
E = p/(γ−1) + ½ρ(u² + v²)
```

### 6.2 Two flux vectors
Per cell, compute both `F(U)` (x-direction) and `G(U)`:
```
F = [ρu, ρu²+p, ρuv, u(E+p)]
G = [ρv, ρuv, ρv²+p, v(E+p)]
```
The divergence operator becomes:
```
∂U/∂t + (F̃_{i+½,j} − F̃_{i−½,j})/Δx + (G̃_{i,j+½} − G̃_{i,j−½})/Δy = 0
```

Apply the 1D Riemann solver **dimension by dimension** (dimensional splitting / unsplit FV) — solve the x-Riemann problem using only `u` as the normal velocity, then the y-Riemann problem using `v`. The HLLC/Rusanov formulas don't change; just feed them the correct normal/tangential velocities.

### 6.3 An extra loop everywhere
Anywhere the 1D code loops `for i in 0..N`, the 2D code adds `for j in 0..M`. Reconstruction, flux assembly, time update, BC application, output — all gain a y-loop. Memory allocations gain a dimension.

A practical strategy: write 2D as the general case (with `M=1` collapsing to 1D), or share a `core/` module of dimension-agnostic primitives (Riemann solver, limiter, EOS) between two thin 1D and 2D drivers. The book takes the former approach in its 2D extension hint (Section 5.4).

---

## 7. Test cases

Every change to numerics or boundary handling must reproduce these before being merged.

### 1D verification

| Case             | Initial state (left / right)                       | t_end  | What to check |
|------------------|----------------------------------------------------|--------|---------------|
| **Sod**          | (1, 0, 1) / (0.125, 0, 0.1)                        | 0.20   | Three-wave structure: leftward rarefaction, contact, rightward shock. Compare against the exact Riemann solver (Toro Ch. 4) for ρ, u, p, e. |
| **Lax**          | (0.445, 0.698, 3.528) / (0.5, 0, 0.571)            | 0.16   | Tests robustness with non-zero initial velocity. |
| **Toro Test 1**  | Sod variant; standard regression                   | 0.25   | Sanity. |
| **Toro Test 2 (123)** | (1, −2, 0.4) / (1, 2, 0.4)                    | 0.15   | Near-vacuum forms at center — any positivity-failing scheme will crash here. |
| **Toro Test 3**  | (1, 0, 1000) / (1, 0, 0.01)                        | 0.012  | Very strong shock; tests CFL robustness and limiter behavior at extreme gradients. |
| **Toro Test 4 (Woodward–Colella)** | (1, 0, 0.01) / (1, 0, 100)       | 0.035  | Left-moving blast wave. |

The exact Sod solution is available analytically — implementing a quick exact Riemann solver (Toro Ch. 4) as a verification utility is recommended.

### 2D verification

| Case | Setup | What to check |
|------|-------|---------------|
| **Sod-radial** | 2D extension of Sod with circular diaphragm at `r = 0.4`, axisymmetric initial state | Radial symmetry should be preserved; any anisotropy is a numerical artifact. |
| **Double Mach reflection** (Woodward–Colella) | Mach-10 shock at 60° to a reflecting wall, domain `[0,4]×[0,1]`, `t_end = 0.2` | Classic 2D shock-capturing benchmark; tests strong oblique shock reflection and slip-line resolution. |
| **Forward-facing step** | Mach-3 wind tunnel with a step; `t_end = 4` | Tests handling of geometry-induced expansion fans and shock reflections. |
| **Kelvin–Helmholtz instability** | Shear layer with small perturbation | Sensitive to numerical dissipation; good scheme-comparison case. |

### Convergence study
Run Sod (or any smooth manufactured solution) at N = 100, 200, 400, 800 cells and verify the L1 error in density drops at the expected order: ~1 for piecewise-constant + Euler-time, ~2 for MUSCL + Euler-time, near 2 for MUSCL + RK2. **Order recovery is the strongest evidence the solver is correct.**

---

## 8. Common pitfalls

These have all been hit by people writing this kind of solver. If you suspect a bug, check these first.

1. **MUSCL without a limiter diverges** on any case with a discontinuity. This is Godunov's theorem in action, not a bug in your code.
2. **CFL > 0.5 with explicit Euler time + nonlinear flux** is unsafe. The book shows divergence at CFL ≈ 0.58 for Sod with MUSCL. If you need larger steps, upgrade to SSP-RK2/3.
3. **Smoothness indicator division-by-zero** on flat regions. Always include the `ε = 1e−8` floor.
4. **Negative pressure / density** can occur in star regions of HLLC for very strong shocks (e.g. Toro Test 3 with γ=1.4). Detection + fallback to HLL or first-order is a standard mitigation.
5. **Wrong ghost-cell sign for reflective walls.** Tangential momentum is **copied** (positive); only the **normal** component negates. Getting this wrong produces flow that magically tunnels through walls or develops spurious vorticity.
6. **Forgetting to recompute Δt each step.** Wavespeeds change as the flow evolves; a fixed Δt set from initial conditions will violate CFL once shocks form.
7. **In 2D, swapping x and y indices** in nested loops, especially in flux assembly. Output looks plausible but wrong. Always test 2D against a known 1D solution by collapsing one direction first.
8. **`auto` in C++ inferring `int` for a flux that should be `double`.** Especially with integer literals in expressions. Annotate types where doubts exist.

---

## 9. Output and post-processing

- **1D:** CSV format, one file per timestep: `solution_{N}_{step:06d}.csv` with header `x,rho,u,p`. Plot with the included Jupyter notebook or any CSV-capable tool.
- **2D:** VTK structured grid (`.vts`) or Tecplot ASCII for ParaView/VisIt. CSV becomes unwieldy past a few thousand cells.
- Filenames must include grid size and step count to allow side-by-side comparisons across resolutions.
- The status line during the run prints time, end time, current step, all overwriting on the same console line (`\r`), so a long run produces one updating line, not thousands.

---

## 10. Style and conventions

- **Variable naming:** physical quantities use their physical symbols (`rho`, `u`, `v`, `p`, `E`, `gamma`, `cfl`). Loop indices: `i, j, k`. Wave-speed estimates: `SL, SR, sStar`. Stick to these — readers familiar with Toro or Teschner should be able to scan the code.
- **Comments:** prefer comments that point at the *equation reference* (`// Eq. 2.70: time update`) over comments that restate what the code does.
- **Magic numbers:** every numerical constant gets a name. `gamma = 1.4`, `cfl = 0.4`, `eps = 1e-8`. No bare floats in the middle of an expression.
- **C++:** C++17 minimum. STL only (no Eigen, no third-party libraries) unless a clear case is made. Memory through `std::vector` / `std::array`; no raw `new`/`delete`.
- **Build (C++):** CMake. Always build with optimization (`-O2` / `Release`); the difference between debug and release Euler solvers can be 10×+.

---

## 11. Extension roadmap

Listed roughly in order of effort, lowest first. Pick one at a time.

- [ ] Add Superbee, van Albada, MC limiters to the `LIMITER` enum.
- [ ] Implement HLL Riemann solver. Verify it matches Rusanov to within expected accuracy on Sod.
- [ ] Implement HLLC. Verify the contact discontinuity in Sod is now resolved sharply (this is the canonical demonstration of why HLLC matters).
- [ ] SSP-RK2 time integration. Verify second-order temporal convergence on a smooth problem.
- [ ] WENO5 reconstruction. Significantly more involved; needs lower-order fallback for the two boundary-adjacent cells, not just one. See Scholarpedia's WENO article or Shu's lecture notes.
- [ ] 2D extension following Section 6 above.
- [ ] Adaptive mesh refinement (AMR). Major rewrite; consider only after the above are solid.
- [ ] Add viscous terms (Navier–Stokes). Requires a separate diffusive flux assembly; the CFL condition becomes more restrictive.

When implementing any of these, add the corresponding test case from Section 7 before merging.

---

## 12. When in doubt

The book is the reference for the framework, the equations, and the conventions. If a request seems to conflict with the book's recommendations, surface that conflict rather than silently going either way. Examples:

- "Add a 4th-order scheme without a limiter." → flag this; Godunov's theorem.
- "Use CFL = 1.5 to speed things up." → flag this; explicit schemes have a hard stability limit.
- "Implement a finite-element solver in this file." → flag this; this codebase is FVM by design.

The goal is a clear, correct, extensible solver — not the most features, not the most clever code. Optimize for the reader.