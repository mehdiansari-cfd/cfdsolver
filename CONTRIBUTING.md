# Contributing

Thanks for your interest in the Euler Laboratory. Contributions — bug reports,
new test cases, additional schemes, or documentation — are all welcome.

## The one hard rule: stay buildless

This project has **no build step and no runtime dependencies**. Open an HTML
file in a browser and it runs. Please keep it that way:

- **No bundler, transpile step, or `npm install` in the runtime path.** React,
  ReactDOM, and Babel-standalone are loaded from a CDN; the math is plain script
  JavaScript (no `import`/`export`, no transpile) in `gas.js`, `solver.js`,
  `nozzle.js`, `facility.js`, `verify.js`, and `solver2d.js`.
- **No frameworks or libraries beyond the CDN React already in use.** No build
  tooling, no TypeScript compile step, no CSS framework.
- The test runner (`node test.js`) and CI use **only Node's built-in modules**
  (`assert`). No test framework, no dev-dependencies.

If you think a dependency is genuinely warranted, open an issue first and make
the case — the default answer is "keep it vanilla."

## Project layout

| File           | Responsibility |
|----------------|----------------|
| `gas.js`       | Gas models in SI (NASA-7 thermo, mixtures, isentropes). No DOM. |
| `solver.js`    | 1D Euler math: EOS adapters, fluxes, MUSCL, SSP-RK2, exact Riemann (ideal + real gas). No DOM. |
| `nozzle.js`    | Quasi-1D nozzle math (area source, reservoir/back-pressure BCs). Builds on `solver.js`. |
| `solver2d.js`  | 2D Cartesian Euler. Self-contained. No DOM. |
| `index.html`   | 1D UI (shock tube + nozzle). React via Babel-standalone. |
| `index2d.html` | 2D UI. |
| `facility.js`  | Shock-tube facility analysis. No DOM. |
| `case.js`, `run.js` | JSON cases, CSV, headless CLI (Node built-ins only). |
| `verify.js`    | Order-of-accuracy studies and GCI. No DOM. |
| `validation/`  | Experimental comparisons. Data must come from the real source — never invent it. |
| `labui.js`     | Shared UI layer: playback engine, canvas plots, controls. No math. |
| `lab.css`      | Shared stylesheet. |
| `test.js`      | Zero-dependency Node verification suite. |

The math layers know nothing about the DOM; the HTML knows nothing about HLLC.
They communicate through plain JS globals (and CommonJS exports for Node).
Browser load order is `gas.js → solver.js → nozzle.js → facility.js → case.js → verify.js` (each uses the globals of the ones before it; in Node each file requires what it needs). After adding a script, load both pages in a browser: top-level `const` names must be unique across all of them.

## Style

- Match the surrounding code: physical quantities use their symbols (`rho`, `u`,
  `p`, `E`, `gamma`, `cfl`); wave speeds are `SL`, `SR`, `sStar`; loop indices
  are `i`, `j`, `k`. Comments point at an equation reference where useful.
- Every numeric constant gets a name (`gamma = 1.4`, `eps = 1e-8`) — no bare
  magic floats mid-expression.
- Keep the time-stepping loop visible and linear. Helper functions are fine;
  deep class hierarchies are not.
- Don't rewrite working solver math to "clean it up." If you change numerics or
  boundary handling, add or update a check in `test.js`.

## Running the tests

```
node test.js
```

It runs with plain `assert` and exits non-zero on failure. CI runs the same
command on every push and pull request. Please make sure it passes before
opening a PR, and add a case for any new numerical behavior. See the
**Validation** section of the README for what the suite currently checks.

## Submitting changes

1. Fork and branch.
2. Make your change; keep it buildless and dependency-free.
3. Run `node test.js` and confirm it passes.
4. Open a PR with a short description of what changed and why. If it's a
   numerical change, mention which test case demonstrates it.
