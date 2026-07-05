# Screenshot / GIF capture guide

The README references a hero image and three animated GIFs that live in this
folder. They are **not** committed yet — this guide says exactly what to record
so the visuals match the captions. Keep each GIF short (5–10 s) and small
(≈ 800–1000 px wide, < 5 MB) so the README stays light.

Any screen recorder works (ScreenToGif on Windows, Kap/QuickTime on macOS,
Peek on Linux). Record the browser content area only — no browser chrome.

Once recorded, drop the files here with these exact names and uncomment the
`<!-- HERO -->` / `<!-- GIFS -->` blocks near the top of `../README.md`.

| File                     | Page          | Preset / setup | What to capture |
|--------------------------|---------------|----------------|-----------------|
| `hero.png`               | `index.html`  | **Sod**, after "run to t_end" | A clean still of the full 1D lab: the tube strip plus the ρ/u/p/e plots with the dashed exact overlay visible. This is the repo's first impression. |
| `sod-shocktube.gif`      | `index.html`  | **Sod**, click **animate flow** | The full animation: rarefaction fan opening left, contact + shock moving right, with the Lagrangian tracer streaks drifting. |
| `nozzle-cd.gif`          | `index.html` → **Nozzle Flow** tab | **Choked CD** (or **Normal Shock**), **animate flow** | The CD nozzle marching to steady state; ideally the **Normal Shock** preset so the shock snaps into the divergent section and the total-pressure plot drops across it. |
| `kelvin-helmholtz.gif`   | `index2d.html` | **Kelvin-Helmholtz**, **animate**, view = **schlieren** or **density** | The shear layer rolling up into the classic KH billows. This is the most visually striking 2D case. |

Tips:
- Use the default resolutions in each preset; they are tuned to stay smooth in
  the browser.
- For the 2D KH GIF, the **schlieren** view shows the vortex roll-up most
  crisply; **density** is also good.
- Trim to the interesting window (skip the first stationary frame).
