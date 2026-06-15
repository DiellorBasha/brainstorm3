# Helmholtz / Vorticity View — differential decomposition of a Dirac source map

**Date:** 2026-06-15
**Status:** Design approved; implementation pending
**Branch:** `feat/helmholtz-view` (brainstorm3)
**Author:** Diellor Basha

## Problem

The Wavelet Designer and Spatial Filter manipulate the source field's **magnitude
spectrum** (scale) and **polarization** (helicity). Neither sees the field's
**rotational structure**: the vortices the user observes — points where `|J|→0`
surrounded by a swirling vector field — are **phase singularities**, a *first-order
differential* feature (curl), not a scale or a polarization. There is no tool that takes
the displayed source vector field and shows its **curl, divergence, and the cores of its
phase singularities**, so the user cannot see *how* a vortex is encoded.

The mathematics is the Helmholtz–Hodge decomposition. On the cortex,
`J = ∇φ + ∇⊥ψ + harmonic`:
- `div J = Δφ` → the **scalar potential φ** (sources/sinks = its extrema);
- `curl J·n̂ = Δψ = ω` (vorticity) → the **stream function ψ** (streamlines = its level
  sets, vortex cores = its critical points).

The first-order intrinsic Dirac `D` (already built and stored in the operator node as
`FirstOrder.Intrinsic`) yields `div` and `curl` in one application; the scalar
Laplace–Beltrami operator (already available via `tess_operators 'Laplace-Beltrami'`)
recovers `φ` and `ψ` by a Poisson solve. So the cores are recoverable *exactly*, not by a
heuristic.

## Goal

A **Helmholtz / Vorticity view**: launched from a Dirac source figure, it opens a
dedicated cortex view + control panel that renders, from the in-view field and following
the time cursor:
- a **switchable signed scalar** colormap — Curl·n̂ (vorticity), Divergence, Stream
  function ψ, Scalar potential φ, or |field|;
- the **field vector quiver**;
- **vortex-core markers** (charge `+1` vortex / `−1` antivortex), with a per-frame
  readout of the core count and net charge.

## Non-goals

- Spectral filtering / wavelets — those are the other two tools.
- Streamline *contours* / line-integral-convolution rendering — ψ as a colormap already
  shows the swirl (level sets); explicit streamlines are a possible later layer.
- In-place override of the source figure — a dedicated figure is used instead (a source
  figure couples its colormap and quiver to one 3-component results source, and there is
  no clean public hook to override the displayed scalar per frame).
- Editing the source data — the view is read-only; it derives and displays.

## Decisions (from brainstorming)

| Question | Decision |
|---|---|
| Display model | **One view, switch the scalar** (radio among curl/div/ψ/φ/|field|), quiver + cores overlaid |
| Render target | **Dedicated cortex figure** (not in-place on the source figure) following the global time cursor |
| Scalars exposed | **Full Hodge**: Curl·n̂, Divergence, Stream ψ, Potential φ, |field| |
| Core detection | **Critical points of ψ** (Poisson-recovered): max/min = `+1` vortex, saddle = `−1` antivortex; handedness from `sign(ω)` |
| Launch | **Source figure popup** → "Helmholtz / vorticity (Dirac)" (3-component surface source guard) |

## Architecture

### The decomposition math — `bst_dirac_helmholtz` (pure, testable core)

`H = bst_dirac_helmholtz(DiracOp, LBO, Surface, J)` where `J` is `[3nV × nT]` (the source
field over time), `DiracOp` is the loaded Dirac operator node, `LBO` the loaded
Laplace–Beltrami operator node, `Surface` the loaded tessellation. Per hemisphere `h`:

1. Embed `J_h` as a pure-imaginary quaternion `ψ_h [4nVh × nT]` (`w=0`, `(x,y,z)=J`).
2. `q_h = DiracOp.FirstOrder.Intrinsic{h} · ψ_h` → `[4nFh × nT]` per face.
3. `div_face = q_h(1:4:end,:)` (w-part); `curl_face = (q_h(2:4:end,:), q_h(3:4:end,:),
   q_h(4:4:end,:))`; `ω_face = curl_face · n̂_face` (per-face normal). Sign conventions
   fixed by the test (a synthetic single vortex → one `+1` core with `ω>0`).
4. Area-weighted **face→vertex** averaging → `div_v`, `ω_v [nVh × nT]`.
5. Poisson solve on vertices with the LBO (stiffness `L_h`, mass `M_h`), in the mean-zero
   subspace (`L` null space = constants): `L_h ψ_h = M_h ω_v`, `L_h φ_h = M_h div_v`.
   Solved via a pre-factorized regularized system (`chol(L + ε·M)`), reused across frames.
6. `|field|_v = sqrt(Jx²+Jy²+Jz²)` per vertex.

Returns `H.Curl, H.Div, H.Psi, H.Phi, H.Fmag` — each `[nV × nT]` (global-vertex indexed) —
and `H.Cores`, a `1×nT` cell; `H.Cores{t}` is a struct array `(iVertex, charge, omega)`
of the critical points of `ψ(:,t)`:

- classify each vertex `v` against its 1-ring neighbours (cyclic order from the faces
  around `v`): all-lower = max, all-higher = min (both `charge=+1`, vortex); ≥4 sign
  changes around the ring = saddle (`charge=-1`, antivortex);
- `omega = ω_v(v,t)` gives the swirl handedness at the core.

This function touches no GUI and is fully unit-testable on synthetic fields/operators.

### Dedicated figure + panel — `view_helmholtz` / `panel_helmholtz`

`view_helmholtz(SrcResultsFile)` (launcher, callable from the popup with the source
figure's results file):

1. Resolve the source `iDS/iResult`; require `nComponents == 3`, surface model.
2. Find-or-create the surface's **Dirac operator** node (`tess_operators(SurfaceFile,
   'Dirac')`) and **LBO** node (`tess_operators(SurfaceFile, 'Laplace-Beltrami')`); load
   both; load the surface (faces, normals, VertConn).
3. Materialize the full field `J [3nV × nT]` (`bst_memory('GetResultsValues', …, 0)`),
   run `bst_dirac_helmholtz` → `H`.
4. Open a **dedicated cortex figure** by displaying the *selected scalar* as a
   **1-component results node** (a temp results in the global default study, like the
   wavelet preview) so the colormap + time-scrub are native. Store `H`, the source
   `iDS/iResult`, and the temp node on the figure appdata.
5. Draw the **field quiver** and **core markers** for the current frame (low-level `line`
   glyphs, the seed-marker pattern), and register a **time-slider listener** that
   refreshes those overlays + the readout whenever the cursor moves.
6. Dock `panel_helmholtz` (the controls), attached to the figure.

`panel_helmholtz`:
- **Scalar** radio: Curl·n̂ / Divergence / Stream ψ / Potential φ / |field|. Switching
  swaps the 1-component results' `ImageGridAmp` to `H.<scalar>` and refreshes the colormap
  (diverging for signed scalars).
- **Show vortex cores** / **Show field quiver** checkboxes.
- **Readout** label: "N+ vortices, N− antivortices, net charge C" at the current frame.
- **Close** → remove the time listener, delete the temp node + figure, undock.

### Time-following

The scalar follows time **natively** (1-component results on a `3DViz` figure; on cursor
move Brainstorm calls `panel_surface('UpdateSurfaceData')`). The **overlays** (quiver +
cores + readout) follow time via a listener added to the time slider on open and removed
on close; each fire reads the current frame index and redraws the precomputed `H` overlays
for that frame (cheap: index + redraw, no recompute). **This listener is implemented and
verified first** (Task 1) as the one novel integration; the fallback, if a slider listener
is unreliable, is to drive both the scalar and the overlays from the same listener.

## Components / files

**Create:**
- `toolbox/math/bst_dirac_helmholtz.m` — the pure decomposition + core detector.
- `toolbox/gui/view_helmholtz.m` — launcher + dedicated figure + overlays + time listener.
- `toolbox/gui/panel_helmholtz.m` — the control panel (scalar radio, toggles, readout).
- `dev/tests/test_dirac_helmholtz.m`, `dev/tests/test_helmholtz_view.m`.

**Modify:**
- `toolbox/gui/figure_3d.m` — add "Helmholtz / vorticity (Dirac)" to the source-figure
  popup (3-component surface source guard), next to "Spatial filter (Dirac)".

**Reuse unchanged:** `tess_operators` (Dirac + LBO find-or-create), `bst_memory`
(`GetDataSetResult`, `GetResultsValues`), `view_surface_data` / `panel_surface` (native
scalar display + time), the low-level `line`-glyph marker pattern from
`view_wavelet_designer`, and `db_add`/teardown patterns from the wavelet/spatial-filter
tools for the temp scalar node.

## Data flow

1. Dirac-dSPM → Display on cortex → right-click → "Helmholtz / vorticity (Dirac)".
2. `view_helmholtz`: load operators + surface, materialize `J`, `bst_dirac_helmholtz` → `H`.
3. Open dedicated figure showing `H.Curl` (default) as a 1-comp results; draw quiver +
   cores; add the time listener; dock the panel.
4. Scrub time → scalar updates natively; the listener redraws quiver + cores + readout.
5. Switch scalar (radio) → swap the results' `ImageGridAmp` to the chosen `H.<scalar>`,
   refresh colormap.
6. Close → remove listener, delete temp node + figure.

## Error handling

- **Not a 3-component surface source:** the popup item is hidden; `view_helmholtz` aborts
  with a message if called otherwise.
- **Operators unavailable (nxr):** `tess_operators` find-or-create reports via
  `bst_progress`/`bst_error`; abort cleanly if they cannot be produced.
- **Vertex/size mismatch** (`size(J,1) ≠ 3·nVert`): abort with a clear message.
- **Poisson singular / all-zero frame:** the mean-zero regularized solve handles the
  constant null space; a zero field yields zero scalars and no cores (no error).
- **Teardown always cleans up:** Close / figure `CloseRequestFcn` remove the time listener
  and the temp results node idempotently.

## Testing

- **`test_dirac_helmholtz`** (pure, synthetic): on a small synthetic surface + Dirac/LBO
  operators, a constructed **single vortex** field yields exactly one `+1` core at the
  expected vertex with `ω>0`; a **pure source** (radial outflow) yields high divergence,
  zero curl, no cores; `curl` of a gradient field ≈ 0 and `div` of a rotational field ≈ 0
  (Hodge orthogonality), within discretization tolerance.
- **`test_helmholtz_view`** (live figure, synthetic 3-comp source on `cortex_20484V`):
  `view_helmholtz` opens a figure + panel; switching the scalar swaps the displayed series
  (`ImageGridAmp` equals `H.<scalar>`); the core overlay count matches
  `numel(H.Cores{iTime})`; Close removes the panel, the temp node, and the listener.

## Build order

1. **Time-slider listener spike** — a minimal dedicated figure that follows the time
   cursor and redraws a marker; confirm the listener fires and is cleanly removed. (De-risk.)
2. **`bst_dirac_helmholtz`** — div/curl from the first-order Dirac + LBO Poisson + core
   classification, with the synthetic test (sign conventions fixed here).
3. **`view_helmholtz`** — load + materialize + dedicated scalar figure + overlays + listener.
4. **`panel_helmholtz`** — scalar radio + toggles + readout, wired to the figure.
5. **`figure_3d` popup item** + end-to-end live test.
