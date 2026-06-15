# Helmholtz components view — show the three Hodge vector fields, not just scalars

**Date:** 2026-06-15
**Status:** Design approved; implementation pending
**Branch:** `feat/helmholtz-view` (brainstorm3)
**Author:** Diellor Basha

## Problem

The Helmholtz view currently colors the cortex by a *scalar* (curl·n̂, divergence, ψ, φ,
|J|) while always drawing the **original** source field `J` as the vector quiver. But a
Helmholtz–Hodge decomposition produces **three vector fields**, and the scalars we show are
either the Poisson *sources* (div, curl·n̂) or the recovered *potentials* (φ, ψ) — never the
decomposition's actual outputs. Two consequences:

1. The selected scalar (e.g. ψ, a *solenoidal* quantity) is shown underneath the *full*
   field's arrows — the colormap and the quiver describe different things.
2. The vortex-core markers (ψ extrema) are centers of the **solenoidal** field `∇⊥ψ`, but
   they're displayed against `J = ∇φ + ∇⊥ψ + h`, so they don't visually sit on the swirls.

We never materialize `∇φ`, `∇⊥ψ`, or `h`.

## Goal

A **component-selected** Helmholtz view. The panel chooses which part of the decomposition
to display, and each state shows **that component's own vector field (quiver) together with
its scalar potential (colormap)**:

| State | Quiver (vector field) | Colormap (scalar) | Markers |
|---|---|---|---|
| **Total** `J` | `J` | `|J|` (one-sided) | none |
| **Irrotational** | `∇φ` | `φ` potential (signed) | sources / sinks = φ extrema (± by div) |
| **Solenoidal** | `∇⊥ψ` | `ψ` stream (signed) | vortex cores = ψ extrema (± by ω) |
| **Harmonic** | `h` | `|h|` (one-sided) | none |

Default state is **Total** (the native source view). The decomposition runs on the
**active frame only** (cached Cholesky factor; recomputed per cursor move), as today.

## Non-goals

- Spectral filtering / wavelets (the other tools).
- A true cohomology/harmonic *basis* — the harmonic part is the **exact residual**
  `h = J − ∇φ − ∇⊥ψ` (see Decisions).
- Editing source data — read-only derivation and display.
- Fixing the core *detection* quality (discrete 1-ring extrema can over/under-count); that
  is tracked separately. This redesign only puts the markers on the correct component.

## Decisions (from brainstorming)

| Question | Decision |
|---|---|
| Panel states | **Total + 3 components**: Total `|J|`, Irrotational `∇φ`, Solenoidal `∇⊥ψ`, Harmonic `h`. Each switches both the quiver and the colormap. |
| Feature markers | **Component-aware**: Solenoidal → vortex cores (ψ extrema, ± by vorticity); Irrotational → sources/sinks (φ extrema, ± by divergence); Total/Harmonic → none. |
| Harmonic definition | **Exact residual** `h = J − ∇φ − ∇⊥ψ` (reconstruction is identically exact); report the harmonic **energy fraction** so we can see how much lands there. |

## Architecture

### Math layer — `bst_dirac_helmholtz` (extend Prepare + Frame)

We already have, per hemisphere `h`: the first-order Dirac `D` (→ per-face divergence
`div` and vorticity `ω = curl·n̂`), the cotan stiffness/mass `K`,`M` with a cached Cholesky
factor, the per-face normals `Nf`, and the area-weighted face→vertex map `Wfv`. We add a
**gradient operator** so we can turn the recovered potentials into vector fields.

**Prepare (new):** build, per hemisphere, the per-face FEM gradient of a per-vertex scalar
as three sparse matrices `Gx,Gy,Gz [nF × nV]`. For a face `(i,j,k)` with area `A` and
normal `n̂`, the gradient coefficient of vertex `i` is `(n̂ × e_i)/(2A)` where `e_i` is the
opposite edge `v_k − v_j` (cyclically). So `∇f|_face = [Gx f, Gy f, Gz f]`. (The same
gradient is obtainable from the first-order Dirac by embedding the scalar in the real/`w`
part and reading the imaginary part on faces; the implementation will use whichever it
validates against the FEM gradient — they must agree up to sign/scale. The FEM form is the
reference.)

**Frame (extended):** for the active frame `Jt [3nV×1]`, per hemisphere:
1. `div`, `ω` from `D` (existing); `φ`,`ψ` from the cached Poisson solves (existing).
2. `∇φ|_face = [Gx φ, Gy φ, Gz φ]` → **irrotational** field; average to vertices with `Wfv`
   → `Virr [nV×3]`.
3. `∇ψ|_face = [Gx ψ, Gy ψ, Gz ψ]`; **skew-gradient** `∇⊥ψ|_face = n̂ × ∇ψ`; average → `Vsol`.
4. `Vtot = reshape(Jt,3,·)'` (per vertex); **harmonic residual** `Vharm = Vtot − Virr − Vsol`.
5. Scalars: `|J|`, `φ`, `ψ`, `|h| = ‖Vharm‖`.
6. Markers (global vertices): vortex cores `FindCores(ψ, VertConn, ω)`; sources/sinks
   `FindCores(φ, VertConn, div)` (same extremum finder, sign taken from div).
7. Harmonic energy fraction `‖Vharm‖²_M / ‖Vtot‖²_M` (mass-weighted, accumulated over hemis).

Return a per-component struct so the view picks one cleanly:

```
Ht.Total = struct('Vec',Vtot[nV×3], 'Scal',|J|[nV×1], 'Signed',false, 'Markers',[],      'Kind','total')
Ht.Irrot = struct('Vec',Virr,       'Scal',φ,          'Signed',true,  'Markers',srcsink, 'Kind','source')
Ht.Solen = struct('Vec',Vsol,       'Scal',ψ,          'Signed',true,  'Markers',cores,   'Kind','vortex')
Ht.Harm  = struct('Vec',Vharm,      'Scal',|h|,        'Signed',false, 'Markers',[],      'Kind','harm')
Ht.HarmFrac = <scalar 0..1>
```

`Markers` is a struct array `(iVertex, charge)` (charge `+1`/`−1`). The whole-series
`Decompose` (the batch primitive) keeps looping `Frame`.

### View layer — `view_helmholtz` (uniform per-component override)

Same shell as now: open the **native** unconstrained-source figure (real vectors + norm),
`Prepare` the operator once, ride the `CustomOverlayFcn` time hook. The per-frame update
becomes **uniform across components** (no Norm special-case):

`UpdateFrame(hFig)`: get the current frame, `Frame`→`Ht` (cached), pick `comp = Ht.(State)`:
- **Colormap/scalar:** `TessInfo(iTess).Data = comp.Scal`; `DataMinMax` symmetric if
  `comp.Signed` else one-sided; figure colormap `stat2` (signed) or `source` (one-sided)
  via `bst_colormaps('AddColormapToFigure', …)`; `panel_surface('UpdateSurfaceColormap')`.
- **Quiver:** `setappdata(hFig,'QuiverVectorOverride', comp.Vec)` then
  `figure_3d('PlotSourceVectors', hFig, iTess)` — the native cones now show *this
  component's* field, gated by its own magnitude.
- **Markers:** if shown, draw `comp.Markers` (red `+` / blue `−`); readout depends on
  `comp.Kind`: vortices/antivortices, sources/sinks, harmonic-energy %, or blank (total).

`SetComponent(hFig, name)` sets the state and calls `UpdateFrame`. `SetVectors`/`SetMarkers`
toggle the two overlays.

### Panel — `panel_helmholtz`

- **Component** radio: `Total field |J|` / `Irrotational (∇φ)` / `Solenoidal (∇⊥ψ)` /
  `Harmonic (h)` — default Total.
- **Show vectors** checkbox (default on) — the component quiver.
- **Show singular points** checkbox (default on) — component-aware markers.
- **Readout** label: per component — e.g. "12 vortices, 9 antivortices (net +3)",
  "8 sources, 5 sinks", "harmonic energy: 4.2% of |J|²", or blank for Total.
- **Close**.

## Components / files

**Modify:**
- `toolbox/math/bst_dirac_helmholtz.m` — Prepare: add `Gx/Gy/Gz`; Frame: compute the three
  component fields + `|h|` + `HarmFrac`, return the per-component struct; reuse `FindCores`
  for both vortex cores and sources/sinks.
- `toolbox/gui/view_helmholtz.m` — component-based uniform override (quiver = component
  field, colormap = component scalar), component-aware markers + readout.
- `toolbox/gui/panel_helmholtz.m` — 4-state component radio + Show vectors + Show singular
  points + per-component readout.

**Reuse unchanged:** the Prepare/Frame on-demand model + cached Cholesky, the
`QuiverVectorOverride`/`PlotSourceVectors` native-vector path, the `CustomOverlayFcn` time
hook, `bst_colormaps` (`source`/`stat2`), `LoadResultsFileFull` launch.

**Tests:**
- `dev/tests/test_dirac_helmholtz.m` — add **decomposition correctness**: exact
  reconstruction `Virr+Vsol+Vharm == Vtot`; `curl(∇φ)` small vs `curl(J)` and `div(∇⊥ψ)`
  small vs `div(J)`; Hodge orthogonality `⟨Virr,Vsol⟩_M ≈ 0`; sources/sinks vs vortex cores
  come from φ vs ψ respectively.
- `dev/tests/test_helmholtz_view.m` — component states: each sets `QuiverVectorOverride` to
  the component field, `TessInfo.Data` to the component scalar with the right colormap
  (`source`/`stat2`); markers are component-aware (vortices on Solenoidal, sources/sinks on
  Irrotational, none on Total/Harmonic); cursor move recomputes; close cleans up.

## Data flow

1. Tree node → "Helmholtz / vorticity (Dirac)" → `view_helmholtz` loads source, `Prepare`,
   opens native figure (Total state: `J` + `|J|`).
2. Pick **Solenoidal** → quiver swaps to `∇⊥ψ`, cortex colored by `ψ` (diverging), vortex
   cores marked; readout shows vortex/antivortex counts.
3. Pick **Irrotational** → quiver `∇φ`, color `φ`, sources/sinks marked.
4. Pick **Harmonic** → quiver `h`, color `|h|`, readout shows harmonic energy %.
5. Scrub time → the active frame is re-decomposed (cached) and the chosen component redraws.

## Error handling

- Non-3-component / unloadable source, bare shared kernel: as today (clear messages).
- Degenerate frame (`J≈0`): components are ~0, `HarmFrac` guarded against /0, no markers.
- Stale figure handle: dispatch/panel guards no-op (as today).
- Reconstruction is exact by construction; the test asserts it to machine precision so a
  broken gradient/averaging is caught immediately.

## Open caveat (documented, not blocking)

How energy splits between the two potentials and the harmonic residual depends on the
Poisson **boundary conditions** (each hemisphere is a disk with a medial-wall boundary) and
on discretization. With the exact-residual choice, reconstruction is always exact, but `h`
absorbs both genuine topological/global circulation **and** boundary/discretization leakage.
The harmonic-energy readout is the diagnostic: if it's large, the boundary handling needs
revisiting (a later refinement, e.g. an explicit harmonic basis).

## Build order

1. `bst_dirac_helmholtz` Prepare gradient + Frame components + `FindCores` reuse, with the
   decomposition-correctness test (reconstruction, curl/div-free, orthogonality).
2. `view_helmholtz` uniform per-component override + component-aware markers/readout.
3. `panel_helmholtz` 4-state radio + two toggles + readout; live component test.
