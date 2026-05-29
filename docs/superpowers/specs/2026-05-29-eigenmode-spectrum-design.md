# Eigenspectrum of Source Activations: Modal Power-Spectrum Viewer — Design

**Date:** 2026-05-29
**Status:** Approved design pending spec review (revised after scope clarification)
**Follows:** `2026-05-28-eigenmode-perhemisphere-viewer-r2-design.md` (eigenmode compute + viewer, shipped on `feature/eigenmode-context-menus`)

## Background

Two eigenmode-based capabilities are being built in parallel and must stay **strictly
isolated**:

1. **Eigenspectrum analysis (this increment)** — *general and source-agnostic*. Given
   **any vertex-mapped scalar field** `S [nVert × nTime]`, it views how that field's
   energy distributes across the cortex's spatial LBO eigenmodes (a spatial spectrum,
   the analogue of the temporal FFT/PSD). It owns **no kernel** and does **not** modify
   the results structure. It simply transforms a realized vertex field into eigenmode
   coefficients and plots the modal power spectrum.

2. **Manifold-harmonics source mapping (separate, already implemented)** — a *new kind
   of inverse* producing results with an `[nEigenmodes × nChannels]` kernel. That work
   is tested and viewed later, on its own. **It is out of scope here.**

This round delivers (1), scoped to its first and most useful input: **source
activations** (dSPM / sLORETA / MNE / etc.). Brainstorm's existing "Cortical
activations → Display on cortex" path already realizes the vertex-mapped scalar field
and provides the global time cursor we step through; this feature adds the eigenmode
transform of that same field and a dedicated spectrum view. Other vertex-mapped scalars
(PET projected to the surface, anatomical overlays, …) are deliberately deferred — the
projection helper is written source-agnostic so they slot in later with no rework.

### Reused building blocks (confirmed against the codebase)

- **Realized vertex field.** `in_bst_results(ResultsFile, 1)` returns `ImageGridAmp`
  by applying the shared kernel on the fly (`ImagingKernel · F(GoodChannel,:)`), exactly
  as cortical display does. For unconstrained models (`nComponents == 3`),
  `bst_source_orient([], 3, [], ImageGridAmp, 'rms')` collapses orientations to a scalar
  `[nVert × nTime]` — the same magnitude shown on the cortex.
- **The eigenmode projection.** `bst_eigenmodes_project(Eigenmodes, Data, MassMatrix)`
  returns `Coeffs [K × nTime] = Φ'·(M·Data)` — the M-weighted projection. Already
  unit-tested (`test_eigenmodes_project_pure.m`).
- **The mass matrix.** `tess_laplacian(Vertices, Faces, 'MassType', MassType)` returns
  `[L, M]`; we rebuild `M` consistent with `Eigenmodes.MassType` at view time.
- **Eigenmodes.** `in_tess_eigenmodes(SurfaceFile)` → `Φ`, `Values (λ)`, `Component`
  (1 = Left, 2 = Right), `CompRank`, `MassType`, `isComputed`.

## Goals

1. A **pure transform** of a realized vertex scalar field into eigenmode coefficients,
   handling both constrained (signed scalar) and unconstrained (norm-collapsed) sources.
2. A dedicated **modal power-spectrum viewer**: power per mode vs **eigenvalue** (toggle
   to spatial wavelength), **left/right hemispheres as separate curves**, **time
   stepping** driven by Brainstorm's global time cursor, **and** a **window-averaged
   overlay**.
3. Maximal reuse of existing machinery; **no** results-structure change, **no** kernel,
   no breakage of any existing workflow.

## Design

### A. Acquire the realized vertex scalar field (reuse)

The viewer obtains the vertex field the same way the cortex display does:

1. `ResultsMat = in_bst_results(ResultsFile, 1)` → `ImageGridAmp`, `nComponents`,
   `SurfaceFile`, `Time`.
2. If `nComponents == 3` (unconstrained): `S = bst_source_orient([], 3, [],
   ImageGridAmp, 'rms')` → `[nVert × nTime]` scalar (non-negative magnitude). If
   `nComponents == 1` (constrained): `S = ImageGridAmp` (signed `[nVert × nTime]`).
3. `[Eig, isComputed] = in_tess_eigenmodes(SurfaceFile)`; **guard `isComputed`**
   (`bst_error('Run "Compute eigenmodes" on this surface first.')` and return).
4. Load the surface (`Vertices`, `Faces`), build `M` via
   `tess_laplacian(Vertices, Faces, 'MassType', Eig.MassType)` (the `M` output only).

### B. Transform to eigenmode coefficients (reuse + thin wrapper)

```
Theta = bst_eigenmodes_project(Eig, S, M);     % [K x nTime] = Phi' * (M * S)
```

A thin wrapper `GetActivationCoeffs(ResultsFile)` performs A+B and returns
`Theta [K × nTime]` plus the metadata the viewer needs (`Values`, `Component`,
`CompRank`, `Time`). This keeps the GUI viewer thin and the data path headlessly
testable. No coefficients are persisted (transient view).

### C. The viewer — `toolbox/gui/view_eigenmode_spectrum.m`

`hFig = view_eigenmode_spectrum(ResultsFile)`:

1. Compute `Theta [K × nTime]` and metadata via `GetActivationCoeffs` (§A+B). On any
   failed guard, `bst_error` and return `[]`.
2. Create a Brainstorm-managed figure **registered with the results/recordings dataset**,
   so the **global time cursor drives it** (same time-sync mechanism source/time-series
   figures use; cf. how `figure_timeseries`/`figure_3d` receive global time updates). On
   each time change, redraw the spectrum at the current sample — no recomputation
   (`Theta` is precomputed once).
3. Render, per current time `t`:
   - power per mode `p_k = |θ_k(t)|²`, split by `Component` into **Left** (Component 1)
     and **Right** (Component 2) curves;
   - x-axis = eigenvalue `λ_k` (`Values`), with a toggle to spatial wavelength
     `≈ 2π/√λ` (`λ ≤ 0` → omitted / "n/a");
   - a **window-averaged overlay** `mean_t |θ_k|²` (dashed), defaulting to the full
     loaded time range and following the time panel's current selection when one is
     active.
4. Legend: current time, averaging window, per-hemisphere mode counts; figure title
   hints at ←/→ stepping and the eigenvalue/wavelength toggle. Guard all-zero /
   degenerate axis scaling.

**Pure helpers (headlessly testable):**

- `GetModalPower(ThetaCol, Component)` → `struct(.left, .right)` power vectors (the
  component split + magnitude-squared).
- `GetSpectrumAxis(Values, mode)` → `struct(.x, .label)` for
  `mode ∈ {'eigenvalue','wavelength'}` (`λ ≤ 0` handled).
- `GetWindowAverage(Theta, iWin)` → mean power per mode over the sample window.

### D. Integration (entry point) — `tree_callbacks.m`

On a **source-results** node (a `'results'`/source node with an `ImagingKernel` or
`ImageGridAmp`), add a **"View eigenspectrum"** item (always enabled; works read-only)
that calls `view_eigenmode_spectrum(ResultsFile)`. If the surface has no eigenmodes, the
viewer shows the friendly error (so the menu build does not read the surface on every
right-click). The node keeps all its existing menus unchanged; this is purely additive.

## Data flow

```
[source results node]  (dSPM/sLORETA/... ; kernel-link or full)
        │  right-click → "View eigenspectrum"
        ▼
GetActivationCoeffs(ResultsFile):
   S = in_bst_results(LoadFull=1).ImageGridAmp
       (nComponents==3 → bst_source_orient 'rms' → [nVert×nTime] scalar;
        nComponents==1 → signed [nVert×nTime])
   [Eig,isComputed] = in_tess_eigenmodes(SurfaceFile)        (isComputed guard)
   [~,M] = tess_laplacian(V, F, 'MassType', Eig.MassType)
   Theta = bst_eigenmodes_project(Eig, S, M)                 [K × nTime]
        ▼
view_eigenmode_spectrum:
   register figure with dataset (global time cursor drives it)
   at time t:  p_k = |θ_k(t)|² ; x = λ_k (toggle wavelength)
               curves: Left (Component 1) / Right (Component 2)
               + dashed window-average overlay mean_t |θ_k|²
   ←/→ step time (reuses Brainstorm time navigation)
```

## Error handling

- Surface has no eigenmodes → `bst_error` → "Run Compute eigenmodes on this surface first."
- Results file not loadable as a vertex field / no `SurfaceFile` → `bst_error`.
- Degenerate (all-zero) power at a time/window → guarded axis scaling.
- `λ ≤ 0` for the wavelength axis → that mode omitted / shown as "n/a".

## Testing strategy

Repo idiom: `dev/tests/*.m` printing `ALL TESTS PASSED`, run via the MATLAB MCP
`evaluate_matlab_code` (not `run_matlab_test_file`); prefix `rng('default')` if the
session is in legacy-RNG mode.

1. **`test_view_eigenmode_spectrum_pure.m`** (new): `GetModalPower` L/R split and
   magnitude-squared; `GetSpectrumAxis` eigenvalue vs wavelength (incl. `λ ≤ 0`
   handling); `GetWindowAverage` equals the manual mean over a window.
2. **`test_eigenmode_spectrum_acquire.m`** (new): on a small synthetic mesh with computed
   eigenmodes + mass matrix, assert (a) a **constrained** scalar field that *is* a pure
   eigenmode projects to a unit coefficient at that mode and ≈0 elsewhere; (b) an
   **unconstrained** `[3·nVert × nTime]` field is norm-collapsed then projected with the
   correct `[K × nTime]` shape and matches the manual `Φ'·M·rms(·)`. Validates the
   reused acquire+project composition end-to-end (no GUI).
3. **Downstream regression:** re-run the existing eigenmode suite
   (`test_eigenmodes_*_pure`, `test_io_eigenmodes_roundtrip`,
   `test_eigenmodes_perhemisphere`, `test_eigenmodes_manifold_gate`, and the
   `process_eigenmodes_*` option tests) — confirm green (this increment adds no changes
   to those paths, so they must remain untouched).
4. **Interactive (user):** right-click a source result → "View eigenspectrum" opens; L/R
   curves update live as the global time cursor steps; the window-averaged overlay; the
   eigenvalue/wavelength toggle; the no-eigenmodes friendly error; both a constrained and
   an unconstrained source model.

## Downstream impact & migration

None. No file format, results structure, or existing process is modified. The feature is
a transient viewer plus one additive menu item. Source maps without computed surface
eigenmodes simply show a friendly error pointing to Compute.

## Out of scope (YAGNI)

- **Manifold-harmonics source mapping** (the `[nEigenmodes × nChannels]` kernel / new
  inverse type) — a separate, already-implemented method, tested and viewed later.
- **Other vertex-mapped scalars** (PET-on-surface, anatomical overlays) — the projection
  helper is source-agnostic, so these are a later menu-surface extension, not a redesign.
- **The 2D joint map (mode × temporal frequency)** — deferred; the same `Theta [K×nTime]`
  coefficient stream can later feed the existing PSD/FFT path to produce a `[K × nFreqs]`
  `timefreq` rendered by `figure_timefreq`.
- **Persisting coefficients** — the view is transient; persistence is added only if/when
  the 2D map or batch analysis needs it.

## Files touched

- `toolbox/gui/view_eigenmode_spectrum.m` — **new** viewer, `GetActivationCoeffs`
  wrapper, and pure render helpers.
- `toolbox/tree/tree_callbacks.m` — "View eigenspectrum" item on source-results nodes.
- `dev/tests/test_view_eigenmode_spectrum_pure.m` (new);
  `dev/tests/test_eigenmode_spectrum_acquire.m` (new).
