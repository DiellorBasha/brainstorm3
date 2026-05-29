# Eigenmode Spectrum: On-the-fly Modal Power Spectrum via an Additive ModeKernel — Design

**Date:** 2026-05-29
**Status:** Approved design pending spec review
**Follows:** `2026-05-28-eigenmode-perhemisphere-viewer-r2-design.md` (eigenmode compute + viewer, shipped on `feature/eigenmode-context-menus`)

## Background

The scalar LBO eigenmode pipeline can now compute per-hemisphere eigenmodes on a
cortex and view them. The next scaffolding step is an **eigenmode spectrum**: take
Brainstorm's *default* source mapping (min-norm / dSPM / sLORETA — explicitly **not**
the team's eigenmode source-mapping work), express the source activity in the
eigenbasis, and view how the source map's energy distributes across spatial
eigenmodes. This is the spatial analogue of the existing FFT/PSD spectrum and the
foundation for later, richer eigenmode-based analysis (a 2D mode × temporal-frequency
map, and eventually eigenmode source mapping).

### Key facts that shaped the design (from codebase exploration)

- **Canonical on-the-fly source mapping.** A *kernel-link* results file stores
  `ImagingKernel [nVertices·nComp × nChannels]` plus a `DataFile` pointer to the
  recordings; source maps are reconstructed lazily as
  `ImageGridAmp = ImagingKernel · F(GoodChannel,:)` (`in_bst_results`, LoadFull=1).
  The global time cursor and all stepping/analysis machinery operate on this.
- **The eigenmode projection is M-weighted.** For an M-orthonormal eigenbasis
  `Φ [nVert × K]`, the coefficients of a vertex field `S` are `Θ = Φ'·M·S`
  (`bst_eigenmodes_project`). The existing eigenmode processes already build the
  mass matrix `M` consistent with `Eigenmodes.MassType`.
- **Linearity holds only for constrained sources.** For `nComponents == 1`
  (one scalar per vertex), `K_vertex` is `[nVert × nCh]` and a single linear
  transform `ModeKernel = Φ'·M·K_vertex [K × nCh]` is exact. For unconstrained
  sources the orientation collapse (a norm) is nonlinear and cannot be folded into
  one matrix — **deferred** to a later increment (planned: a globally consistent
  tangent reference frame via the vector heat method).
- **The results structure is extensible.** `db_template('resultsmat')` defines the
  results fields; downstream code reads specific fields and ignores unknown ones, so
  an additive optional field is safe.

## Goals

1. A **one-shot transform** that turns an existing constrained source kernel into an
   eigenmode kernel `ModeKernel [K × nChannels]`, stored **additively** on the
   results structure (leaving `ImagingKernel` untouched).
2. A dedicated **modal power-spectrum viewer**: power per mode vs **eigenvalue**
   (toggle to spatial wavelength), **left/right hemispheres as separate curves**,
   with **time stepping** (driven by Brainstorm's global time cursor) **and** a
   **window-averaged overlay**.
3. Maximal reuse of the canonical on-the-fly source-mapping machinery; no breakage of
   any downstream consumer of `ImagingKernel`.

## Design

### A. Storage — additive `ModeKernel` on the results structure

The eigenmode kernel is stored **alongside** the standard imaging kernel in the same
(shared kernel-link) results file:

- `ImagingKernel [nVert × nChannels]` — **unchanged**; continues to power all existing
  vertex-space source visualizations, scouts, stats, exports.
- `ModeKernel [K × nChannels]` — **new**; `ModeKernel = Φ'·M·ImagingKernel`. Because it
  is built from the same `ImagingKernel`, `Θ = ModeKernel·F` is provably the M-weighted
  projection of the default source map `S = ImagingKernel·F` (faithfulness guarantee).
- `ModeInfo` — **new** companion snapshot captured at transform time so the viewer is
  self-contained and immune to a later surface recompute changing `K`/ordering:
  - `.Values [K×1]` — eigenvalues `λ_k` (ascending within each component)
  - `.Component [K×1]` — hemisphere/component id per mode (1 = Left, 2 = Right)
  - `.CompRank [K×1]` — within-component rank
  - `.nComponents` — number of **mesh** connected components (hemispheres), typically 2
    (distinct from the results file's source-orientation `nComponents`, which the
    transform requires to be 1)
  - `.SurfaceFile`, `.MassType` — provenance
  - `.SourceMethod` — the source method that produced `ImagingKernel` (provenance)

`ModeKernel` and `ModeInfo` are registered as optional fields (default `[]`) in
`db_template('resultsmat')` so they survive any re-save and are documented. Column
alignment of `ModeKernel` matches `ImagingKernel` (same `GoodChannel`), so the same
`F(GoodChannel,:)` selection applies.

The node remains an ordinary source-results node (it still has `ImagingKernel`), so it
keeps all its existing menus and **additionally** offers "Eigenmode spectrum" when
`ModeKernel` is present. No tree gating is needed to suppress cortical views.

### B. The transform — `process_eigenmodes_kernel.m` + pure core

**Pure core** `toolbox/math/bst_eigenmode_kernel.m`:

```
ModeKernel = bst_eigenmode_kernel(Phi, M, Kvertex)
  Phi      [nVert × K]   M-orthonormal eigenmodes
  M        [nVert × nVert] sparse mass matrix
  Kvertex  [nVert × nCh] constrained imaging kernel
  ->
  ModeKernel [K × nCh] = Phi' * (M * Kvertex)
```

Headlessly testable; no I/O.

**`ComputeInteractive(iStudy, ResultsFile)`** (mirrors the R2 idiom):

1. Load the results link (`ImagingKernel`, `GoodChannel`, `SurfaceFile`, `DataFile`,
   `nComponents`, source method).
2. **Guard `nComponents == 1`** (constrained). Otherwise `bst_error` naming the planned
   unconstrained support and return (no-op).
3. Load eigenmodes from `SurfaceFile` (`in_tess_eigenmodes` → `Φ`, `Values`,
   `Component`, `CompRank`, `nComponents_mesh`, `MassType`). **Guard `isComputed`**;
   otherwise `bst_error('Run "Compute eigenmodes" first.')` and return.
4. Build the mass matrix `M` from the surface, consistent with `MassType`, via the
   same mass-matrix routine the existing eigenmode processes use. *(Implementation
   note for the plan: if that routine is currently inline in `tess_eigenmodes`, expose
   a small reusable helper rather than duplicating it.)*
5. `ModeKernel = bst_eigenmode_kernel(Φ, M, ImagingKernel)`.
6. Write `ModeKernel` + `ModeInfo` into the results file (additive), save, and
   `db_reload_studies` so the new capability appears.

If `ModeKernel` already exists on the file, confirm overwrite (decline → no-op).

### C. The viewer — `toolbox/gui/view_eigenmode_spectrum.m`

`hFig = view_eigenmode_spectrum(ResultsFile)`:

1. Load `ModeKernel` + `ModeInfo` from the results file; if absent → `bst_error`
   directing to run the transform.
2. Load the linked recordings once and **precompute** `Θ = ModeKernel · F(GoodChannel,:)`
   `[K × nTime]` (mirrors how kernel-link source figures realize `ImageGridAmp` on
   open — responsive stepping with no per-step matmul).
3. Create a Brainstorm-managed figure **registered with the recordings dataset**, so
   the **global time cursor drives it** (the same time-sync mechanism source/time-series
   figures use; cf. how `figure_timeseries`/`figure_3d` receive global time updates).
   On each time change, redraw the spectrum at the current sample.
4. Render, per current time `t`:
   - power per mode `p_k = |θ_k(t)|²`, split by `ModeInfo.Component` into **Left** and
     **Right** curves;
   - x-axis = eigenvalue `λ_k` (`ModeInfo.Values`), with a toggle to spatial wavelength
     `≈ 2π/√λ` (`λ ≤ 0` → omitted/"n/a");
   - a **window-averaged overlay** `mean_t |θ_k|²` over the selected time window
     (dashed), updated when the window changes (defaults to the full loaded time range;
     follows the time panel's current selection when one is active).
5. Legend: current time, window, per-hemisphere mode counts; figure title hints at
   ←/→ stepping and the eigenvalue/wavelength toggle. Guard all-zero/degenerate scaling.

**Pure helpers (headlessly testable):**

- `GetModalPower(ThetaCol, Component)` → `struct(.left, .right)` power vectors (the
  component split + magnitude-squared).
- `GetSpectrumAxis(Values, mode)` → `struct(.x, .label)` for `mode ∈ {'eigenvalue','wavelength'}`.
- `GetWindowAverage(Theta, iWin)` → mean power per mode over the sample window.

### D. Integration (entry points) — `tree_callbacks.m`

On a **source-results** node that has an `ImagingKernel` and a constrained model
(`nComponents == 1`), add (gated by `~bst_get('ReadOnly')` for the transform path):

- **"Eigenmode spectrum"** — if the file already has `ModeKernel`, call
  `view_eigenmode_spectrum(ResultsFile)`; otherwise call
  `process_eigenmodes_kernel('ComputeInteractive', iStudy, ResultsFile)` (which persists
  `ModeKernel`) and then open the viewer. A single consolidated item, transform-if-needed
  then view (mirrors R2's compute-then-auto-open).

## Data flow

```
[default CONSTRAINED source kernel link]   (existing; ImagingKernel [nVert×nCh])
        │  right-click → "Eigenmode spectrum"
        ▼
ComputeInteractive (if ModeKernel absent):
   load ImagingKernel, GoodChannel, SurfaceFile, DataFile, nComponents(==1 guard)
   load Φ, Values, Component, CompRank from SurfaceFile (isComputed guard)
   build M per MassType
   ModeKernel = Φ' · (M · ImagingKernel)          [K × nCh]
   write ModeKernel + ModeInfo into results file; db_reload_studies
        ▼
view_eigenmode_spectrum:
   precompute Θ = ModeKernel · F(GoodChannel,:)    [K × nTime]
   register figure with recordings dataset (global time cursor drives it)
   at time t:  p_k = |θ_k(t)|² ; x = λ_k (toggle wavelength)
               curves: Left (Component 1) / Right (Component 2)
               + dashed window-average overlay mean_t |θ_k|²
   ←/→ step time (reuses Brainstorm time navigation)
```

## Error handling

- No eigenmodes on the surface → `bst_error` → "Run Compute eigenmodes first."
- Unconstrained kernel (`nComponents ≠ 1`) → `bst_error` naming the planned unconstrained
  (vector-heat tangent-frame) support; no-op.
- Results file without `ImagingKernel`/`DataFile` (not a kernel link) → `bst_error`.
- `ModeKernel` already present → overwrite confirm; decline → no-op.
- Viewer on a file lacking `ModeKernel` → `bst_error` directing to the transform.
- Degenerate (all-zero) power at a time/window → guarded axis scaling.

## Testing strategy

Repo idiom: `dev/tests/*.m` printing `ALL TESTS PASSED`, run via the MATLAB MCP
`evaluate_matlab_code` (not `run_matlab_test_file`); prefix `rng('default')` if the
session is in legacy-RNG mode.

1. **`test_eigenmode_kernel_pure.m`** (new): on a small synthetic mesh + eigenbasis +
   random `Kvertex`, assert (a) dimensions `[K × nCh]`; (b) associativity
   `ModeKernel·F == Φ'·M·(Kvertex·F)`; (c) M-orthonormality identity — feeding a
   `Kvertex` whose single column is a pure eigenmode field yields a unit coefficient at
   that mode and ≈0 elsewhere.
2. **`test_view_eigenmode_spectrum_pure.m`** (new): `GetModalPower` L/R split and
   magnitude-squared; `GetSpectrumAxis` eigenvalue vs wavelength (incl. `λ ≤ 0`
   handling); `GetWindowAverage` equals the manual mean over a window.
3. **`test_eigenmode_kernel_transform.m`** (new, in-memory integration): build a tiny
   synthetic results struct (`ImagingKernel`, `GoodChannel`, `nComponents=1`) + a
   synthetic recordings matrix; run the transform core; assert `Θ = ModeKernel·F` equals
   the projected default source map and that `ModeInfo` carries correct `Values`/`Component`.
4. **Downstream regression:** re-run the existing eigenmode suite
   (`test_eigenmodes_*_pure`, `test_io_eigenmodes_roundtrip`,
   `test_eigenmodes_perhemisphere`, `test_eigenmodes_manifold_gate`, and the
   `process_eigenmodes_*` option tests) — confirm green.
5. **Interactive (user):** right-click constrained source result → "Eigenmode spectrum"
   transforms (first time) and opens; L/R curves update live as the time cursor steps;
   the window-averaged overlay; the eigenvalue/wavelength toggle; the unconstrained and
   no-eigenmodes friendly errors.

## Downstream impact & migration

`ImagingKernel` is untouched, so existing source workflows are unaffected. `ModeKernel`
is purely additive and small (`K·nCh`, e.g. 600×300 ≈ 1.4 MB vs the ~36 MB imaging
kernel). Older files without `ModeKernel` simply lack the "Eigenmode spectrum" capability
until the transform is run.

## Out of scope (YAGNI)

- **Unconstrained sources** — deferred to a later increment (globally consistent tangent
  frame via the vector heat method).
- **The 2D joint map (mode × temporal frequency)** — deferred, but kept nearly free:
  the same `ModeKernel·F` coefficient stream can later feed the existing PSD/FFT path to
  produce a `[K × nFreqs]` `timefreq` rendered by `figure_timefreq`.
- **Eigenmode → cortex reconstruction (`Φ·Θ`) as a first-class inverse type** — the
  Option-B "ModeKernel-in-ImagingKernel + on-the-fly inversion" design; revisited when
  the eigenmode source-mapping work matures (`ModeKernel` is forward-compatible with it).
- **Surface/source types beyond constrained cortex** for the menu.

## Files touched

- `toolbox/math/bst_eigenmode_kernel.m` — **new** pure transform core.
- `toolbox/process/functions/process_eigenmodes_kernel.m` — **new** `ComputeInteractive`
  + thin wrapper around the core; writes `ModeKernel`/`ModeInfo`.
- `toolbox/gui/view_eigenmode_spectrum.m` — **new** viewer + pure helpers.
- `toolbox/db/db_template.m` — add optional `ModeKernel`/`ModeInfo` to `resultsmat`.
- `toolbox/tree/tree_callbacks.m` — "Eigenmode spectrum" item on constrained source-results nodes.
- `dev/tests/test_eigenmode_kernel_pure.m` (new); `dev/tests/test_view_eigenmode_spectrum_pure.m`
  (new); `dev/tests/test_eigenmode_kernel_transform.m` (new).
