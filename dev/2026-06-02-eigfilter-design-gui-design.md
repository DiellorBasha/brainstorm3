# Eigenmode filter design: extend panel_eigenmodes with spectral kernels

**Date:** 2026-06-02
**Author:** Diellor Basha
**Status:** Design approved (v1 scope), pending implementation plan

> Supersedes the earlier standalone-tool sketch. The existing `panel_eigenmodes`
> "EigenModes" panel already provides spectral weighting of the cortical field
> with a live 3D render and a source-filtering path; we extend it rather than
> build a parallel tool.

## Motivation

The `bst_eigfilter` library defines spectral kernels `g(λ)` on the cortical LBO
eigenmodes, but a kernel's *spatial* effect is hard to read from its analytic
form. The existing **EigenModes panel** (`panel_eigenmodes.m`) is already a live
"spectral-weighting of the cortical field" control: a per-mode weight `W` is
chosen (box / tapered / Gaussian window over a mode-index band) and rendered on
the cortex, and — when active on a source figure — it filters the source map via
`bst_eigenmodes_filter(Eig, u, M, 'custom', 'TransferFn', …)`. We **add spectral
kernels** as a new weighting source, plus a `g(λ)` preview and an impulse
(point-spread) display, so the user can design a kernel and see both its
frequency response and its spatial response in real time. Because the panel
already filters source maps, kernels immediately work as an analysis knob on real
data too.

## What already exists (reused verbatim)

- **State:** `GlobalData.UserModes` (SurfaceFile, nModes, Weights, WindowShape,
  Band, isActive, **CacheEig**, **CacheMass**).
- **Cache:** `EnsureCache(SurfaceFile)` loads `Eig` (Vectors/Values/CompRank) +
  mass matrix `M`.
- **Filter-apply path:** `ApplyToColumn(SurfaceFile, u)` →
  `bst_eigenmodes_filter(Eig, u, M, 'custom', 'TransferFn', @(l) wRaw(:))`.
- **Live-render loop:** control change → `RecomputeWeights` → `NotifyChanged` →
  `bst_figures('FireModesChanged')` → for an eigenmode-view figure
  `view_eigenmodes('ModesChangedCallback', hFig)`, else
  `panel_surface('UpdateSurfaceData', hFig)`.
- **Launch:** surface tree menu "View eigenmodes" → `view_eigenmodes(SurfaceFile)`
  → opens the cortex figure + shows the EigenModes panel.

## Scope (v1)

In scope:
1. **Kernel weighting** in the panel — a "Kernel" option (alongside Box/Taper/
   Gauss) with a kernel dropdown + parameter sliders; the weight/transfer function
   becomes `g(λ)` from `bst_eigfilter_kernel`.
2. **g(λ) companion figure** — a small live plot of `g(λ)` over the eigenvalue
   axis.
3. **Delta point-spread display** — in a "View eigenmodes" figure, show the kernel
   applied to an impulse at a chosen vertex (the spatial point-spread), with a
   default/index/pick vertex selector.
4. **Analysis knob (free):** kernels apply to a real source map when the panel
   lever is Active on a source figure (existing `ApplyToColumn` path).

Out of scope (future): filterbank preview (vector-scale kernels); explicit
UI-range hints in kernel metadata; saving a kernel-filtered source as a new node.

## Design decisions (resolved in brainstorming)

1. **Extend `panel_eigenmodes`**, not a standalone tool.
2. **Add "Kernel" alongside** the existing window shapes (additive).
3. **Delta point-spread** is the design preview display (impulse → point-spread).
4. **g(λ)** lives in a **small companion MATLAB figure** (not embedded in the
   Swing panel — embedding an axes in Swing is fragile).

## Component design

### A. panel_eigenmodes — new controls

Added to the existing Swing panel (keep all current controls):
- **"Kernel" radio** in the shape `ButtonGroup` → `{Box, Taper, Gauss, Kernel}`.
- **Kernel dropdown** populated from `bst_eigfilter_kernel('list')` (display names
  from each kernel's `'meta'`), enabled only when "Kernel" is selected.
- **Two generic parameter sliders + readouts** (`Param1`, `Param2`). Every current
  kernel has ≤ 2 scalar parameters (e.g. heat `t`; matern `kappa`,`nu`; dog
  `t1`,`t2`; `ideal`'s 2-element `band` → the two sliders). On kernel change, read
  `meta.params`, label/range/enable the sliders from it (log-scaled bounds
  heuristic `[default/100, default·100]`, clamped to the meta `range`; defaults
  from `meta.params.<p>.default`). Unused sliders are disabled.
- **Vertex row** (for the delta display): an index field + a "Pick on surface"
  toggle (default vertex = nearest to surface centroid; pick reuses the
  scout-style click→nearest-vertex on the cortex figure).
- **Display toggle** (eigenmode-view figures only): "Synthesis" vs "Delta
  point-spread".

### B. Weight / transfer-function computation

Extend `RecomputeWeights`:
- **Window shapes** (box/tapered/gauss): unchanged (`BuildWeights`).
- **Kernel:** build `g = bst_eigfilter_kernel(name, params)`. Store the handle in
  state (`GlobalData.UserModes.KernelFn = g`) and the paired-rank display weights
  `W(k) = g(λ_paired_k)` for the synthesis view (λ_paired_k = the paired rank's
  representative eigenvalue; hemispheres are near-symmetric so the pair shares
  ≈ one λ).

Extend `ApplyToColumn`: when in kernel mode, pass `TransferFn = KernelFn`
(evaluated by `bst_eigenmodes_filter` at each **raw** mode's λ — exact); otherwise
keep `@(l) wRaw(:)`. This makes the source-map analysis knob exact per mode.

### C. Delta point-spread display (view_eigenmodes)

Extend `view_eigenmodes`'s `EigenView` appdata with `DisplayMode`
('synthesis' | 'delta') and `DeltaVertex`. In `ModesChangedCallback`:
- **synthesis** (existing): `col = SynthColumn(PairedGrid, W)`.
- **delta:** `e = zeros(nVert,1); e(DeltaVertex)=1;
  col = bst_eigenmodes_filter(Eig, e, M, 'custom', 'TransferFn', g)` — the kernel's
  point-spread. Auto-scale the color limits to the response range (symmetric for
  sign-changing band-pass kernels so the diverging colormap centers at zero).
The render path (write `ImageGridAmp`, `panel_surface('UpdateSurfaceData')`) is
unchanged.

### D. g(λ) companion figure

A lightweight MATLAB figure (created lazily, reused) with one axes: `g(λ)` over
`λgrid = linspace(0, max(Eig.Values), N)`, the real `λ_k` drawn as faint x-axis
ticks, title = kernel display name. `RecomputeWeights` (kernel mode) updates the
line via `set(hLine, 'YData', …)`. Closed/ignored in window-shape mode. A small
helper `view_eigfilter_response(...)` (or a subfunction) owns create/update.

## Data flow (kernel mode, delta display)

```
slider drag
  -> panel: read kernelName + params
  -> g = bst_eigfilter_kernel(name, params)
  -> update g(λ) companion figure (set YData)
  -> store KernelFn; RecomputeWeights; NotifyChanged
  -> bst_figures('FireModesChanged')
       -> view_eigenmodes('ModesChangedCallback', hFig)   [delta mode]
            col = bst_eigenmodes_filter(Eig, δ_vertex, M, 'custom', g)
            write ImageGridAmp; panel_surface('UpdateSurfaceData')   [cortex repaints]
```

## Testing

- **Pure** (headless):
  - kernel→weights: for a kernel `g`, the panel's paired weights equal
    `g(λ_paired)` and `ApplyToColumn`'s transfer function reproduces
    `bst_eigenmodes_filter(..., 'custom', g)`.
  - delta point-spread character on a small synthetic 2-component mesh: `flat` →
    concentrated at the vertex (band-limited spike, not exact); `heat`/`tikhonov` →
    smooth non-negative blob centered at the vertex, broadening with the scale;
    `mexhat`/`dog` → sign-changing center-surround.
  - slider-range heuristic maps a param meta to sensible bounds.
- **Smoke** (needs a surface with eigenmodes): the panel builds with the new
  controls; selecting a kernel + dragging a param updates the cortex `Data` and the
  g(λ) figure; switching Synthesis↔Delta works; toggling back to a window shape
  restores prior behavior.
- **Manual:** drag across all 10 kernels in a "View eigenmodes" figure; pick
  vertices; verify the point-spread and g(λ) update; then make the lever Active on
  a real source figure and confirm the kernel filters it.

## Files

**Modified:**
- `toolbox/gui/panel_eigenmodes.m` — Kernel radio + kernel dropdown + 2 param
  sliders + vertex row + display toggle; `RecomputeWeights` kernel branch (+
  `KernelFn` state); `ApplyToColumn` kernel branch; drive the g(λ) figure.
- `toolbox/gui/view_eigenmodes.m` — `DisplayMode`/`DeltaVertex` in `EigenView`;
  `ModesChangedCallback` delta branch (point-spread).

**New:**
- `toolbox/gui/view_eigfilter_response.m` — the small g(λ) companion figure
  (create/update), or an equivalent subfunction in `panel_eigenmodes`.
- `dev/tests/test_eigfilter_design_pure.m` — kernel→weights, point-spread
  character, slider-range heuristic.
- `dev/tests/test_eigfilter_design_smoke.m` — build-and-drive the extended panel
  (skips if no eigenmode surface).

No new tree menu needed — reached via the existing "View eigenmodes" launch.

## Future work (documented, not built)

- Save a kernel-filtered source map as a new results node (explicit analysis
  output, vs the live lever).
- Filterbank preview (vector-scale kernels → multiple g(λ) curves / scale-space).
- Explicit `uirange`/`uiscale` hints in kernel metadata (replace the GUI heuristic).
- Compose two kernels (`bst_eigfilter_compose`) from the panel.
