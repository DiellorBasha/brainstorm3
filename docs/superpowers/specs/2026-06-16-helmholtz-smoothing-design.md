# Helmholtz + smoothing — one panel: band-limit the field, then decompose

**Date:** 2026-06-16
**Status:** Design approved; implementation pending
**Branch:** `feat/helmholtz-view` (brainstorm3)
**Author:** Diellor Basha

## Problem

The Helmholtz components view detects vortex cores / sources as raw 1-ring extrema of the
stream function ψ / potential φ. On a noisy single-frame MEG reconstruction this over-fires
badly (≈59 vortex cores, ≈940 sources on one real frame) because the detector is local,
threshold-free, and single-scale — it marks every noise-scale wiggle, not the handful of
real vortices. `div`/`curl` are derivatives, so they amplify high-frequency noise (φ is much
rougher than ψ, hence the 940 vs 59).

The principled fix is to **band-limit the source field to a chosen spatial scale before
decomposing** — exactly what the Dirac eigenmode filter does. A band-limited field has a
bounded number of critical points, so the cores thin out *by construction*, and the field
you see (the quiver) matches the field the cores come from. The standalone Spatial Filter
already implements this filtering; it should live inside the Helmholtz tool so you can tune
the smoothing scale and watch the cores settle, instead of being a separate round-trip.

## Goal

Fold the Spatial Filter into the Helmholtz panel as a **Smoothing** section, so one tool:
1. low-passes / band-selects the **active frame's** field in the Dirac eigenbasis
   (kernel + scale, the shared `bst_eigfilter_panel` section),
2. decomposes the **smoothed** field (Total / Irrotational / Solenoidal / Harmonic),
3. detects component-aware markers on the smoothed field, further pruned by a **magnitude
   gate** (drop near-zero `|ω|` cores / `|div|` sources),
4. can **save** the smoothed source as a new results node (whole-series export).

The smoothing scale + the gate are the two knobs for optimizing vortex-core tracking.

## Non-goals

- Spectral *design* (wavelets) — that's the Wavelet Designer.
- A topological-persistence core detector — the scale filter + magnitude gate are the first
  pass; persistence is a later upgrade if still needed.
- Changing the Dirac decomposition math (`bst_dirac_helmholtz` unchanged — the Dirac-based
  version the user kept).
- Keeping the standalone Spatial Filter — it is removed (folded in).

## Decisions (from brainstorming)

| Question | Decision |
|---|---|
| Standalone Spatial Filter | **Fold into Helmholtz**: remove the figure-popup item and `panel_spatial_filter` (+ its test); the smoothing lives in the Helmholtz panel. `process_dirac_filter` (batch) is untouched. |
| Smoothing scope | **Active frame, on-demand**: filter the displayed frame's field before decomposing; slider re-filters + re-decomposes the current frame; the Total quiver shows the smoothed field. A **Save smoothed file** button exports the whole filtered series. |
| Marker pruning | **Scale + magnitude gate**: the smoothing band-limits div/curl; a threshold slider additionally drops extrema with `|ω|`/`|div|` below a fraction of the frame max. |

## Architecture

### Loading (view_helmholtz launch)

In addition to today's Dirac **operator** + LBO operator (for div/curl + Poisson), load the
Dirac **eigenbasis** for the filter (as the Spatial Filter does):
`EigenMat = tess_eigen(SurfaceFile, 'Dirac')`, `OpMat = load(EigenMat.OperatorFile)`,
`Mass = OpMat.Mass`, `Lambda = EigenMat.Lambda{1}`. Store `EigenMat, Mass, Lambda` on the
figure state.

### Per-frame pipeline (UpdateFrame)

1. Fetch the active frame `Jt [3nV×1]` (as today).
2. **Smooth** (if on): `g = bst_eigfilter_kernel(name, params)`;
   `Jt = real(bst_dirac_eigenmodes_filter(EigenMat, Mass, Jt, 'custom', 'TransferFn', g))`.
3. **Decompose** the smoothed `Jt`: `Ht = bst_dirac_helmholtz('Frame', Op, Jt)` (cached per
   frame *at the current smoothing setting*).
4. Pick the component (Total/Irrot/Solen/Harm) → override the quiver (`QuiverVectorOverride`
   = the smoothed component field) and the cortex colormap (the component scalar) as today.
5. **Markers**: `comp.Markers` (ψ extrema for Solen, φ extrema for Irrot), then **gate**:
   keep markers with `|omega| ≥ τ · max(|omega|)` over that frame's markers, τ from the
   threshold slider (`τ=0` keeps all). Draw + readout the gated counts.

The decomposition cache (`St.Cache`, keyed by frame index) is **invalidated when the
smoothing changes** (new kernel/scale → the cached frames were decomposed on the old field);
it is *not* invalidated by component / gate / vector-toggle changes (those are
post-decomposition).

### Panel (`panel_helmholtz`)

Two stacked sections + controls:
- **Smoothing** (the shared `bst_eigfilter_panel` section): kernel dropdown (default
  **Heat / low-pass**) + auto scale sliders, plus a **Smoothing on** checkbox (default
  **off** so the raw baseline shows first). Slider-settle / kernel-change →
  `view_helmholtz('SetSmoothing', hFig, isOn, name, params)`.
- **Component** radio: Total / Irrotational / Solenoidal / Harmonic.
- **Show vectors**, **Show singular points** checkboxes.
- **Marker threshold** slider (`τ`, 0→aggressive) → `view_helmholtz('SetGate', hFig, τ)`.
- **Readout**: component-aware counts (post-gate) + harmonic %.
- **Save smoothed file** button → `view_helmholtz('SaveSmoothed', hFig)`.
- **Close**.

The panel needs `Lambda` to build the scale sliders, passed in at create
(`panel_helmholtz('CreatePanel', hFig, Lambda)`).

### Save smoothed file

Materialize the full series `J [3nV×nT]` from the source, filter it with the current `g`,
and save a new results node in the source's study — reusing `panel_spatial_filter`'s
`SaveFiltered` logic (`in_bst_results` → replace `ImageGridAmp`, `Comment = '<src> | dirac
filter(<info>)'`, `bst_history` + `db_add_data`).

### State (`HelmholtzState`) additions

`EigenMat, Mass, Lambda, Smooth = struct('on',bool,'name',kernelKey,'params',struct),
GateFrac` — alongside the existing `Op, srcDS, srcResult, Component, ShowVectors,
ShowMarkers, iTess, nV, Cache`.

## Components / files

**Modify:**
- `toolbox/gui/view_helmholtz.m` — load eigenbasis+mass; smooth the active frame before
  decomposing; gate markers; dispatch `SetSmoothing`/`SetGate`/`SaveSmoothed`; cache
  invalidation on smoothing change.
- `toolbox/gui/panel_helmholtz.m` — add the shared **Smoothing** kernel section + on/off,
  the **Marker threshold** slider, the **Save smoothed file** button; wire callbacks; take
  `Lambda` at create.
- `toolbox/gui/figure_3d.m` — remove the "Spatial filter (Dirac)" popup item.

**Remove:**
- `toolbox/gui/panel_spatial_filter.m`, `dev/tests/test_spatial_filter.m` (folded in).

**Reuse unchanged:** `bst_eigfilter_panel`, `bst_eigfilter_kernel`,
`bst_dirac_eigenmodes_filter`, `tess_eigen`, `bst_dirac_helmholtz` (Dirac Frame), the
`QuiverVectorOverride` native-vector path, `bst_colormaps`.

**Test:**
- `dev/tests/test_helmholtz_view.m` — extend: turning smoothing on (heat low-pass) reduces
  the vortex-core count vs raw; raising the marker threshold reduces drawn markers
  monotonically; the decomposition cache clears on a smoothing change; `SaveSmoothed`
  writes a 3-component results node into the source study.

## Data flow

1. Tree node → "Helmholtz / vorticity (Dirac)" → load operators **+ Dirac eigenbasis**,
   open native figure (Total, smoothing off).
2. Tick **Smoothing on**, pick **Heat**, drag the **scale** slider → the active frame is
   re-filtered + re-decomposed; the quiver shows the smoothed field; cores thin out.
3. Pick **Solenoidal** → smoothed `∇⊥ψ` + vortex cores (now far fewer); nudge **Marker
   threshold** to drop the last weak ones.
4. Scrub time → each frame is filtered + decomposed at the current setting.
5. **Save smoothed file** → filtered source saved as a new node.

## Error handling

- No Dirac eigenbasis / nxr unavailable: `tess_eigen('Dirac')` find-or-creates (progress
  bar); abort cleanly if it can't be produced (as the Spatial Filter does).
- Smoothing-off path is identical to today (no filter applied).
- Field/basis vertex mismatch: abort with a clear message (mirrors the Spatial Filter).
- Gate with no markers, or all-equal `|omega|`: keep all (guard `max(|omega|)=0`).
- Teardown unchanged (delete figure, drop overlays; no temp node).

## Testing

- **Smoothing thins cores** (`test_helmholtz_view`, real or synthetic): core count with a
  heat low-pass on < count with smoothing off.
- **Gate monotonic**: drawn-marker count is non-increasing as `τ` rises; `τ=0` shows all.
- **Cache invalidation**: changing the smoothing scale clears `St.Cache` (recomputes).
- **Save smoothed**: the saved node loads to a 3-component series in the source study.
- **Regression**: the component states / colormaps / vectors / time-following / close +
  stale-guard checks still pass.

## Build order

1. view_helmholtz: load eigenbasis + smooth the active frame (+ cache invalidation);
   `SetSmoothing` dispatch.
2. Marker magnitude gate (`SetGate`) + gated readout.
3. panel_helmholtz: Smoothing kernel section + on/off + threshold slider + Save smoothed;
   wire callbacks; extend the test.
4. Remove the Spatial filter popup item + `panel_spatial_filter` + its test.
