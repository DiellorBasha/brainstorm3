# Spatial Filter — live in-view filtering of a Dirac source map

**Date:** 2026-06-15
**Status:** Design approved; implementation pending
**Branch:** `feat/spatial-filter` (brainstorm3)
**Author:** Diellor Basha

## Problem

The Wavelet Designer is *localized* (a vertex + direction + kernel atom). The other
half of the Dirac filterbank is a pure **filter**: no localization, it just reshapes an
existing cortical source vector field across the Dirac eigenvalue spectrum to isolate a
spatial scale (smooth / detail / band). The batch version exists
(`process_dirac_filter`: filter a whole results file → new file), but there is no
**interactive** way to filter the source map you are already looking at.

The target workflow: compute **Dirac-dSPM** unconstrained sources, **Display on cortex**
(a 3-D view that shows the source vector field at each time step), then open a small
**Spatial Filter** control and watch the displayed field get filtered — live, in place,
following the time cursor.

## Goal

A **Spatial Filter** control panel that:

1. Attaches to an open Dirac source figure (launched from its popup).
2. Filters the displayed source vector field with a kernel + scale (the Wavelet
   Designer's *Filter kernel* section, reused).
3. Applies **in place** (the figure's vectors + colormap become the filtered field),
   **non-destructively** (the original is restored on off/close), and **follows time**
   (each frame shows its filtered spatial field).
4. Can **save** the filtered result as a new results node.

## Non-goals

- Localization / direction / chirality / tiling — those belong to the Wavelet Designer.
- A persistent docked tool tab — the panel is launched per source figure.
- Temporal filtering — the filter is purely spatial (independent per time frame).
- Changing the Dirac filter math (`bst_dirac_eigenmodes_filter`) — reused unchanged.

## Decisions (from brainstorming)

| Question | Decision |
|---|---|
| Display model | **In-place toggle** on the same figure (non-destructive; re-filters on time scrub) |
| Launch point | **Source figure popup** → "Spatial filter (Dirac)"; panel attaches to that figure |
| Save | **Live view + a "Save filtered file" button** (new results node; batch already exists in `process_dirac_filter`) |
| Kernel section | **Factor out a shared helper** from `panel_wavelet_designer` used by both panels |
| Compute model | **Whole-series swap**: filter every time step's spatial field up front, swap the in-memory `ImageGridAmp`; scrubbing is then instant. Each frame's spatial field is filtered independently — exactly the per-time-step result, with no per-scrub recompute |

## Architecture

### Shared kernel-section helper (refactor)

Extract the *Filter kernel* UI + logic currently inside `panel_wavelet_designer` into a
reusable module so both panels share one implementation:

`toolbox/gui/bst_eigfilter_panel.m` (dispatched like a panel helper):
- `kctrl = bst_eigfilter_panel('Create', jParent, Lambda, onChange)` — builds the kernel
  dropdown (curated display names via `bst_eigfilter_kernel`) + the auto mode-index scale
  sliders (rebuilt per kernel from kernel `meta.params`), wires kernel-change (rebuild
  sliders) and slider-settle to the caller's `onChange` handle, and returns a handle
  bundle `kctrl` (combobox, the sliders' parent panel, `Lambda`).
- `[kernelName, params] = bst_eigfilter_panel('Read', kctrl)` — reads the current kernel
  name and its params (mode index → param via the shared `1/lambda_k` / `lambda_k` map),
  ordering `t1 < t2` for dog.

`panel_wavelet_designer` is refactored to build its Section 2 via this helper (and its
`BuildParamWidgets`/`ReadParams`/`i_param_value`/`i_param_label` move into the helper).
The wavelet designer test suite is re-run to confirm no behaviour change.

### Spatial Filter panel

`toolbox/gui/panel_spatial_filter.m`:

- `panel_spatial_filter('Start', hFig)` — the launcher: resolve the figure's
  `iDS`/`iResult` (`bst_memory('GetDataSetResult')` on the displayed results), verify it
  is an unconstrained (3-component) surface source, find-or-create the surface's Dirac
  eigenbasis (`tess_eigen('Dirac')`) + operator mass, **back up** the original
  `ImageGridAmp`/`ImagingKernel`, build the panel (shared kernel section + Filter on/off
  + Save filtered file + Close) attached to `hFig`, and dock it (`gui_show`,
  `BrainstormTab`).
- Controls: the shared kernel section (default kernel = **Heat / low-pass**), a **Filter
  on/off** checkbox (**default OFF** — the figure shows the raw source until the user
  ticks it on), **Save filtered file**, **Close**.
- State (on the figure appdata, keyed `SpatialFilterState`): `iDS, iResult, EigenMat,
  Mass, OrigImageGridAmp, OrigImagingKernel, isOn`.

### Filter mechanism (whole-series, non-destructive)

`Apply(panelName)` (called on filter-on and on any kernel/scale change while on):
1. Materialize the full field `J [3nV × nT]` from the ORIGINAL data
   (`bst_memory('GetResultsValues', iDS, iResult, [], 'all', 0)` → reshape to `3nV × nT`;
   reconstructs from the imaging kernel if the source is kernel-only).
2. `[name, params] = bst_eigfilter_panel('Read', kctrl)`;
   `Jf = real(bst_dirac_eigenmodes_filter(EigenMat, Mass, J, 'custom', 'TransferFn',
   bst_eigfilter_kernel(name, params)))`.
3. Swap into the figure's in-memory results:
   `GlobalData.DataSet(iDS).Results(iResult).ImageGridAmp = Jf; .ImagingKernel = [];`
   reset `TessInfo.DataMinMax = []`, then `panel_surface('UpdateSurfaceData', hFig)` +
   `('UpdateSurfaceColormap', hFig)` + `figure_3d('SetShowSourceVectors', hFig, iTess, 1)`.

`Restore(panelName)` (filter-off and teardown): write the backed-up
`OrigImageGridAmp`/`OrigImagingKernel` back, reset `DataMinMax`, refresh.

Because the whole series is swapped, the normal time pipeline shows the filtered field
at every frame with no per-frame hook and no flicker. Recompute happens only on
filter-on and on slider-settle (a few seconds for a long series; instant scrubbing
after).

### Save filtered file

`SaveFiltered(panelName)`: build a results struct = a copy of the original results
metadata (`SurfaceFile, Time, nComponents=3, HeadModelType='surface', GoodChannel,
ColormapType, Function, History, …`) with `ImageGridAmp = Jf` (the current filtered
series) and `ImagingKernel = []`, `Comment = '<orig> | dirac filter(<info>)'`. Save into
the SAME study as the source via the standard results-save path (`db_add` / `bst_save`
+ `db_reload_studies`), reusing `process_dirac_filter`'s comment/labeling.

## Components / files

**Create:**
- `toolbox/gui/bst_eigfilter_panel.m` — shared kernel-section UI helper (Create/Read).
- `toolbox/gui/panel_spatial_filter.m` — the Spatial Filter panel + Start/Apply/Restore/
  SaveFiltered/teardown.
- `dev/tests/test_eigfilter_panel.m`, `dev/tests/test_spatial_filter.m`.

**Modify:**
- `toolbox/gui/panel_wavelet_designer.m` — build Section 2 via `bst_eigfilter_panel`
  (remove the now-shared subfunctions).
- `toolbox/gui/figure_3d.m` — add "Spatial filter (Dirac)" to the figure popup
  (`DisplayFigurePopup`), guarded to a 3-component surface Dirac source figure; calls
  `panel_spatial_filter('Start', hFig)`.

**Reuse unchanged:** `bst_dirac_eigenmodes_filter`, `bst_eigfilter_kernel`, `tess_eigen`,
`bst_memory` (`GetDataSetResult`, `GetResultsValues`), `panel_surface`
(`UpdateSurfaceData`/`UpdateSurfaceColormap`), `figure_3d('SetShowSourceVectors')`,
`process_dirac_filter` (labeling/args reference for save).

## Data flow

1. User: Dirac-dSPM → Display on cortex → right-click figure → "Spatial filter (Dirac)".
2. `Start(hFig)`: resolve iDS/iResult, check 3-comp, load eigenbasis+mass, back up
   original data, dock the panel.
3. Choose kernel + scale, tick **Filter on** → `Apply`: materialize, filter, swap,
   refresh. Scrub time → filtered frames shown instantly.
4. Adjust scale → `Apply` recomputes. Untick **Filter on** → `Restore`.
5. **Save filtered file** → new results node in the study.
6. Close panel / close figure → `Restore` (original data intact).

## Error handling

- **Not a 3-component surface Dirac source:** the popup item is hidden; `Start` aborts
  with a message if called anyway.
- **No Dirac eigenbasis / nxr unavailable:** `tess_eigen('Dirac')` find-or-creates it
  (progress bar); if it cannot be produced, `Start` aborts cleanly.
- **Vertex/size mismatch** (field rows ≠ `3·nVert` of the basis): abort with a clear
  message (mirrors `process_dirac_filter`).
- **Teardown always restores:** `Restore` runs from a single idempotent path used by
  filter-off, panel close, and the figure `CloseRequestFcn`, so the in-memory source is
  never left filtered. If the same results is shown in another figure it filters there
  too (shared in-memory data) and reverts together — acceptable and documented.

## Testing

- **Shared helper** (`test_eigfilter_panel`): `Create` builds a kernel section whose
  `Read` returns the selected kernel + params for a synthetic `Lambda`; switching kernel
  rebuilds the sliders; mode-index → param mapping matches the documented map
  (`t = 1/lambda_k`, `beta = lambda_k`, dog `t1 < t2`). Pure/headless where possible.
- **Wavelet designer regression:** the existing wavelet suite passes unchanged after the
  refactor (the kernel section still drives synthesis).
- **Spatial filter apply/restore** (`test_spatial_filter`, live figure): open a Dirac
  source figure, `Start`, `Apply` with a low-pass → the figure's `ImageGridAmp` equals
  `bst_dirac_eigenmodes_filter` of the original (assert on a sample frame); `Restore`
  → `ImageGridAmp` bit-identical to the backup; teardown restores and removes the panel.
- **Save filtered file:** the saved node loads to the filtered series, sits in the source
  study, and has `nComponents == 3`.

## Build order

1. `bst_eigfilter_panel` shared helper + its test.
2. Refactor `panel_wavelet_designer` Section 2 onto the helper; re-run the wavelet suite.
3. `panel_spatial_filter`: `Start` + backup + dock + the shared kernel section.
4. `Apply`/`Restore` (whole-series swap) + the on/off toggle + apply/restore test.
5. `Save filtered file` + its test.
6. `figure_3d` popup item (guarded) + end-to-end live test.
