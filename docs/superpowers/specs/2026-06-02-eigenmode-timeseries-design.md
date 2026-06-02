# Eigenmode Time Series Viewer — Design

**Date:** 2026-06-02
**Author:** Diellor Basha (with Claude)
**Status:** Approved design — ready for implementation plan
**Builds on:** `2026-06-01-eigenmode-viewer-lever-design.md` (the lever + paired-rank
foundation, merged), `2026-05-27-eigenmode-transform-spectrum-design.md`
(`bst_eigenmodes_transform`)

## Goal

Add an **eigenmode time series viewer**: one trace per eigenmode showing its sensor→mode
coefficient `θₖ(t)`, displayed in the standard butterfly/column time series layout used by
scout time series. The viewer tracks the `EigenModes` panel band — moving the slider or
changing the width adds/removes traces live — so it becomes the *temporal* counterpart to
`view_eigenmodes` (which shows the selection's *spatial* pattern `Φ·w` on the cortex).

## Conceptual model

The eigenmode selector is to **spatial frequency** what `panel_time` is to **time**: one
shared selection state, many clients. This spec adds a new client:

- **`view_eigenmodes`** (built) — synthesizes the selection `Φ·w` on the cortex (spatial).
- **Eigenmode time series** (this spec) — shows how the *real data's* coefficients on the
  selected modes evolve over time: `θₖ(t)` for each mode `k` in the band (temporal).

Where scout time series answers "what is this *region's* signal over time," the eigenmode
time series answers "what is this *spatial-frequency mode's* coefficient over time." Mode ≈
scout; the coefficient `θₖ(t)` ≈ the scout's extracted signal.

## Data: how θₖ(t) is computed

Sensor→eigenmode transform (unregularized), via the existing
`bst_eigenmodes_transform(Gain, Phi, Data)`:

```
L̃     = Gain * Phi            [nChan × nRawModes]   compressed lead field
Kernel = pinv(L̃)  (via SVD)   [nRawModes × nChan]
Theta  = Kernel * Data         [nRawModes × nTime]   per-raw-mode coefficients
```

- **`Gain`** — lead field from the study's **default head model** (auto-selected; no prompt).
- **`Phi`** — eigenvectors from the displayed surface's eigenmodes (`in_tess_eigenmodes`).
- **`Data`** — sensor recordings `[nChan × nTime]` from the associated data file, good
  channels for the head model's modality only.

`Theta` is computed **once** at launch over **all raw modes** and cached in the figure.
Band changes never recompute — they reselect rows (see live tracking).

## Traces: band → rows (paired rank, two traces per rank)

The panel band lives in **paired-rank** space (`[kLo, kHi]` over `1..K_paired`,
`K_paired = max(Eig.CompRank)`). Each paired rank `k` maps to its left and right raw
columns and yields **two traces**:

- For paired rank `k`: `iL = find(CompRank==k & Component==1)`,
  `iR = find(CompRank==k & Component==2)`. Emit `Theta(iL,:)` and `Theta(iR,:)`.
- **Labels:** `"Mode k L"`, `"Mode k R"`.
- **Colors:** keyed by hemisphere — left and right each get a fixed hue, so a butterfly
  plot reads as two interleaved hemisphere families. (Exact RGB chosen at implementation;
  left warm / right cool is the intent.)
- A rank missing one hemisphere (asymmetric mesh / backward-compat single component) simply
  emits the one trace it has.

**Raw coefficients, not weighted.** The band selects *which* modes appear; the panel's
window weights (Box / Taper / Gauss) are a *spatial-filter* concept and are **not**
multiplied into the time series. Every shown trace is the raw `θₖ(t)`.

## Display: reuse `view_timeseries_matrix`

Hand the assembled rows to the existing generic displayer so butterfly/column, legend,
time cursor, montage menus, and resize all come for free:

```matlab
view_timeseries_matrix( ...
    BaseFile, ...                 % the data file the figure is associated with
    {Theta_band}, ...             % F: cell with one [nTraces × nTime] matrix
    TimeVector, ...
    'eigenmodes', ...             % Modality tag
    {'Eigenmodes'}, ...           % AxesLabels
    LineLabels, ...               % {'Mode 1 L','Mode 1 R', ...}
    LineColors, ...               % [nTraces × 3] by hemisphere
    hFig);                        % reuse on refresh, [] on first open
```

Default `TsInfo.DisplayMode = 'butterfly'`; the existing Butterfly/Column menu toggle
works unchanged.

## Architecture (4 units)

1. **`view_eigenmodes_timeseries(DataFile, hFig)`** (new, `toolbox/gui/`). Entry point.
   Resolves context (channel data, default head model → `Gain`, surface eigenmodes → `Phi`,
   recordings → `Data`); runs `bst_eigenmodes_transform` once; caches `Theta`, `Eig`
   metadata (`Component`, `CompRank`), `TimeVector`, and `SurfaceFile` in figure appdata;
   calls `RefreshTraces` to do the first plot.

2. **Band → traces (`RefreshTraces`, private).** Reads the current panel band via
   `panel_eigenmodes('GetState')`, expands paired ranks to raw L/R columns, slices the
   cached `Theta`, builds labels + hemisphere colors, and calls `view_timeseries_matrix`
   (reusing `hFig`). Pure row-selection over cached data — cheap.

3. **Live tracking (`ModesChangedCallback`, private + dispatch wiring).** The panel's
   existing `bst_figures('FireModesChanged')` broadcast dispatches by client (the pattern
   `view_eigenmodes` already uses). Add an `EigenTimeSeries` figure-tag branch that calls
   `view_eigenmodes_timeseries('ModesChangedCallback', hFig)` → `RefreshTraces`. Moving the
   slider/width updates the traces in place without recompute.

4. **Launch (figure popup menu).** Add an "Eigenmode time series" item to the **3D source
   figure** popup (`figure_3d.m`) — that figure already has the cortical surface in context,
   so "which surface's eigenmodes" is unambiguous, and it carries the associated data file
   and head model. Enabled only when (a) the displayed surface has computed eigenmodes and
   (b) a head model is available in the study. Clicking resolves the `DataFile` from the
   figure and calls `view_eigenmodes_timeseries(DataFile, [])`. (A recordings/2D figure is
   *not* a launch site: it has no surface context.)

## Data flow (one selection change)

```
user edits band -> panel_eigenmodes SetBand/SetCurrentMode -> NotifyChanged
   -> RefreshControls (panel) + FireModesChanged
       -> per figure on the lever's surface, dispatch by client:
            EigenTimeSeries -> view_eigenmodes_timeseries('ModesChangedCallback', hFig)
                                  -> RefreshTraces (reselect cached Theta rows, replot)
            EigenView       -> view_eigenmodes('ModesChangedCallback', hFig)   (existing)
            source-map      -> panel_surface('UpdateSurfaceData', hFig)        (existing)
```

## Edge cases & errors

- **No eigenmodes on surface** — menu item disabled; if reached programmatically, error
  dialog "Compute eigenmodes first."
- **No head model in study** — menu item disabled; same guard in the entry point.
- **Empty band** (width 0 around a single center) — one rank → up to two traces; valid.
- **Channel/lead-field mismatch** — restrict `Data` and `Gain` to the head model's good
  channels for its modality (reuse the channel-selection the inverse path already does).
- **Backward-compat (no `Component`/`CompRank`)** — `in_tess_eigenmodes` already fills
  these; if a single component, each rank emits one trace.
- **Figure closed** — listener deregistered on close (mirror `view_eigenmodes` cleanup).

## Testing

Following the project's pure-function + smoke convention (`dev/tests/`):

- **`test_view_eigenmodes_timeseries_pure.m`** — given a synthetic `Theta`, `Component`,
  `CompRank`, and a band, assert the row-selection helper returns the correct raw-column
  indices, ordering (L then R per rank), labels, and trace count.
- **Reuse** `test_eigenmodes_transform_pure.m` for the `Theta` computation (already covered).
- **Smoke** (manual / e2e harness like `test_eigenmode_viewer_e2e.m`): open on a real data
  file, confirm a figure with the expected number of traces, toggle butterfly/column, move
  the band and confirm traces add/remove without recompute.

## Out of scope (YAGNI)

- Saving `Theta` as a DB result file (the `process_eigenmodes_transform` path already does
  that for batch use; this is a live, interactive viewer).
- Weighting traces by the panel window, combining L+R into one trace, or showing all modes
  independent of the band — all considered and rejected during brainstorming.
- Projecting source maps (`Φ'·M·u`) instead of the sensor transform — the sensor transform
  was chosen as the data source.

## File-level summary

| Unit | File | Change |
|------|------|--------|
| Entry point + refresh + callback | `toolbox/gui/view_eigenmodes_timeseries.m` | **new** |
| Multi-client dispatch | `bst_figures.m` (`FireModesChanged`) | add `EigenTimeSeries` branch |
| Launch menu | 3D source figure popup (`figure_3d.m`) | add menu item + enable guard |
| Pure test | `dev/tests/test_view_eigenmodes_timeseries_pure.m` | **new** |
