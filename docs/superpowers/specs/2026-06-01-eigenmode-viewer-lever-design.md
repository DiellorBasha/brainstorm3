# Eigenmode Viewer ↔ Lever — Design

**Date:** 2026-06-01
**Author:** Diellor Basha (with Claude)
**Status:** Approved design — ready for implementation plan
**Builds on:** `2026-06-01-eigenmode-scale-lever-design.md` (the lever foundation, merged)

## Goal

Make `view_eigenmodes` a client of the eigenmode selection lever (the `EigenModes`
panel) instead of hijacking the time stepper, and reindex the lever to **paired rank**
so selections are hemisphere-symmetric. The viewer becomes the lever's primary,
immediate-feedback display: select a single mode or superpose a range and see `Φ·w`
on the cortex live.

## Conceptual model

The eigenmode selector is to **spatial frequency** what `panel_time` is to **time**: one
shared selection state, many clients, each reacting to a change broadcast in its own way.

- **`view_eigenmodes`** (this spec) — the lever's *primary* partner: *synthesizes* the
  selection as `Φ·w` on the cortex (immediate visual feedback of the selection itself).
- **Cortical activations / source maps** (already built) — *filters* real source data
  (`Φ·diag(w)·Φᵀ·M·u`); needs an input (a vertex-delta impulse or a source activation).
- **`view_filter_design`** (future) — would render the transfer function / kernel.

The selection *is* the picture in the viewer; filtering is a separate, data-dependent
*application* of the same selection. The broadcast already exists
(`bst_figures('FireModesChanged')`); this spec makes it dispatch by client (the shape of
`FireCurrentTimeChanged`) so future clients slot in by adding a handler.

## Why paired rank (the asymmetry fix)

`tess_eigenmodes` stores modes as concatenated per-component blocks: raw columns
`1..K_L` are the left hemisphere's modes (rank 1…K_L), then `K_L+1..K_L+K_R` the right's;
each mode has disjoint support across components. The merged lever indexes these **raw
columns**, so a band `[1,15]` selects the left hemisphere's first 15 modes **and none of
the right's** — smoothing/showing one hemisphere and blanking the other. The natural,
correct unit is **paired rank**: rank k = both hemispheres' k-th mode (same spatial
frequency on each). This spec moves the lever to paired-rank space, fixing the asymmetry
for *every* client.

## Architecture (4 units)

1. **Paired-rank reindex (`panel_eigenmodes`).** State works in paired-rank space:
   `nModes = max(Eig.CompRank) = K_paired`; band/center/weights over `[1, K_paired]`. One
   expansion serves every client: `w_raw = W(Eig.CompRank)` (length = raw `Eig.nModes`).
2. **Viewer as synthesis client (`view_eigenmodes`).** Displays `PairedGrid·W` (where
   `PairedGrid` is `BuildPairedGrid`'s `[nV × K_paired]`) as a **single-frame** transient
   source result (no `Time = 1:K` hijack). Tags its figure as an eigenmode view; caches
   `PairedGrid` at open; owns a `ModesChangedCallback` that re-synthesizes on selection
   change. Arrow keys route to `panel_eigenmodes('SetCurrentMode', …)`.
3. **Multi-client dispatch (`bst_figures('FireModesChanged')`).** Iterates figures on the
   lever's surface and dispatches by client (mirrors `FireCurrentTimeChanged`):
   eigenmode-view → `view_eigenmodes('ModesChangedCallback', hFig)`; source-map →
   `panel_surface('UpdateSurfaceData', hFig)` (the existing filter path).
4. **Lifecycle & Active-gating (`panel_eigenmodes`).** `UpdatePanel` runs on figure
   open/close/current-3D-change; the panel goes inert when no eligible view is in front;
   the filtering **Active** toggle is enabled only in source-map context.

## Paired-rank mechanics

- `K_paired = max(Eig.CompRank)`. Backward-compat: if `CompRank` absent/empty, treat each
  column as its own rank (`CompRank = (1:nModes)'`, so `K_paired = nModes`).
- **Expansion (one place, used by all clients):** `w_raw = W(CompRank(:))` — each raw
  column takes its rank's weight. `Eig.CompRank` is already cached with `Eig`.
- **Filter client (`ApplyToColumn`):** transfer fn changes from `@(l) W(:)` to
  `@(l) W(Eig.CompRank(:))`; the `numel(W)` guard changes from `Eig.nModes` to
  `max(Eig.CompRank)` (= `K_paired`). No other change — fixes the asymmetry.
- **Synth client (`view_eigenmodes`):** display column `= PairedGrid * W(:)`. Single mode
  (`single` shape, center k) → the rank-k paired column; band → the weighted sum.
- **`UpdatePanel`:** `K = max(Eig.CompRank)`; slider maxima = `K_paired`; state resets when
  `SurfaceFile` or `K_paired` changes.

## Data flow (one selection change)

```
user edits selection -> SetBand/SetCurrentMode -> NotifyChanged
   -> RefreshControls (panel) + FireModesChanged
       -> per figure on the lever's surface, dispatch by client:
            eigenmode-view fig -> view_eigenmodes('ModesChangedCallback', hFig)
                  -> col = PairedGrid * W(:) -> repaint (colormap UI intact)
            source-map fig     -> panel_surface('UpdateSurfaceData', hFig)
                  -> ApplyToColumn filters with W(CompRank)
```

The viewer's figure carries `setappdata(hFig, 'EigenView', struct(...))` holding its cached
`PairedGrid`, the per-rank eigenvalues, and the `SurfaceFile`. `ModesChangedCallback`
recomputes the displayed column from that cache and repaints via the standard surface
colormap path (so thresholding / colormap UI keep working).

## Panel lifecycle & Active-gating

The panel stays registered (present, like Surface/Scout) but is **inert when no eligible
view is in front**. `UpdatePanel(hFig)` is driven on figure open / close / current-3D-change
(the events `panel_surface('UpdatePanel')` already responds to) and classifies the front
figure:

| Front figure | Selection controls (band/center/shape) | Filtering **Active** toggle |
|---|---|---|
| Eigenmode view | enabled (always live → drives `Φ·w`) | hidden/disabled — N/A |
| Source map with eigenmodes | enabled | enabled (toggles filtering) |
| Anything else / none | disabled; readout "no eigenmode view"; `isActive→0` | disabled |

Eligibility: a `3DViz` figure that is either tagged `EigenView`, or carries a source
overlay (`DataSource.FileName`) on a surface with eigenmodes.

## Error handling

- **No eligible figure** → disabled panel, forced `isActive = 0` (no stale filtering).
- **Surface without eigenmodes** → disabled (existing).
- **Viewer closed** → `UpdatePanel` on close drives the panel inert; the viewer's transient
  result is cleaned up by the existing `CleanupResult` (`DeleteFcn`).
- **Source-map ↔ eigenmode-view switch** (possibly different surfaces) → `UpdatePanel`
  repopulates; state resets when `SurfaceFile` or `K_paired` changes.
- **Single-component / legacy eigenmodes** (no `CompRank`) → `K_paired = nModes`; the
  expansion is the identity, so all clients behave as before.

## Testing

- **Pure (headless):** `w_raw = W(CompRank)` expansion and `K_paired = max(CompRank)` on a
  synthetic `CompRank` (incl. the no-`CompRank` identity fallback).
- **Asymmetry-fix (key, headless):** build a **2-component** mesh (two disjoint spheres,
  concatenated with offset faces so `tess_eigenmodes` returns two components). Assert a low
  paired-band reconstruction has energy on **both** components (the old raw-column band put
  energy on only one). Assert symmetry: components' kept-energy ratio ≈ matched.
- **Synthesis (headless):** `PairedGrid·W` equals the rank-k paired column for `single`, and
  the weighted sum for a band.
- **Lifecycle (headless-limited):** `UpdatePanel` disables controls + forces `isActive=0`
  when the front figure is ineligible; the Active toggle is gated by context.
- **e2e (live):** open View Eigenmodes → step modes with arrows *and* the panel band;
  superpose a range; confirm the cortex updates and the filtering Active toggle is inactive
  in viewer context; close the view → panel goes inert.

## Out of scope (roadmap)

- `view_filter_design` client (transfer-function / kernel view) — the dispatch is designed
  to accept it, but it is not built here.
- Filtering UX that supplies an impulse/source activation as a standalone evaluation input.
- The cache-invalidation hook and other follow-ups already tracked in the lever spec.

## Tracked follow-ups (from final review — non-blocking)

Implemented and verified (e2e ran live on a single-component cortex; 2-component
asymmetry fix unit-covered by `test_eigenmode_lever_paired`):

- **Live 2-component coverage.** Run `test_eigenmode_viewer_e2e` on a genuine bilateral
  cortex to exercise paired synthesis + repaint on two components live (the headless
  paired test covers the math; the live viewer path ran only single-component).
- **Pre-tag broadcasts.** `ViewFigure` open fires `FireModesChanged` from
  `SetWindowShape`/`SetCurrentMode` before the figure is tagged `EigenView` — harmless
  (a redundant `UpdateSurfaceData` on any same-surface source map), optionally suppressed
  by seeding state without broadcast.
