# Eigenmode Scale Lever — Design

**Date:** 2026-06-01
**Author:** Diellor Basha (with Claude)
**Status:** Approved design — ready for implementation plan

## Goal

A stateful, surface-scoped, live-broadcast **selection lever in eigenmode space** — the
spatial-frequency / *spatial scale* control for cortical maps. Changing the selection
live-coarsens/smooths a displayed source map on the cortex (band-limited reconstruction),
with no file written and no destructive edit.

## Motivation & conceptual frame

Brainstorm already has selection "levers" in three of four quadrants of a
{temporal, spatial} × {localized, spectral} matrix. The eigenmode lever fills the empty
cell. Every lever has the same shape — **a center + an extent (window)**:

| Axis | Lever | Center | Extent (window) | Status |
|---|---|---|---|---|
| Temporal — localized | time stepper (`panel_time`) | current time | time window | exists |
| Temporal — spectral | frequency smoother | center frequency | band | deferred (later spec) |
| Spatial — localized | scouts / brush (`panel_scout`, + future nxr heat-distance) | seed vertex | brush radius | exists / improvable |
| **Spatial — spectral** | **eigenmode lever (this spec)** | **center mode** | **mode range + window** | **building now** |

North star: at any instant the full state is the 4-tuple
`(time window @ t, spatial brush @ vertex, mode window @ k, freq band @ f)`, with a filter
layer on top for differential/spectral analysis. **This spec builds only the eigenmode
quadrant + its live coarsen/smooth consumer.** The temporal-frequency smoother, the
nxr heat-distance brush upgrade, and the differential/spectral filter library are roadmap,
not this build — but the lever is designed so its pattern (center + extent + window,
stateful, live-broadcast, surface-scoped) is clean enough to clone for the frequency lever.

**True template:** the eigenmode lever is the spatial-spectral sibling of `panel_freq`
(temporal frequency) for *selection*, and of `GlobalData.VisualizationFilters` /
`bst_memory:FilterLoadedData` for the *display-time transform* mechanism. It is **stateful
and live-broadcast** like `panel_time`/`panel_freq` (consumers track it continuously and
re-render), **not** act-on-demand like scouts. It is **surface-scoped** like the eigenbasis
itself (each surface has its own `K`, `λ`, `Φ`).

## Eigenmode ordering primer

Modes are ordered by eigenvalue `λ` (≈ spatial frequency): low index = smooth/global
patterns, high index = fine/local detail. A low-pass band `[1, k_hi]` therefore coarsens
(smooths); widening the band toward `K` sharpens. Edge truncation of a hard band causes
Gibbs-like ringing on the surface — hence the tapered/gain window shapes.

## Architecture & components

Four units, each with one responsibility:

1. **State** (`GlobalData.UserModes`) — the lever's source of truth. Surface-scoped,
   refreshed when the active surface changes.
2. **`panel_eigenmodes`** (new GUI panel) — the control surface (center slider + band +
   window-shape selector + active toggle). Reads/writes State, fires one broadcast. No math,
   no rendering.
3. **Live-filter consumer** — on broadcast, applies a display-time band-limited
   reconstruction to each open source figure on the matching surface and repaints via the
   existing `figure_3d('UpdateSurfaceColor')` path.
4. **Core reuse (no new math)** — `bst_eigenmodes_project` (projection/reconstruction),
   `bst_eigenmodes_filter_gain` (gain/window shapes), `in_tess_eigenmodes` (`Φ`, `λ`, mass
   type). The lever contributes *state + wiring*, not algorithms.

Broadcast/subscribe mirrors `bst_figures('FireCurrentTimeChanged')`; the new event is
`bst_figures('FireModesChanged')` so future consumers (deferred TS viewer, single-mode 3D
browser) subscribe identically.

## State structure

```matlab
GlobalData.UserModes = struct(...
    'SurfaceFile',  '',       ... % owning surface (eigenmodes are per-surface)
    'nModes',       0,        ... % K available on this surface
    'iCurrentMode', 1,        ... % CENTER (stepper cursor; band center — COUPLED)
    'Weights',      [],       ... % [1 x K] window/taper w(k) in [0,1] — CANONICAL state
    'WindowShape',  'single', ... % 'single'|'box'|'tapered'|'gain' — how UI authors Weights
    'Band',         [1 1],    ... % [k_lo k_hi] extent for box/tapered shapes
    'isActive',     0);           % filtering on/off (off => raw map shown, lever inert)
```

`Weights` is canonical (consumers read only this). `WindowShape`/`Band`/`iCurrentMode` are
how the panel authors `Weights`:

- `single`  → `w = e_k` at `iCurrentMode` (zero-width band at center)
- `box`     → `w = 1` on `Band`, else 0
- `tapered` → Tukey/Hann taper across `Band` edges (anti-ringing); `1` in the interior
- `gain`    → smooth curve from `bst_eigenmodes_filter_gain` (e.g. heat-kernel rolloff)

**Center↔band coupling: COUPLED.** `iCurrentMode` is the band center; moving the center
slides the whole window. One "scale" knob. `single` is the zero-width degenerate case.
Decoupling (independent inspection cursor) is deferred until the single-mode 3D browser needs it.

## Panel API (`panel_eigenmodes`)

Thin verbs mirroring `panel_freq`/`panel_time`; each setter recomputes `Weights`, stores
it, then fires `bst_figures('FireModesChanged')`:

```matlab
panel_eigenmodes('SetCurrentMode', k)        % move center (coupled: slides band); step with arrows
panel_eigenmodes('SetBand', kLo, kHi)        % set extent (recenters iCurrentMode to band center)
panel_eigenmodes('SetWindowShape', shape)    % 'single'|'box'|'tapered'|'gain'
panel_eigenmodes('SetActive', 0|1)           % bypass vs apply (instant, non-destructive)
W = panel_eigenmodes('GetWeights')           % canonical read for consumers
panel_eigenmodes('UpdatePanel', SurfaceFile) % repopulate on active-surface change
```

## Data flow (one lever change → repaint)

```
user drags band -> panel_eigenmodes('SetBand', ...)
   -> recompute Weights -> store in GlobalData.UserModes -> FireModesChanged
       -> consumer, for each open source figure on this SurfaceFile:
            u   = current-time source column           (bst_memory, already loaded)
            c   = bst_eigenmodes_project(u, Phi, M)     [K x 1]   (= Phi' * M * u)
            uF  = Phi * (Weights' .* c)                 [nV x 1]  band-limited reconstruction
            push uF as displayed values -> figure_3d('UpdateSurfaceColor')
```

- For display we transform **only the current-time column**, so cost is two mat-vec products
  with `Phi` (~15k×200) ≈ a few million flops — sub-millisecond; smooth to drag. No debounce
  required; debounce-on-release is the safety valve only if a pathological surface appears.
- Composes with the time stepper: each new time column flows through the same `c -> uF`
  transform; stepping time and dragging the lever both feed `UpdateSurfaceColor`.
- `isActive = 0` ⇒ consumer is a no-op and the raw map shows. Toggling is instant.

## Tree / GUI wiring (entry point)

- **Panel lifecycle:** `panel_eigenmodes` shows/enables when the active figure is a source
  map whose `SurfaceFile` has stored eigenmodes (`in_tess_eigenmodes` → `isComputed`);
  hidden/greyed otherwise. Same lifecycle as `panel_freq` appearing only for time-freq data.
- **Activation:** a "Spatial scale (eigenmodes)" toggle in the source figure context calls
  `SetActive(1)` (reveal panel + start live filtering) / `SetActive(0)` (restore raw map).
- **No new anatomy-node menu** — the lever acts on a *displayed source map*, not on the
  surface file itself.

## Error handling

- **No eigenmodes on surface** → panel disabled; toggle shows the standard "Run Compute
  eigenmodes first" message (reuse `view_eigenmodes` wording).
- **Band exceeds stored K** → clamp `Band`/`iCurrentMode` to `[1, K]`.
- **Active-surface change / figure close** → `UpdatePanel` repopulates from the new surface
  (or clears state if none qualifies); consumer skips figures whose `SurfaceFile` ≠ state's.
- **Mass matrix** for projection comes from the eigenmodes' stored `MassType` via the core;
  if absent, fall back to `bst_eigenmodes_project`'s default — never a silent wrong-metric
  reconstruction.
- **Non-destructive guarantee:** consumer writes only *displayed* values; the stored results
  file is never modified. `SetActive(0)` always restores the original.

## Testing

- **Pure-math (no GUI):** weight-shape builders — `single`→delta; `box`→boxcar on `Band`;
  `tapered`→edge taper (values ≤1, monotone shoulders, interior = 1); `gain`→matches
  `bst_eigenmodes_filter_gain`. Round-trip identity: full band (`w = 1`) ⇒ `uF == u` to
  machine precision.
- **State logic:** `SetBand`/`SetCurrentMode` clamp to `[1,K]`; coupled center slides band;
  `isActive = 0` ⇒ consumer no-op.
- **Integration (MATLAB MCP, headless):** load a real surface with eigenmodes + a source map,
  drive `SetBand`, assert displayed values equal the analytic band-limited reconstruction and
  that a low band ⇒ smoother map (lower surface gradient / high-mode energy ≈ 0).
- **Performance smoke:** one lever change recomputes a single column + repaints in a few ms
  for K ≈ 200.

## Out of scope (roadmap)

- Temporal-frequency smoother (`panel_freq` + `bst_eigenmodes_wavelet` time coarsening).
- nxr heat-distance scout brushing.
- `view_eigenmode_timeseries` (on-the-fly sensor→mode time series viewer) and the
  single-mode 3D browser — both future consumers of `FireModesChanged`.
- Differential/spectral analysis filter library.
- Persisting lever state to disk (state is session-scoped, like `panel_time`/`panel_freq`).
```
