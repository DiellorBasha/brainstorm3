# Dynamics "Focus" State — Bidirectional Time & Frequency

**Date:** 2026-06-25
**Status:** Design (approved, pre-plan)
**Component:** `panel_bst_dynamics` / `figure_timeseries` / `figure_spectrum` / Dynamics Atoms system

## 1. Motivation & Scope

The Dynamics panel implements a three-state workflow — **navigate → focus → save** —
centered on *Atoms* (joint spatial-temporal-frequency-scale localization markers).
Navigation already piggybacks on Brainstorm's global time stepper. The **save**
state commits atoms. The intermediate **focus** state — the ephemeral "I am
currently looking at *this* time window and *this* frequency band" — has no
principled UI feedback yet.

This spec implements the **focus** state for the **time** and **frequency** axes
only, by binding the panel's active selection to the native, user-editable
*selection box* on the corresponding data-measurement figure, in **two-way sync**.

**Explicitly deferred to their own specs:**
- Source/scale focus (the geodesic Region tool is the existing source analogue).
- Downstream differential analyses (phase via `process_evt_refphase`,
  divergence / curl / potential / stream critical points).

The focus state is the skeleton an *Atom* is built from: later it adapts to
`view_eigenspectrum` and to source-map filters.

## 2. Core Concept — Focus ↔ Native Selection Box

Each axis binds the panel's active selection to the native selection idiom of one
measurement figure. Both idioms already exist, are draggable, and render a
"Selection: […]" label for free.

| Axis | Figure (type) | Native selection (existing) | Units | Panel fields |
|------|---------------|-----------------------------|-------|--------------|
| Time | `figure_timeseries` (recording, `DataTimeSeries`) | `GraphSelection = [tStart, tEnd]` | seconds | active atom-group window |
| Freq | `figure_spectrum` (PSD butterfly, `Spectrum`) | `GraphSelection = [fLo, fHi]` | Hz | `jFreqC` / `jFreqW` |

**Decision:** reuse the *native* selection rather than painting a static overlay
patch. A static patch cannot be dragged, which would break the bidirectional
"adjust the overlay manually" requirement. The native selection already supports
mouse drag, redraw, linked propagation across figures, and the text label.

### Verified mechanics (from source investigation)

**Time series** (`toolbox/gui/figure_timeseries.m`):
- `GraphSelection` appdata, single `1x2` `[min max]` interval (single window only).
- Programmatic set: `figure_timeseries('SetTimeSelectionManual', hFig, [tStart tEnd])`
  (snaps to samples) → `SetTimeSelectionLinked` (propagates to all
  `DataTimeSeries`/`ResultsTimeSeries` figures) → `DrawTimeSelection` (patch tag
  `TimeSelectionPatch`, label tag `TextTimeSel`: "Selection: […] Duration: […]").

**Spectrum / PSD** (`toolbox/gui/figure_spectrum.m`):
- `view_spectrum(TimefreqFile, 'Spectrum')` opens PSD in butterfly mode, X axis =
  **frequency in real Hz**, axes tag `AxesGraph`, figure `FigureId.Type='Spectrum'`.
- `GraphSelection = [fLo, fHi]` **in Hz**; `DrawSelection(hFig)` draws the patch
  (tag `SelectionPatch`) and a "Selection: [%.2f Hz - %.2f Hz]" label.
- Find an open one: `bst_figures('GetFiguresByType', 'Spectrum')`.

**Dynamics panel** (`toolbox/gui/panel_bst_dynamics.m`):
- Band presets in `i_bands()`: δ`[2 4]`, θ`[4 8]`, α`[8 13]`, β`[13 30]`, γ`[30 60]`.
- Preset selection: `OnFreqPreset` → `i_freq_preset` fills `jFreqC`/`jFreqW` →
  `OnAxisChange('freq')` → `i_drive('freq', loc)` → `panel_filter('SetFilters', …)`
  (already applies a live bandpass to the time-series figure).
- Atoms navigation: `OccurSel_Callback` already highlights the marker and jumps
  the time cursor (`panel_time('SetCurrentTime', …)`).
- Detect: `OnDetect` runs `process_evt_refphase('Compute', …)` and stages
  **extended** event groups (band-strong windows) plus point phase markers.

## 3. Data Flow (bidirectional, both axes)

### Panel → View (drive)

**Frequency** (on band preset / `i_drive('freq')`):
1. Keep the existing `panel_filter` bandpass on the time-series view (unchanged).
2. **Ensure a PSD spectrum figure exists** (§4.4 helper):
   find open `Spectrum` figure for this recording → else open precomputed PSD via
   `view_spectrum` → else auto-compute `process_psd` (averaged, **magnitude**,
   freq range covering 0–60 Hz) then open.
3. Set that figure's `GraphSelection = [fLo, fHi]` and redraw (`DrawSelection`).

**Time** (on detect / Atoms-list navigation):
- On detect: set the recording figure's Time Selection to the **first** atom
  group's `[onset, offset]` via `SetTimeSelectionManual`.
- On Atoms-list selection (`OccurSel_Callback`): re-point the single box to the
  selected group's `[onset, offset]` (in addition to the existing time jump).
- Only **extended** atom groups map to a window; point phase markers
  (peak/trough/rising/falling) are time-delta functions and do **not** set a box.

### View → Panel (sync back)

A hook fires when the user edits a selection on either figure:
- **freq drag** → write `jFreqC=(lo+hi)/2`, `jFreqW=(hi-lo)/2`; set band combo to
  the matching preset, else `custom`; re-apply the time-series bandpass —
  **without** re-driving the figure selection (loop guard, §4.3).
- **time drag** → update the *staged/active* focus window the panel reports. In
  the pre-save focus state this adjusts the currently-selected (staged) atom
  group's `[onset, offset]`; it does not retroactively mutate already-saved atoms
  (post-save editing is out of scope for this spec).

## 4. Components (unit boundaries)

### 4.1 `panel_bst_dynamics('NotifySelection', hFig, axis, range)`
Single new panel entry point. Dispatcher the figures call when a native selection
changes. `axis ∈ {'time','freq'}`, `range = [lo hi]`. Updates the panel fields per
§3 (View → Panel). No-op if no Dynamics target owns `hFig`, or if the
re-entrancy guard is set.

### 4.2 Sync hooks in the figures
Minimal additions:
- `figure_timeseries.SetTimeSelectionLinked` (or `DrawTimeSelection`): after
  updating `GraphSelection`, call `panel_bst_dynamics('NotifySelection', hFig,
  'time', sort(GraphSelection))` **iff** a Dynamics target owns the figure.
- `figure_spectrum.DrawSelection`: same call with `'freq'`.
Guarded so non-Dynamics usage is unaffected (cheap appdata check first).

### 4.3 Re-entrancy guard
Module-global flag (e.g. persistent `bDrivingSelection` in the panel) set while
the panel programmatically drives a figure selection; checked at the top of
`NotifySelection` to suppress the echo and prevent infinite ping-pong.

### 4.4 PSD lifecycle helper
`i_ensure_psd(DataFile)` (panel-local): find-or-open-or-compute the PSD spectrum
figure for the current recording and return a tracked handle. Encodes the
"always ensure a TF view exists" rule. Auto-compute uses `process_psd` with
averaged output, `Measure='magnitude'`, window/overlap defaults, frequency range
covering 0–60 Hz.

## 5. Persistence Trade-off (decided)

Native selections clear on a plain (no-drag) click. **Decision: accept native
semantics** — clicking to clear = defocus — rather than forcing a separate
persistent overlay patch. The user can re-assert focus by re-selecting a preset
or an atom group, and the panel re-asserts the selection when an atom group is
(re)selected. This keeps both axes fully editable and the code minimal.

## 6. Validation

Brainstorm has no unit tests; validate manually (per the run-tutorial idiom) on
the **alpha example block** (`S01_AEF_01_notch`, band 7–13 Hz, right
parieto-occipital ~10.55 Hz burst at ~22.6 s):

1. Open the recording + Dynamics panel.
2. Pick the **alpha** preset → time-series bandpass applies **and** the PSD
   spectrum opens (auto-compute if absent) with a band strip at 8–13 Hz.
3. **Drag the PSD band strip** → panel `jFreqC`/`jFreqW` update, band combo flips
   to `custom`, time-series bandpass follows. No feedback loop.
4. **Detect** → Time Selection box lands on the first extended atom group.
5. **Navigate** the Atoms list → the single box moves to each group's window.
6. **Drag the Time Selection box** → panel's active window updates. No loop.

## 7. Risks / Notes

- **Loop safety** is the main risk; the re-entrancy guard (§4.3) is mandatory and
  must be verified in both directions.
- The time selection is **single-window** by design; multi-burst navigation is
  served by the Atoms list, not by drawing N boxes.
- `process_psd` auto-compute on continuous/raw data has a cost; it runs only when
  no PSD figure or precomputed PSD file is available, and only on band-preset
  selection (not on every panel interaction).
- Hooks in `figure_timeseries`/`figure_spectrum` must be inert for all non-Dynamics
  users (guard on a Dynamics-target appdata check before any work).

## 8. Validation (2026-06-25)

Implemented via subagent-driven development (7 tasks, per-task spec+quality review).
Controller end-to-end live validation on the raw `S01_AEF_01_notch` alpha recording
(real GUI Dynamics session, PSD auto-computed) — **all paths pass, zero runtime errors**:

- **Freq drive:** alpha preset → PSD spectrum opens, band strip exactly `[8 13]`,
  X-axis `[0 60]`, `curBandName=alpha`.
- **Freq sync-back:** simulated drag to `[18 24]` → panel fields center/half = 21/3,
  combo flips to `custom`; the `i_driving` guard suppresses the echo (no loop).
- **Time drive:** Detect found 62 alpha windows; the Time Selection box lands on the
  first window `[0.453 0.890]`.
- **Staged navigation:** the Atoms tree shows `detwinroot` + 62 per-window `detwin`
  leaves; selecting a later window re-points the box (tree branch `i_jump`s first to
  load the raw page, then focuses).
- **Time sync-back:** dragging the box over a staged window rewrites that window's
  `[onset offset]`; `focusTime` is recorded; the recording's `isModified` flag is
  preserved (render-only, never dirties the file); the guard suppresses the echo.
- **Inertness:** with no Dynamics session, the figure hooks no-op cleanly.

**Known raw-recording nuance (minor):** `i_focus_time` snaps to the currently-loaded
raw page; focusing a window in another page requires loading that page first. The tree
selection branches do this (`i_jump` precedes `i_focus_time`); `OnDetect`'s first-window
focus does not, which is harmless when the first window is on the initial page. Non-raw
(imported/averaged) recordings are single-page and unaffected.
