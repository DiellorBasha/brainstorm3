# Dynamics panel — substrate cleanup (sub-project A) — design

**Date:** 2026-06-30
**Status:** design (approved; to be implemented in a fresh session)
**Scope:** `toolbox/gui/panel_bst_dynamics.m` only (behavior-preserving refactor)
**Part of:** the larger Dynamics-portal program. Decomposed into four sub-projects:
- **A — Substrate cleanup** (this doc): prune dead/retired code + consolidate the realise/reduce duplication.
- **B — Frame core:** turn the atom list into a real frame (live frame bounds via `bst_eigenwavelet('Bounds')`, "Design tight frame", cached source projection `C` for instant re-synthesis).
- **C — Analyze / reconstruct:** whole-bank frame transform → scalogram / per-member energy + synthesis residual, `JTVAnalysis → JTVAtoms` wiring, likely a `process_` for the whole series.
- **D — Preview completion:** operator/source compatibility gate + Dirac-eigenbasis sensor forward (Task 3).

Order: A → B → C → D. This doc covers **A only**.

---

## 1. Motivation

`panel_bst_dynamics.m` is 1861 lines spanning **three generations** of UI bolted onto one file:

1. **Truly unreachable** code that references Swing controls `CreatePanel` no longer builds (the old tree + 4-axis "navigator" + PSD/freq-overlay machinery). It cannot execute, but it obscures the live code and references functions that look alive.
2. A **reachable older mission** — band-power detection (`process_evt_refphase`), Helmholtz differential maps (Measure), and recording/extrema capture (Record) — that embodies a *different* notion of "atom" (a thresholded marker / field extremum) than the current filterbank model (an atom = a localized eigenfilter generator).
3. The **current filterbank** portal (Create atom → edit params → Localize → Apply impulse/Preview), plus the realise/reduce helpers — which are **duplicated three ways**.

Goal: leave a clean, legible filterbank portal (plus the differential Measure overlay, which is kept) as the substrate for sub-project B. No behavior change to anything retained.

## 2. Decisions (from brainstorming)

- **Stratum 1 (unreachable):** delete.
- **Stratum 2 — refphase detection + record/capture:** delete. **Differential Measure (Div/Curl/Potential/Stream):** keep.
- **Phase-display UI:** keep as latent capability (Show-phases submenu, Show-all toggle, the `showPhase` argument threaded through `view_dynamics('Redraw')`). It currently has **no producer** (refphase, its only previous producer, is removed); this is an explicit choice to leave the door open for a later producer (e.g. the frame / `JTVAtoms` work in sub-project C).
- **Consolidation depth (#6):** Approach 1 — collapse the magnitude-reduction duplication into one local helper `i_paintable_scalar`; keep realise as two entry points (impulse vs apply) sharing the reduction tail; single file; no extraction to a new module (the reduce/normalize logic is single-caller — `view_atom_designer`'s realiser does not reduce/normalize — so module-locality says keep it in the panel).
- **Verification:** MCP-driven launch + click-through, preceded by a `checkcode` pre-pass.

## 3. Invariant

The following must behave **identically** after the refactor:

- Filterbank flow: Create atom → select in list → edit a kernel slider → Localize (click cortex to re-seed) → Apply toggle (OFF = impulse-response Design paint; ON = filtered-source Preview paint).
- **Measure** differential overlay (None / Divergence / Curl / Potential / Stream) via `view_dynamics`'s per-frame Helmholtz decomposition + `PickScalar`.
- Phase-display controls (inert but present).
- File menu (Open / Save / Save as), Atoms menu group management (Add / Rename / Delete / Set color / Set operator / Sort), Threshold, Save filterbank, Close-session teardown.

No edits to `bst_dynamics.m`, `bst_eigen*.m`, `view_dynamics.m`, `view_atom_designer.m`, `bst_geodesic_tool.m`, or `panel_eigenfilter_design.m`.

## 4. Deletion inventory (verified against the current file)

### 4.1 Stratum 1 — unreachable

**Tree / occurrence list** (reference `ctrl.jTree`, `ctrl.jListOccur`, never built):
- `BuildTree` — reduce to a thin wrapper `function BuildTree(); UpdateAtomList(); end` (2 callers: `SetTarget`, `i_apply`); delete its ~55-line unreachable body.
- `TreeSel_Callback`, `OccurSel_Callback`, `i_window_atoms`, `i_group_atoms`.

**4-axis navigator** (reference `ctrl.jFreqC/jFreqW/jFreqBand`, `jSrcC/jSrcW`, `jScaleC/jScaleW`, never built):
- `i_axis_block`, `i_read_block`, `OnAxisChange`, `OnFreqPreset`, `i_drive`, `i_freq_preset`, `i_freq_name`, `i_fill_block`, `OnLoadAtom`, `i_band_match`, `i_bands`.

**Frequency / PSD overlay + time/freq sync-back** (Navigator-era):
- `NotifySelection`, `i_rec_figure`, `i_owns_rec`, `i_owns_spec`, `i_ensure_psd`, `i_fix_spec_xlim`, `i_find_psd_file`, `i_compute_psd`, `i_freq_overlay`, `i_freq_overlay_clear`, `i_sync_freq`, `i_focus_time`, `i_sync_time`, `i_driving`.

### 4.2 Stratum 2a — refphase detection

- Functions: `OnDetect`, `OnSaveDetection`, `OnClearDetection`, `i_detect_events`, `i_remove_band`, `i_load_meg`, `i_has_staged_detection`, `OnSave`.
- Toolbar buttons: **Detect windows** and **Clear** (the two `gui_component('ToolbarButton'…@OnDetect / @OnClearDetection)` lines and their separators as appropriate).

### 4.3 Stratum 2b — record / capture

- Functions: `OnRecord`, `OnSaveCursor`, `OnCaptureRegion`, `OnRegionTool`, `ctrl_region_state`.
- Atoms-menu item: **"Record at cursor"** (`@OnRecord`).
- Now-orphaned helpers (verified single-use by the above): `i_peaks`, `i_find_group`, `i_op_color`, `i_disp_band`, `i_scale_name`, `i_first_results`.

### 4.4 State-struct prune (`SetTarget`)

`SetTarget` currently seeds a struct with many fields used only by deleted paths. Drop:
`nav`, `occMap`, `nodeList`, `nodeInfo`, `Lambda`, `hSpec`, `focusTime`, `detSel`, `curGroup`, `curBand`, `curBandName`, `curScale`.

Retain: `hFig`, `T`, `file`, `curAtom`, `atomSeed`, `atomThreshold`, `showPhase`, `curOp` (read by the retained Measure path), plus `atomAx`/`atomBounds` populated lazily by `i_atom_ensure_axes`.

> Verification cross-check: after the prune, `grep` each dropped field name in the file — there must be zero remaining reads.

## 5. Retained code (explicit keep-list)

- **Filterbank:** `CreatePanel` (minus deleted toolbar/menu items), `OnCreateAtom`, `UpdateAtomList`, `SetSelectedAtom`, `AtomsListValueChanged_Callback`, `i_select_atom_load`, `i_atom_writeback`, `OnSaveFilterbank`, `OnThresholdMenu`, `OnSetOperator`, `i_select_op_radio`, `OnKernelChange`, `OnParamSettle`, `OnLocalize`, `i_atom_ensure_axes`, `i_atom_axes`, `i_atom_surface`, `i_atom_default_bounds`, `i_atom_current_kernel`, `i_default_atom`, `i_atom_detail`, `i_launch_operator`, `i_atom_realise`, `i_seed_block`, `i_atom_normalize`, `i_atom_preview`, `i_atom_preview_impulse`, `i_atom_op`, `i_atom_bounds`, `i_atom_filter_field`, `i_cursor_window`/`_core`/`_test`, `OnApply`, `i_atom_apply`, `i_overlay_nv`.
- **Localize feedback:** `SyncSource` (invoked by `bst_geodesic_tool.m:168`).
- **Measure / differential:** `OnMeasureMenu`, `OnMeasurement`, the **Measure** toolbar button.
- **Phase display (latent):** `OnShowAll`, `OnTogglePhase`, `i_phase_index`, `i_phase_type`, `jShow`, `jPhaseItems`, the Show-phases submenu, and the `showPhase` argument in `i_apply`→`view_dynamics('Redraw', …, showPhase)`.
- **Infrastructure:** `SetTarget`, `OnCloseSession`, `OnFigureDeleted`, `i_close_panel`, `i_cs`, `i_apply`, `i_field`, `i_firsttime`, `i_jump`, `i_str`, File menu (`FileOpen`/`FileSave`/`FileSaveAs`), Atoms group management (`AtomAddGroup`/`AtomRenameGroup`/`AtomDeleteGroup`/`AtomSetColor`/`AtomSort`/`i_selected`).

## 6. Consolidation (#6)

Introduce one local helper:

```matlab
% Reduce a real/complex/vector field [k*nRows x nT] to a per-row magnitude [nRows x nT]
% (scalar passes through). k = 3 (Dirac/vector) or 4 (quaternion) -> per-row L2 norm;
% complex (Connection Laplacian) -> abs. nRows is the divisor for THIS call: the basis-support
% count nGv at the impulse site, the full-surface count nV at the apply site.
function s = i_paintable_scalar(F, nRows)
    if ~isreal(F), F = abs(F); end
    if size(F,1) == nRows, s = F; return; end
    if mod(size(F,1), nRows) == 0
        nc = size(F,1) / nRows;
        s = reshape(sqrt(sum(reshape(F, nc, nRows, []).^2, 1)), nRows, []);
    else
        s = F;            % unexpected shape -> caller guards
    end
end
```

The three reduction sites it replaces (verified against the current file):
1. `i_atom_realise` — the `~isreal(W)` + vector-reshape block (current `:1563-1567`); call `i_paintable_scalar(W, nGv)` (reduces on the **basis-support** row count `nGv`, not surface `nV`).
2. `i_atom_apply` — the pre-filter reduction `Fr = i_vec2scalar(F, nV)` (current `:1700`) → `i_paintable_scalar(F, nV)`.
3. `i_atom_apply` — the post-filter reduction `Ffilt = i_vec2scalar(Ffilt, nV)` (current `:1710`) → `i_paintable_scalar(Ffilt, nV)`.

Then **delete `i_vec2scalar`** (its only callers are the two in (2)–(3)).

`i_atom_normalize` is **left as-is** — it is a *separate concern* (sign-class density-vs-peak choice driving sequential-vs-diverging colormap), not a shape reduction. Realise keeps two entry points (`i_atom_realise` for the impulse response via `bst_eigenfilter('Atom')`; `i_atom_filter_field` for the filtered real source via `bst_eigenfilter('Analysis')` / `bst_eigenwavelet('JTVAnalysis')`) because they invoke different backend verbs; they share only the reduction tail.

## 7. Verification plan

1. **Static:** run MATLAB `checkcode` (M-Lint) on `panel_bst_dynamics.m`. Acceptance: no "undefined function or variable" and no reference to any deleted control field (`jTree`, `jFreqC`, `jFreqW`, `jFreqBand`, `jSrcC`, `jSrcW`, `jScaleC`, `jScaleW`, `jListOccur`, `jRegionTool`, `jPeaks`). Filter out the standard Brainstorm idioms (`%#ok` suppressions, `bst_call` patterns).
2. **Grep cross-check:** for each deleted function and each pruned state field, confirm zero remaining references in `panel_bst_dynamics.m` (and no external caller in `toolbox/` other than `bst_geodesic_tool`→`SyncSource`, which is retained).
3. **Live (MCP):** launch Brainstorm; open a Dirac-dSPM result; `view_dynamics('FromResult', <result>)`; then exercise, capturing a screenshot at each step:
   - Create atom → atom appears in the list, impulse paint on cortex.
   - Edit a kernel slider → paint updates.
   - Localize → click a cortex vertex → seed moves, paint re-centres.
   - Apply ON → filtered real-source paint; Apply OFF → back to impulse.
   - Measure → Curl → differential overlay renders.
4. **Acceptance:** panel opens with no console errors; all five live steps behave as before the refactor.

## 8. Risks & mitigations

- **Missed live caller of a "dead" helper.** Mitigation: the grep cross-check (step 2) + `checkcode` (step 1) catch undefined references before launch; the live run catches runtime breakage. Fix by restoring the minimal helper.
- **A retained function transitively calls a deleted one.** Mitigation: the keep-list in §5 was derived by tracing each retained function's calls; the grep cross-check confirms.
- **Phase UI with no producer looks broken to a user.** Accepted by explicit decision; revisit in C.

## 9. Out of scope

- Any behavior change, new feature, or new file.
- The frame model, cached projection, scalogram, operator gate, Dirac sensor forward (sub-projects B–D).
- `view_atom_designer` consolidation (its realiser does not duplicate the reduce/normalize logic).
- Moving the Measure/differential workflow elsewhere (kept in place, unchanged).
