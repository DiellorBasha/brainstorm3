# Dynamics Panel Substrate Cleanup — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Prune three generations of dead/retired code from `toolbox/gui/panel_bst_dynamics.m` and collapse the duplicated magnitude-reduction logic into one helper, leaving a clean filterbank portal (plus the differential Measure overlay) with identical behavior.

**Architecture:** Single-file, behavior-preserving refactor. Delete in dependency order so each task is self-contained — **refphase → record/capture → navigator-last** (the navigator task owns the four functions `i_focus_time`/`i_rec_figure`/`i_driving`/`i_read_block` that the earlier clusters also call, since by then all their callers are gone). Then prune the session-state struct (repointing the one retained reader from `curGroup` to the live `curAtom`) and consolidate three magnitude-reduction sites into `i_paintable_scalar`. Each task is verified by `checkcode` + targeted `grep` cross-checks; the final task is an MCP-driven live smoke-test of the panel.

**Note on `NotifySelection`:** it is a public hook still called by `figure_spectrum.m` and `figure_timeseries.m` (external files, out of scope to edit). It is therefore **not deleted** — its body (which referenced now-deleted sync helpers) is reduced to a no-op stub, preserving the external contract.

**Tech Stack:** MATLAB (Brainstorm fork). Verification via the brainstorm-dev MATLAB MCP (`check_matlab_code` for M-Lint, `evaluate_matlab_code` / `run_matlab_file` for live checks) and `git` + `grep` from bash.

## Global Constraints

- Edit **only** `toolbox/gui/panel_bst_dynamics.m`. No changes to `bst_dynamics.m`, `bst_eigen*.m`, `view_dynamics.m`, `view_atom_designer.m`, `bst_geodesic_tool.m`, `panel_eigenfilter_design.m`.
- **Behavior-preserving** for all retained features (filterbank Create→edit→Localize→Apply, Measure differential overlay, phase-display UI, File/Atoms menus, close/teardown). The sole intentional deviation is the `i_selected` `curGroup`→`curAtom` repoint (Task 4), which fixes an already-broken menu path.
- Do **not** add new files, new features, or new dependencies.
- The phase-display UI (`OnShowAll`, `OnTogglePhase`, `i_phase_index`, `i_phase_type`, `jShow`, `jPhaseItems`, Show-phases submenu, and the `showPhase` arg to `view_dynamics('Redraw')`) is **retained** even though its only former producer (refphase) is removed.
- After every deletion, two gates must pass: (1) `check_matlab_code` reports no new "undefined function/variable" and no reference to a deleted symbol; (2) `grep` shows zero remaining references to each deleted symbol in the file (and no external caller in `toolbox/`).
- Commit after each task. Commit messages: `refactor(dynamics): <summary>` with the standard session trailers.
- Branch: work on `development` (current). Do not push unless asked.

**Reference — controls that `CreatePanel` builds (anything else is a dead-control reference):**
`jListAtoms, jMenuFile, jMenuAtoms, jKernel, jAtomParams, jLocalize, jAtomInfo, jApply, jOpItems, opVariants, atomKeys, jShow, jPhaseItems`.

---

### Task 1: Delete refphase detection + its toolbar buttons

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m`

**Interfaces:**
- Consumes: nothing.
- Produces: toolbar with the **Detect** and **Clear** buttons removed; **Measure** and **Show-all** buttons retained.

**Note:** the deleted `OnDetect` calls `i_focus_time` (and the navigator chain calls `i_rec_figure`/`i_driving`/`i_read_block`). Those four functions are **NOT** deleted here — they are owned by Task 3 (navigator-last), where their last callers disappear. Deleting only the refphase functions in this task leaves no dangling references because the four shared functions still exist.

- [ ] **Step 1: Delete the refphase subfunctions**

Delete in full: `OnDetect`, `OnSaveDetection`, `OnClearDetection`, `i_detect_events`, `i_remove_band`, `i_load_meg`, `i_has_staged_detection`, `OnSave`.

- [ ] **Step 2: Delete the Detect and Clear toolbar buttons in `CreatePanel`**

In `CreatePanel`, delete exactly these two `gui_component` lines (currently near lines 124 and 127):
```matlab
gui_component('ToolbarButton', jToolbar2, [], '', {IconLoader.ICON_EVT_TYPE_ADD, TB_DIM}, 'Detect windows: run the band-power detector on the selected band (preview events; not saved)', @(h,e)bst_call(@OnDetect));
```
```matlab
gui_component('ToolbarButton', jToolbar2, [], '', {IconLoader.ICON_EVT_TYPE_DEL, TB_DIM}, 'Clear: discard the staged detection preview (no save)', @(h,e)bst_call(@OnClearDetection));
```
Keep the `jShow` (Show-all) toggle line between them. After removal, collapse any now-doubled `jToolbar2.addSeparator();` so there are not two separators in a row.

- [ ] **Step 3: Run M-Lint**

Call `check_matlab_code` on the file. Expected: no undefined-reference errors. (`i_focus_time` etc. are still defined, so calls to them from any not-yet-deleted code are fine.)

- [ ] **Step 4: Grep cross-check**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
for fn in OnDetect OnSaveDetection OnClearDetection i_detect_events i_remove_band \
  i_load_meg i_has_staged_detection OnSave; do
    n=$(grep -c "\b$fn\b" toolbox/gui/panel_bst_dynamics.m)
    if [ "$n" -ne 0 ]; then echo "STILL REFERENCED: $fn ($n)"; fi
done
echo "done"
```
Expected: `done`, no "STILL REFERENCED". (Note: `\bOnSave\b` must not match `OnSaveDetection`/`OnSaveFilterbank`/`OnSaveCursor` — the `\b` word boundaries handle this; `OnSaveFilterbank` and `OnSaveCursor` are different tokens and are not deleted in this task.)

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m
git commit -m "refactor(dynamics): remove refphase detection (Detect/Clear)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```

---

### Task 2: Delete record/capture + its menu item + orphaned helpers

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m`

**Interfaces:**
- Consumes: nothing.
- Produces: Atoms menu with **Record at cursor** removed. `SyncSource` (the Localize-tool feedback called by `bst_geodesic_tool.m:168`) is **retained**.

**Note:** the deleted `OnSaveCursor` calls `i_read_block('time')`. `i_read_block` is **NOT** deleted here — it is owned by Task 3 (its other caller, `OnAxisChange`, is in the navigator cluster). Deleting only the record/capture functions leaves no dangling reference because `i_read_block` still exists.

- [ ] **Step 1: Delete the record/capture subfunctions**

Delete in full: `OnRecord`, `OnSaveCursor`, `OnCaptureRegion`, `OnRegionTool`, `ctrl_region_state`.

Delete the now-orphaned helpers (verified single-use by the above): `i_peaks`, `i_find_group`, `i_op_color`, `i_disp_band`, `i_scale_name`, `i_first_results`.

**Do NOT delete `SyncSource`.**

- [ ] **Step 2: Delete the "Record at cursor" Atoms-menu item in `CreatePanel`**

Delete this line (currently near line 89):
```matlab
gui_component('MenuItem', jMenuAtoms, [], 'Record at cursor', IconLoader.ICON_EVT_TYPE_ADD, [], @(h,e)bst_call(@OnRecord));
```
If this leaves a dangling `jMenuAtoms.addSeparator();` with a separator immediately adjacent on both sides, remove one separator so the menu has no doubled divider.

- [ ] **Step 3: Run M-Lint**

Call `check_matlab_code` on the file. Expected: no undefined-reference errors.

- [ ] **Step 4: Grep cross-check (internal + external)**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
for fn in OnRecord OnSaveCursor OnCaptureRegion OnRegionTool ctrl_region_state \
  i_peaks i_find_group i_op_color i_disp_band i_scale_name i_first_results; do
    n=$(grep -c "\b$fn\b" toolbox/gui/panel_bst_dynamics.m)
    if [ "$n" -ne 0 ]; then echo "STILL REFERENCED: $fn ($n)"; fi
done
echo "--- SyncSource must REMAIN (expect >=1 here and an external caller) ---"
grep -c "\bSyncSource\b" toolbox/gui/panel_bst_dynamics.m
grep -rn "SyncSource" toolbox/dynamics/bst_geodesic_tool.m
echo "done"
```
Expected: no "STILL REFERENCED"; `SyncSource` count ≥ 1; the `bst_geodesic_tool.m` line is present.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m
git commit -m "refactor(dynamics): remove record/capture (keep SyncSource localize hook)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```

---

### Task 3: Delete the unreachable navigator / tree / PSD-overlay cluster (navigator-last)

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m`

**Interfaces:**
- Consumes: refphase (Task 1) and record/capture (Task 2) already deleted — so `OnDetect`/`OnSaveCursor` no longer call the shared functions.
- Produces: `BuildTree()` reduced to a thin wrapper `function BuildTree(); UpdateAtomList(); end` (still called by `SetTarget` and `i_apply`); `NotifySelection` reduced to a no-op stub (retained for its external callers).

- [ ] **Step 1: Delete the navigator / tree / PSD subfunctions**

Delete these subfunctions in full from `panel_bst_dynamics.m`:

Tree / occurrence list:
- `TreeSel_Callback`, `OccurSel_Callback`, `i_window_atoms`, `i_group_atoms`

4-axis navigator (includes the four shared functions whose other callers were removed in Tasks 1-2):
- `i_axis_block`, `i_read_block`, `OnAxisChange`, `OnFreqPreset`, `i_drive`, `i_freq_preset`, `i_freq_name`, `i_fill_block`, `OnLoadAtom`, `i_band_match`, `i_bands`

Frequency / PSD overlay + time/freq sync:
- `i_rec_figure`, `i_owns_rec`, `i_owns_spec`, `i_ensure_psd`, `i_fix_spec_xlim`, `i_find_psd_file`, `i_compute_psd`, `i_freq_overlay`, `i_freq_overlay_clear`, `i_sync_freq`, `i_focus_time`, `i_sync_time`, `i_driving`

- [ ] **Step 2: Reduce `BuildTree` to a wrapper**

`BuildTree` is called by retained code (`SetTarget`, `i_apply`). Do **not** delete it; replace its entire body so it becomes exactly:
```matlab
function BuildTree()
    UpdateAtomList();
end
```

- [ ] **Step 3: Reduce `NotifySelection` to a no-op stub**

`NotifySelection` is a public hook still called by `figure_spectrum.m` and `figure_timeseries.m` (out-of-scope files). Do **not** delete it. Replace its entire body (which references the now-deleted `i_owns_spec`/`i_owns_rec`/`i_sync_freq`/`i_sync_time`/`i_driving`) so it becomes exactly:
```matlab
function NotifySelection(hFig, axis, range) %#ok<DEFNU,INUSD>
    % Retired no-op. The freq/time selection-sync was removed with the Navigator
    % strip; figure_spectrum / figure_timeseries still call this hook, so the entry
    % point is retained as a no-op to preserve their contract.
end
```

- [ ] **Step 4: Run M-Lint on the file**

Call the MCP tool `check_matlab_code` on `toolbox/gui/panel_bst_dynamics.m`.
Expected: no "undefined function or variable" referencing any deleted name; no errors. (Pre-existing `%#ok` style notices are fine.)

- [ ] **Step 5: Grep cross-check — deleted symbols have zero references**

Run:
```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
for fn in TreeSel_Callback OccurSel_Callback i_window_atoms i_group_atoms \
  i_axis_block i_read_block OnAxisChange OnFreqPreset i_drive i_freq_preset \
  i_freq_name i_fill_block OnLoadAtom i_band_match i_bands \
  i_rec_figure i_owns_rec i_owns_spec i_ensure_psd i_fix_spec_xlim \
  i_find_psd_file i_compute_psd i_freq_overlay i_freq_overlay_clear i_sync_freq \
  i_focus_time i_sync_time i_driving; do
    n=$(grep -c "\b$fn\b" toolbox/gui/panel_bst_dynamics.m)
    if [ "$n" -ne 0 ]; then echo "STILL REFERENCED: $fn ($n)"; fi
done
echo "--- NotifySelection RETAINED as a stub (expect 1 def + its external callers) ---"
grep -c "\bNotifySelection\b" toolbox/gui/panel_bst_dynamics.m
echo "done"
```
Expected output: no "STILL REFERENCED" lines; `NotifySelection` count = 1 (its definition only — its former internal call has been removed with the body); `done`.

- [ ] **Step 6: Grep cross-check — `NotifySelection` external callers still resolve**

Run:
```bash
grep -rln "NotifySelection" toolbox/ | grep -v panel_bst_dynamics.m
```
Expected: `figure_spectrum.m` and `figure_timeseries.m` appear (their `panel_bst_dynamics('NotifySelection', ...)` calls now hit the retained no-op stub — this is correct, not a defect).

- [ ] **Step 7: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m
git commit -m "refactor(dynamics): delete unreachable navigator/tree/PSD; stub NotifySelection

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```

---

### Task 4: Prune the session-state struct + repoint `i_selected` to `curAtom`

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m`

**Interfaces:**
- Consumes: `st.curAtom` (the live list selection, set by `AtomsListValueChanged_Callback`/`SetSelectedAtom`).
- Produces: `i_selected()` returns `g = st.curAtom`; the `SetTarget` state struct carries only the retained fields.

- [ ] **Step 1: Repoint `i_selected` to the live selection**

In `i_selected` (currently near line 1384), change:
```matlab
    g = st.curGroup;
    if g < 1, java_dialog('warning', 'Select a band atom in the tree first.', 'Atoms'); end
```
to:
```matlab
    g = i_field(st, 'curAtom', 0);
    if g < 1, java_dialog('warning', 'Select an atom in the list first.', 'Atoms'); end
```

- [ ] **Step 2: Update the two `curGroup` resets**

In `AtomDeleteGroup` (currently near line 1362), change `st.curGroup = 0;` to `st.curAtom = 0;`.
In `AtomSort` (currently near line 1378), change `st.curGroup = 0;` to `st.curAtom = 0;`.

- [ ] **Step 3: Prune the `SetTarget` state struct**

In `SetTarget` (currently near line 1118), replace the `struct(...)` literal so it seeds only the retained fields. Replace:
```matlab
    setappdata(0, 'DynamicsTarget', struct('hFig',hFig, 'T',T, 'file',file, 'curGroup',0, ...
        'nodeList',{ {} }, 'nodeInfo',[], 'occMap',[], 'Lambda',[], 'showPhase',[1 1 1 1], ...
        'hSpec',[], 'focusTime',[], 'detSel',[], ...
        'nav', bst_dynamics('NewGroup', 'cursor')));
```
with:
```matlab
    setappdata(0, 'DynamicsTarget', struct('hFig',hFig, 'T',T, 'file',file, ...
        'curAtom',0, 'atomSeed',[], 'showPhase',[1 1 1 1], 'curOp','none'));
```
(`atomThreshold`, `atomAx`, `atomBounds` are added lazily by existing code via `i_field`/`i_atom_ensure_axes`, so they need no seed.)

- [ ] **Step 4: Run M-Lint**

Call `check_matlab_code` on the file. Expected: no undefined-reference errors.

- [ ] **Step 5: Grep cross-check — pruned fields have zero remaining reads**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
for f in nav occMap nodeList nodeInfo Lambda hSpec focusTime detSel curGroup curBand curBandName curScale; do
    n=$(grep -c "\b$f\b" toolbox/gui/panel_bst_dynamics.m)
    if [ "$n" -ne 0 ]; then echo "STILL PRESENT: $f ($n)"; fi
done
echo "done"
```
Expected: `done`, no "STILL PRESENT". (`Lambda` as a struct field must be gone; `ax.Lambda` was already removed with no such token because the grep is `\bLambda\b` — if `ax.Lambda` usages in retained code match, that is expected and acceptable: verify any hits are `ax.Lambda`, not `st.Lambda`. If the only hits are `ax.Lambda`, treat as PASS.)

> Note for the implementer: `\bLambda\b` will match the retained `ax.Lambda` usages. Re-run `grep -n "\bLambda\b" toolbox/gui/panel_bst_dynamics.m` and confirm **every** hit is `ax.Lambda` (axes struct), not `st.Lambda` or a bare `Lambda` field. That is the real acceptance for `Lambda`.

- [ ] **Step 6: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m
git commit -m "refactor(dynamics): prune dead state fields; i_selected uses live curAtom

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```

---

### Task 5: Consolidate magnitude-reduction into `i_paintable_scalar`

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m`

**Interfaces:**
- Consumes: nothing.
- Produces: `function s = i_paintable_scalar(F, nRows)` — reduces a real/complex/vector field `[k*nRows x nT]` to per-row magnitude `[nRows x nT]`; scalar passes through. Used by `i_atom_realise` (with `nGv`) and `i_atom_apply` (with `nV`). `i_vec2scalar` is deleted.

- [ ] **Step 1: Add the helper**

Add this subfunction (place it next to the current `i_vec2scalar`, near line 1720):
```matlab
% Reduce a real/complex/vector field [k*nRows x nT] to a per-row magnitude [nRows x nT]
% (scalar passes through). nRows is the divisor for THIS call: the basis-support count nGv
% at the impulse site, the full-surface count nV at the apply site.
function s = i_paintable_scalar(F, nRows)
    if ~isreal(F), F = abs(F); end
    if size(F,1) == nRows, s = F; return; end
    if mod(size(F,1), nRows) == 0
        nc = size(F,1) / nRows;
        s = reshape(sqrt(sum(reshape(F, nc, nRows, []).^2, 1)), nRows, []);
    else
        s = F;                                                     % unexpected shape -> caller guards
    end
end
```

- [ ] **Step 2: Replace the reduction block in `i_atom_realise`**

In `i_atom_realise` (currently near lines 1563-1567), replace:
```matlab
    if ~isreal(W), W = abs(W); end                                  % complex tangent (Connection Laplacian) -> magnitude
    if (size(W,1) > nGv) && (mod(size(W,1), nGv) == 0)              % vector/quaternion basis (Dirac k=4) -> per-vertex magnitude
        nc = size(W,1) / nGv;
        W  = reshape(sqrt(sum(reshape(W, nc, nGv, []).^2, 1)), nGv, []);
    end
    if size(W,1) ~= nGv, W = [];  return; end                       % unexpected shape -> not paintable (guarded)
```
with:
```matlab
    W = i_paintable_scalar(W, nGv);
    if size(W,1) ~= nGv, W = [];  return; end                       % unexpected shape -> not paintable (guarded)
```

- [ ] **Step 3: Replace the two `i_vec2scalar` calls in `i_atom_apply`**

In `i_atom_apply`, change the pre-filter reduction (currently near line 1700):
```matlab
        Fr = i_vec2scalar(F, nV);                                   % scalar operator: per-vertex magnitude
```
to:
```matlab
        Fr = i_paintable_scalar(F, nV);                             % scalar operator: per-vertex magnitude
```
and the post-filter reduction (currently near line 1710):
```matlab
    Ffilt = i_vec2scalar(Ffilt, nV);
```
to:
```matlab
    Ffilt = i_paintable_scalar(Ffilt, nV);
```

- [ ] **Step 4: Delete `i_vec2scalar`**

Delete the `i_vec2scalar` subfunction (currently near lines 1719-1728) in full.

- [ ] **Step 5: Run M-Lint**

Call `check_matlab_code` on the file. Expected: no undefined-reference errors; no remaining reference to `i_vec2scalar`.

- [ ] **Step 6: Grep cross-check**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
echo "i_vec2scalar (expect 0):"; grep -c "\bi_vec2scalar\b" toolbox/gui/panel_bst_dynamics.m
echo "i_paintable_scalar (expect >=3: def + 3 calls):"; grep -c "\bi_paintable_scalar\b" toolbox/gui/panel_bst_dynamics.m
```
Expected: `i_vec2scalar` = 0; `i_paintable_scalar` ≥ 4 (1 definition + 3 call sites).

- [ ] **Step 7: MCP math sanity-check of the helper**

Precondition: Brainstorm path loaded in MATLAB (`brainstorm` started or its `toolbox/` on the path so `panel_bst_dynamics` dispatches via `macro_method`).

Call `evaluate_matlab_code` with:
```matlab
a = panel_bst_dynamics('i_paintable_scalar', [3;4;0], 1);      % 3-vector -> 5
b = panel_bst_dynamics('i_paintable_scalar', [1;2;2;0], 1);    % quaternion -> 3
c = panel_bst_dynamics('i_paintable_scalar', [3+4i;0], 2);     % complex -> [5;0]
d = panel_bst_dynamics('i_paintable_scalar', [1 2;3 4], 2);    % scalar passthrough
assert(abs(a-5)<1e-12 && abs(b-3)<1e-12 && isequal(c,[5;0]) && isequal(d,[1 2;3 4]));
disp('i_paintable_scalar OK');
```
Expected: prints `i_paintable_scalar OK` with no assertion error.

- [ ] **Step 7b: If the dispatch call errors** (e.g. `macro_method` does not expose the subfunction in this build), fall back to pasting the function body inline into `evaluate_matlab_code` as a local anonymous re-implementation and asserting the same four cases. Record which path was used.

- [ ] **Step 8: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m
git commit -m "refactor(dynamics): single i_paintable_scalar replaces 3 reduction sites

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```

---

### Task 6: MCP live smoke-test (acceptance)

**Files:**
- No code change expected. If a regression is found, fix it in `panel_bst_dynamics.m` and re-run.

**Interfaces:**
- Consumes: the cleaned panel from Tasks 1-5.
- Produces: confirmation (screenshots + console-clean) that all retained features behave as before.

- [ ] **Step 1: Launch Brainstorm**

Use the `brainstorm-dev:brainstorm-start` skill (or `evaluate_matlab_code` running `brainstorm` if a desktop session is already up). Confirm the GUI is up with no startup errors.

- [ ] **Step 2: Open a Dirac-dSPM source result and the Dynamics session**

In a protocol containing an unconstrained Dirac result (e.g. `results_DiracEig_KERNEL_*` on `S01_AEF_01_notch` per the project's test data), run via `evaluate_matlab_code`:
```matlab
% Resolve a Dirac results file path from the current protocol, then:
view_dynamics('FromResult', ResultsFile);
```
Expected: a `figure_3d` cortex opens and the docked **Dynamics** panel appears, no console error.

- [ ] **Step 3: Exercise the filterbank — Create + edit + Localize**

Programmatically or via screenshot-guided clicks:
1. `panel_bst_dynamics('OnCreateAtom')` → an atom appears in the list; impulse paint on the cortex.
2. Change a kernel slider (or call `panel_bst_dynamics('OnParamSettle')` after setting a value) → paint updates.
3. `panel_bst_dynamics('OnLocalize')` to arm, click a cortex vertex → seed moves, paint re-centres.

Capture a screenshot after each. Expected: no error; paint tracks the edits.

- [ ] **Step 4: Exercise Apply (Design ↔ Preview)**

Toggle Apply ON (`panel_bst_dynamics('OnApply')` with the toggle selected) → filtered real-source paint appears; toggle OFF → impulse paint returns. Screenshot both. Expected: both paints render; no error in the `i_paintable_scalar` path (this exercises the 3-vector and scalar reductions on real data).

- [ ] **Step 5: Exercise Measure (differential overlay)**

Open the Measure menu and select **Curl** (`panel_bst_dynamics('OnMeasurement','Curl')`). Expected: the differential overlay renders (purple Curl map); selecting **None** restores the native source paint. Screenshot.

- [ ] **Step 6: Exercise the repointed Atoms menu**

With an atom selected in the list, invoke `panel_bst_dynamics('AtomSetColor')` (pick a color) and `panel_bst_dynamics('AtomRenameGroup')` — confirm they now act on the selected atom (no "Select an atom in the list first" warning). Expected: color/label of the selected atom changes. This validates the `curGroup`→`curAtom` repoint.

- [ ] **Step 7: Record acceptance**

Confirm all of Steps 2-6 succeeded with no console errors and behavior matching pre-refactor. If any step regressed, debug and fix in `panel_bst_dynamics.m`, then re-run the affected step and commit:
```bash
git add toolbox/gui/panel_bst_dynamics.m
git commit -m "fix(dynamics): <regression fix from smoke-test>

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```

---

## Acceptance criteria (whole plan)

- `panel_bst_dynamics.m` reduced from ~1861 to ~1050 lines; one file.
- `check_matlab_code` clean (no undefined refs, no dead-control refs).
- All deleted symbols: zero references in the file; no broken external callers (only `SyncSource` is called externally, and it is retained).
- The five live smoke-test flows (open session, Create/edit/Localize, Apply on/off, Measure Curl, repointed Atoms menu) all work with no console errors.
- Filterbank, Measure overlay, and phase-display UI behave as before; Atoms-menu group management now works against the live list selection.
