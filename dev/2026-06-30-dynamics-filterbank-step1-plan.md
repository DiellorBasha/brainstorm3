# Dynamics filterbank (step 1) Implementation Plan

> **For agentic workers:** Execute INLINE (superpowers:executing-plans) with live checkpoints — Brainstorm GUI (gui_component / BstClusterList, NOT raw Swing) + a live dynamics-session gate. Task 1 (data-model helpers) is headless-TDD'd; Task 2 (panel layout + wiring) is built + validated live, following panel_scout conventions.

**Goal:** Reorient `panel_bst_dynamics` toward an eigenwavelet filterbank: a "+ Create atom" button adds a default Diffusion **filter** atom (no thresholding) that previews on the source map; the atom list becomes a `BstClusterList`, the Atom section is parameters-only (the selected atom's properties), and Localize/Threshold/Save move to the east toolbar.

**Architecture (panel_scout-aligned):** An atom = a generator `atomgroup` (`KernelName/KernelParams/vertices`, threshold/Scout/Event unset). The CENTER tree|list split is **replaced by a single `BstClusterList`** (colored dot + `atomN` per row, like scouts/events). The SOUTH **Atom section is the properties** — selecting an atom loads its kernel/params there and previews via the existing `i_atom_realise`/`SetAtomField`. Refresh uses scout's `UpdateAtomList` / `UpdateAtomProps` / `SetSelectedAtom` chain with the **callback-suppression idiom**. Create/Localize/Threshold/Save are east-toolbar buttons.

**Tech Stack:** MATLAB R2023b, Brainstorm dev fork. Components: `gui_component`, `org.brainstorm.list.BstClusterList` + `BstClusterListRenderer`, `gui_river`, `IconLoader`, `java_scaled`, `bst_call`, `BstPanel`. Reuse: `bst_dynamics`, `bst_eigfilter_controls`, `panel_eigenfilter_design`, `bst_eigenfilter('Atom')`, `view_dynamics('SetAtomField')`, `bst_geodesic_tool`.

## Global Constraints

- No new dependencies; **Brainstorm components only** (no raw Swing layout). Live validation in the Brainstorm session (preventad; restart `brainstorm` if it drops). Source-map preview ONLY (no sensor projection). Legacy detection toolbar actions (Detect/Show/Clear/Measure) stay untouched.
- An atom is a FILTER: created atoms carry their generator and leave `Threshold`/`region`/`times` UNSET. `AtomFromKernel` is no longer on the create path.
- panel_scout conventions to mirror (verified by analysis): `BstClusterList` + `BstClusterListRenderer('I', fontSize)` for the list; `BstListItem(itemType, [], itemText, i, R,G,B)` rows; the callback-suppression idiom in `SetSelectedAtom`; `bst_call`-wrapped callbacks; `java_scaled`/`TB_DIM` sizing; `gui_river(gap, insets, 'Title')` titled sub-panels; the Create button = `ToolbarButton`/`ICON_SCOUT_NEW` + HTML-bold tooltip.
- `lint` every edited `.m`; commit after each task with `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

## File structure
- Modify: `toolbox/gui/panel_bst_dynamics.m` — data-model helpers (Task 1) + layout/wiring (Task 2).
- Test: `dev/tests/test_dynamics_filter_atom.m` (Task 1).

---

### Task 1: default-filter-atom builder + detail rows (headless)

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m`; Test `dev/tests/test_dynamics_filter_atom.m`

**Interfaces (produces):**
- `G = i_default_atom(kernelName, kp, seed, surfaceFile, label)` — builds a generator `atomgroup`: `NewGroup(label)`, `KernelName/KernelParams/vertices=seed/SurfaceFile`; `Threshold`/`region`/`times` left UNSET. The atom-as-filter.
- `s = i_atom_detail(G)` — a one-line summary string for the Atom-section readout: `'<kernel> · vtx <seed> · <param>=<val> …'` (skip `lmax`).

- [ ] **Step 1: failing test** — `dev/tests/test_dynamics_filter_atom.m`:

```matlab
% test_dynamics_filter_atom - a created atom is a FILTER (generator set, threshold/scout/event unset)
kp = struct('lmax',40,'tau',0.3);
G  = panel_bst_dynamics('i_default_atom', 'diffusion', kp, 13, 'surf.mat', 'atom1');
assert(strcmp(G.label,'atom1') && strcmp(G.KernelName,'diffusion'), 'label + kernel');
assert(isequal(G.KernelParams, kp) && isequal(G.vertices,13), 'generator carried');
assert(isempty(G.Threshold) && isempty(G.region) && isempty(G.times), 'NOT thresholded (filter, not marker)');
s = panel_bst_dynamics('i_atom_detail', G);
assert(ischar(s) && contains(s,'diffusion') && contains(s,'13'), 'detail shows kernel + seed');
disp('OK');
```

- [ ] **Step 2: run → fail** (function undefined).
- [ ] **Step 3: implement** `i_default_atom` + `i_atom_detail` as local functions (reached via `eval(macro_method)`). `i_default_atom` uses `NewGroup` + sets the generator fields. `i_atom_detail` formats kernel/seed + each `KernelParams` field except `lmax`.
- [ ] **Step 4: run → pass** (`OK`).
- [ ] **Step 5: lint + commit** `feat(dynamics): default-filter-atom builder + detail (atom = generator, no threshold)`.

---

### Task 2: panel layout (BstClusterList + toolbar) + create/select/save wiring (live)

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m`. **Live-validated, panel_scout-aligned.**

**Consumes:** Task 1's `i_default_atom`/`i_atom_detail`; existing `i_atom_ensure_axes`/`i_atom_realise`/`i_atom_preview`, `bst_eigfilter_controls`, `panel_eigenfilter_design`, `bst_geodesic_tool`.

- [ ] **Step 1: replace the split with a `BstClusterList`.** In `CreatePanel` (the `jTree`/`jScrollTree`/`jListOccur`/`jSplit` block, ~lines 89-129), remove the split and the two old widgets; build the atom list like panel_scout (L154-160):

```matlab
jListAtoms = java_create('org.brainstorm.list.BstClusterList');
jListAtoms.setCellRenderer(org.brainstorm.list.BstClusterListRenderer('I', fontSize));
java_setcb(jListAtoms, 'ValueChangedCallback', @(h,ev)bst_call(@AtomsListValueChanged_Callback,h,ev));
jScrollList = JScrollPane(jListAtoms);  jScrollList.setBorder(java_scaled('titledborder',''));
jPanelMain.add(jScrollList, BorderLayout.CENTER);
```
(`fontSize = round(11 * bst_get('InterfaceScaling')/100)`.) The SOUTH Atom section stays; it is the per-atom properties.

- [ ] **Step 2: east toolbar.** In `jToolbar2` (~L117-124): add **Create atom** (top), **Localize**, **Threshold**, repoint **Save**; keep Detect/Show/Clear/Measure. All sized with `TB_DIM = java_scaled('dimension',25,25)`:

```matlab
gui_component('ToolbarButton', jToolbar2, [], '', {IconLoader.ICON_SCOUT_NEW, TB_DIM}, ...
    '<HTML><B>Create atom</B>:<BR><BLOCKQUOTE> - Click to add a default (diffusion) atom to the filterbank<BR> - Edit its parameters in the Atom section below</BLOCKQUOTE>', @(h,e)bst_call(@OnCreateAtom));
gui_component('ToolbarButton', jToolbar2, [], '', {IconLoader.ICON_SAVE, TB_DIM}, 'Save the filterbank (atom table) to disk', @(h,e)bst_call(@OnSaveFilterbank));
jToolbar2.addSeparator();
jLocalize = gui_component('ToolbarToggle', jToolbar2, [], '', {IconLoader.ICON_SCOUT_NEW, TB_DIM}, 'Localize: click a cortex vertex to re-seed the selected atom', @(h,e)bst_call(@OnLocalize));
gui_component('ToolbarButton', jToolbar2, [], '', {IconLoader.ICON_PROPERTIES, TB_DIM}, 'Threshold: set the level-set threshold for the optional Scout+Event export', @(h,e)bst_call(@OnThresholdMenu));
jToolbar2.addSeparator();
```
(Existing Detect/Show/Clear/Measure follow; remove the old standalone Save line L118.)

- [ ] **Step 3: Atom section = params only + a readout.** In the `jAtom`/`jRowA` block (~L143-146), DELETE `jLocalize`/`jThresh`/`jStore`. Add a read-only detail label at the top of the Atom section: `jAtomInfo = gui_component('label', jAtom, 'br', '');`. Update the `ctrl` struct (~L156): replace `jTree/jListOccur/jStore` with `jListAtoms`, `jLocalize` (toolbar), `jAtomInfo`; keep `jKernel/jAtomParams/jConn/jGeom`. Store `st.atomThreshold` (default 0.5) on the target instead of a `jThresh` widget.

- [ ] **Step 4: `OnCreateAtom`.** Default-diffusion filter atom → append → select (which previews):

```matlab
function OnCreateAtom() %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    st = i_atom_ensure_axes(st);  if isempty(i_field(st,'atomAx',[])), return; end
    ax = st.atomAx;  lmax = max(ax.Lambda{1}(:));  seed = ax.GlobalVertices{1}(1);
    vals = panel_eigenfilter_design('ReadAtomVals', ctrl.jAtomParams);
    kp = bst_eigfilter_controls('ToKernel', 'diffusion', vals, lmax);
    G  = i_default_atom('diffusion', kp, seed, ax.SurfaceFile, sprintf('atom%d', numel(st.T.Groups)+1));
    st.T = bst_dynamics('AddGroup', st.T, G);  setappdata(0,'DynamicsTarget',st);
    UpdateAtomList();  SetSelectedAtom(numel(st.T.Groups));
end
```

- [ ] **Step 5: list build + select + props (scout idioms).**
  - `UpdateAtomList()` — rebuild the `BstClusterList` model from `st.T.Groups`: for each atom a `BstListItem('', [], G.label, i, R,G,B)` (assign a per-atom colour, e.g. `panel_scout`-style or a fixed palette); `jListAtoms.setModel(listModel)`.
  - `SetSelectedAtom(iAtom)` — **callback-suppressed** (save `ValueChangedCallback`, `setSelectedIndices(iAtom-1)`, restore), set `st.curAtom=iAtom`, then `i_select_atom_load(iAtom)`.
  - `AtomsListValueChanged_Callback(h,ev)` — `if ev.getValueIsAdjusting, return; end`; `i_select_atom_load(jListAtoms.getSelectedIndices()+1)`.
  - `i_select_atom_load(iAtom)` — `G = st.T.Groups(iAtom)`; set the Filter combobox to `G.KernelName` + `RebuildSliders` + `SetAtomVals` from `G.KernelParams`; set operator toggles from `G` (or current); `st.atomSeed = G.vertices`; `ctrl.jAtomInfo.setText(i_atom_detail(G))`; `i_atom_preview()`.

- [ ] **Step 6: Save / Localize / Threshold / param-edit writeback.**
  - `OnSaveFilterbank()` — `if ~isempty(st.file), bst_dynamics('Save', st.file, st.T); else <prompt path>`.
  - `OnLocalize()` (toolbar toggle) — `bst_geodesic_tool('Toggle', state)`; `SyncSource` already sets `st.atomSeed`+previews; on a settled seed also write `st.T.Groups(st.curAtom).vertices = st.atomSeed`.
  - `OnThresholdMenu()` — `java_dialog('input',…)` → `st.atomThreshold` (future export only).
  - `OnParamSettle()` (existing) — also write the edited kp back: `st.T.Groups(st.curAtom).KernelParams = bst_eigfilter_controls('ToKernel', CurrentKernel, ReadAtomVals, lmax)` + refresh `jAtomInfo`, before `i_atom_preview`.

- [ ] **Step 7: retire `OnStore`** (Atom-section Store path); repoint `BuildTree` callers to `UpdateAtomList` (or delete `BuildTree` if now unused). `AtomFromKernel` stays for the future export.

- [ ] **Step 8: live gate.** Open a dynamics session → **+ Create atom** adds `atom1` (a coloured row) and the diffusion pattern paints on the source map; the Atom section shows its params + a readout line; selecting/editing re-previews + updates the row/readout; **Localize** re-seeds; **Save** writes the file; the CENTER is a single atom list (no split); Detect/Measure still present. Fix anything that misbehaves.

- [ ] **Step 9: commit** `feat(dynamics): filterbank atom list (BstClusterList) + Create/Localize/Save toolbar, params-only Atom section`.

---

## Done criteria
- `+ Create atom` adds a default Diffusion filter atom (generator set, no threshold) that previews on the source map; the CENTER is a `BstClusterList` of atoms; the Atom section is parameters-only (loads on selection, writes back on edit); Localize/Threshold/Save are on the toolbar; `test_dynamics_filter_atom` passes; the live gate passes.

## Risks / notes
- Built with Brainstorm components to panel_scout conventions (BstClusterList, callback-suppression, java_scaled, gui_river titled panels) — NOT raw Swing.
- Reconcile prior atom-tool controls (`jStore`/`jThresh` in the Atom section, `BuildTree`, `OnStore`) — removed/repurposed; update every `ctrl.jStore/jThresh/jTree/jListOccur` reference.
- The legacy detection wiring referenced `jTree`/`jListOccur` (BuildTree, TreeSel/OccurSel callbacks); with the split removed, repoint or stub those so the legacy Detect/Show/Clear path doesn't error (it's out of scope but must not crash).
