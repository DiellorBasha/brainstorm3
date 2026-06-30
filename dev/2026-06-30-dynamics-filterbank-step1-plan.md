# Dynamics filterbank (step 1) Implementation Plan

> **For agentic workers:** Execute INLINE (superpowers:executing-plans) with live checkpoints — Swing GUI + a live dynamics-session gate. Task 1 (data-model helpers) is headless-TDD'd; Task 2 (panel layout + wiring) is built + validated live by the controller.

**Goal:** Reorient `panel_bst_dynamics` toward an eigenwavelet filterbank: a "+ Create atom" button adds a default Diffusion **filter** atom (no thresholding) that previews on the source map; atom controls/actions are reorganized (stacked list, params-only Atom section, Localize/Threshold/Save on the toolbar).

**Architecture:** An atom = a generator `atomgroup` (`KernelName/KernelParams/vertices`, threshold/Scout/Event unset). Create builds a default-diffusion atom and previews via the existing `i_atom_realise`/`i_atom_preview`/`SetAtomField`. The tree|list split flips to vertical (atom list / detail); Localize/Threshold/Save move to the east toolbar; the Atom section keeps only Operator+Filter+sliders.

**Tech Stack:** MATLAB R2023b, Brainstorm dev fork (Java/Swing via `gui_component`). Reuse: `bst_dynamics`, `bst_eigfilter_controls`, `panel_eigenfilter_design`, `bst_eigenfilter('Atom')`, `view_dynamics('SetAtomField')`, `bst_geodesic_tool`.

## Global Constraints

- No new dependencies. Live validation in the Brainstorm session (preventad; the session has been unstable — restart `brainstorm` if it drops). Source-map preview ONLY this step (no sensor projection). Legacy detection toolbar actions (Detect/Show/Clear/Measure) stay untouched.
- An atom is a FILTER: created atoms carry their generator and leave `Threshold`/`region`/`times` UNSET. `AtomFromKernel` is no longer on the create path.
- The "+ Create atom" button copies `panel_scout`'s Create-scout button design (a `ToolbarButton` with `ICON_SCOUT_NEW` + an HTML-bold tooltip), adapted to click-to-create-default-atom.
- `lint` every edited `.m`; commit after each task with `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

## File structure
- Modify: `toolbox/gui/panel_bst_dynamics.m` — the data-model helpers (Task 1) + the layout/wiring (Task 2).
- Test: `dev/tests/test_dynamics_filter_atom.m` (Task 1).

---

### Task 1: default-filter-atom builder + detail rows (headless)

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m`; Test `dev/tests/test_dynamics_filter_atom.m`

**Interfaces (produces):**
- `G = i_default_atom(kernelName, kp, seed, surfaceFile, label)` — builds a generator `atomgroup`: `NewGroup(label)`, `KernelName=kernelName`, `KernelParams=kp`, `vertices=seed`, `SurfaceFile=surfaceFile`; `Threshold`/`region`/`times` left at template defaults (unset). The atom-as-filter, no thresholding.
- `rows = i_atom_detail(G)` — a cellstr of property rows for the bottom pane: `{'kernel: <name>', 'seed: vtx <n>', '<param>: <val>', …}` from the generator.

- [ ] **Step 1: failing test** — `dev/tests/test_dynamics_filter_atom.m`:

```matlab
% test_dynamics_filter_atom - a created atom is a FILTER (generator set, threshold/scout/event unset)
kp = struct('lmax',40,'tau',0.3);
G  = panel_bst_dynamics('i_default_atom', 'diffusion', kp, 13, 'surf.mat', 'atom1');
assert(strcmp(G.label,'atom1') && strcmp(G.KernelName,'diffusion'), 'label + kernel');
assert(isequal(G.KernelParams, kp) && isequal(G.vertices,13), 'generator carried');
assert(isempty(G.Threshold) && isempty(G.region) && isempty(G.times), 'NOT thresholded (filter, not marker)');
rows = panel_bst_dynamics('i_atom_detail', G);
assert(iscellstr(rows) && any(contains(rows,'diffusion')) && any(contains(rows,'13')), 'detail rows show kernel + seed');
disp('OK');
```

- [ ] **Step 2: run → fail** (function undefined).
- [ ] **Step 3: implement** `i_default_atom` + `i_atom_detail` as local functions (reachable via the panel's `eval(macro_method)`). `i_default_atom` uses `NewGroup` + sets the generator fields. `i_atom_detail` formats `kernel`/`seed`/each `KernelParams` field (skip `lmax`).
- [ ] **Step 4: run → pass** (`OK`).
- [ ] **Step 5: lint + commit** `feat(dynamics): default-filter-atom builder + detail rows (atom = generator, no threshold)`.

---

### Task 2: panel layout + create/select/save wiring (live)

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m`. **Live-validated.**

**Consumes:** Task 1's `i_default_atom`/`i_atom_detail`; existing `i_atom_ensure_axes`/`i_atom_realise`/`i_atom_preview`, `bst_eigfilter_controls`, `panel_eigenfilter_design`, `bst_geodesic_tool`.

- [ ] **Step 1: stacked split.** `panel_bst_dynamics.m:103` — change `JSplitPane(JSplitPane.HORIZONTAL_SPLIT, jScrollTree, jScrollOccur)` to `JSplitPane(JSplitPane.VERTICAL_SPLIT, jScrollTree, jScrollOccur)` (top = atom list `jTree`; bottom = detail `jListOccur`).

- [ ] **Step 2: toolbar buttons.** In the `jToolbar2` block (lines ~117-124): add a **Create-atom** button at the top (copied from `panel_scout`, adapted) and **Localize** + **Threshold** buttons; keep Detect/Show/Clear/Measure; repoint **Save** to `OnSaveFilterbank`:

```matlab
gui_component('ToolbarButton', jToolbar2, [], '', IconLoader.ICON_SCOUT_NEW, ...
    '<HTML><B>Create atom</B>:<BR><BLOCKQUOTE> - Click to add a default (diffusion) atom to the filterbank<BR> - Edit its parameters in the Atom section below</BLOCKQUOTE>', @(h,e)bst_call(@OnCreateAtom));
gui_component('ToolbarButton', jToolbar2, [], '', IconLoader.ICON_SAVE,     'Save the filterbank (atom table) to disk', @(h,e)bst_call(@OnSaveFilterbank));
jToolbar2.addSeparator();
jLocalize = gui_component('ToolbarToggle', jToolbar2, [], '', IconLoader.ICON_SCOUT_NEW, 'Localize: click a cortex vertex to re-seed the selected atom', @(h,e)bst_call(@OnLocalize)); %#ok<NASGU>
gui_component('ToolbarButton', jToolbar2, [], '', IconLoader.ICON_PROPERTIES, 'Threshold: set the level-set threshold for the optional Scout+Event export', @(h,e)bst_call(@OnThresholdMenu));
```
(Keep the existing Detect/Show/Clear/Measure buttons after these.) Remove the old standalone Save button line (118) — replaced above.

- [ ] **Step 3: Atom section = params only.** In the `jAtom`/`jRowA` block (lines ~143-146), DELETE the `jLocalize`/`jThresh`/`jStore` controls (they move to the toolbar). The Atom section keeps the operator toggles + filter combobox + `jAtomParams`. Update the `ctrl` struct (line ~156): drop `jStore`; `jLocalize`/`jThresh` are now toolbar handles — store the toolbar `jLocalize` + a panel field `jThresh` value holder (keep a `jThresh` text in `ctrl` with a default `'0.5'` value, set by `OnThresholdMenu`).

- [ ] **Step 4: `OnCreateAtom`.** Append a default-diffusion filter atom, select it, preview:

```matlab
function OnCreateAtom() %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    st = i_atom_ensure_axes(st);  if isempty(i_field(st,'atomAx',[])), return; end
    ax = st.atomAx;  lmax = max(ax.Lambda{1}(:));
    seed = ax.GlobalVertices{1}(1);                                  % default seed (first support vertex)
    vals = panel_eigenfilter_design('ReadAtomVals', ctrl.jAtomParams);  % current slider values if diffusion, else defaults
    kp = bst_eigfilter_controls('ToKernel', 'diffusion', vals, lmax);
    n  = numel(st.T.Groups) + 1;
    G  = i_default_atom('diffusion', kp, seed, ax.SurfaceFile, sprintf('atom%d', n));
    st.T = bst_dynamics('AddGroup', st.T, G);
    st.curAtom = numel(st.T.Groups);
    setappdata(0,'DynamicsTarget', st);
    BuildTree();  i_select_atom(st.curAtom);                         % select -> load params + preview + detail
end
```

- [ ] **Step 5: select/load/preview + detail.** `i_select_atom(idx)`: set `st.curAtom`, set the Filter combobox + sliders from `G.KernelName/KernelParams` (reuse `RebuildSliders` + `SetAtomVals`), set `st.atomSeed = G.vertices`, `i_atom_preview()`, and fill the bottom pane via `i_atom_detail(G)`. Wire the `jTree` selection callback to `i_select_atom`. `BuildTree` shows one flat node per `st.T.Groups` atom (label).

- [ ] **Step 6: Save / Localize / Threshold.** `OnSaveFilterbank`: `if ~isempty(st.file), bst_dynamics('Save', st.file, st.T)` (else prompt a path). `OnLocalize` (toolbar toggle): `bst_geodesic_tool('Toggle', state)` → `SyncSource` already sets `st.atomSeed` + previews; on a settled seed also write it into `st.T.Groups(st.curAtom).vertices`. `OnThresholdMenu`: a `java_dialog('input', …)` to set `ctrl.jThresh` value (stored on `st` as `st.atomThreshold`; used only by the future export).

- [ ] **Step 7: retire OnStore.** Remove `OnStore` (the Atom-section Store/threshold→Scout+Event path) — it is no longer the create path. (`AtomFromKernel` stays for the future export.)

- [ ] **Step 8: live gate.** Open a dynamics session → **+ Create atom** adds `atom1`, the diffusion pattern paints on the source map, the bottom pane shows kernel/seed/params; editing sliders re-previews; **Localize** re-seeds; **Save** writes the file; the split is stacked; the Atom section holds only parameters; Detect/Measure still present. Fix anything that misbehaves.

- [ ] **Step 9: commit** `feat(dynamics): + Create-atom filterbank flow, stacked list, toolbar actions, params-only Atom section`.

---

## Done criteria
- `+ Create atom` adds a default Diffusion filter atom (generator set, no threshold) that previews on the source map; the atom list is stacked above its detail; Localize/Threshold/Save are on the toolbar; the Atom section is parameters-only; `test_dynamics_filter_atom` passes; the live gate passes.

## Risks / notes
- Swing layout built live; the `+` button is copied from `panel_scout`.
- Reconcile the prior atom-tool controls (`jStore`/`jLocalize`/`jThresh` in the Atom section) — moved to the toolbar; update every `ctrl.jStore/jThresh` reference (OnStore removed; `OnThresholdMenu` owns the threshold value).
- `BuildTree` currently renders the old detection nest; repoint it to render flat filter atoms (one node per `st.T.Groups`).
