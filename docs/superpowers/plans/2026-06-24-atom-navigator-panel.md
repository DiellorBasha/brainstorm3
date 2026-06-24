# Atom Navigator Panel Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `panel_bst_dynamics`'s Frequency/Space/Record control sections with four symmetric `(center, extent)` axis blocks (Time/Frequency/Source/Scale) + a Measurement row, driven through `bst_atom` into a transient cursor atom — the live Navigate state.

**Architecture:** One `i_axis_block` builder produces the symmetric row (numeric center + window on the left, the axis selector in the right slot) for all four axes. Editing a field or operating the selector calls `OnAxisChange(axis)`, which writes the axis Localization into a transient cursor group `st.nav` via `bst_atom('Set', …)` and dispatches to a thin per-axis driver that reuses the existing engines (`panel_time`/`panel_filter`/`bst_geodesic_tool`/`view_helmholtz`). The drivers also keep the legacy coordinate fields (`curBand`/`curBandName`/`curScale`/`curOp`) populated so Detect/Record/Capture stay unchanged.

**Tech Stack:** MATLAB, Brainstorm GUI (`gui_component`, `gui_river`, `BstPanel`), `bst_atom` (Phase 1), `bst_geodesic_tool` (Phase 2); tested headless under `brainstorm nogui`.

## Global Constraints

- Phase 3 of `docs/superpowers/specs/2026-06-24-atom-tensor-architecture-analysis.md`; spec `docs/superpowers/specs/2026-06-24-atom-navigator-panel-design.md`.
- **Symmetric layout:** every block is `center [ ] · window [± ]` on the LEFT and the **selector/preset in the same right-hand slot**. Source uses center (seed vertex) + window (radius) like the rest; its Region-tool toggle sits in the right slot (where the band combobox sits for Frequency).
- **Navigate-only:** the blocks drive the linked viewers live; nothing persists. Detect/Record/Capture remain the only writers, unchanged — so the drivers MUST keep `st.curBand`/`st.curBandName` (freq), `st.curOp` (measurement), `st.curScale` (scale) populated, since OnDetect/OnRecord read them.
- **Basic Scale:** the Scale `window` drives a low-pass heat eigenfilter (`view_helmholtz('SetSmoothing', hFig, 1, 'heat', struct('t', window))`); the Scale `center` field is present for symmetry but **reserved** (Phase 5). Kernel fixed to `'heat'`; no kernel/slider UI.
- **Measurement** (Φ/Ψ/|J|) is a descriptor row, not an axis — `view_helmholtz('SetComponent', …)`, mutually exclusive, all-off = Total.
- Freq band combobox = the frequency-atlas presets `i_bands()` (delta/theta/alpha/beta/gamma) + 'custom'; selecting a band fills the freq center/window fields.
- The Atoms table (tree/list/File+Atoms menus, Show-phases filter) and the Detect/Record/Capture/Peaks actions are KEPT; the Region-tool toggle relocates from the Record row into the Source block.
- Do not start implementation on `development`; this work is on branch `feature/atom-navigator-panel` (spec already committed there).
- Tests run headless with Brainstorm live in **`brainstorm nogui` (GuiLevel 0)** — NOT `brainstorm server` (-1), under which `bst_get('PanelControls','Dynamics')` is `[]`. Do not `clear`/restart Brainstorm or close all figures; `rehash` and re-run.

---

### Task 1: Navigator rebuild — four blocks, drivers, cursor atom

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (rebuild the `jCtrl` control area in `CreatePanel`; add `i_axis_block`/`OnAxisChange`/`i_read_block`/`i_drive`/`OnMeasurement`/`i_freq_preset`; fold `OnBand`→freq driver, `OnSpaceComp`→`OnMeasurement`; drop `OnSpaceSmooth`/`OnSpaceKernel`/`SetupSpace` slider UI; add `st.nav` in `SetTarget`)
- Modify: `dev/test_dynamics_atoms.m` (retarget T4's `jBands(3)`/`jSpaceStr` to the new handles)
- Test: `dev/test_nav_panel.m` (new)

**Interfaces:**
- Consumes: `bst_atom('Set'/'Get'/'NewLoc', …)` (Phase 1); `bst_dynamics('NewGroup', …)`; `bst_geodesic_tool('Toggle', …)` (Phase 2); `panel_filter`/`panel_time`/`view_helmholtz`; `i_bands()` (existing local).
- Produces:
  - Controls struct handles: `jTimeC, jTimeW, jFreqC, jFreqW, jFreqBand, jSrcC, jSrcW, jRegionTool, jScaleC, jScaleW, jMeasPot, jMeasStr` (plus kept `jTree, jListOccur, jMenuFile, jMenuAtoms, jPhaseItems, jPeaks`).
  - `st.nav` — a transient single-occurrence cursor group on the `DynamicsTarget` appdata.
  - `panel_bst_dynamics('OnAxisChange', axis)` / `('OnMeasurement', which)` verbs.

- [ ] **Step 1: Write the failing test**

Create `dev/test_nav_panel.m`:

```matlab
function test_nav_panel()
% TEST_NAV_PANEL: the four-axis (center,extent) navigator drives bst_atom + the engines.
%
% USAGE:  test_nav_panel   % Brainstorm running in nogui (GuiLevel 0)
%
% Authors: Diellor Basha, 2026

    PF = {'FAIL','PASS'};  pass = true;

    % open the panel via view_dynamics on a Dirac result (reuses the atom-suite fixture)
    [linkFile, relData] = i_find_kernel_nav();
    if isempty(linkFile)
        fprintf('SKIPPED (no unconstrained kernel link)\n');
        fprintf('\n==== SUITE: %s ====\n', PF{pass+1});  return;
    end
    sStudy = bst_get('DataFile', relData);
    R = '';  for j=1:numel(sStudy.Result), if ~isempty(regexp(sStudy.Result(j).Comment,'MN: MEG\(Unconstr\)','once')) && ~isempty(regexp(sStudy.Result(j).FileName,'KERNEL','once')), R = sStudy.Result(j).FileName; break; end; end
    hFig = view_dynamics('FromResult', R);  drawnow;
    ctrl = bst_get('PanelControls', 'Dynamics');

    % ---------- T1: controls struct has the new block handles, not the old ones ----------
    ok1 = all(isfield(ctrl, {'jFreqC','jFreqW','jFreqBand','jTimeC','jTimeW','jSrcC','jSrcW','jScaleC','jScaleW','jMeasPot','jMeasStr','jRegionTool'})) ...
       && ~isfield(ctrl,'jBands') && ~isfield(ctrl,'jSpaceStr');
    fprintf('T1 handles: newBlocks=%d oldGone=%d => %s\n', all(isfield(ctrl,{'jFreqC','jMeasStr'})), ~isfield(ctrl,'jBands'), PF{ok1+1});
    pass = pass && ok1;

    % ---------- T2: freq block -> st.nav freq Localization + display filter ----------
    ctrl.jFreqC.setText('10');  ctrl.jFreqW.setText('2');
    panel_bst_dynamics('OnAxisChange', 'freq');  drawnow;
    st = getappdata(0,'DynamicsTarget');
    lf = bst_atom('Get', st.nav, 'freq');
    ok2 = (abs(lf.center-10)<1e-9) && (abs(lf.extent-2)<1e-9) && isequal(st.curBand,[8 12]);
    fprintf('T2 freq block: center=%g extent=%g curBand=%s => %s\n', lf.center, lf.extent, mat2str(st.curBand), PF{ok2+1});
    pass = pass && ok2;

    % ---------- T3: band combobox preset fills the freq fields (alpha -> 10.5 / 2.5) ----------
    ctrl.jFreqBand.setSelectedItem('alpha');  panel_bst_dynamics('OnFreqPreset');  drawnow;
    c = str2double(char(ctrl.jFreqC.getText()));  w = str2double(char(ctrl.jFreqW.getText()));
    ok3 = (abs(c-10.5)<1e-6) && (abs(w-2.5)<1e-6);
    fprintf('T3 band preset: center=%g window=%g => %s\n', c, w, PF{ok3+1});
    pass = pass && ok3;

    % ---------- T4: measurement Psi -> view_helmholtz component + curOp ----------
    panel_bst_dynamics('OnMeasurement', 'Solen');  drawnow;
    st = getappdata(0,'DynamicsTarget');
    St = getappdata(st.hFig, 'HelmholtzState');
    ok4 = strcmp(st.curOp,'Solen') && ~isempty(St) && strcmpi(St.Component,'Solen');
    fprintf('T4 measurement: curOp=%s figComp=%s => %s\n', st.curOp, St.Component, PF{ok4+1});
    pass = pass && ok4;

    if ishandle(hFig), close(hFig); end
    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end

function [linkFile, relData] = i_find_kernel_nav()
    linkFile = '';
    relData = 'Subject01/S01_AEF_20131218_01_notch/data_block001_02.mat';
    [sStudy, ~] = bst_get('DataFile', relData);
    if isempty(sStudy), return; end
    comments = {sStudy.Result.Comment};  fnames = {sStudy.Result.FileName};
    isMN = ~cellfun(@isempty, regexp(comments, 'MN: MEG\(Unconstr\)', 'once')) & ...
           ~cellfun(@isempty, regexp(fnames,   'KERNEL', 'once'));
    for j = find(isMN)
        try
            r = in_bst_results(fnames{j}, 0, 'nComponents','ImagingKernel');
            if (r.nComponents==3) && ~isempty(r.ImagingKernel), linkFile = ['link|' fnames{j} '|' relData];  return; end
        catch
        end
    end
end
```

(Note: `St.Component` is the field `view_helmholtz` stores the current component under — verify the exact field name in `view_helmholtz('SetComponent')` while implementing; if it differs, adjust the T4 assertion to the actual field.)

- [ ] **Step 2: Run test to verify it fails**

Run (MATLAB, Brainstorm in `nogui`):
```matlab
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3/dev'); rehash; test_nav_panel
```
Expected: FAIL — T1 fails because the panel still has `jBands`/`jSpaceStr` and lacks `jFreqC` etc.; `panel_bst_dynamics('OnAxisChange', …)` is an unknown verb.

- [ ] **Step 3: Add the block-builder + cursor/driver functions**

In `toolbox/gui/panel_bst_dynamics.m`, add these functions (place them near the other helpers, e.g. after `i_bands`):

```matlab
%% ===== UNIFORM AXIS BLOCK BUILDER =====
% Symmetric row: [center field] [window field] [right selector slot].
%   axis      'time'|'freq'|'source'|'scale'
%   cLabel    left label for center; wLabel for window
%   rightSel  the axis selector component already created (combobox/toggle) or [] (none)
% Returns the center/window text fields.
function [jC, jW] = i_axis_block(jCtrl, axis, title, cLabel, wLabel, rightSel)
    import java.awt.*;
    BH = java_scaled('value', 22);  FW = java_scaled('value', 52);
    jB = gui_river([2 2], [0 7 2 7], title);
    gui_component('label', jB, '', cLabel, [], [], [], []);
    jC = gui_component('text', jB, '', '', {Dimension(FW,BH)}, ['Center (' axis ')'], []);
    gui_component('label', jB, 'tab', wLabel, [], [], [], []);
    jW = gui_component('text', jB, '', '', {Dimension(FW,BH)}, ['Window/extent (' axis ')'], []);
    if ~isempty(rightSel), jB.add('tab', rightSel); end
    java_setcb(jC, 'ActionPerformedCallback', @(h,e)bst_call(@()OnAxisChange(axis)));
    java_setcb(jW, 'ActionPerformedCallback', @(h,e)bst_call(@()OnAxisChange(axis)));
    jCtrl.add(jB);
end


%% ===== READ a block's (center, extent) into a Localization =====
function loc = i_read_block(ctrl, axis)
    loc = bst_atom('NewLoc', axis);
    switch axis
        case 'time',   jC = ctrl.jTimeC;  jW = ctrl.jTimeW;
        case 'freq',   jC = ctrl.jFreqC;  jW = ctrl.jFreqW;
        case 'source', jC = ctrl.jSrcC;   jW = ctrl.jSrcW;
        case 'scale',  jC = ctrl.jScaleC; jW = ctrl.jScaleW;
        otherwise, return;
    end
    c = str2double(char(jC.getText()));  w = str2double(char(jW.getText()));
    if ~isnan(c), loc.center = c; end
    if ~isnan(w), loc.extent = abs(w); else, loc.extent = 0; end
    if isfinite(loc.center), if loc.extent>0, loc.state='window'; else, loc.state='point'; end; end
end


%% ===== AXIS CHANGE: write the cursor atom + drive the engine =====
function OnAxisChange(axis) %#ok<DEFNU>
    [ctrl, st] = i_cs();
    if isempty(ctrl) || isempty(st), return; end
    loc = i_read_block(ctrl, axis);
    st.nav = bst_atom('Set', st.nav, axis, 1, loc);
    setappdata(0, 'DynamicsTarget', st);
    i_drive(axis, loc);
end


%% ===== FREQ PRESET: the band combobox fills the freq fields, THEN drives =====
% Only the combobox calls this (not the field edits), so typing a custom value is never
% overwritten by the selected band.
function OnFreqPreset() %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'Dynamics');  if isempty(ctrl), return; end
    i_freq_preset(ctrl);
    OnAxisChange('freq');
end


%% ===== per-axis engine driver (reuses existing engines; keeps legacy coords) =====
function i_drive(axis, loc)
    [ctrl, st] = i_cs();  if isempty(st), return; end
    switch axis
        case 'time'
            if isfinite(loc.center), try, panel_time('SetCurrentTime', loc.center); catch, end; end %#ok<CTCH>
        case 'freq'
            if isfinite(loc.center) && (loc.extent>0)
                lo = loc.center - loc.extent;  hi = loc.center + loc.extent;
                panel_filter('SetFilters', 1, hi, 1, lo, 0, [], 0, 1);
                st.curBand = [lo hi];  st.curBandName = i_freq_name(ctrl);
            else
                panel_filter('SetFilters', 0, [], 0, [], 0, [], 0, 0);
                st.curBand = [];  st.curBandName = '';
            end
        case 'source'
            % center/window are populated by the Region tool (Task 2 syncs them); nothing to drive here
        case 'scale'
            if ~isempty(st.hFig) && ishandle(st.hFig) && ~isempty(st.Lambda) && (loc.extent>0)
                params = struct('t', loc.extent);
                try, view_helmholtz('SetSmoothing', st.hFig, 1, 'heat', params); catch, end %#ok<CTCH>
                st.curScale = struct('on',1,'name','heat','params',params);
            elseif ~isempty(st.hFig) && ishandle(st.hFig)
                try, view_helmholtz('SetSmoothing', st.hFig, 0, 'heat', struct('t',1)); catch, end %#ok<CTCH>
                st.curScale = struct('on',0,'name','heat','params',[]);
            end
    end
    setappdata(0, 'DynamicsTarget', st);
end


%% ===== MEASUREMENT (operator descriptor; not an axis) =====
function OnMeasurement(which) %#ok<DEFNU>
    [ctrl, st] = i_cs();
    if isempty(ctrl) || isempty(st) || ~ishandle(st.hFig), return; end
    if strcmp(which, 'Irrot')
        if ctrl.jMeasPot.isSelected(), ctrl.jMeasStr.setSelected(false); name = 'Irrot'; else, name = 'Total'; end
    else
        if ctrl.jMeasStr.isSelected(), ctrl.jMeasPot.setSelected(false); name = 'Solen'; else, name = 'Total'; end
    end
    view_helmholtz('SetComponent', st.hFig, name);
    st.curOp = name;  setappdata(0, 'DynamicsTarget', st);
end


%% ===== frequency-atlas preset: a chosen band fills the freq center/window fields =====
function i_freq_preset(ctrl)
    if ~isfield(ctrl,'jFreqBand') || isempty(ctrl.jFreqBand), return; end
    sel = char(ctrl.jFreqBand.getSelectedItem());
    b = i_bands();  k = find(strcmpi(b(:,1), sel), 1);
    if isempty(k), return; end                                 % 'custom' -> leave fields as typed
    lo = b{k,2}(1);  hi = b{k,2}(2);
    ctrl.jFreqC.setText(num2str((lo+hi)/2));  ctrl.jFreqW.setText(num2str((hi-lo)/2));
end
function nm = i_freq_name(ctrl)
    nm = '';
    if isfield(ctrl,'jFreqBand') && ~isempty(ctrl.jFreqBand), nm = char(ctrl.jFreqBand.getSelectedItem()); end
    if strcmpi(nm,'custom'), nm = ''; end
end
```

- [ ] **Step 4: Rebuild the `jCtrl` control area in `CreatePanel`**

In `toolbox/gui/panel_bst_dynamics.m`, `CreatePanel`, replace the whole control area — from the line `% ===== CONTROL area (above Atoms): the atom-coordinate selectors =====` down through `jCtrl.add(jRec);` (the Frequency + Space + Record sections) — with:

```matlab
    % ===== CONTROL area: the 4-axis (center,extent) navigator =====
    jCtrl = JPanel();  jCtrl.setLayout(BoxLayout(jCtrl, BoxLayout.Y_AXIS));
    BW = java_scaled('value', 30);  BH = java_scaled('value', 22);

    % TIME block (no preset yet)
    [jTimeC, jTimeW] = i_axis_block(jCtrl, 'time', 'Time', 'center', char(177), []);

    % FREQUENCY block + band-atlas preset combobox (right slot)
    bnames = i_bands();  bandItems = [bnames(:,1); {'custom'}];
    jFreqBand = gui_component('combobox', [], [], [], {bandItems}, [], [], []);
    jFreqBand.setSelectedItem('custom');
    java_setcb(jFreqBand, 'ActionPerformedCallback', @(h,e)bst_call(@OnFreqPreset));
    [jFreqC, jFreqW] = i_axis_block(jCtrl, 'freq', 'Frequency', 'center', char(177), jFreqBand);

    % SOURCE block + Region tool (right slot) -- the seed/radius picker
    jRegionTool = gui_component('toggle', [], '', 'Region', {Insets(0,0,0,0), Dimension(java_scaled('value',54),BH)}, 'Heat-disk tool: click a cortex vertex to seed (center), scroll to grow the radius (window)', @(h,e)bst_call(@()bst_geodesic_tool('Toggle', ctrl_region_state())));
    [jSrcC, jSrcW] = i_axis_block(jCtrl, 'source', 'Source', 'center', 'radius', jRegionTool);

    % SCALE block (basic: window -> heat smoothing; center reserved for Phase 5)
    [jScaleC, jScaleW] = i_axis_block(jCtrl, 'scale', 'Scale', 'center', char(177), []);

    % MEASUREMENT row (descriptor, not an axis) + actions
    jMeas = gui_river([2 2], [0 7 2 7], 'Measurement');
    jMeasPot = gui_component('toggle', jMeas, '', char(934), {Insets(0,0,0,0), Dimension(BW,BH)}, 'Potential \Phi (divergence: sources / sinks)', @(h,e)bst_call(@()OnMeasurement('Irrot')));
    jMeasStr = gui_component('toggle', jMeas, '', char(936), {Insets(0,0,0,0), Dimension(BW,BH)}, 'Stream \Psi (curl: vortices)', @(h,e)bst_call(@()OnMeasurement('Solen')));
    gui_component('label', jMeas, 'tab', '  Peaks:', [], [], [], []);
    jPeaks = gui_component('text', jMeas, '', '3', {Dimension(java_scaled('value',26), BH)}, 'Extrema kept per sign', []);
    jCtrl.add(jMeas);

    % ACTIONS row (kept: Detect / Record / Capture)
    jAct = gui_river([2 2], [0 7 2 7], 'Actions');
    gui_component('button', jAct, 'hfill', 'Detect windows', [], 'Run the band-power detector (refphase) on the selected band: writes the band-window stack + phase markers', @(h,e)bst_call(@OnDetect));
    gui_component('button', jAct, 'br hfill', 'Record at cursor', [], 'Store the shaped field''s extrema at the cursor as atoms', @(h,e)bst_call(@OnRecord));
    gui_component('button', jAct, 'br hfill', 'Capture region -> active atom', [], 'Snapshot the Region tool''s heat-disk into the selected atom', @(h,e)bst_call(@OnCaptureRegion));
    jCtrl.add(jAct);
```

Then replace the `BstPanel(...)` controls struct (the final `struct(...)` argument) with the new handles:

```matlab
    bstPanelNew = BstPanel(panelName, jPanelNew, struct( ...
        'jTree',jTree, 'jListOccur',jListOccur, 'jMenuFile',jMenuFile, 'jMenuAtoms',jMenuAtoms, ...
        'jTimeC',jTimeC, 'jTimeW',jTimeW, 'jFreqC',jFreqC, 'jFreqW',jFreqW, 'jFreqBand',jFreqBand, ...
        'jSrcC',jSrcC, 'jSrcW',jSrcW, 'jRegionTool',jRegionTool, 'jScaleC',jScaleC, 'jScaleW',jScaleW, ...
        'jMeasPot',jMeasPot, 'jMeasStr',jMeasStr, 'jPeaks',jPeaks, 'jPhaseItems',jPhaseItems));
```

- [ ] **Step 5: Initialize `st.nav`; remove dead callbacks**

In `SetTarget`, add `nav` to the `DynamicsTarget` struct (a transient cursor group). Find the `setappdata(0, 'DynamicsTarget', struct(...))` call and add the field:
```matlab
        'showPhase',[1 1 1 1], 'nav', bst_dynamics('NewGroup', 'cursor'));
```

Delete the now-dead functions: `OnBand`, `OnSpaceSmooth`, `OnSpaceKernel`, `OnSpaceComp`, and `SetupSpace` (its eigenfilter-slider construction is gone; the scale driver builds heat params directly). Keep `i_cs`, `i_bands`, `i_field`, `ctrl_region_state`, and the `i_space_enable` helper only if still referenced (grep; if unreferenced after the deletions, remove it too). Keep `st.Lambda` read in `SetTarget` (the scale driver needs it) — if it was set inside `SetupSpace`, move that one read (`St = getappdata(hFig,'HelmholtzState'); if ~isempty(St)&&isfield(St,'Lambda'), st.Lambda = St.Lambda; end`) into `SetTarget` directly.

- [ ] **Step 6: Retarget the existing atom-suite steps that used the old handles**

In `dev/test_dynamics_atoms.m`, T4 (lines ~122-123) currently does:
```matlab
    ctrl.jBands(3).doClick();  drawnow;          % alpha band
    ctrl.jSpaceStr.doClick();  drawnow;          % Psi (stream / curl) -> signed
```
Replace with the new-handle equivalents (select alpha via the freq band combobox; Psi via the measurement toggle):
```matlab
    ctrl.jFreqBand.setSelectedItem('alpha');  panel_bst_dynamics('OnAxisChange','freq');  drawnow;   % alpha band
    ctrl.jMeasStr.setSelected(true);  panel_bst_dynamics('OnMeasurement','Solen');  drawnow;          % Psi (stream)
```
(T4 line ~144 `haveBand = ~isempty(st0.curBand)` is unchanged — the freq driver still sets `curBand`.)

- [ ] **Step 7: Run the tests**

Run:
```matlab
rehash; test_nav_panel
test_dynamics_atoms
```
Expected: `test_nav_panel` ends `==== SUITE: PASS ====` (T1 handles, T2 freq block, T3 band preset, T4 measurement); `test_dynamics_atoms` ends `==== SUITE: PASS ====` (T1–T8, with T4 retargeted).

- [ ] **Step 8: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m dev/test_nav_panel.m dev/test_dynamics_atoms.m
git commit -m "feat(dynamics): 4-axis (center,extent) navigator panel via bst_atom + thin drivers"
```

---

### Task 2: Source block live-sync from the geodesic tool

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (sync the Source center/window fields from `bst_geodesic_tool('GetState')`)
- Modify: `toolbox/dynamics/bst_geodesic_tool.m` (notify the panel after a seed/grow, if a Dynamics panel is open)
- Test: `dev/test_nav_panel.m` (add T5)

**Interfaces:**
- Consumes: `bst_geodesic_tool('GetState')` → `struct(seed,pos,radius,vertices,SurfaceFile)`.
- Produces: `panel_bst_dynamics('SyncSource')` — reads `GetState` and writes the Source center (seed vertex id) + window (radius in mm) fields, and updates `st.nav` source Localization.

- [ ] **Step 1: Write the failing test**

In `dev/test_nav_panel.m`, before the final `if ishandle(hFig)` cleanup, add T5:

```matlab
    % ---------- T5: Source block syncs from the geodesic tool ----------
    st = getappdata(0,'DynamicsTarget');  surf = st.T.SurfaceFile;
    if isempty(surf), rs = in_bst_results(R,0,'SurfaceFile');  surf = rs.SurfaceFile; end
    vi = round(size(in_tess_bst(surf,0).Vertices,1)/3);
    bst_geodesic_tool('Seed', surf, vi);
    panel_bst_dynamics('SyncSource');  drawnow;
    st = getappdata(0,'DynamicsTarget');
    cv = str2double(char(ctrl.jSrcC.getText()));
    rw = str2double(char(ctrl.jSrcW.getText()));
    ls = bst_atom('Get', st.nav, 'source');
    ok5 = (cv==vi) && (rw>0) && (ls.center==vi);
    fprintf('T5 source sync: center=%g radius=%g navSeed=%g => %s\n', cv, rw, ls.center, PF{ok5+1});
    pass = pass && ok5;
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```matlab
rehash; test_nav_panel
```
Expected: FAIL at T5 — `panel_bst_dynamics('SyncSource')` is an unknown verb; the Source fields stay empty.

- [ ] **Step 3: Add `SyncSource` to the panel**

In `toolbox/gui/panel_bst_dynamics.m`, add:

```matlab
%% ===== SYNC the Source block fields from the geodesic tool state =====
function SyncSource() %#ok<DEFNU>
    [ctrl, st] = i_cs();
    if isempty(ctrl) || isempty(st) || ~isfield(ctrl,'jSrcC'), return; end
    gs = bst_geodesic_tool('GetState');
    if isempty(gs), return; end
    ctrl.jSrcC.setText(num2str(double(gs.seed)));
    ctrl.jSrcW.setText(num2str(round(gs.radius*1000)));      % radius in mm
    loc = bst_atom('NewLoc', 'source');
    loc.center = double(gs.seed);  loc.extent = gs.radius;  loc.pos = gs.pos;  loc.state = 'window';
    st.nav = bst_atom('Set', st.nav, 'source', 1, loc);
    setappdata(0, 'DynamicsTarget', st);
end
```

- [ ] **Step 4: Notify the panel from the geodesic tool after a seed/grow**

In `toolbox/dynamics/bst_geodesic_tool.m`, at the end of `Draw` (the overlay is redrawn on every seed/grow), notify the Dynamics panel if it is loaded — so the Source fields track the live disk without polling. Add after the `patch(...)` call in `Draw`:

```matlab
    % keep the Dynamics Source block fields in sync with the live disk (if the panel is open)
    if ~isempty(bst_get('PanelControls', 'Dynamics'))
        try, panel_bst_dynamics('SyncSource'); catch, end %#ok<CTCH>
    end
```

- [ ] **Step 5: Run the tests**

Run:
```matlab
rehash; test_nav_panel
test_dynamics_atoms
```
Expected: `test_nav_panel` `==== SUITE: PASS ====` (T1–T5); `test_dynamics_atoms` still `==== SUITE: PASS ====` (the `SyncSource` notify in `Draw` is a no-op for the atom suite's region rendering since it guards on the Dynamics panel + GetState).

- [ ] **Step 6: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m toolbox/dynamics/bst_geodesic_tool.m dev/test_nav_panel.m
git commit -m "feat(dynamics): Source navigator block live-syncs seed/radius from bst_geodesic_tool"
```

---

## Self-Review

**1. Spec coverage:**
- Four symmetric blocks (center/window left, selector right) → Task 1 `i_axis_block` + CreatePanel (Step 4).
- `bst_atom` cursor write per axis → `OnAxisChange` (Step 3).
- Drivers reuse engines + keep legacy `curBand`/`curOp`/`curScale` → `i_drive` + `OnMeasurement` (Step 3).
- Freq band-combobox preset fills fields → `i_freq_preset` + `jFreqBand` (Steps 3-4).
- Region tool in the Source right slot → Step 4; live seed/radius sync → Task 2.
- Basic Scale (window→heat smoothing, center reserved) → `i_drive` scale case.
- Measurement descriptor row → `OnMeasurement` + `jMeas*` (Steps 3-4).
- Actions kept (Detect/Record/Capture/Peaks) → Step 4 Actions row.
- T1–T8 retarget → Step 6; new `test_nav_panel` → Steps 1, Task 2 Step 1.
- Out of scope (Navigate/Detect/Save contract, full Scale, non-freq atlases, load-atom-into-blocks) — not in any task. Correct.

**2. Placeholder scan:** none — every code step is complete; every run step has an exact command + expected line. The one verify-while-implementing note (T4's `St.Component` field name) is an explicit instruction to confirm an existing field, not a placeholder.

**3. Type consistency:** the controls handles (`jTimeC/jTimeW/jFreqC/jFreqW/jFreqBand/jSrcC/jSrcW/jRegionTool/jScaleC/jScaleW/jMeasPot/jMeasStr`) are identical in `i_axis_block` returns, the `BstPanel` struct, `i_read_block`, the tests, and `SyncSource`. Verbs `OnAxisChange`/`OnMeasurement`/`SyncSource` match between definition, callbacks, and tests. `st.nav` is a `bst_dynamics('NewGroup')` group written via `bst_atom('Set', st.nav, axis, 1, loc)` and read via `bst_atom('Get', st.nav, axis)` consistently. The legacy `st.curBand`/`curBandName`/`curOp`/`curScale` written by the drivers are the same fields `OnDetect`/`OnRecord` read (unchanged).
