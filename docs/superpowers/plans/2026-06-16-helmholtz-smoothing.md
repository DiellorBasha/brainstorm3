# Helmholtz + Smoothing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Fold the Dirac eigenmode smoothing (Spatial Filter) into the Helmholtz panel as an active-frame control, plus a marker magnitude gate, so vortex-core tracking is tunable. Remove the standalone Spatial Filter.

**Architecture:** `view_helmholtz` loads the Dirac eigenbasis at launch and low-passes the active frame's field (`bst_dirac_eigenmodes_filter`) before decomposing; the decomposition cache is keyed by frame and cleared on any smoothing change. Markers are gated by a fraction of the per-frame `max|omega|`. `panel_helmholtz` gains the shared `bst_eigfilter_panel` smoothing section + an on/off + a threshold slider.

**Tech Stack:** MATLAB (Brainstorm); `tess_eigen('Dirac')` eigenbasis + `bst_dirac_eigenmodes_filter` + `bst_eigfilter_kernel` + the shared `bst_eigfilter_panel` UI; Java-Swing. Tests via the MATLAB MCP against Subject01 `cortex_pial_low`.

**Repo / branch:** `brainstorm3` on `feat/helmholtz-view`. **Spec:** `docs/superpowers/specs/2026-06-16-helmholtz-smoothing-design.md`.

---

## File Structure

**Modify:**
- `toolbox/gui/view_helmholtz.m` — load eigenbasis; smooth the active frame; gate markers; `SetSmoothing`/`SetGate` dispatch; cache invalidation; pass `Lambda` to the panel.
- `toolbox/gui/panel_helmholtz.m` — smoothing kernel section + on/off + marker-threshold slider; `CreatePanel(hFig, Lambda)`.
- `toolbox/gui/figure_3d.m` — remove the "Spatial filter (Dirac)" popup item.
- `dev/tests/test_helmholtz_view.m` — smoothing/gate/cache checks.

**Remove:**
- `toolbox/gui/panel_spatial_filter.m`, `dev/tests/test_spatial_filter.m`.

---

## Task 1: view_helmholtz — load eigenbasis, smooth the active frame, SetSmoothing/SetGate

**Files:** Modify `toolbox/gui/view_helmholtz.m`

- [ ] **Step 1: Extend the dispatch list (add SetSmoothing/SetGate)**

Replace the dispatch block at the top of `view_helmholtz` with:

```matlab
    if (nargin >= 1) && ischar(SrcResultsFile) && any(strcmp(SrcResultsFile, {'SetComponent','SetVectors','SetMarkers','SetSmoothing','SetGate','Close','UpdateFrame'}))
        if any(strcmp(SrcResultsFile, {'SetComponent','SetVectors','SetMarkers','SetSmoothing','SetGate','UpdateFrame'})) && ...
                (isempty(varargin) || isempty(varargin{1}) || ~all(ishandle(varargin{1})))
            return;
        end
        feval(SrcResultsFile, varargin{:});
        return;
    end
```

- [ ] **Step 2: Load the Dirac eigenbasis at launch + seed smoothing/gate state**

In the launch path, replace the block from `Op = bst_dirac_helmholtz('Prepare', ...)` through the `St = struct(...)` / `setappdata(hFig, 'HelmholtzState', St)` with:

```matlab
    Op = bst_dirac_helmholtz('Prepare', Dirac, LBO, Surf);
    bst_progress('text', 'Loading Dirac eigenbasis...');
    EigenMat = tess_eigen(SurfaceFile, 'Dirac');
    OpMat    = load(file_fullpath(EigenMat.OperatorFile));
    Lambda   = double(EigenMat.Lambda{1}(:));
    bst_progress('stop');

    [hFig, iDSf] = view_surface_data(SurfaceFile, SrcResultsFile, [], 'NewFigure');
    if isempty(hFig); return; end
    iTess = i_find_tess(hFig);

    St = struct('Op',Op, 'srcDS',iDSf, 'srcResult',iResult, 'Component','Total', ...
                'ShowVectors',true, 'ShowMarkers',true, 'iTess',iTess, 'nV',nV, ...
                'EigenMat',EigenMat, 'Mass',{OpMat.Mass}, 'Lambda',Lambda, ...
                'Smooth',struct('on',false,'name','heat','params',struct()), 'GateFrac',0, ...
                'Cache',containers.Map('KeyType','double','ValueType','any'));
    setappdata(hFig, 'HelmholtzState', St);
```

> `'Mass',{OpMat.Mass}` cell-wraps so `St.Mass` is the operator's per-hemisphere mass cell (the form `bst_dirac_eigenmodes_filter` expects), exactly as `panel_spatial_filter` stores it.

- [ ] **Step 3: Pass Lambda when docking the panel**

Change the panel create call (end of launch) from `panel_helmholtz('CreatePanel', hFig)` to:

```matlab
    bstPanel = panel_helmholtz('CreatePanel', hFig, St.Lambda);
```

- [ ] **Step 4: Smooth the active frame inside UpdateFrame**

In `UpdateFrame`, immediately after the `if size(Jt,1) ~= 3*St.nV; return; end` line and before the cache/`Frame` block, insert:

```matlab
    % low-pass / band-limit the active frame in the Dirac eigenbasis before decomposing
    if St.Smooth.on
        g  = bst_eigfilter_kernel(St.Smooth.name, St.Smooth.params);
        Jt = real(bst_dirac_eigenmodes_filter(St.EigenMat, St.Mass, Jt, 'custom', 'TransferFn', g));
    end
```

- [ ] **Step 5: Gate the markers in UpdateFrame**

Replace the markers block in `UpdateFrame` (from `delete(findobj(hAx,'Tag','HelmholtzCore'));` through the `i_readout(comp, Ht);` call) with:

```matlab
    % component markers, pruned by the magnitude gate (fraction of the frame's max |omega|)
    mk = comp.Markers;
    if ~isempty(mk) && St.GateFrac > 0
        om = abs([mk.omega]);  mx = max(om);
        if mx > 0; mk = mk(om >= St.GateFrac * mx); end
    end
    delete(findobj(hAx,'Tag','HelmholtzCore'));
    if St.ShowMarkers && ~isempty(mk)
        V = get(TessInfo(St.iTess).hPatch, 'Vertices');
        for k = 1:numel(mk)
            v = mk(k).iVertex; col = [1 0 0]; if mk(k).charge < 0; col = [0 0 1]; end
            line('Parent',hAx,'XData',V(v,1),'YData',V(v,2),'ZData',V(v,3), 'Marker','o', ...
                'MarkerSize',9,'MarkerFaceColor',col,'MarkerEdgeColor','k','LineStyle','none', ...
                'Tag','HelmholtzCore','Clipping','off');
        end
    end
    i_readout(comp.Kind, mk, Ht);
```

- [ ] **Step 6: Add SetSmoothing/SetGate; update i_readout signature**

Add these two functions after `SetMarkers`:

```matlab
function SetSmoothing(hFig, isOn, name, params) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.Smooth = struct('on',logical(isOn), 'name',name, 'params',params);
    St.Cache  = containers.Map('KeyType','double','ValueType','any');   % decompositions now stale
    setappdata(hFig, 'HelmholtzState', St);  UpdateFrame(hFig);
end
function SetGate(hFig, frac) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.GateFrac = max(0, min(1, frac));  setappdata(hFig, 'HelmholtzState', St);  UpdateFrame(hFig);
end
```

Replace `i_readout` with the gated-markers version:

```matlab
function i_readout(kind, mk, Ht)
    switch kind
        case 'vortex'
            if isempty(mk); txt = '0 vortices, 0 antivortices';
            else; np=sum([mk.charge]>0); nn=sum([mk.charge]<0); txt=sprintf('%d vortices (+), %d antivortices (-), net %+d', np, nn, np-nn); end
        case 'source'
            if isempty(mk); txt = '0 sources, 0 sinks';
            else; np=sum([mk.charge]>0); nn=sum([mk.charge]<0); txt=sprintf('%d sources (+), %d sinks (-), net %+d', np, nn, np-nn); end
        case 'harm'
            txt = sprintf('harmonic energy: %.1f%% of |J|^2', 100*Ht.HarmFrac);
        otherwise
            txt = 'total field |J|';
    end
    try, panel_helmholtz('SetReadout', txt); catch, end %#ok<CTCH>
end
```

- [ ] **Step 7: Parse check**

Run (MCP `evaluate_matlab_code`): `rehash; exist('view_helmholtz','file')` → expect `2`.

- [ ] **Step 8: Commit**

```bash
git add toolbox/gui/view_helmholtz.m
git commit -m "feat(helmholtz): smooth the active frame in the Dirac eigenbasis + marker gate"
```

---

## Task 2: panel_helmholtz — smoothing section + threshold slider

**Files:** Modify `toolbox/gui/panel_helmholtz.m`

- [ ] **Step 1: Replace CreatePanel with the smoothing+gate version**

Replace `CreatePanel` with:

```matlab
function bstPanelNew = CreatePanel(hFig, Lambda) %#ok<DEFNU>
    import javax.swing.*;
    panelName = 'Helmholtz';
    jPanelNew = gui_component('Panel');
    jOpt = JPanel(); jOpt.setLayout(BoxLayout(jOpt, BoxLayout.Y_AXIS));
    jSec = gui_river([2 2], [2 8 3 6], 'Helmholtz / Hodge components');

    % --- Smoothing (Dirac eigenmodes) ---
    gui_component('label', jSec, 'br', 'Smoothing (Dirac eigenmodes):');
    [keys, displays] = bst_eigfilter_panel('Kernels');
    jKernel = gui_component('combobox', jSec, 'br hfill', [], {displays}, [], [], []);
    iHeat = find(strcmp(keys,'heat'),1); if ~isempty(iHeat); jKernel.setSelectedIndex(iHeat-1); end
    jParams = gui_river([2 2], [0 2 0 2]);  jSec.add('br hfill', jParams);
    jSmoothOn = gui_component('checkbox', jSec, 'br', 'Smoothing on');
    bst_eigfilter_panel('BuildSliders', jParams, bst_eigfilter_panel('CurrentKernel', jKernel, keys), Lambda, @() OnSmooth(panelName));
    java_setcb(jKernel,   'ActionPerformedCallback', @(h,e) OnKernel(panelName));
    java_setcb(jSmoothOn, 'ActionPerformedCallback', @(h,e) OnSmooth(panelName));

    % --- Component ---
    gui_component('label', jSec, 'br', 'Component:');
    names  = {'Total','Irrot','Solen','Harm'};
    labels = {'Total field |J|','Irrotational (grad phi)','Solenoidal (curl psi)','Harmonic (h)'};
    grp = ButtonGroup(); jRadio = javaArray('javax.swing.JRadioButton', numel(names));
    for i = 1:numel(names)
        jRadio(i) = gui_component('radio', jSec, 'br', labels{i});
        grp.add(jRadio(i));
        java_setcb(jRadio(i), 'ActionPerformedCallback', @(h,e) OnComponent(panelName, names{i}));
    end
    jRadio(1).setSelected(true);
    jVec  = gui_component('checkbox', jSec, 'br', 'Show vectors');           jVec.setSelected(true);
    jMark = gui_component('checkbox', jSec, 'br', 'Show singular points');   jMark.setSelected(true);
    java_setcb(jVec,  'ActionPerformedCallback', @(h,e) OnVectors(panelName));
    java_setcb(jMark, 'ActionPerformedCallback', @(h,e) OnMarkers(panelName));

    % --- Marker threshold (magnitude gate) ---
    gui_component('label', jSec, 'br', 'Marker threshold:');
    jThresh = JSlider(0, 100, 0);  jThresh.setPreferredSize(java_scaled('dimension', 120, 22));
    jSec.add('br hfill', jThresh);
    java_setcb(jThresh, 'StateChangedCallback', @(h,e) OnGate(panelName));

    jReadout = gui_component('label', jSec, 'br', '');
    jClose  = gui_component('button', jSec, 'br', 'Close');
    java_setcb(jClose, 'ActionPerformedCallback', @(h,e) OnClose(panelName));

    jOpt.add(jSec); jPanelNew.add(jOpt, java.awt.BorderLayout.NORTH);
    ctrl = struct('hFig',hFig, 'jVec',jVec, 'jMark',jMark, 'jReadout',jReadout, ...
                  'jKernel',jKernel, 'KernelKeys',{keys}, 'jParams',jParams, ...
                  'jSmoothOn',jSmoothOn, 'Lambda',Lambda, 'jThresh',jThresh);
    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);
end
```

- [ ] **Step 2: Add the smoothing/gate callbacks**

Add after `OnMarkers`:

```matlab
function OnKernel(panelName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    key = bst_eigfilter_panel('CurrentKernel', ctrl.jKernel, ctrl.KernelKeys);
    bst_eigfilter_panel('BuildSliders', ctrl.jParams, key, ctrl.Lambda, @() OnSmooth(panelName));
    OnSmooth(panelName);
end
function OnSmooth(panelName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    name   = bst_eigfilter_panel('CurrentKernel', ctrl.jKernel, ctrl.KernelKeys);
    params = bst_eigfilter_panel('ReadParams', ctrl.jParams, ctrl.Lambda);
    view_helmholtz('SetSmoothing', ctrl.hFig, ctrl.jSmoothOn.isSelected(), name, params);
end
function OnGate(panelName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    if ctrl.jThresh.getValueIsAdjusting(); return; end
    view_helmholtz('SetGate', ctrl.hFig, double(ctrl.jThresh.getValue())/100);
end
```

Update the file header comment to mention the smoothing section + threshold. `SetReadout`, `OnClose`, `i_valid` are unchanged.

- [ ] **Step 3: Parse check**

Run: `rehash; exist('panel_helmholtz','file')` → expect `2`.

- [ ] **Step 4: Commit**

```bash
git add toolbox/gui/panel_helmholtz.m
git commit -m "feat(helmholtz): panel smoothing section (Dirac eigenmodes) + marker-threshold slider"
```

---

## Task 3: Live test — smoothing thins cores, gate is monotonic, cache invalidates

**Files:** Modify `dev/tests/test_helmholtz_view.m`

- [ ] **Step 1: Add the smoothing/gate checks**

In `test_helmholtz_view`, just before the `% close + stale guard` block, insert:

```matlab
    % --- smoothing: a heat low-pass on the active frame thins the vortex cores ---
    view_helmholtz('SetComponent', hFig, 'Solen'); view_helmholtz('SetMarkers', hFig, true);
    view_helmholtz('SetGate', hFig, 0); drawnow;
    nRaw = numel(findobj(hAx,'Tag','HelmholtzCore'));
    St = getappdata(hFig,'HelmholtzState');  Lam = St.Lambda;
    tt = 1 / Lam(max(1, round(numel(Lam)/3)));               % heat scale at a low-ish mode
    view_helmholtz('SetSmoothing', hFig, true, 'heat', struct('t',tt)); drawnow;
    nSmooth = numel(findobj(hAx,'Tag','HelmholtzCore'));
    nFail = nFail + chk('smoothing thins vortex cores', nSmooth < nRaw);
    St2 = getappdata(hFig,'HelmholtzState');
    nFail = nFail + chk('smoothing change cleared the cache', St2.Cache.Count <= 1);

    % --- magnitude gate is monotonic (more pruning -> fewer markers) ---
    view_helmholtz('SetGate', hFig, 0);   drawnow; g0 = numel(findobj(hAx,'Tag','HelmholtzCore'));
    view_helmholtz('SetGate', hFig, 0.5); drawnow; g5 = numel(findobj(hAx,'Tag','HelmholtzCore'));
    view_helmholtz('SetGate', hFig, 0.95);drawnow; g9 = numel(findobj(hAx,'Tag','HelmholtzCore'));
    nFail = nFail + chk('gate monotonic (g0>=g5>=g9)', (g0>=g5) && (g5>=g9));
    view_helmholtz('SetSmoothing', hFig, false, 'heat', struct('t',tt)); view_helmholtz('SetGate', hFig, 0); drawnow;
```

- [ ] **Step 2: Run the test**

Run: `dev/tests/test_helmholtz_view.m`
Expected: PASS — all prior checks plus the three new ones. If `smoothing thins vortex cores` is flaky on the synthetic random field, relax `<` to `<=`; the intent is "not more cores after low-pass".

- [ ] **Step 3: Commit**

```bash
git add dev/tests/test_helmholtz_view.m
git commit -m "test(helmholtz): smoothing thins cores, gate monotonic, cache invalidation"
```

---

## Task 4: Remove the standalone Spatial Filter

**Files:** Modify `toolbox/gui/figure_3d.m`; remove `toolbox/gui/panel_spatial_filter.m`, `dev/tests/test_spatial_filter.m`

- [ ] **Step 1: Remove the popup item**

In `figure_3d.m` `DisplayFigurePopup`, delete the whole `% === SPATIAL FILTER (unconstrained Dirac source) ===` block (the `if ~isempty(ResultsFile) ... 'Spatial filter (Dirac)' ... end` guarding it).

- [ ] **Step 2: Delete the folded-in files**

```bash
git rm toolbox/gui/panel_spatial_filter.m dev/tests/test_spatial_filter.m
```

- [ ] **Step 3: Verify nothing else references them**

Run: `grep -rn "panel_spatial_filter" toolbox/ dev/` → expect no matches.

- [ ] **Step 4: End-to-end live check (real Dirac source link)**

Run (MCP `evaluate_matlab_code`): open `view_helmholtz` on the real link; assert it opens; turn smoothing on (heat) and confirm Solenoidal vortex-core count drops vs off; close. Expected: prints a lower core count with smoothing on.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/figure_3d.m
git commit -m "refactor(helmholtz): remove the standalone Spatial Filter (folded into Helmholtz)"
```

---

## Self-Review

**Spec coverage:** smoothing section folded into Helmholtz (Task 2) filtering the active frame (Task 1 Step 4) ✓; active-frame-only, cache invalidation on smoothing change (Task 1 Steps 4/6) ✓; marker magnitude gate (Task 1 Step 5 + Task 2 slider) ✓; standalone Spatial Filter removed (Task 4) ✓; no whole-series / save ✓.

**Placeholder scan:** none — every step has complete code.

**Type consistency:** `St.Smooth = struct('on','name','params')` set in launch + `SetSmoothing`, read in `UpdateFrame`. `St.GateFrac` set in launch + `SetGate`, read in `UpdateFrame`. `St.EigenMat/Mass/Lambda` set in launch, read in `UpdateFrame` (filter) + passed to the panel. Dispatch names `{'SetComponent','SetVectors','SetMarkers','SetSmoothing','SetGate','Close','UpdateFrame'}` match the panel callbacks (`OnSmooth`→`SetSmoothing`, `OnGate`→`SetGate`). `i_readout(kind, mk, Ht)` new signature matches its one call site. `panel_helmholtz('CreatePanel', hFig, Lambda)` matches the new 2-arg signature. ✓

---

## Build order

1. Task 1 — view: eigenbasis load + active-frame smoothing + gate + dispatch.
2. Task 2 — panel: smoothing section + threshold slider.
3. Task 3 — live test.
4. Task 4 — remove the standalone Spatial Filter.
