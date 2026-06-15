# Spatial Filter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** An interactive **Spatial Filter** panel that filters the Dirac source vector field shown on a cortex figure — non-destructive in-place toggle, whole-series swap so time-scrubbing is instant — plus a Save-filtered-file button.

**Architecture:** Factor the Wavelet Designer's *Filter kernel* section (kernel dropdown + mode-index scale sliders) into a shared `bst_eigfilter_panel` helper used by both panels. The new `panel_spatial_filter` is launched from a Dirac source figure's popup, finds-or-creates the surface's Dirac eigenbasis, and on toggle-on filters every time step's spatial field with `bst_dirac_eigenmodes_filter`, swapping the figure's in-memory `ImageGridAmp` (backing up the original; restoring on off/close).

**Tech Stack:** MATLAB (Brainstorm toolbox), Java-Swing GUI (`gui_component`/`gui_show`), the `bst_eigfilter_kernel` registry, `bst_dirac_eigenmodes_filter`, `tess_eigen`, `bst_memory` (`GetDataSetResult`/`GetResultsValues`), `panel_surface` refresh. Tests via the MATLAB MCP against the live dev protocol (Subject01 `cortex_20484V` has a Dirac eigen node).

**Repo / branch:** `brainstorm3` on `feat/spatial-filter`.

**Conventions:** Tests live in `dev/tests/`, are plain functions that `error()` on failure, run with the MATLAB MCP `run_matlab_file`. Brainstorm panel functions dispatch subfunctions via `eval(macro_method)`. Commit after each task with the message in its final step.

---

## File Structure

**Create:**
- `toolbox/gui/bst_eigfilter_panel.m` — shared kernel-section helper (`Kernels`/`CurrentKernel`/`BuildSliders`/`ReadParams`).
- `toolbox/gui/panel_spatial_filter.m` — the Spatial Filter panel (`Start`/`Apply`/`Restore`/`OnToggle`/`OnKernelChanged`/`SaveFiltered`/`Close`).
- `dev/tests/test_eigfilter_panel.m`, `dev/tests/test_spatial_filter.m`.

**Modify:**
- `toolbox/gui/panel_wavelet_designer.m` — build Section 2 via `bst_eigfilter_panel`; remove the now-shared subfunctions.
- `toolbox/gui/figure_3d.m` — add the "Spatial filter (Dirac)" popup item (guarded to a 3-component surface source).

---

## Task 1: Shared kernel-section helper (bst_eigfilter_panel)

**Files:**
- Create: `toolbox/gui/bst_eigfilter_panel.m`
- Test: `dev/tests/test_eigfilter_panel.m`

- [ ] **Step 1: Write the failing test**

```matlab
function test_eigfilter_panel()
% Headless test of the shared kernel-section helper: build sliders into a JPanel,
% read back the kernel name + params for a synthetic Lambda.
% Authors: Diellor Basha, 2026
    nFail = 0;
    Lambda = sort(rand(400,1) * 3e-5);     % synthetic Dirac eigenvalues
    [keys, displays] = bst_eigfilter_panel('Kernels');
    nFail = nFail + chk('lists curated kernels', numel(keys) >= 4 && numel(keys)==numel(displays));
    nFail = nFail + chk('heat is in the list', any(strcmp(keys,'heat')));

    % build a heat section into a river panel, read params back
    jP = gui_river([2 2], [0 2 0 2]);
    bst_eigfilter_panel('BuildSliders', jP, 'heat', Lambda, @() []);
    names = bst_eigfilter_panel('ParamNames', jP);
    nFail = nFail + chk('heat has one param (t)', numel(names)==1 && strcmp(names{1},'t'));
    p = bst_eigfilter_panel('ReadParams', jP, Lambda);
    js = jP.getClientProperty('slider_t'); k = double(js.getValue());
    nFail = nFail + chk('t = 1/lambda_k', abs(p.t - 1/max(Lambda(k),eps)) < 1e-9);

    % dog: two params, ordered t1<t2 regardless of slider positions
    bst_eigfilter_panel('BuildSliders', jP, 'dog', Lambda, @() []);
    names = bst_eigfilter_panel('ParamNames', jP);
    nFail = nFail + chk('dog has t1,t2', numel(names)==2);
    jP.getClientProperty('slider_t1').setValue(300);   % force t1 mode > t2 mode
    jP.getClientProperty('slider_t2').setValue(50);
    p = bst_eigfilter_panel('ReadParams', jP, Lambda);
    nFail = nFail + chk('dog t1 < t2 enforced', p.t1 < p.t2);

    fprintf('\n==== test_eigfilter_panel: %d failed ====\n', nFail);
    if nFail > 0, error('test_eigfilter_panel FAILED'); end
end

function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
```

- [ ] **Step 2: Run it to verify it fails**

Run (MCP `run_matlab_file`): `dev/tests/test_eigfilter_panel.m`
Expected: FAIL — `Undefined function 'bst_eigfilter_panel'`.

- [ ] **Step 3: Implement the helper**

Create `toolbox/gui/bst_eigfilter_panel.m`:

```matlab
function varargout = bst_eigfilter_panel(varargin)
% BST_EIGFILTER_PANEL: Shared "Filter kernel" UI section (kernel dropdown + mode-index
% scale sliders), reused by panel_wavelet_designer and panel_spatial_filter.
%
% API (dispatched via macro_method):
%   [keys, displays] = bst_eigfilter_panel('Kernels')
%   key   = bst_eigfilter_panel('CurrentKernel', jKernel, keys)
%   bst_eigfilter_panel('BuildSliders', jParams, kernelKey, Lambda, onSettle)
%   names = bst_eigfilter_panel('ParamNames', jParams)
%   params = bst_eigfilter_panel('ReadParams', jParams, Lambda)
%
% The scale sliders live in mode-index space (1..K): mode k maps to eigenvalue
% Lambda(k), then to the kernel's scale parameter (t = 1/lambda_k for heat/mexhat/dog,
% beta = lambda_k for tikhonov). dog's t1<t2 is enforced in ReadParams. onSettle is a
% no-arg handle called when a slider drag settles (the owner's live recompute).
%
% Authors: Diellor Basha, 2026
    eval(macro_method);
end

function [keys, displays] = Kernels() %#ok<DEFNU>
    keys = {'mexhat','dog','heat','inverse_heat','tikhonov'};
    displays = cell(1, numel(keys));
    for i = 1:numel(keys)
        try
            m = bst_eigfilter_kernel('info', keys{i});  displays{i} = m.display;
        catch
            displays{i} = keys{i};
        end
    end
end

function key = CurrentKernel(jKernel, keys) %#ok<DEFNU>
    idx = max(1, min(numel(keys), jKernel.getSelectedIndex() + 1));
    key = keys{idx};
end

function BuildSliders(jParams, kernelKey, Lambda, onSettle) %#ok<DEFNU>
    meta = bst_eigfilter_kernel('info', kernelKey);
    K = numel(Lambda);
    jParams.removeAll();
    pf = fieldnames(meta.params);
    nP = numel(pf);
    for i = 1:nP
        nm = pf{i};
        % stagger defaults so multi-scale kernels (dog: t1,t2) start at distinct modes
        defMode = max(1, min(K, round(K * i/(nP+1))));
        [js, jTitle] = i_labeled_slider(jParams, ...
            sprintf('%s: mode %d', i_param_label(nm), defMode), 'coarse', 'fine', 1, K, defMode);
        jParams.putClientProperty(['slider_' nm], js);
        jParams.putClientProperty(['title_'  nm], jTitle);
        java_setcb(js, 'StateChangedCallback', @(h,e) i_slider_changed(jParams, nm, Lambda, h, onSettle));
    end
    jParams.putClientProperty('ParamNames', strjoin(pf(:).', ','));
    jParams.revalidate(); jParams.repaint();
end

function names = ParamNames(jParams) %#ok<DEFNU>
    s = jParams.getClientProperty('ParamNames');
    if isempty(s); names = {}; else; names = strsplit(char(s), ','); end
end

function params = ReadParams(jParams, Lambda) %#ok<DEFNU>
    params = struct();
    K = numel(Lambda);
    names = ParamNames(jParams);
    for i = 1:numel(names)
        nm = names{i};
        js = jParams.getClientProperty(['slider_' nm]);
        if isempty(js); continue; end
        k = max(1, min(K, double(js.getValue())));
        params.(nm) = i_param_value(nm, Lambda(k));
    end
    % dog requires t1 < t2 (the two sliders are independent): order + separate
    if isfield(params,'t1') && isfield(params,'t2')
        lo = min(params.t1, params.t2);  hi = max(params.t1, params.t2);
        if hi <= lo * (1 + 1e-3); hi = lo * 1.5; end
        params.t1 = lo;  params.t2 = hi;
    end
end

%% ===== internal =====
function i_slider_changed(jParams, name, Lambda, js, onSettle)
    jt = jParams.getClientProperty(['title_' name]);
    if ~isempty(jt)
        k = max(1, min(numel(Lambda), double(js.getValue())));
        jt.setText(sprintf('%s: mode %d', i_param_label(name), k));
    end
    if ~js.getValueIsAdjusting() && ~isempty(onSettle); onSettle(); end
end

function v = i_param_value(name, lamk)
    if strcmpi(name, 'beta'); v = max(lamk, eps); else; v = 1 ./ max(lamk, eps); end
end

function lab = i_param_label(name)
    switch lower(name)
        case 't',    lab = 'Scale';
        case 't1',   lab = 'Band edge 1';
        case 't2',   lab = 'Band edge 2';
        case 'beta', lab = 'Scale';
        otherwise,   lab = name;
    end
end

function js = i_slider(jParent, constraints, mn, mx, val)
    import javax.swing.*;
    js = JSlider(mn, mx, val);
    js.setPreferredSize(java_scaled('dimension', 40, 22));
    jParent.add(constraints, js);
end

function [js, jTitle] = i_labeled_slider(jParent, titleText, loLabel, hiLabel, mn, mx, val)
    jTitle = gui_component('label', jParent, 'br', titleText);
    gui_component('label', jParent, 'br', loLabel);
    js = i_slider(jParent, 'hfill', mn, mx, val);
    gui_component('label', jParent, '', hiLabel);
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `dev/tests/test_eigfilter_panel.m`
Expected: PASS — all `PASS`, `0 failed`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/bst_eigfilter_panel.m dev/tests/test_eigfilter_panel.m
git commit -m "feat(gui): bst_eigfilter_panel shared kernel-section helper"
```

---

## Task 2: Refactor the Wavelet Designer onto the shared helper

**Files:**
- Modify: `toolbox/gui/panel_wavelet_designer.m`

The wavelet panel currently has its own kernel-section subfunctions. Replace Section 2's build + reads with the helper, and delete the duplicated subfunctions. Keep `i_slider`/`i_labeled_slider` (still used by the direction/tiles sliders).

- [ ] **Step 1: Build Section 2 via the helper in CreatePanel**

In `panel_wavelet_designer.m` `CreatePanel`, replace the Section 2 block:

```matlab
    % ===== SECTION 2: FILTER KERNEL =====
    [keys, displays] = i_kernel_list();
    jSec2 = gui_river([2 2], [2 8 3 6], '2. Filter kernel');
    gui_component('label', jSec2, 'br', 'Kernel:');
    jKernel = gui_component('combobox', jSec2, 'br hfill', [], {displays}, [], [], []);
    jParams = gui_river([2 2], [0 2 0 2]);        % scale sliders, rebuilt per kernel
    jSec2.add('br hfill', jParams);
    jOpt.add(jSec2);
```
with:
```matlab
    % ===== SECTION 2: FILTER KERNEL (shared helper) =====
    [keys, displays] = bst_eigfilter_panel('Kernels');
    jSec2 = gui_river([2 2], [2 8 3 6], '2. Filter kernel');
    gui_component('label', jSec2, 'br', 'Kernel:');
    jKernel = gui_component('combobox', jSec2, 'br hfill', [], {displays}, [], [], []);
    jParams = gui_river([2 2], [0 2 0 2]);        % scale sliders, rebuilt per kernel
    jSec2.add('br hfill', jParams);
    jOpt.add(jSec2);
```
(Only the `i_kernel_list()` → `bst_eigfilter_panel('Kernels')` line changes here.)

- [ ] **Step 2: Replace the initial slider build + kernel callback wiring**

Find the line that builds the initial sliders (currently `BuildParamWidgets(ctrl, hFig);`) and replace with the helper, wiring slider-settle and kernel-change to the panel's `Refresh`:

```matlab
    % --- build the initial scale sliders for the default kernel (shared helper) ---
    bst_eigfilter_panel('BuildSliders', jParams, keys{1}, S.Lambda, @() Refresh('WaveletDesigner'));
```
And in the callback-wiring block, replace the kernel `ActionPerformedCallback` so it rebuilds via the helper:
```matlab
    java_setcb(jKernel,   'ActionPerformedCallback', @(h,e) OnKernelChanged(panelName));
```
(unchanged line — keep it; `OnKernelChanged` is rewritten in Step 4.)

- [ ] **Step 3: Drop `ParamNames` from the panel state**

In the `S = struct(...)` initializer in `CreatePanel`, remove the `'ParamNames',{{}}` field (the helper now tracks param names on `jParams`). The line:
```matlab
               'SeedCoeffs',[], 'ActiveTile',1, 'iVertex',[], 'Tiles',[], 'ParamNames',{{}});
```
becomes:
```matlab
               'SeedCoeffs',[], 'ActiveTile',1, 'iVertex',[], 'Tiles',[]);
```

- [ ] **Step 4: Rewrite `OnKernelChanged`, `BuildWavelet`'s reads; delete the shared subfunctions**

Replace `OnKernelChanged`:
```matlab
function OnKernelChanged(panelName)
    ctrl = bst_get('PanelControls', panelName);
    if isempty(ctrl); return; end
    S = getappdata(ctrl.hFig, 'WaveletDesignerState');
    key = bst_eigfilter_panel('CurrentKernel', ctrl.jKernel, ctrl.KernelKeys);
    bst_eigfilter_panel('BuildSliders', ctrl.jParams, key, S.Lambda, @() Refresh('WaveletDesigner'));
    Refresh(panelName);
end
```

In `BuildWavelet`, replace the kernel/params reads:
```matlab
function wavelet = BuildWavelet(S, ctrl)
    wavelet = struct('Kernel', bst_eigfilter_panel('CurrentKernel', ctrl.jKernel, ctrl.KernelKeys), ...
                     'Params', bst_eigfilter_panel('ReadParams', ctrl.jParams, S.Lambda), ...
                     'Direction', SeedDirection(S, ctrl), ...
                     'Chirality', 0, 'Axis', [0 0 1]);
    if S.isDirac
        if ctrl.jChirPlus.isSelected();      wavelet.Chirality = +1;
        elseif ctrl.jChirMinus.isSelected(); wavelet.Chirality = -1;
        else;                                wavelet.Chirality = 0;  end
    end
end
```
In `LoadBank`, replace `bst_eigfilter_panel` usages similarly: the kernel set via `ctrl.jKernel.setSelectedIndex`, then `bst_eigfilter_panel('BuildSliders', ctrl.jParams, base.Kernel, S.Lambda, @() Refresh('WaveletDesigner'));`, then set each slider by closest mode (the existing loop, but read param names via `bst_eigfilter_panel('ParamNames', ctrl.jParams)` instead of `S.ParamNames`).

Then **delete** these now-shared subfunctions from `panel_wavelet_designer.m` (they live in the helper):
`i_kernel_list`, `i_current_kernel`, `BuildParamWidgets`, `ReadParams`, `i_param_value`, `i_param_label`, `i_param_label_update`, `OnParamSlider`. **Keep** `i_slider` and `i_labeled_slider` (still used by the direction + tiles sliders).

- [ ] **Step 5: Run the wavelet suite (regression — no behaviour change)**

Run each (MCP `run_matlab_file`):
`dev/tests/test_wavelet_direction.m`, `dev/tests/test_wavelet_designer_session.m`.
Then a live check (MCP `evaluate_matlab_code`):
```matlab
rehash;
try, gui_hide('WaveletDesigner'); catch, end
hOld = findall(0,'Type','figure','Tag','3DViz'); for k=1:numel(hOld); try, bst_figures('DeleteFigure',hOld(k),[]); catch, end; end
EigenFile = bst_get('Subject',1).Surface(5).Eigen(1).FileName;
hFig = view_wavelet_designer(EigenFile); drawnow;
ctrl = bst_get('PanelControls','WaveletDesigner');
panel_wavelet_designer('SetSeedVertex','WaveletDesigner',100); drawnow;
% kernel switch + scale read still work via the helper
for k=1:numel(ctrl.KernelKeys); if strcmp(ctrl.KernelKeys{k},'dog'); ctrl.jKernel.setSelectedIndex(k-1); break; end; end
panel_wavelet_designer('OnKernelChanged','WaveletDesigner'); drawnow;
S = panel_wavelet_designer('GetState','WaveletDesigner');
fprintf('dog params via helper: %s\n', strjoin(bst_eigfilter_panel('ParamNames', ctrl.jParams), ','));
panel_wavelet_designer('OnCancel','WaveletDesigner'); drawnow;
```
Expected: both tests `0 failed`; the live check prints `dog params via helper: t1,t2` with no error.

- [ ] **Step 6: Commit**

```bash
git add toolbox/gui/panel_wavelet_designer.m
git commit -m "refactor(gui): wavelet designer kernel section uses bst_eigfilter_panel"
```

---

## Task 3: panel_spatial_filter — Start (attach, eigenbasis, backup, dock)

**Files:**
- Create: `toolbox/gui/panel_spatial_filter.m`

- [ ] **Step 1: Create the panel skeleton with Start + CreatePanel + GetState**

Create `toolbox/gui/panel_spatial_filter.m`:

```matlab
function varargout = panel_spatial_filter(varargin)
% PANEL_SPATIAL_FILTER: Live in-place spatial filter of a Dirac source map shown on a
% cortex figure. Launched from the figure popup; filters every time step's spatial
% field with bst_dirac_eigenmodes_filter and swaps the in-memory ImageGridAmp
% (non-destructive: the original is restored on toggle-off / close). Controls are the
% shared Filter-kernel section (bst_eigfilter_panel) + a Filter on/off toggle +
% Save filtered file + Close.
%
% Dispatched: Start(hFig), GetState(panelName), Apply/Restore/OnToggle/OnKernelChanged/
%             SaveFiltered/Close.
% Authors: Diellor Basha, 2026
    eval(macro_method);
end

%% ===== LAUNCH (from the source figure popup) =====
function Start(hFig) %#ok<DEFNU>
    global GlobalData;
    panelName = 'SpatialFilter';
    % resolve the displayed source results
    TessInfo = getappdata(hFig, 'Surface');
    iTess = find(arrayfun(@(t) ~isempty(t.DataSource) && strcmpi(t.DataSource.Type,'Source'), TessInfo), 1);
    if isempty(iTess)
        bst_error('No source map on this figure.', 'Spatial filter', 0); return;
    end
    [iDS, iResult] = bst_memory('GetDataSetResult', TessInfo(iTess).DataSource.FileName);
    if isempty(iResult)
        bst_error('Could not resolve the source results.', 'Spatial filter', 0); return;
    end
    R = GlobalData.DataSet(iDS).Results(iResult);
    if isempty(R.nComponents) || (R.nComponents ~= 3)
        bst_error('Spatial filter requires an unconstrained (3-component) source.', 'Spatial filter', 0); return;
    end
    SurfaceFile = R.SurfaceFile;

    % Dirac eigenbasis + operator mass (find-or-create)
    bst_progress('start', 'Spatial filter', 'Loading Dirac eigenbasis...');
    EigenMat = tess_eigen(SurfaceFile, 'Dirac');
    OpMat    = load(file_fullpath(EigenMat.OperatorFile));
    bst_progress('stop');
    nVert = double(max(cellfun(@(x) max(x(:)), EigenMat.GlobalVertices)));

    % materialize the full displayed field [3nVert x nT] and back it up
    J0 = i_full_field(iDS, iResult, nVert);
    if isempty(J0)
        bst_error(sprintf('Field size mismatch (expected 3*%d rows).', nVert), 'Spatial filter', 0); return;
    end

    % session state on the figure
    St = struct('iDS',iDS, 'iResult',iResult, 'hFig',hFig, 'iTess',iTess, ...
                'EigenMat',EigenMat, 'Mass',{OpMat.Mass}, 'Lambda',double(EigenMat.Lambda{1}(:)), ...
                'Orig',J0, 'OrigKernel',R.ImagingKernel, 'isOn',false);
    setappdata(hFig, 'SpatialFilterState', St);

    % build + dock the panel
    gui_hide(panelName);
    bstPanel = CreatePanel(St);
    gui_show(bstPanel, 'BrainstormTab', 'tools');
    try, gui_brainstorm('SetSelectedTab', panelName, 0); catch, end %#ok<CTCH>
    % restore the original if the figure closes
    set(hFig, 'DeleteFcn', @(h,e) Close(panelName));
end

%% ===== PANEL =====
function bstPanelNew = CreatePanel(St)
    import javax.swing.*;
    panelName = 'SpatialFilter';
    jPanelNew = gui_component('Panel');
    jOpt = JPanel(); jOpt.setLayout(BoxLayout(jOpt, BoxLayout.Y_AXIS));

    jSec = gui_river([2 2], [2 8 3 6], 'Spatial filter (Dirac)');
    [keys, displays] = bst_eigfilter_panel('Kernels');
    gui_component('label', jSec, 'br', 'Kernel:');
    jKernel = gui_component('combobox', jSec, 'br hfill', [], {displays}, [], [], []);
    % default to heat (low-pass)
    iHeat = find(strcmp(keys,'heat'),1); if ~isempty(iHeat); jKernel.setSelectedIndex(iHeat-1); end
    jParams = gui_river([2 2], [0 2 0 2]);
    jSec.add('br hfill', jParams);
    jFilterOn = gui_component('checkbox', jSec, 'br', 'Filter on');
    jSave   = gui_component('button', jSec, 'br', 'Save filtered file');
    jClose  = gui_component('button', jSec, '', 'Close');
    jOpt.add(jSec);
    jPanelNew.add(jOpt, java.awt.BorderLayout.NORTH);

    ctrl = struct('jKernel',jKernel, 'KernelKeys',{keys}, 'jParams',jParams, ...
                  'jFilterOn',jFilterOn, 'jSave',jSave, 'jClose',jClose, 'hFig',St.hFig);

    bst_eigfilter_panel('BuildSliders', jParams, char(jKernel.getSelectedItem()), St.Lambda, @() OnKernelOrScale(panelName));
    java_setcb(jKernel,   'ActionPerformedCallback', @(h,e) OnKernelChanged(panelName));
    java_setcb(jFilterOn, 'ActionPerformedCallback', @(h,e) OnToggle(panelName));
    java_setcb(jSave,     'ActionPerformedCallback', @(h,e) SaveFiltered(panelName));
    java_setcb(jClose,    'ActionPerformedCallback', @(h,e) Close(panelName));

    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);
end

function [St, ctrl] = GetState(panelName) %#ok<DEFNU>
    if (nargin < 1) || isempty(panelName); panelName = 'SpatialFilter'; end
    ctrl = bst_get('PanelControls', panelName);
    St = [];
    if isempty(ctrl) || ~isfield(ctrl,'hFig') || ~ishandle(ctrl.hFig); return; end
    St = getappdata(ctrl.hFig, 'SpatialFilterState');
end

%% ===== materialize the full displayed field as [3nVert x nT] =====
function J = i_full_field(iDS, iResult, nVert)
    global GlobalData;
    R = GlobalData.DataSet(iDS).Results(iResult);
    if ~isempty(R.ImageGridAmp)
        J = double(R.ImageGridAmp);
    else
        J = double(bst_memory('GetResultsValues', iDS, iResult, [], [], 0));   % [3nVert x nT], no orient
    end
    if isempty(J) || (size(J,1) ~= 3*nVert); J = []; end
end
```

- [ ] **Step 2: Smoke-check Start opens the panel on a displayed Dirac source**

This needs a displayed unconstrained Dirac source figure. Run (MCP `evaluate_matlab_code`):
```matlab
rehash;
% find a displayed source figure, else skip with a message
hSrc = findall(0,'Type','figure','Tag','3DViz');
hasSrc = false;
for k=1:numel(hSrc)
    TI = getappdata(hSrc(k),'Surface');
    if ~isempty(TI) && any(arrayfun(@(t) ~isempty(t.DataSource) && strcmpi(t.DataSource.Type,'Source'), TI))
        hasSrc = true; hF = hSrc(k); break;
    end
end
if ~hasSrc
    fprintf('SKIP: open a Dirac-dSPM source on cortex first (no source figure found)\n');
else
    panel_spatial_filter('Start', hF); drawnow;
    St = panel_spatial_filter('GetState','SpatialFilter');
    fprintf('panel up: %d ; backed up field [%s] ; isOn=%d\n', ...
        ~isempty(bst_get('PanelControls','SpatialFilter')), num2str(size(St.Orig)), St.isOn);
    panel_spatial_filter('Close','SpatialFilter'); drawnow;
end
```
Expected: with a source figure open, `panel up: 1`, a `[3nVert x nT]` backup, `isOn=0`. (If no source figure is open, it prints SKIP — open one and re-run.)

- [ ] **Step 3: Commit**

```bash
git add toolbox/gui/panel_spatial_filter.m
git commit -m "feat(gui): panel_spatial_filter Start + backup + dock (no apply yet)"
```

---

## Task 4: Apply / Restore (whole-series swap) + on/off toggle

**Files:**
- Modify: `toolbox/gui/panel_spatial_filter.m`
- Test: `dev/tests/test_spatial_filter.m`

- [ ] **Step 1: Write the failing apply/restore test (synthetic dataset)**

```matlab
function test_spatial_filter()
% Apply/Restore the spatial filter on a synthetic in-memory unconstrained source,
% checking the swapped ImageGridAmp equals the Dirac filter and Restore is exact.
% Requires Brainstorm running with a Dirac eigen node (Subject01 surface 5).
% Authors: Diellor Basha, 2026
    nFail = 0;
    EigenFile = bst_get('Subject',1).Surface(5).Eigen(1).FileName;
    EigenMat  = tess_eigen(EigenFile, 'Dirac');
    Op = load(file_fullpath(EigenMat.OperatorFile));
    nVert = double(max(cellfun(@(x) max(x(:)), EigenMat.GlobalVertices)));
    rng(0); J0 = randn(3*nVert, 4);          % synthetic [3nVert x 4 time]

    % reference: filter via the core (heat)
    g = bst_eigfilter_kernel('heat', struct('t', 1/median(EigenMat.Lambda{1})));
    Jref = real(bst_dirac_eigenmodes_filter(EigenMat, Op.Mass, J0, 'custom', 'TransferFn', g));

    % exercise the panel's pure ComputeFiltered against the reference
    St = struct('EigenMat',EigenMat, 'Mass',{Op.Mass}, 'Orig',J0);
    Jf = panel_spatial_filter('ComputeFiltered', St, 'heat', struct('t', 1/median(EigenMat.Lambda{1})));
    nFail = nFail + chk('ComputeFiltered == core filter', max(abs(Jf(:)-Jref(:))) < 1e-9);
    nFail = nFail + chk('filtered differs from original', max(abs(Jf(:)-J0(:))) > 0);

    fprintf('\n==== test_spatial_filter: %d failed ====\n', nFail);
    if nFail > 0, error('test_spatial_filter FAILED'); end
end

function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `dev/tests/test_spatial_filter.m`
Expected: FAIL — `ComputeFiltered` is not a dispatched subfunction yet.

- [ ] **Step 3: Add ComputeFiltered + Apply + Restore + OnToggle + OnKernelChanged**

In `panel_spatial_filter.m`, add:

```matlab
%% ===== PURE: filter the whole series with a kernel =====
function Jf = ComputeFiltered(St, kernelName, params) %#ok<DEFNU>
    g  = bst_eigfilter_kernel(kernelName, params);
    Jf = real(bst_dirac_eigenmodes_filter(St.EigenMat, St.Mass, St.Orig, 'custom', 'TransferFn', g));
end

%% ===== APPLY (swap in the filtered series) =====
function Apply(panelName)
    global GlobalData;
    [St, ctrl] = GetState(panelName);
    if isempty(St); return; end
    name   = bst_eigfilter_panel('CurrentKernel', ctrl.jKernel, ctrl.KernelKeys);
    params = bst_eigfilter_panel('ReadParams', ctrl.jParams, St.Lambda);
    Jf = ComputeFiltered(St, name, params);
    GlobalData.DataSet(St.iDS).Results(St.iResult).ImageGridAmp  = Jf;
    GlobalData.DataSet(St.iDS).Results(St.iResult).ImagingKernel = [];
    St.isOn = true; setappdata(ctrl.hFig, 'SpatialFilterState', St);
    i_refresh(ctrl.hFig, St.iTess);
end

%% ===== RESTORE (put the original back) =====
function Restore(panelName)
    global GlobalData;
    [St, ctrl] = GetState(panelName);
    if isempty(St) || ~ishandle(ctrl.hFig); return; end
    GlobalData.DataSet(St.iDS).Results(St.iResult).ImageGridAmp  = St.Orig;
    GlobalData.DataSet(St.iDS).Results(St.iResult).ImagingKernel = St.OrigKernel;
    St.isOn = false; setappdata(ctrl.hFig, 'SpatialFilterState', St);
    i_refresh(ctrl.hFig, St.iTess);
end

function OnToggle(panelName) %#ok<DEFNU>
    [St, ctrl] = GetState(panelName);
    if isempty(St); return; end
    if ctrl.jFilterOn.isSelected(); Apply(panelName); else; Restore(panelName); end
end

function OnKernelChanged(panelName) %#ok<DEFNU>
    [St, ctrl] = GetState(panelName); %#ok<ASGLU>
    if isempty(ctrl); return; end
    key = bst_eigfilter_panel('CurrentKernel', ctrl.jKernel, ctrl.KernelKeys);
    bst_eigfilter_panel('BuildSliders', ctrl.jParams, key, St.Lambda, @() OnKernelOrScale(panelName));
    OnKernelOrScale(panelName);
end

% kernel/scale changed: re-apply only if the filter is currently on
function OnKernelOrScale(panelName)
    [St, ctrl] = GetState(panelName); %#ok<ASGLU>
    if ~isempty(ctrl) && ctrl.jFilterOn.isSelected(); Apply(panelName); end
end

%% ===== refresh the figure display =====
function i_refresh(hFig, iTess)
    TessInfo = getappdata(hFig, 'Surface');
    for k = 1:numel(TessInfo); TessInfo(k).DataMinMax = []; end
    setappdata(hFig, 'Surface', TessInfo);
    panel_surface('UpdateSurfaceData', hFig);
    panel_surface('UpdateSurfaceColormap', hFig);
    try, figure_3d('SetShowSourceVectors', hFig, iTess, 1); catch, end %#ok<CTCH>
end
```

- [ ] **Step 4: Add Close (restore + hide)**

```matlab
%% ===== CLOSE (restore original, undock) =====
function Close(panelName) %#ok<DEFNU>
    [St, ctrl] = GetState(panelName);
    if ~isempty(St) && ~isempty(ctrl) && ishandle(ctrl.hFig)
        if St.isOn; Restore(panelName); end
        try, set(ctrl.hFig, 'DeleteFcn', ''); catch, end %#ok<CTCH>
        try, rmappdata(ctrl.hFig, 'SpatialFilterState'); catch, end %#ok<CTCH>
    end
    gui_hide(panelName);
end
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `dev/tests/test_spatial_filter.m`
Expected: PASS — both `PASS`, `0 failed`.

- [ ] **Step 6: Commit**

```bash
git add toolbox/gui/panel_spatial_filter.m dev/tests/test_spatial_filter.m
git commit -m "feat(gui): spatial filter apply/restore (whole-series swap) + toggle"
```

---

## Task 5: Save filtered file

**Files:**
- Modify: `toolbox/gui/panel_spatial_filter.m`

- [ ] **Step 1: Add SaveFiltered**

In `panel_spatial_filter.m`, add:

```matlab
%% ===== SAVE the filtered series as a new results node =====
function SaveFiltered(panelName) %#ok<DEFNU>
    global GlobalData;
    [St, ctrl] = GetState(panelName);
    if isempty(St); return; end
    name   = bst_eigfilter_panel('CurrentKernel', ctrl.jKernel, ctrl.KernelKeys);
    params = bst_eigfilter_panel('ReadParams', ctrl.jParams, St.Lambda);
    Jf = ComputeFiltered(St, name, params);

    % copy the source results metadata, replace the data with the filtered series
    srcFile = GlobalData.DataSet(St.iDS).Results(St.iResult).FileName;
    R = in_bst_results(srcFile, 0);                 % full struct (no kernel reconstruction)
    R.ImageGridAmp  = Jf;
    R.ImagingKernel = [];
    R.Comment       = [R.Comment ' | dirac filter(' name ')'];
    if isfield(R,'History'); R.History(end+1,1:3) = {datestr(now,'dd-mmm-yyyy'), 'filter', ['Spatial filter: dirac ' name]}; end

    % save into the same study as the source
    [sStudy, iStudy] = bst_get('AnyFile', srcFile);
    ProtocolInfo = bst_get('ProtocolInfo');
    c = clock; strTime = sprintf('%02.0f%02.0f%02.0f_%02.0f%02.0f', c(1)-2000, c(2:5));
    OutFile = bst_fullfile(ProtocolInfo.STUDIES, bst_fileparts(sStudy.FileName), ['results_diracfilt_' strTime '.mat']);
    OutFile = file_unique(OutFile);
    bst_save(OutFile, R, 'v6');
    db_add_data(iStudy, OutFile, R);
    panel_protocols('UpdateNode', 'Study', iStudy);
    bst_progress('text', 'Saved filtered source map.');
end
```

- [ ] **Step 2: Verify the save round-trips (live, if a source figure is open)**

Run (MCP `evaluate_matlab_code`):
```matlab
rehash;
hSrc = findall(0,'Type','figure','Tag','3DViz'); hF=[];
for k=1:numel(hSrc); TI=getappdata(hSrc(k),'Surface'); if ~isempty(TI) && any(arrayfun(@(t) ~isempty(t.DataSource)&&strcmpi(t.DataSource.Type,'Source'),TI)); hF=hSrc(k); break; end; end
if isempty(hF); fprintf('SKIP: no source figure open\n'); else
  panel_spatial_filter('Start', hF); drawnow;
  ctrl = bst_get('PanelControls','SpatialFilter');
  [~,iStudy] = bst_get('AnyFile', getappdata(hF,'ResultsFile'));
  n0 = numel(bst_get('Study',iStudy).Result);
  panel_spatial_filter('SaveFiltered','SpatialFilter'); drawnow;
  n1 = numel(bst_get('Study',iStudy).Result);
  fprintf('results nodes %d -> %d (expect +1)\n', n0, n1);
  panel_spatial_filter('Close','SpatialFilter'); drawnow;
end
```
Expected: `results nodes N -> N+1` (a new filtered node), or SKIP if no source figure.

- [ ] **Step 3: Commit**

```bash
git add toolbox/gui/panel_spatial_filter.m
git commit -m "feat(gui): spatial filter 'Save filtered file' -> new results node"
```

---

## Task 6: figure_3d popup item (guarded) + end-to-end

**Files:**
- Modify: `toolbox/gui/figure_3d.m`

- [ ] **Step 1: Add the guarded popup item**

In `figure_3d.m` `DisplayFigurePopup`, after the `ResultsFile`/`TessInfo` are fetched (around the source-results menu region, ~line 1685 where "View sources" lives), add:

```matlab
        % === SPATIAL FILTER (unconstrained Dirac source) ===
        if ~isempty(ResultsFile)
            iTessSrc = find(arrayfun(@(t) ~isempty(t.DataSource) && strcmpi(t.DataSource.Type,'Source'), TessInfo), 1);
            if ~isempty(iTessSrc)
                [iDSr, iResr] = bst_memory('GetDataSetResult', TessInfo(iTessSrc).DataSource.FileName);
                if ~isempty(iResr) && isequal(GlobalData.DataSet(iDSr).Results(iResr).nComponents, 3)
                    gui_component('MenuItem', jPopup, [], 'Spatial filter (Dirac)', IconLoader.ICON_RESULTS, [], @(h,ev)bst_call(@panel_spatial_filter, 'Start', hFig));
                end
            end
        end
```
(`GlobalData` is already global in this function; `ResultsFile`/`TessInfo` already defined above.)

- [ ] **Step 2: End-to-end live check (if a source figure is open)**

Run (MCP `evaluate_matlab_code`):
```matlab
rehash;
hSrc = findall(0,'Type','figure','Tag','3DViz'); hF=[];
for k=1:numel(hSrc); TI=getappdata(hSrc(k),'Surface'); if ~isempty(TI) && any(arrayfun(@(t) ~isempty(t.DataSource)&&strcmpi(t.DataSource.Type,'Source'),TI)); hF=hSrc(k); break; end; end
if isempty(hF); fprintf('SKIP: no source figure open\n'); else
  global GlobalData;
  panel_spatial_filter('Start', hF); drawnow;
  ctrl = bst_get('PanelControls','SpatialFilter');
  [iDS,iResult] = bst_memory('GetDataSetResult', getappdata(hF,'ResultsFile'));
  before = GlobalData.DataSet(iDS).Results(iResult).ImageGridAmp(:,1);
  ctrl.jFilterOn.setSelected(true); panel_spatial_filter('OnToggle','SpatialFilter'); drawnow;
  afterF = GlobalData.DataSet(iDS).Results(iResult).ImageGridAmp(:,1);
  fprintf('filter ON changed frame1: %d\n', max(abs(before-afterF))>0);
  ctrl.jFilterOn.setSelected(false); panel_spatial_filter('OnToggle','SpatialFilter'); drawnow;
  afterR = GlobalData.DataSet(iDS).Results(iResult).ImageGridAmp(:,1);
  fprintf('filter OFF restored frame1 exactly: %d\n', isequal(before, afterR));
  panel_spatial_filter('Close','SpatialFilter'); drawnow;
end
```
Expected: `filter ON changed frame1: 1` and `filter OFF restored frame1 exactly: 1` (or SKIP).

- [ ] **Step 3: Commit**

```bash
git add toolbox/gui/figure_3d.m
git commit -m "feat(gui): 'Spatial filter (Dirac)' figure popup item (3-comp source)"
```

---

## Self-Review

**Spec coverage:**
- Shared kernel-section helper → Task 1; wavelet refactor onto it → Task 2. ✓
- Launch from source figure popup, attach to figure → Tasks 3 (Start) + 6 (popup). ✓
- In-place toggle, non-destructive, whole-series swap, follows time → Task 4 (Apply/Restore swap `ImageGridAmp`; time pipeline reads the swapped series). ✓
- Controls = shared kernel section + Filter on/off + Save + Close, default heat / off → Task 3 (CreatePanel) + Task 4 (toggle). ✓
- Save filtered file → Task 5. ✓
- Requires unconstrained 3-comp Dirac source; eigenbasis find-or-create → Task 3 (Start checks + `tess_eigen`); popup guard → Task 6. ✓
- Teardown restores → Task 4 (`Close`/`Restore`) + Task 3 (`DeleteFcn`). ✓

**Placeholder scan:** No "TBD/TODO". Every code step shows full code; the wavelet refactor (Task 2) lists exact replacements and the exact subfunctions to delete vs keep. ✓

**Type consistency:** `bst_eigfilter_panel` API (`Kernels`/`CurrentKernel`/`BuildSliders`/`ParamNames`/`ReadParams`) is defined in Task 1 and used identically in Tasks 2, 3, 4, 5. `panel_spatial_filter` state fields (`iDS,iResult,hFig,iTess,EigenMat,Mass,Lambda,Orig,OrigKernel,isOn`) set in Task 3 `Start` and read in Tasks 4/5. `ComputeFiltered(St,name,params)` defined in Task 4, used in Tasks 4/5 and the test. `i_full_field`/`i_refresh` defined where first used. `St.Mass` stored as `{OpMat.Mass}` (a 1x2 cell) matches `bst_dirac_eigenmodes_filter`'s `MassCell` arg. ✓

---

## Build order summary

1. Task 1 — `bst_eigfilter_panel` shared helper (+ test)
2. Task 2 — refactor the Wavelet Designer onto it (+ regression)
3. Task 3 — `panel_spatial_filter` Start (attach, eigenbasis, backup, dock)
4. Task 4 — Apply/Restore whole-series swap + toggle (+ test)
5. Task 5 — Save filtered file
6. Task 6 — figure popup item + end-to-end

Tasks 1–2 are the shared-helper refactor (verified by the helper test + the wavelet regression). Tasks 3–6 build the new panel incrementally, each verified live against an open Dirac source figure (the live checks SKIP gracefully if none is open).
