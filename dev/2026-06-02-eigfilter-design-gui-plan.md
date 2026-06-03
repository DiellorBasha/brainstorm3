# Eigenmode Filter-Design GUI Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the existing EigenModes panel so the user can pick a spectral kernel `g(λ)` from `bst_eigfilter`, adjust its parameters live, see `g(λ)` in a companion plot, and watch its spatial point-spread (impulse response) on the cortex — reusing the panel's existing live-render and source-filter plumbing.

**Architecture:** Add a "Kernel" weighting mode to `panel_eigenmodes.m` (alongside Box/Taper/Gauss). Kernel weights come from `bst_eigfilter_kernel`. The panel's existing `RecomputeWeights → NotifyChanged → bst_figures('FireModesChanged') → render` loop is reused. A new "Delta point-spread" display in `view_eigenmodes.m` shows the kernel applied to an impulse. A small companion figure plots `g(λ)`. Because the panel's `ApplyToColumn` already filters source maps, kernels become an analysis knob on real data for free.

**Tech Stack:** MATLAB, Brainstorm Java-Swing `BstPanel` GUI, `bst_eigfilter` library, MATLAB MCP for tests.

---

## Background references (read before starting)

- Spec: `dev/2026-06-02-eigfilter-design-gui-design.md`
- `toolbox/gui/panel_eigenmodes.m` (578 lines) — the panel. Key anchors:
  - `GetState()` default struct (`:376-394`) — add kernel/display fields here.
  - `RecomputeWeights()` (`:411-416`) — add the kernel branch.
  - `BuildWeights(shape,…)` (`:53-92`) — existing pure window weights (do not change).
  - `CreatePanel()` (`:96-146`) — add controls + register in `ctrl`.
  - `RefreshControls()` (`:304-346`), `SetSelectEnabled()` (`:293-300`).
  - `ApplyToColumn(SurfaceFile,u)` (`:551-578`) — add the kernel branch.
- `toolbox/gui/view_eigenmodes.m` (221 lines):
  - `SynthColumn(PairedGrid, W)` (`:66-72`), `ModesChangedCallback(hFig)` (`:189-220`),
    `EigenView` appdata set at `:138`.
- Library: `bst_eigfilter_kernel('list')` / `('info',name)` (meta `.params.<p>.default/.range`,
  `.display`); `bst_eigfilter_evaluate(g, lambdas)` → `[K x 1]`.
  `bst_eigenmodes_filter(Eig, u, M, 'custom', 'TransferFn', g)` → `Φ·diag(g(λ))·Φ'·M·u`.
- Panel subfunctions are callable via the `eval(macro_method)` dispatch, e.g.
  `panel_eigenmodes('KernelPairedWeights', Eig, name, params)`. Mark new pure
  subfunctions `%#ok<DEFNU>`.

**Run a test (MATLAB MCP):** load `mcp__plugin_brainstorm-dev_MATLAB__run_matlab_file`, call on the absolute `.m` path. Tests `addpath('/Users/diellorbasha/workspace/research/code/brainstorm3')` and the `eigfilter` subfolder if needed; a MATLAB+Brainstorm session is running. Success prints `ALL TESTS PASSED`.

**Branch:** before Task 1:
```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git checkout development && git checkout -b feature/eigfilter-design-gui
```

---

## File structure

| File | Change |
|------|--------|
| `toolbox/gui/panel_eigenmodes.m` (modify) | kernel state + pure `KernelPairedWeights`/`KernelSliderRange`/`DeltaPointSpread` subfns; `RecomputeWeights` + `ApplyToColumn` + `GetDisplayColumn` kernel/delta branches; UI controls + callbacks |
| `toolbox/gui/view_eigenmodes.m` (modify) | `ModesChangedCallback` calls `GetDisplayColumn` (synthesis or delta) |
| `toolbox/gui/view_eigfilter_response.m` (new) | small companion figure plotting `g(λ)` |
| `dev/tests/test_eigfilter_design_pure.m` (new) | pure: kernel weights, slider range, delta point-spread character |
| `dev/tests/test_eigfilter_design_smoke.m` (new) | build/drive the extended panel (skips w/o eigenmode surface) |

---

## Task 1: Kernel weighting — state + pure logic

**Files:** Modify `toolbox/gui/panel_eigenmodes.m`. Test: `dev/tests/test_eigfilter_design_pure.m`.

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_eigfilter_design_pure.m`:
```matlab
function test_eigfilter_design_pure
% Pure logic for the eigenmode filter-design panel: kernel->paired-weights,
% slider-range heuristic, and the delta point-spread (impulse response).
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3');
if ~brainstorm('status'); brainstorm nogui; end
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3/toolbox/math/eigfilter');

% Synthetic 2-component eigenmodes (LH verts 1:4, RH verts 5:8), 4 paired ranks.
nV = 8; m = 4;
Vectors = zeros(nV, 2*m);
for k = 1:m
    for n = 1:4; Vectors(n, k)     = cos(pi/4*(n-0.5)*(k-1)); end
    for n = 1:4; Vectors(4+n, m+k) = cos(pi/4*(n-0.5)*(k-1)); end
end
Values   = [1;4;9;16;  1;4;9;16];          % both components share the same spectrum
Eig = struct('Vectors',Vectors,'Values',Values,'nModes',2*m, ...
    'Component',[ones(m,1);2*ones(m,1)],'CompRank',[(1:m)';(1:m)'],'nComponents',2);

% --- KernelPairedWeights: W(k) = g(lambda_paired_k) ---
[W, g] = panel_eigenmodes('KernelPairedWeights', Eig, 'heat', struct('t',0.1));
assert(isequal(size(W),[1 m]), 'W must be 1 x Kpaired.');
lamPaired = [1;4;9;16];
assert(max(abs(W(:) - exp(-0.1*lamPaired))) < 1e-12, 'paired weights must equal g(lambda_paired).');
assert(isa(g,'function_handle'), 'must return the kernel handle.');

% --- KernelSliderRange: log heuristic around the default, clamped to range ---
[smin, smax] = panel_eigenmodes('KernelSliderRange', 0.01, [0 Inf]);
assert(smin > 0 && smin < 0.01 && smax > 0.01, 'range must bracket the default.');

% --- DeltaPointSpread: heat point-spread is smooth, concentrated, finite ---
M = speye(nV);                                   % identity mass for the test
gf = panel_eigenmodes('KernelPairedWeights', Eig, 'heat', struct('t',0.05));
ps = panel_eigenmodes('DeltaPointSpread', Eig, M, gf, 3);   % impulse at vertex 3 (LH)
assert(isequal(size(ps),[nV 1]), 'point-spread must be [nVert x 1].');
assert(all(isfinite(ps)), 'point-spread must be finite.');
assert(abs(ps(3)) == max(abs(ps)), 'heat point-spread must peak at the impulse vertex.');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test, confirm it FAILS** (`Undefined ... KernelPairedWeights`).

- [ ] **Step 3: Add kernel/display fields to `GetState()`**

In `panel_eigenmodes.m` `GetState()` default struct (`:380-391`), add fields before the closing `);`:
```matlab
            'CacheMass',    [], ...
            'WeightMode',   'window', ...    % 'window' | 'kernel'
            'KernelName',   'heat', ...
            'KernelParams', struct('t',0.01), ...
            'KernelFn',     [], ...
            'DisplayMode',  'synthesis', ... % 'synthesis' | 'delta'
            'DeltaVertex',  1);
```
(Move the `);` to the new last line; the previous last field `'CacheMass', []` gets a trailing `, ...`.)

- [ ] **Step 4: Add the pure subfunctions**

Add near `BuildWeights` (after line 92):
```matlab
%% ===== PURE: representative eigenvalue per paired rank =====
function lam = PairedLambda(Eig) %#ok<DEFNU>
    cr  = Eig.CompRank(:);
    val = Eig.Values(:);
    K   = max(cr);
    lam = zeros(K,1);
    for k = 1:K
        lam(k) = mean(val(cr == k));   % L/R near-symmetric -> shared representative lambda
    end
end

%% ===== PURE: paired weights from a spectral kernel =====
function [W, g] = KernelPairedWeights(Eig, name, params) %#ok<DEFNU>
    g   = bst_eigfilter_kernel(name, params);
    lam = PairedLambda(Eig);
    W   = bst_eigfilter_evaluate(g, lam)';    % 1 x Kpaired
end

%% ===== PURE: slider bounds for a kernel parameter (log heuristic) =====
function [smin, smax] = KernelSliderRange(pdef, prange) %#ok<DEFNU>
    if isempty(pdef) || ~isfinite(pdef) || pdef <= 0; pdef = 1; end
    smin = pdef / 100;
    smax = pdef * 100;
    if numel(prange) == 2
        if isfinite(prange(1)); smin = max(smin, prange(1) + eps); end
        if isfinite(prange(2)); smax = min(smax, prange(2)); end
    end
    if smax <= smin; smax = smin * 100; end
end

%% ===== PURE: delta (impulse) point-spread = Phi*diag(g)*Phi'*M*e_i =====
function ps = DeltaPointSpread(Eig, MassMatrix, g, iVertex) %#ok<DEFNU>
    nV = size(Eig.Vectors, 1);
    iVertex = min(max(round(iVertex), 1), nV);
    e = zeros(nV, 1); e(iVertex) = 1;
    ps = bst_eigenmodes_filter(Eig, e, MassMatrix, 'custom', 'TransferFn', g);
end
```

- [ ] **Step 5: Add the kernel branch to `RecomputeWeights()`**

Replace `RecomputeWeights()` (`:411-416`) with:
```matlab
function RecomputeWeights()
    global GlobalData;
    st = GlobalData.UserModes;
    if strcmpi(st.WeightMode, 'kernel')
        % Kernel weighting: W(k)=g(lambda_paired_k); store the raw-lambda handle too.
        if ~EnsureCache(st.SurfaceFile)
            return;
        end
        Eig = GlobalData.UserModes.CacheEig;
        [W, g] = KernelPairedWeights(Eig, st.KernelName, st.KernelParams);
        GlobalData.UserModes.Weights  = W;
        GlobalData.UserModes.KernelFn = g;
    else
        GlobalData.UserModes.Weights = BuildWeights(st.WindowShape, ...
            st.Band(1), st.Band(2), st.iCurrentMode, st.nModes);
        GlobalData.UserModes.KernelFn = [];
    end
end
```

- [ ] **Step 6: Run test, confirm PASS** (`ALL TESTS PASSED`).

- [ ] **Step 7: Lint** `panel_eigenmodes.m` (new issues only).

- [ ] **Step 8: Commit**
```bash
git add toolbox/gui/panel_eigenmodes.m dev/tests/test_eigfilter_design_pure.m
git commit -m "feat(eigfilter-gui): kernel weighting state + pure logic in panel_eigenmodes

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Filtering + display-column dispatch

**Files:** Modify `toolbox/gui/panel_eigenmodes.m`. Test: extend `dev/tests/test_eigfilter_design_pure.m`.

- [ ] **Step 1: Extend the test (failing first)**

Append before `disp('ALL TESTS PASSED')` in `test_eigfilter_design_pure.m`:
```matlab
% --- GetDisplayColumn: synthesis = PairedGrid*W; delta = point-spread ---
% Set up panel state for the synthetic surface (no real DB needed: use a fake file).
global GlobalData;
panel_eigenmodes('ResetState', 'fake_surf.mat', m);
GlobalData.UserModes.CacheSurfaceFile = 'fake_surf.mat';
GlobalData.UserModes.CacheEig  = Eig;
GlobalData.UserModes.CacheMass = speye(nV);
GlobalData.UserModes.WeightMode = 'kernel';
GlobalData.UserModes.KernelName = 'heat';
GlobalData.UserModes.KernelParams = struct('t',0.05);
panel_eigenmodes('RecomputeWeights');   % fills Weights + KernelFn

PairedGrid = zeros(nV, m);
for k = 1:m; PairedGrid(:,k) = sum(Eig.Vectors(:, Eig.CompRank==k), 2); end

GlobalData.UserModes.DisplayMode = 'synthesis';
colS = panel_eigenmodes('GetDisplayColumn', 'fake_surf.mat', PairedGrid);
assert(isequal(size(colS),[nV 1]), 'synthesis column shape.');
assert(max(abs(colS - PairedGrid*GlobalData.UserModes.Weights(:))) < 1e-12, 'synthesis = PairedGrid*W.');

GlobalData.UserModes.DisplayMode = 'delta';
GlobalData.UserModes.DeltaVertex = 6;   % RH vertex
colD = panel_eigenmodes('GetDisplayColumn', 'fake_surf.mat', PairedGrid);
assert(isequal(size(colD),[nV 1]) && all(isfinite(colD)), 'delta column shape/finite.');
assert(abs(colD(6)) == max(abs(colD)), 'delta point-spread peaks at the impulse vertex.');
```

- [ ] **Step 2: Run, confirm FAIL** (`Undefined ... GetDisplayColumn`).

- [ ] **Step 3: Add the kernel branch to `ApplyToColumn()`**

In `ApplyToColumn` (`:551-578`), replace the final two lines (`wRaw = …; uF = bst_eigenmodes_filter(…)`) with:
```matlab
    if strcmpi(st.WeightMode, 'kernel') && ~isempty(GlobalData.UserModes.KernelFn)
        % Kernel: apply g at each RAW mode's lambda (exact), not per paired rank.
        uF = bst_eigenmodes_filter(Eig, u, M, 'custom', 'TransferFn', GlobalData.UserModes.KernelFn);
    else
        wRaw = W(CompRank);                    % expand paired -> raw columns
        uF = bst_eigenmodes_filter(Eig, u, M, 'custom', 'TransferFn', @(l) wRaw(:));
    end
```

- [ ] **Step 4: Add `GetDisplayColumn()`**

Add after `ApplyToColumn` (end of file):
```matlab
%% ===== DISPLAY: column for an eigenmode-view figure (synthesis or delta) =====
function col = GetDisplayColumn(SurfaceFile, PairedGrid) %#ok<DEFNU>
    global GlobalData;
    st = GetState();
    W  = st.Weights;
    if strcmpi(st.DisplayMode, 'delta')
        if ~EnsureCache(SurfaceFile)
            col = PairedGrid * W(:); return;   % fallback if cache missing
        end
        Eig = GlobalData.UserModes.CacheEig;
        M   = GlobalData.UserModes.CacheMass;
        if strcmpi(st.WeightMode, 'kernel') && ~isempty(st.KernelFn)
            tf = st.KernelFn;
        else
            wRaw = W(Eig.CompRank(:));
            tf = @(l) wRaw(:);
        end
        col = DeltaPointSpread(Eig, M, tf, st.DeltaVertex);
    else
        col = PairedGrid * W(:);
    end
end
```

- [ ] **Step 5: Run, confirm PASS.**
- [ ] **Step 6: Lint.**
- [ ] **Step 7: Commit**
```bash
git add toolbox/gui/panel_eigenmodes.m dev/tests/test_eigfilter_design_pure.m
git commit -m "feat(eigfilter-gui): kernel ApplyToColumn + synthesis/delta display dispatch

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: g(λ) companion figure

**Files:** Create `toolbox/gui/view_eigfilter_response.m`.

- [ ] **Step 1: Create the figure helper**

`toolbox/gui/view_eigfilter_response.m` (prepend the standard Brainstorm license header, copied from `toolbox/gui/view_eigenmodes.m`):
```matlab
function hFig = view_eigfilter_response(g, lambdas, titleStr)
% VIEW_EIGFILTER_RESPONSE: Small companion plot of a spectral kernel g(lambda).
% USAGE:  view_eigfilter_response(g, lambdas, titleStr)   % create/update (reused by tag)
%         view_eigfilter_response('close')                % close it
% g        : kernel function handle @(l)
% lambdas  : eigenvalues to mark on the x-axis (real lambda_k)
% titleStr : kernel display name
TAG = 'EigfilterResponse';
if (nargin >= 1) && ischar(g) && strcmpi(g, 'close')
    delete(findall(0, 'Type','figure', 'Tag', TAG));
    hFig = [];
    return;
end
hFig = findall(0, 'Type','figure', 'Tag', TAG);
if isempty(hFig)
    hFig = figure('Tag', TAG, 'Name', 'Kernel response g(\lambda)', ...
        'NumberTitle','off', 'Color','w', 'MenuBar','none', 'ToolBar','none');
    ax = axes('Parent', hFig); hold(ax,'on');
    xlabel(ax, '\lambda (eigenvalue)'); ylabel(ax, 'g(\lambda)'); grid(ax,'on');
    setappdata(hFig, 'Axes', ax);
else
    hFig = hFig(1);
end
ax = getappdata(hFig, 'Axes');
lambdas = double(lambdas(:));
lmax = max([lambdas; eps]);
xg = linspace(0, lmax, 400)';
yg = g(xg);
cla(ax); hold(ax,'on');
% faint ticks at the real eigenvalues
yl = [min([0; yg(:)]), max([1; yg(:)])];
plot(ax, [lambdas lambdas]', repmat(yl(:), 1, numel(lambdas)), '-', 'Color',[.85 .85 .85]);
plot(ax, xg, yg, 'b-', 'LineWidth', 2);
xlim(ax, [0 lmax]); ylim(ax, yl + [0 max(eps, 0.05*range(yl))]);
xlabel(ax, '\lambda (eigenvalue)'); ylabel(ax, 'g(\lambda)'); grid(ax,'on');
if nargin >= 3 && ~isempty(titleStr); title(ax, titleStr, 'Interpreter','none'); end
end
```

- [ ] **Step 2: Smoke-check it builds (MATLAB MCP, throwaway, not committed)**
```matlab
if ~brainstorm('status'); brainstorm nogui; end
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3/toolbox/math/eigfilter');
g = bst_eigfilter_kernel('heat', struct('t',0.05));
h = view_eigfilter_response(g, [1 4 9 16], 'Heat');
disp(class(h)); view_eigfilter_response('close'); disp('RESP_OK');
```
Expected: prints a figure class, then `RESP_OK`, no error.

- [ ] **Step 3: Lint + Commit**
```bash
git add toolbox/gui/view_eigfilter_response.m
git commit -m "feat(eigfilter-gui): g(lambda) companion response figure

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Panel UI controls + callbacks

**Files:** Modify `toolbox/gui/panel_eigenmodes.m`. Test: `dev/tests/test_eigfilter_design_smoke.m`.

This task adds Swing controls following the existing `CreatePanel`/callback patterns in the file. No headless unit test for Swing; verified by a build smoke test (Step 7) + manual.

- [ ] **Step 1: Add controls in `CreatePanel()`**

After the shape radios block (`:131`, before the `jRadioBox.setSelected(1)` line is fine to keep), add a "Kernel" radio to the same `jGroup`, then the kernel controls. Insert after line 131:
```matlab
    jRadioKernel = gui_component('Radio', jPanelNew, '',  'Kernel', jGroup, '', @(h,ev)Shape_Callback('kernel'));
    % Kernel dropdown (populated from the eigfilter registry)
    jLabelKernel = gui_component('Label', jPanelNew, 'br', 'Kernel:');
    kernelNames  = bst_eigfilter_kernel('list');
    jComboKernel = gui_component('ComboBox', jPanelNew, 'hfill', [], {kernelNames}, '', @(h,ev)KernelChanged_Callback());
    % NOTE: confirm the dropdown widget before coding this line. Check
    % toolbox/gui/gui_component.m for the combo/dropdown type, and grep an
    % existing usage (e.g. `grep -rn "'ComboBox'\|JComboBox" toolbox/gui/panel_*.m`)
    % to match the exact arg pattern for the items list + ActionPerformedCallback.
    % If gui_component has no combo type, create it directly:
    %   jComboKernel = javax.swing.JComboBox(kernelNames);
    %   java_setcb(jComboKernel, 'ActionPerformedCallback', @(h,ev)KernelChanged_Callback());
    %   jPanelNew.add('hfill', jComboKernel);
    % Two generic parameter sliders (relabeled per kernel; integer 0..1000 -> log value)
    jLabelP1 = gui_component('Label', jPanelNew, 'br', 'p1');
    jSliderP1 = JSlider(0, 1000, 500);
    java_setcb(jSliderP1, 'StateChangedCallback', @(h,ev)ParamPreview_Callback(1), 'MouseReleasedCallback', @(h,ev)Param_Callback(1));
    jPanelNew.add('hfill', jSliderP1);
    jReadoutP1 = gui_component('Label', jPanelNew, '', '');
    jLabelP2 = gui_component('Label', jPanelNew, 'br', 'p2');
    jSliderP2 = JSlider(0, 1000, 500);
    java_setcb(jSliderP2, 'StateChangedCallback', @(h,ev)ParamPreview_Callback(2), 'MouseReleasedCallback', @(h,ev)Param_Callback(2));
    jPanelNew.add('hfill', jSliderP2);
    jReadoutP2 = gui_component('Label', jPanelNew, '', '');
    % Display mode + delta vertex (eigenmode-view figures)
    jGroupDisp = ButtonGroup();
    jLabelDisp = gui_component('Label', jPanelNew, 'br', 'Show:');
    jRadioSynth = gui_component('Radio', jPanelNew, '', 'Synthesis', jGroupDisp, '', @(h,ev)Display_Callback('synthesis'));
    jRadioDelta = gui_component('Radio', jPanelNew, '', 'Delta',     jGroupDisp, '', @(h,ev)Display_Callback('delta'));
    jRadioSynth.setSelected(1);
    jLabelVtx = gui_component('Label', jPanelNew, 'br', 'Delta vertex:');
    jTextVtx  = gui_component('Text', jPanelNew, '', '1', [], '', @(h,ev)Vertex_Callback());
```
Then register all new handles in the `ctrl = struct(...)` (`:134-144`): add
`'jRadioKernel', jRadioKernel, 'jComboKernel', jComboKernel, 'jLabelP1', jLabelP1, 'jSliderP1', jSliderP1, 'jReadoutP1', jReadoutP1, 'jLabelP2', jLabelP2, 'jSliderP2', jSliderP2, 'jReadoutP2', jReadoutP2, 'jRadioSynth', jRadioSynth, 'jRadioDelta', jRadioDelta, 'jTextVtx', jTextVtx`.

- [ ] **Step 2: Slider value <-> log-parameter mapping helpers**

Add subfunctions (near `KernelSliderRange`):
```matlab
%% ===== map JSlider int (0..1000) to a log-scaled parameter value =====
function v = SliderToValue(islide, smin, smax)
    t = min(max(islide,0),1000) / 1000;
    v = smin * (smax/smin)^t;       % geometric (log) interpolation
end
function islide = ValueToSlider(v, smin, smax)
    v = min(max(v, smin), smax);
    islide = round(1000 * log(v/smin) / log(smax/smin));
end
```

- [ ] **Step 3: Kernel-mode state verbs + callbacks**

Add state verbs (near `SetWindowShape`):
```matlab
%% ===== STATE: set weighting mode (window vs kernel) =====
function SetWeightMode(mode) %#ok<DEFNU>
    global GlobalData;
    GetState();
    GlobalData.UserModes.WeightMode = lower(mode);
    RecomputeWeights();
    NotifyChanged();
end
%% ===== STATE: set kernel name (reset params to meta defaults) =====
function SetKernelName(name) %#ok<DEFNU>
    global GlobalData;
    GetState();
    meta = bst_eigfilter_kernel('info', name);
    p = struct();
    fn = fieldnames(meta.params);
    for i = 1:numel(fn); p.(fn{i}) = meta.params.(fn{i}).default; end
    GlobalData.UserModes.KernelName   = name;
    GlobalData.UserModes.KernelParams = p;
    GlobalData.UserModes.WeightMode   = 'kernel';
    RecomputeWeights();
    NotifyChanged();
end
%% ===== STATE: set one kernel parameter by index (1..2) =====
function SetKernelParam(idx, value) %#ok<DEFNU>
    global GlobalData;
    GetState();
    fn = fieldnames(GlobalData.UserModes.KernelParams);
    if idx <= numel(fn)
        GlobalData.UserModes.KernelParams.(fn{idx}) = value;
        RecomputeWeights();
        NotifyChanged();
    end
end
%% ===== STATE: display mode + delta vertex =====
function SetDisplayMode(mode) %#ok<DEFNU>
    global GlobalData; GetState();
    GlobalData.UserModes.DisplayMode = lower(mode);
    NotifyChanged();
end
function SetDeltaVertex(v) %#ok<DEFNU>
    global GlobalData; GetState();
    GlobalData.UserModes.DeltaVertex = max(round(v),1);
    NotifyChanged();
end
```
Add UI callbacks (near `Shape_Callback`). Extend `Shape_Callback` so 'kernel' sets the weight mode; add the rest:
```matlab
function KernelChanged_Callback()
    ctrl = bst_get('PanelControls', 'EigenModes');
    name = char(ctrl.jComboKernel.getSelectedItem());
    SetKernelName(name);
    RefreshKernelControls();          % relabel/range sliders from meta
    UpdateResponsePlot();
end
function Param_Callback(idx)
    ctrl = bst_get('PanelControls', 'EigenModes');
    [v, ~] = ReadParamSlider(ctrl, idx);
    SetKernelParam(idx, v);
    UpdateResponsePlot();
end
function ParamPreview_Callback(idx)
    ctrl = bst_get('PanelControls', 'EigenModes');
    sl = ctrl.(sprintf('jSliderP%d', idx));
    if ~sl.getValueIsAdjusting(); return; end
    [v, ro] = ReadParamSlider(ctrl, idx);
    ro2 = ctrl.(sprintf('jReadoutP%d', idx));
    ro2.setText(sprintf('%.4g', v)); %#ok<NASGU>
end
function Display_Callback(mode)
    SetDisplayMode(mode);
end
function Vertex_Callback()
    ctrl = bst_get('PanelControls', 'EigenModes');
    v = str2double(char(ctrl.jTextVtx.getText()));
    if ~isnan(v); SetDeltaVertex(v); end
end
```
In `Shape_Callback(shape)` (`:210-218`) add a leading branch:
```matlab
    if strcmpi(shape, 'kernel')
        SetWeightMode('kernel');
        return;
    end
    SetWeightMode('window');     % any window shape leaves kernel mode
```
(insert before the existing width-based logic).

- [ ] **Step 4: Param-slider read + relabel helpers**

```matlab
function [v, ro] = ReadParamSlider(ctrl, idx)
    st = GetState();
    meta = bst_eigfilter_kernel('info', st.KernelName);
    fn = fieldnames(meta.params);
    if idx > numel(fn); v = []; ro = ''; return; end
    pinfo = meta.params.(fn{idx});
    [smin, smax] = KernelSliderRange(pinfo.default, pinfo.range);
    sl = ctrl.(sprintf('jSliderP%d', idx));
    v = SliderToValue(sl.getValue(), smin, smax);
    ro = sprintf('%.4g', v);
end
function RefreshKernelControls()
    ctrl = bst_get('PanelControls', 'EigenModes');
    if isempty(ctrl) || ~isfield(ctrl,'jComboKernel'); return; end
    st = GetState();
    meta = bst_eigfilter_kernel('info', st.KernelName);
    fn = fieldnames(meta.params);
    for idx = 1:2
        L  = ctrl.(sprintf('jLabelP%d', idx));
        S  = ctrl.(sprintf('jSliderP%d', idx));
        RO = ctrl.(sprintf('jReadoutP%d', idx));
        on = (idx <= numel(fn));
        L.setEnabled(on); S.setEnabled(on); RO.setEnabled(on);
        if on
            pinfo = meta.params.(fn{idx});
            [smin, smax] = KernelSliderRange(pinfo.default, pinfo.range);
            L.setText(fn{idx});
            S.setValue(ValueToSlider(st.KernelParams.(fn{idx}), smin, smax));
            RO.setText(sprintf('%.4g', st.KernelParams.(fn{idx})));
        else
            L.setText(sprintf('p%d', idx)); RO.setText('');
        end
    end
end
function UpdateResponsePlot()
    st = GetState();
    if ~strcmpi(st.WeightMode,'kernel') || isempty(st.CacheEig); return; end
    g = st.KernelFn;
    if isempty(g); g = bst_eigfilter_kernel(st.KernelName, st.KernelParams); end
    meta = bst_eigfilter_kernel('info', st.KernelName);
    view_eigfilter_response(g, st.CacheEig.Values(:), meta.display);
end
```

- [ ] **Step 5: Extend `SetSelectEnabled` and `RefreshControls`**

In `SetSelectEnabled` (`:293-300`) add the kernel controls to the `sel` list:
`'jRadioKernel','jComboKernel','jSliderP1','jSliderP2','jRadioSynth','jRadioDelta','jTextVtx'`.
In `RefreshControls` (`:304-346`), after the existing shape-radio sync, call
`RefreshKernelControls();` and set the kernel/display radios from state:
```matlab
    if isfield(ctrl,'jRadioKernel')
        ctrl.jRadioKernel.setSelected(strcmpi(st.WeightMode,'kernel'));
        ctrl.jRadioSynth.setSelected(strcmpi(st.DisplayMode,'synthesis'));
        ctrl.jRadioDelta.setSelected(strcmpi(st.DisplayMode,'delta'));
        ctrl.jTextVtx.setText(num2str(st.DeltaVertex));
    end
```

- [ ] **Step 6: Build smoke test**

Create `dev/tests/test_eigfilter_design_smoke.m`:
```matlab
function test_eigfilter_design_smoke
% Smoke: the extended EigenModes panel builds, and kernel-mode state verbs drive
% weights without error. Skips if no eigenmode surface is loaded.
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3');
if ~brainstorm('status'); brainstorm nogui; end

% Panel builds with the new controls
p = panel_eigenmodes('CreatePanel');
assert(~isempty(p), 'panel must build.');

% Find a surface with eigenmodes (else skip)
sSubjects = bst_get('ProtocolSubjects');
SurfaceFile = '';
if ~isempty(sSubjects)
    for iS = 1:numel(sSubjects.Subject)
        for iC = 1:numel(sSubjects.Subject(iS).Surface)
            sf = sSubjects.Subject(iS).Surface(iC).FileName;
            [~, ok] = in_tess_eigenmodes(sf);
            if ok; SurfaceFile = sf; break; end
        end
        if ~isempty(SurfaceFile); break; end
    end
end
if isempty(SurfaceFile); disp('SKIP: no eigenmode surface.'); return; end

[Eig,~] = in_tess_eigenmodes(SurfaceFile);
Kp = double(max(Eig.CompRank));
panel_eigenmodes('ResetState', SurfaceFile, Kp);
global GlobalData;
GlobalData.UserModes.CacheSurfaceFile = SurfaceFile;
GlobalData.UserModes.CacheEig  = Eig;
GlobalData.UserModes.CacheMass = Eig.MassMatrix;
panel_eigenmodes('SetKernelName', 'heat');       % -> kernel mode, recompute
W = panel_eigenmodes('GetWeights');
assert(numel(W) == Kp && all(isfinite(W)), 'kernel weights must be finite, length Kpaired.');
panel_eigenmodes('SetKernelParam', 1, 0.2);      % adjust t
assert(all(isfinite(panel_eigenmodes('GetWeights'))), 'weights finite after param change.');
% delta display column
PairedGrid = zeros(size(Eig.Vectors,1), Kp);
for k=1:Kp; PairedGrid(:,k)=sum(Eig.Vectors(:,Eig.CompRank==k),2); end
panel_eigenmodes('SetDisplayMode','delta');
panel_eigenmodes('SetDeltaVertex', 1);
col = panel_eigenmodes('GetDisplayColumn', SurfaceFile, PairedGrid);
assert(isequal(size(col),[size(Eig.Vectors,1) 1]) && all(isfinite(col)), 'delta column ok.');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 7: Run the smoke test** (`ALL TESTS PASSED` or `SKIP: …`). Also re-run `dev/tests/test_eigfilter_design_pure.m` (`ALL TESTS PASSED`).
- [ ] **Step 8: Lint + Commit**
```bash
git add toolbox/gui/panel_eigenmodes.m dev/tests/test_eigfilter_design_smoke.m
git commit -m "feat(eigfilter-gui): kernel/param/display controls in EigenModes panel

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Wire the delta display into view_eigenmodes

**Files:** Modify `toolbox/gui/view_eigenmodes.m`.

- [ ] **Step 1: Route the displayed column through the panel**

In `ModesChangedCallback` (`:189-220`), replace the column computation
```matlab
    W = panel_eigenmodes('GetWeights');
    if isempty(W) || (numel(W) ~= size(ev.PairedGrid,2)), return; end
    col = SynthColumn(ev.PairedGrid, W);
```
with:
```matlab
    W = panel_eigenmodes('GetWeights');
    if isempty(W) || (numel(W) ~= size(ev.PairedGrid,2)), return; end
    % The panel decides synthesis vs delta point-spread (and window vs kernel).
    col = panel_eigenmodes('GetDisplayColumn', ev.SurfaceFile, ev.PairedGrid);
```
Leave the rest (write `ImageGridAmp`, `panel_surface('UpdateSurfaceData')`, legend) unchanged. The legend's "Modes/lambda" text is still meaningful in window mode; in kernel/delta mode it harmlessly shows the weight-based summary.

- [ ] **Step 2: Manual/scripted check (MATLAB MCP)** on a loaded eigenmode surface:
```matlab
if ~brainstorm('status'); brainstorm nogui; end
% (use the SurfaceFile from the smoke test if available)
% hFig = view_eigenmodes(SurfaceFile);
% panel_eigenmodes('SetKernelName','heat'); panel_eigenmodes('SetDisplayMode','delta');
% panel_eigenmodes('SetDeltaVertex', 5000);   % cortex should show a heat blob at vertex 5000
```
Expected (manual, with GUI): the cortex shows the kernel's point-spread; dragging `t` broadens it.

- [ ] **Step 3: Lint + Commit**
```bash
git add toolbox/gui/view_eigenmodes.m
git commit -m "feat(eigfilter-gui): view_eigenmodes renders panel's display column (synthesis/delta)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: "Pick on surface" vertex selection + final integration

**Files:** Modify `toolbox/gui/panel_eigenmodes.m`.

- [ ] **Step 1: Investigate the existing vertex-pick mechanism**

Read how scouts get a clicked vertex (search `panel_scout.m` and `figure_3d.m` for a click→nearest-vertex helper, e.g. a `GetNearestVertex`/`SelectVertex`/`CurrentPoint`-based routine). Report the exact callable. If a reusable helper exists, use it; otherwise implement the fallback in Step 2.

- [ ] **Step 2: Add a "Pick on surface" toggle**

In `CreatePanel` add (after the vertex text field):
```matlab
    jTogglePick = gui_component('Toggle', jPanelNew, '', 'Pick', [], 'Click the cortex to place the delta', @(h,ev)Pick_Callback());
```
Register `'jTogglePick', jTogglePick` in `ctrl`. Add the callback:
```matlab
function Pick_Callback()
    ctrl = bst_get('PanelControls', 'EigenModes');
    hFig = bst_figures('GetCurrentFigure', '3D');
    if isempty(hFig) || ~ishandle(hFig); ctrl.jTogglePick.setSelected(0); return; end
    if ctrl.jTogglePick.isSelected()
        setappdata(hFig, 'EigfilterPickBak', get(hFig, 'WindowButtonDownFcn'));
        set(hFig, 'WindowButtonDownFcn', @(h,e) PickVertex(h));
    else
        if ~isempty(getappdata(hFig,'EigfilterPickBak'))
            set(hFig, 'WindowButtonDownFcn', getappdata(hFig,'EigfilterPickBak'));
        end
    end
end
function PickVertex(hFig)
    % Map the click to the nearest cortex vertex and place the delta there.
    st = GetState();
    if isempty(st.CacheEig); return; end
    TessInfo = getappdata(hFig, 'Surface');
    iTess = 1;                                  % cortex patch
    V = get(TessInfo(iTess).hPatch, 'Vertices');
    cp = get(gca, 'CurrentPoint');              % 2x3 click ray
    p = cp(1,:);
    [~, iv] = min(sum((V - p).^2, 2));          % nearest vertex (ray near-point heuristic)
    ctrl = bst_get('PanelControls', 'EigenModes');
    ctrl.jTextVtx.setText(num2str(iv));
    SetDeltaVertex(iv);
end
```
If Step 1 found a precise Brainstorm vertex-pick helper, replace the `CurrentPoint`/nearest-vertex heuristic in `PickVertex` with it (more accurate); keep the rest.

- [ ] **Step 3: Final verification (MATLAB MCP)** — run both tests:
  - `dev/tests/test_eigfilter_design_pure.m` → `ALL TESTS PASSED`
  - `dev/tests/test_eigfilter_design_smoke.m` → `ALL TESTS PASSED` or `SKIP`
  And re-run the eigfilter library tests to confirm no breakage:
  - `dev/tests/test_eigfilter_pure.m` → `ALL TESTS PASSED`

- [ ] **Step 4: Manual GUI checklist** (launch via `/brainstorm-start`):
  1. Right-click a cortex surface with eigenmodes → "View eigenmodes".
  2. In the EigenModes panel, select **Kernel** → pick `heat` → drag its `t` slider; the g(λ) figure updates and (in Synthesis) the cortex changes.
  3. Switch **Show: Delta**, set a vertex (or use **Pick**) → the cortex shows the point-spread; drag `t` to broaden it; try `mexhat`/`dog` for center-surround.
  4. Switch back to **Box/Taper/Gauss** → original window behavior restored.
  5. Open a real source map, make the lever **Active**, select a kernel → the source map is filtered by g(λ).

- [ ] **Step 5: Lint + Commit**
```bash
git add toolbox/gui/panel_eigenmodes.m
git commit -m "feat(eigfilter-gui): pick-on-surface delta vertex selection

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Final verification

- [ ] `test_eigfilter_design_pure.m`, `test_eigfilter_design_smoke.m` pass (smoke may SKIP without an eigenmode surface).
- [ ] `test_eigfilter_pure.m` (library) still passes.
- [ ] Manual checklist (Task 6 Step 4) complete with a real surface.
- [ ] `git diff development -- toolbox/math/` is empty (the library is untouched; the GUI only consumes it).
- [ ] Use `superpowers:finishing-a-development-branch`.
