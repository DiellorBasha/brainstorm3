# Helmholtz Derivative-Order Selector — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a Field / Velocity / Acceleration selector to `view_helmholtz` so the cortex colormap + quiver (and detection/tracking) show the k-th time derivative of the source field.

**Architecture:** A pure `bst_time_derivative` finite-differences the last (k+1) frames. `UpdateFrame` fetches those frames, forms `Dᵏ J(t)`, and feeds it into the *existing* pipeline unchanged — so decomposition, colormap, quiver, gate, markers, and tracking all describe the derivative field. A panel radio sets the order.

**Tech Stack:** MATLAB (Brainstorm). Tests: plain functions with `chk(label,cond)` run via the MATLAB MCP.

**Commits:** user-managed on `development`; commit steps OPTIONAL/user-gated.

**Preconditions:** MATLAB up, Brainstorm running, `TutorialAuditory` protocol, toolbox + `dev/tests` on path.

---

### Task 1: `bst_time_derivative` (pure finite difference)

**Files:** Create `toolbox/math/bst_time_derivative.m`; Create `dev/tests/test_time_derivative.m`.

- [ ] **Step 1: Failing test** — create `dev/tests/test_time_derivative.m`:

```matlab
function test_time_derivative()
% Unit tests for bst_time_derivative (backward finite differences).
% Author: Diellor Basha, 2026
    nFail = 0;
    a=[1;2]; b=[4;6]; c=[9;12]; dt=0.5;
    nFail = nFail + chk('order 0 = current frame',     isequal(bst_time_derivative([a b c], dt, 0), c));
    nFail = nFail + chk('order 1 = (b-a)/dt',          isequal(bst_time_derivative([a b],   dt, 1), (b-a)/dt));
    nFail = nFail + chk('order 2 = (c-2b+a)/dt^2',     isequal(bst_time_derivative([a b c], dt, 2), (c-2*b+a)/dt^2));
    nFail = nFail + chk('uses newest cols (extra ok)', isequal(bst_time_derivative([a b c], dt, 1), (c-b)/dt));
    ok = false; try, bst_time_derivative(a, dt, 1); catch, ok = true; end
    nFail = nFail + chk('too few frames errors', ok);
    fprintf('\n==== test_time_derivative: %d failed ====\n', nFail);
    if nFail > 0, error('test_time_derivative FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
```

- [ ] **Step 2: Run, expect FAIL** (`Undefined function 'bst_time_derivative'`). Run `dev/tests/test_time_derivative.m`.

- [ ] **Step 3: Implement** — create `toolbox/math/bst_time_derivative.m`:

```matlab
function D = bst_time_derivative(F, dt, order)
% BST_TIME_DERIVATIVE  Backward finite-difference time derivative of a field.
%   F     : [n x m] consecutive frames, columns oldest -> newest (last col = current).
%           Must have at least order+1 columns; the newest order+1 are used.
%   dt    : time step (s)
%   order : 0 (field) | 1 (velocity) | 2 (acceleration)
%   D     : [n x 1]
% Author: Diellor Basha, 2026
    if size(F,2) < order+1
        error('bst_time_derivative:tooFewFrames', 'need %d frames for order %d (got %d)', order+1, order, size(F,2));
    end
    switch order
        case 0, D = F(:,end);
        case 1, D = (F(:,end) - F(:,end-1)) / dt;
        case 2, D = (F(:,end) - 2*F(:,end-1) + F(:,end-2)) / dt^2;
        otherwise, error('bst_time_derivative:badOrder', 'order must be 0, 1, or 2');
    end
end
```

- [ ] **Step 4: Run, expect PASS** (`test_time_derivative: 0 failed`).

- [ ] **Step 5: Lint** `toolbox/math/bst_time_derivative.m`.

- [ ] **Step 6: Commit (optional)** `git add -A && git commit -m "feat(helmholtz): bst_time_derivative finite-difference helper"`

---

### Task 2: Wire the derivative into `view_helmholtz`

**Files:** Modify `toolbox/gui/view_helmholtz.m`.

- [ ] **Step 1: Register `SetDeriv` in the dispatch.** Add `'SetDeriv'` to both the outer command list and the handle-requiring inner list (the same pattern used for `SetTrack`):
```matlab
    if (nargin >= 1) && ischar(SrcResultsFile) && any(strcmp(SrcResultsFile, {'SetComponent','SetVectors','SetMarkers','SetSmoothing','SetGate','SetTrack','SetDeriv','Close','UpdateFrame'}))
        if any(strcmp(SrcResultsFile, {'SetComponent','SetVectors','SetMarkers','SetSmoothing','SetGate','SetTrack','SetDeriv','UpdateFrame'})) && ...
                (isempty(varargin) || isempty(varargin{1}) || ~all(ishandle(varargin{1})))
            return;
        end
```

- [ ] **Step 2: Add `Deriv` to the state struct.** In the `St = struct(...)` constructor, add `'Deriv',0` (anywhere; e.g. next to `'GateFrac',0`):
```matlab
                'Smooth',struct('on',false,'name','heat','params',struct()), 'GateFrac',0, 'Deriv',0, ...
```

- [ ] **Step 3: Add the `SetDeriv` handler** (next to `SetTrack`):
```matlab
function SetDeriv(hFig, order) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.Deriv  = max(0, min(2, round(order)));
    St.Cache  = containers.Map('KeyType','double','ValueType','any');  % decompositions now stale
    St.Tracks = [];  St.LastIT = [];                                    % derivative field changed
    setappdata(hFig, 'HelmholtzState', St);
    UpdateFrame(hFig);
end
```

- [ ] **Step 4: Replace the frame-fetch block in `UpdateFrame`.** Find (near the top of `UpdateFrame`):
```matlab
    [~, iT] = bst_memory('GetTimeVector', St.srcDS, St.srcResult, 'CurrentTimeIndex');
    if isempty(iT) || iT < 1; iT = 1; end
    Jt = double(bst_memory('GetResultsValues', St.srcDS, St.srcResult, [], iT, 0));
    if size(Jt,1) ~= 3*St.nV; return; end
    % low-pass / band-limit the active frame in the Dirac eigenbasis before decomposing
    if St.Smooth.on
        g  = bst_eigfilter_kernel(St.Smooth.name, St.Smooth.params);
        Jt = real(bst_dirac_eigenmodes_filter(St.EigenMat, St.Mass, Jt, 'custom', 'TransferFn', g));
    end
```
Replace with:
```matlab
    [TimeVec, iT] = bst_memory('GetTimeVector', St.srcDS, St.srcResult, 'CurrentTimeIndex');
    if isempty(iT) || iT < 1; iT = 1; end
    k = St.Deriv;
    if iT <= k                                            % not enough history for this order
        i_blank_display(hFig, St, hAx, k);
        St.Tracks = []; St.LastIT = iT; setappdata(hFig, 'HelmholtzState', St);
        return;
    end
    % fetch the (k+1) consecutive frames, oldest -> newest (last col = current)
    F = zeros(3*St.nV, k+1);
    for j = 0:k
        fj = double(bst_memory('GetResultsValues', St.srcDS, St.srcResult, [], iT-(k-j), 0));
        if size(fj,1) ~= 3*St.nV; return; end
        F(:,j+1) = fj;
    end
    if k > 0, dt = TimeVec(iT) - TimeVec(iT-1); else, dt = 1; end
    Jt = bst_time_derivative(F, dt, k);                   % Dᵏ J(t)  (Field / Velocity / Acceleration)
    % spatial eigenmode smoothing (linear -> applying after the difference is equivalent)
    if St.Smooth.on
        g  = bst_eigfilter_kernel(St.Smooth.name, St.Smooth.params);
        Jt = real(bst_dirac_eigenmodes_filter(St.EigenMat, St.Mass, Jt, 'custom', 'TransferFn', g));
    end
```

- [ ] **Step 5: Add the blank-display helper** (near `i_readout` / `i_count_str`):
```matlab
function i_blank_display(hFig, St, hAx, k)
    TessInfo = getappdata(hFig,'Surface');
    TessInfo(St.iTess).Data       = zeros(St.nV,1);
    TessInfo(St.iTess).DataMinMax = [0 1];
    TessInfo(St.iTess).ColormapType = 'source';
    setappdata(hFig,'Surface',TessInfo);
    panel_surface('UpdateSurfaceColormap', hFig);
    try, figure_3d('SetShowSourceVectors', hFig, St.iTess, 0); catch, end %#ok<CTCH>
    delete(findobj(hAx,'Tag','HelmholtzCore'));
    delete(findobj(hAx,'Tag','HelmholtzTrack'));
    delete(findobj(hAx,'Tag','HelmholtzCentroid'));
    nm = 'velocity'; if k >= 2, nm = 'acceleration'; end
    try, panel_helmholtz('SetReadout', sprintf('%s: needs %d earlier frame(s)', nm, k)); catch, end %#ok<CTCH>
end
```

- [ ] **Step 6: Tag the readout with the order.** In the `if needCores` readout branch at the end of `UpdateFrame`, prefix the order name when k>0. Replace:
```matlab
    if needCores
        i_readout(comp.Kind, mk, Ht);
    else
        try, panel_helmholtz('SetReadout', 'singular points hidden'); catch, end %#ok<CTCH>
    end
```
with:
```matlab
    ordName = {'', 'velocity ', 'acceleration '};  pfx = ordName{St.Deriv+1};
    if needCores
        i_readout(comp.Kind, mk, Ht, pfx);
    else
        try, panel_helmholtz('SetReadout', [pfx 'singular points hidden']); catch, end %#ok<CTCH>
    end
```
And update `i_readout` to accept an optional prefix (default ''):
```matlab
function i_readout(kind, mk, Ht, pfx)
    if nargin < 4, pfx = ''; end
```
and prepend `pfx` to the final `txt` before `SetReadout`, i.e. change the final line to:
```matlab
    try, panel_helmholtz('SetReadout', [pfx txt]); catch, end %#ok<CTCH>
```

- [ ] **Step 7: Lint** `toolbox/gui/view_helmholtz.m` (idioms only).

- [ ] **Step 8: Commit (optional)** `git add -A && git commit -m "feat(helmholtz): compute Dᵏ source field per derivative order"`

---

### Task 3: Panel selector

**Files:** Modify `toolbox/gui/panel_helmholtz.m`.

- [ ] **Step 1: Add the Derivative radio.** After the Component radio block (after `jRadio(1).setSelected(true);`), add:
```matlab
    % --- Derivative order ---
    gui_component('label', jSec, 'br', 'Derivative:');
    dnames = {'Field','Velocity','Acceleration'};
    grpD = ButtonGroup(); jDeriv = javaArray('javax.swing.JRadioButton', numel(dnames));
    for i = 1:numel(dnames)
        jDeriv(i) = gui_component('radio', jSec, 'br', dnames{i});
        grpD.add(jDeriv(i));
        java_setcb(jDeriv(i), 'ActionPerformedCallback', @(h,e) OnDeriv(panelName, i-1));
    end
    jDeriv(1).setSelected(true);   % Field (order 0)
```

- [ ] **Step 2: Add `jDeriv` to the controls struct** (the `ctrl = struct(...)` line):
```matlab
                  'jSmoothOn',jSmoothOn, 'Lambda',Lambda, 'jThresh',jThresh, 'jTrack',jTrack, 'jDeriv',jDeriv);
```

- [ ] **Step 3: Add the callback** (near `OnComponent`):
```matlab
function OnDeriv(panelName, order) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    view_helmholtz('SetDeriv', ctrl.hFig, order);
end
```

- [ ] **Step 4: Verify the panel builds** (MCP `evaluate_matlab_code`):
```matlab
bp = panel_helmholtz('CreatePanel', [], (1:50)'); disp(class(bp));
view_helmholtz('SetDeriv', [], 1);   % no-op without handle -> must not error
fprintf('panel + SetDeriv dispatch OK\n');
```
Expected: `BstPanel`, then `panel + SetDeriv dispatch OK`.

- [ ] **Step 5: Commit (optional)** `git add -A && git commit -m "ui(helmholtz): Field/Velocity/Acceleration selector"`

---

### Task 4: Integration + regression

**Files:** none (verification only)

- [ ] **Step 1: Run all suites** — `test_time_derivative`, `test_dirac_helmholtz`, `test_vortex_track`, `test_vortex_persistence`, `test_helmholtz_track`: all `0 failed`.

- [ ] **Step 2: Live check — open, switch orders, screenshot** (MCP `evaluate_matlab_code`):
```matlab
df='Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat'; [sS,~]=bst_get('DataFile',df);
link=['link|' sS.Result(1).FileName '|' df];
Tv=in_bst_data(df,'Time'); Time=Tv.Time; k0=find(Time>=22.60,1);
hFig = view_helmholtz(link);
view_helmholtz('SetComponent', hFig, 'Total');
panel_time('SetCurrentTime', Time(k0));
% Velocity
view_helmholtz('SetDeriv', hFig, 1);
TI=getappdata(hFig,'Surface'); iT=find(arrayfun(@(t)~isempty(t.DataSource)&&strcmpi(t.DataSource.Type,'Source'),TI),1);
dvel = TI(iT).Data; fprintf('velocity: colormap range [%.2g %.2g], nonzero=%d\n', min(dvel), max(dvel), nnz(dvel)>0);
% Acceleration
view_helmholtz('SetDeriv', hFig, 2);
TI=getappdata(hFig,'Surface'); dacc = TI(iT).Data; fprintf('accel: colormap range [%.2g %.2g]\n', min(dacc), max(dacc));
% insufficient history: jump to frame 1 with acceleration
panel_time('SetCurrentTime', Time(1));
TI=getappdata(hFig,'Surface'); fprintf('frame1 accel blank: maxData=%.2g (expect 0)\n', max(TI(iT).Data));
png='/Users/diellorbasha/workspace/research/code/brainstorm3/dev/benchmarks/helmholtz_velocity_view.png';
panel_time('SetCurrentTime', Time(k0)); view_helmholtz('SetDeriv', hFig, 1); saveas(hFig, png); fprintf('saved %s\n', png);
```
Expected: velocity + acceleration colormaps are non-trivial (nonzero range), frame-1 acceleration blanks to 0, PNG saved. Inspect the PNG (velocity magnitude on cortex + quiver).

- [ ] **Step 3: Final lint sweep** of the new/modified files.

- [ ] **Step 4: Commit (optional)** `git add -A && git commit -m "test(helmholtz): derivative-order integration"`

---

## Self-review notes

- **Spec coverage:** pure finite difference (Task 1) ✓; total-field derivative fed to existing pipeline so detection/colormap/quiver/tracking act on it (Task 2 Step 4) ✓; selector UI (Task 3) ✓; insufficient-history blank (Task 2 Step 5) ✓; cache+track invalidation on order change (Task 2 Step 3) ✓; readout tagging (Task 2 Step 6) ✓; tests unit + integration (Tasks 1, 4) ✓.
- **Naming consistency:** `bst_time_derivative(F, dt, order)` signature identical in test, helper, and `UpdateFrame`; `St.Deriv` defined (Step 2) and used (Steps 3,4,6); `i_blank_display(hFig,St,hAx,k)` defined (Step 5) and called (Step 4); `i_readout(...,pfx)` optional arg is backward-compatible (other callers pass 3 args).
- **Placeholder scan:** none.
- **Risk:** `GetResultsValues` per previous frame is one extra kernel*data column per order (≤2 extra) — negligible; the heavier `Frame` decomposition path is unchanged and still skips detection when markers/track are off.
```
