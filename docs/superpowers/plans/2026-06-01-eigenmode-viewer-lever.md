# Eigenmode Viewer ↔ Lever Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Reindex the eigenmode lever to paired rank (hemisphere-symmetric), and make `view_eigenmodes` a client of the lever — synthesizing `Φ·w` (single mode or superposed range) with the panel as its control — instead of hijacking the time stepper. Fix the panel lifecycle so it is inert when no eligible view is in front.

**Architecture:** The lever (`panel_eigenmodes`) selects in paired-rank space; one expansion `w_raw = W(Eig.CompRank)` serves every client. `bst_figures('FireModesChanged')` dispatches by client: eigenmode-view figures → `view_eigenmodes('ModesChangedCallback')` (synthesize), source-map figures → `panel_surface('UpdateSurfaceData')` (filter). `view_eigenmodes` caches its `PairedGrid` in figure appdata and owns its synthesis.

**Tech Stack:** MATLAB, Brainstorm; `eval(macro_method)` dispatch; `GlobalData`. Tests run via MATLAB MCP tool `mcp__plugin_brainstorm-dev_MATLAB__evaluate_matlab_code` (function-style, prints `ALL TESTS PASSED`).

**Spec:** `docs/superpowers/specs/2026-06-01-eigenmode-viewer-lever-design.md`

## Verified facts
- `in_tess_eigenmodes` always returns `Eig.CompRank` (backfilled `(1:nModes)'` for legacy), `Eig.Component`, `Eig.Vectors [nV×nModes]`, `Eig.Values`, `Eig.nModes`, `Eig.MassType`, `Eig.MassMatrix`.
- The merged lever stores weights over **raw columns** today: `ApplyToColumn` builds `@(l) W(:)` and guards `numel(W) ~= Eig.nModes`; `UpdatePanel` sets `K = Eig.nModes`.
- `view_eigenmodes` has a pure `BuildPairedGrid(Eig)` → `[Grid [nV×K_paired], K_paired, Info]` where `K_paired = max(CompRank)` and column k sums each component's rank-k mode (disjoint support).
- `bst_figures('FireModesChanged')` (merged) loops visible `3DViz` figures matching `GlobalData.UserModes.SurfaceFile` and calls `panel_surface('UpdateSurfaceData', hFig)`.
- `[iDS,iResult] = bst_memory('GetDataSetResult', ResultsFile)` resolves a loaded results file; its values live at `GlobalData.DataSet(iDS).Results(iResult).ImageGridAmp`.
- `bst_figures('GetCurrentFigure','3D')` returns the current 3D figure (or `[]`).

---

## Task 1: Reindex the lever to paired rank

Make `ApplyToColumn` expand paired weights to raw columns via `CompRank`, and `UpdatePanel` use `K_paired = max(CompRank)`. This fixes the hemisphere asymmetry for the filter client. The single-component integration test still passes unchanged (identity expansion).

**Files:** Modify `toolbox/gui/panel_eigenmodes.m`; create `dev/tests/test_eigenmode_lever_paired.m`.

- [ ] **Step 1: Write the failing test** — `dev/tests/test_eigenmode_lever_paired.m`:

```matlab
function test_eigenmode_lever_paired
% Paired-rank reindex: weights are over paired rank and expand to raw columns
% via CompRank; a low paired-band keeps BOTH components (hemisphere symmetry).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(fullfile(repoRoot, 'toolbox', 'gui'));
addpath(fullfile(repoRoot, 'toolbox', 'math'));
addpath(fullfile(repoRoot, 'toolbox', 'anatomy'));
global GlobalData; %#ok<GVMIS>

% Build a TWO-component eigenbasis by block-diagonalizing two single spheres.
[V1,F1] = tess_sphere(162);
[E1,~,M1] = tess_eigenmodes(V1, F1, 'nModes', 20, 'MassType','barycentric','RemoveDC',1,'Verbose',0);
[V2,F2] = tess_sphere(162);
[E2,~,M2] = tess_eigenmodes(V2, F2, 'nModes', 20, 'MassType','barycentric','RemoveDC',1,'Verbose',0);
K1 = E1.nModes; K2 = E2.nModes; Kp = min(K1,K2);
% Truncate both to Kp ranks so paired ranks line up 1..Kp on each component
P1 = E1.Vectors(:,1:Kp); P2 = E2.Vectors(:,1:Kp);
nV1 = size(P1,1); nV2 = size(P2,1);
Vectors = [ P1, zeros(nV1,Kp); zeros(nV2,Kp), P2 ];     % disjoint support
Values  = [ E1.Values(1:Kp); E2.Values(1:Kp) ];
Comp    = [ ones(Kp,1); 2*ones(Kp,1) ];
CompRank= [ (1:Kp)'; (1:Kp)' ];
Eig = struct('Vectors',Vectors,'Values',Values,'nModes',2*Kp, ...
             'Component',Comp,'CompRank',CompRank,'MassType','barycentric');
M = blkdiag(M1, M2);

SurfaceFile = '/synthetic/twocomp.mat';
panel_eigenmodes('ResetState', SurfaceFile, Kp);   % K_paired = Kp
panel_eigenmodes('SetCache', SurfaceFile, Eig, M);

% A field with energy on BOTH components
rng(7); u = Vectors * randn(2*Kp, 1);
iC1 = 1:nV1; iC2 = nV1+(1:nV2);

% Low paired-band [1,3], box -> keep rank 1..3 of BOTH components
panel_eigenmodes('SetActive', 1);
panel_eigenmodes('SetWindowShape', 'box');
panel_eigenmodes('SetBand', 1, 3);
uF = panel_eigenmodes('ApplyToColumn', SurfaceFile, u);

e1 = norm(uF(iC1)); e2 = norm(uF(iC2));
assert(e1 > 1e-6 && e2 > 1e-6, 'low paired-band must keep BOTH components (symmetry)');
assert(abs(e1 - e2)/max(e1,e2) < 0.6, 'kept energy should be comparable across components');

% Analytic: w_raw = W(CompRank); uF == Vectors*(w_raw.*(Vectors'*M*u))
W = panel_eigenmodes('GetWeights');           % length Kp (paired)
assert(numel(W) == Kp, 'weights are paired-length');
wRaw = W(CompRank);
analytic = Vectors * (wRaw(:) .* (Vectors' * (M * u)));
assert(max(abs(uF - analytic)) < 1e-9, 'ApplyToColumn must expand via CompRank');

fprintf('ALL TESTS PASSED: test_eigenmode_lever_paired\n');
end
```

- [ ] **Step 2: Run it — verify it FAILS** (current `ApplyToColumn` uses raw `W`, guard `numel(W)~=Eig.nModes` → `numel(W)=Kp ≠ 2*Kp` → returns `u` unchanged → `e2` assertion or analytic assertion fails).
  Run: `run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_paired.m')`

- [ ] **Step 3: Reindex `ApplyToColumn`** in `toolbox/gui/panel_eigenmodes.m`. Replace the weight guard + filter call so weights are paired and expand via `CompRank`:

```matlab
    CompRank = Eig.CompRank(:);
    Kpaired  = max(CompRank);
    W = GlobalData.UserModes.Weights;
    if isempty(W) || (numel(W) ~= Kpaired)
        return;
    end
    wRaw = W(CompRank);                         % expand paired -> raw columns
    % Reconstruct via the core spectral filter (custom transfer = expanded weights)
    uF = bst_eigenmodes_filter(Eig, u, M, 'custom', 'TransferFn', @(l) wRaw(:));
```
(Keep all the earlier guards — inactive / surface mismatch / scalar-field — unchanged; only the weight-length guard and the transfer-fn expansion change.)

- [ ] **Step 4: Reindex `UpdatePanel`** — change `K = Eig.nModes;` to:
```matlab
    K = double(max(Eig.nModes >= 1) * max(Eig.CompRank));   % paired rank count
```
Simpler and clearer — use:
```matlab
    K = double(max(Eig.CompRank));   % K_paired (per-component rank count)
```
Use the second form. Everything else in `UpdatePanel` (slider maxima `= K`, `ResetState` on `nModes ~= K`) stays.

- [ ] **Step 5: Run the paired test + the existing integration test (identity case)** — both must pass:
```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_paired.m')
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_integration.m')
```
Expected: both `ALL TESTS PASSED` (the single-sphere integration test has `CompRank=(1:n)'`, so `K_paired=nModes` and the expansion is the identity — it stays green).

- [ ] **Step 6: Commit**
```bash
git add toolbox/gui/panel_eigenmodes.m dev/tests/test_eigenmode_lever_paired.m
git commit -m "Eigenmode lever: reindex selection to paired rank (hemisphere-symmetric)"
```

---

## Task 2: Multi-client dispatch + view_eigenmodes synthesis

Make `FireModesChanged` dispatch eigenmode-view figures to a synthesis callback, and rewrite `view_eigenmodes` to consume the lever (`PairedGrid·W`), drop the `panel_time` hijack, and route arrows to `SetCurrentMode`.

**Files:** Modify `toolbox/core/bst_figures.m`, `toolbox/gui/view_eigenmodes.m`; create `dev/tests/test_eigenmode_viewer_synth.m`.

- [ ] **Step 1: Write the failing synthesis test** — `dev/tests/test_eigenmode_viewer_synth.m`:

```matlab
function test_eigenmode_viewer_synth
% Pure synthesis: view_eigenmodes('SynthColumn', PairedGrid, W) == PairedGrid*W(:).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(fullfile(repoRoot, 'toolbox', 'gui'));

nV = 50; K = 8;
PairedGrid = reshape(1:(nV*K), nV, K) / 100;

% single mode k=3 -> the k-th paired column
W = zeros(1,K); W(3) = 1;
col = view_eigenmodes('SynthColumn', PairedGrid, W);
assert(isequal(size(col), [nV 1]), 'column is nV x 1');
assert(max(abs(col - PairedGrid(:,3))) < 1e-12, 'single -> the rank-k paired column');

% band [2..4] box -> weighted sum
W = zeros(1,K); W(2:4) = 1;
col = view_eigenmodes('SynthColumn', PairedGrid, W);
assert(max(abs(col - sum(PairedGrid(:,2:4),2))) < 1e-12, 'band -> sum of selected columns');

fprintf('ALL TESTS PASSED: test_eigenmode_viewer_synth\n');
end
```

- [ ] **Step 2: Run it — verify it FAILS** (`SynthColumn` not implemented).
  Run: `run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_viewer_synth.m')`

- [ ] **Step 3: Rewrite `toolbox/gui/view_eigenmodes.m`** to consume the lever. Keep `BuildPairedGrid`. Add the pure `SynthColumn`, dispatch entries for `SynthColumn` and `ModesChangedCallback`, and rewrite `ViewFigure` to: build `PairedGrid`, create a single-frame transient source result (initial column = mode 1), tag the figure `EigenView`, set the lever to this surface (single shape, mode 1), and route arrows to `SetCurrentMode`. Replace the dispatch head and the relevant subfunctions:

```matlab
function varargout = view_eigenmodes(varargin)
% ... (keep the existing header comment) ...
if (nargin >= 1) && ischar(varargin{1}) && any(strcmp(varargin{1}, {'BuildPairedGrid','SynthColumn','ModesChangedCallback'}))
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
[varargout{1:nargout}] = ViewFigure(varargin{:});
end

%% ===== PURE: synthesized display column = PairedGrid * W(:) =====
function col = SynthColumn(PairedGrid, W) %#ok<DEFNU>
    W = W(:);
    if (numel(W) ~= size(PairedGrid,2))
        error('view_eigenmodes:SynthColumn: weight length (%d) must equal K_paired (%d).', numel(W), size(PairedGrid,2));
    end
    col = PairedGrid * W;
end
```

Keep `BuildPairedGrid` as-is. Rewrite `ViewFigure(SurfaceFile, ~)`:
- Load `Eig` via `in_tess_eigenmodes` (error if none, existing wording).
- `[Grid, Kp, Info] = BuildPairedGrid(Eig);`
- Initialise the lever for this surface and pick mode 1:
  `panel_eigenmodes('ResetState', SurfaceFile, Kp);`
  `panel_eigenmodes('SetWindowShape', 'single');`
  `panel_eigenmodes('SetCurrentMode', 1);`
  `W0 = panel_eigenmodes('GetWeights');`
- Build a single-frame Source result: `ResMat.ImageGridAmp = SynthColumn(Grid, W0);` `ResMat.Time = [0 1];` (duplicate the column so it is a valid 2-sample static result: `ResMat.ImageGridAmp = [col col];`), `ResMat.nComponents=1; ResMat.SurfaceFile=SurfaceFile; ResMat.HeadModelType='surface'; ResMat.ColormapType='stat2';` Comment as before. Save transient `results_eigenview` to the intra study and register (reuse the existing save/register/cleanup logic).
- `hFig = view_surface_data(SurfaceFile, file_short(OutputFile));`
- Tag the figure: `setappdata(hFig, 'EigenView', struct('SurfaceFile', SurfaceFile, 'PairedGrid', Grid, 'Info', Info, 'ResultsFile', file_short(OutputFile)));`
- Bring up + populate the panel: `gui_brainstorm('ShowToolTab', 'EigenModes'); panel_eigenmodes('UpdatePanel', hFig);`
- Replace the keyboard handler so arrows drive the lever:
```matlab
    function KeyPress_Callback(h, keyEvent)
        st = []; %#ok<NASGU>
        cur = 1;
        try, cur = panel_eigenmodes('GetCurrentMode'); catch, end
        switch (keyEvent.Key)
            case 'leftarrow',  panel_eigenmodes('SetCurrentMode', cur - 1);
            case 'rightarrow', panel_eigenmodes('SetCurrentMode', cur + 1);
            case 'pageup',     panel_eigenmodes('SetCurrentMode', cur + 10);
            case 'pagedown',   panel_eigenmodes('SetCurrentMode', cur - 10);
            case 'h', java_dialog('msgbox', ['Eigenmode viewer:' 10 '  Left/Right: step mode' 10 '  PgUp/PgDn: +/-10' 10 '  Use the "Spatial scale (eigenmodes)" panel to superpose a range.'], 'Eigenmode viewer');
            otherwise
                if ~isempty(KeyPressFcn_bak), KeyPressFcn_bak(h, keyEvent); end
        end
    end
```
- Keep the `DeleteFcn`/`CleanupResult` transient-result cleanup. Remove the old `SetMode`/`panel_time('SetCurrentTime')` stepping and the `Time = 1:K` grid.

Add the synthesis callback (dispatched from `FireModesChanged`):
```matlab
%% ===== Re-synthesize the viewer's displayed column on a lever change =====
function ModesChangedCallback(hFig) %#ok<DEFNU>
    global GlobalData;
    ev = getappdata(hFig, 'EigenView');
    if isempty(ev), return; end
    W = panel_eigenmodes('GetWeights');
    if isempty(W) || (numel(W) ~= size(ev.PairedGrid,2)), return; end
    col = SynthColumn(ev.PairedGrid, W);
    % Push into the loaded result and repaint via the standard surface path
    [iDS, iResult] = bst_memory('GetDataSetResult', ev.ResultsFile);
    if isempty(iDS), return; end
    GlobalData.DataSet(iDS).Results(iResult).ImageGridAmp = [col col];
    panel_surface('UpdateSurfaceData', hFig);
    figure_3d('UpdateSurfaceColor', hFig);
end
```

Add a tiny accessor used by the arrow keys — in `panel_eigenmodes.m` add:
```matlab
function k = GetCurrentMode() %#ok<DEFNU>
    st = GetState();
    k = st.iCurrentMode;
end
```

- [ ] **Step 4: Dispatch eigenmode-view figures in `FireModesChanged`** (`toolbox/core/bst_figures.m`). In the match branch, replace the single `panel_surface('UpdateSurfaceData', ...)` call with a client dispatch:
```matlab
            if ~isempty(getappdata(sFig.hFigure, 'EigenView'))
                view_eigenmodes('ModesChangedCallback', sFig.hFigure);
            else
                panel_surface('UpdateSurfaceData', sFig.hFigure);
            end
```

- [ ] **Step 5: Run the synthesis test + the lever tests** (no regression):
```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_viewer_synth.m')
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_paired.m')
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_integration.m')
```
Also `checkcode` both modified files — no new errors.

- [ ] **Step 6: Commit**
```bash
git add toolbox/gui/view_eigenmodes.m toolbox/core/bst_figures.m toolbox/gui/panel_eigenmodes.m dev/tests/test_eigenmode_viewer_synth.m
git commit -m "Eigenmode viewer: synthesize Phi*w from the lever; dispatch FireModesChanged by client"
```

---

## Task 3: Panel lifecycle & Active-gating

Make the panel inert when no eligible view is in front, and gate the filtering **Active** toggle to source-map context. Drive `UpdatePanel` on figure close / current-3D-change so it clears when the view goes away.

**Files:** Modify `toolbox/gui/panel_eigenmodes.m`, `toolbox/core/bst_figures.m`; create `dev/tests/test_eigenmode_lever_lifecycle.m`.

- [ ] **Step 1: Write the failing test** — `dev/tests/test_eigenmode_lever_lifecycle.m` (headless: drive the classification helper, which must not require Java):

```matlab
function test_eigenmode_lever_lifecycle
% Classify the front-figure context for the panel (pure, no Java).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(fullfile(repoRoot, 'toolbox', 'gui'));

% none
c = panel_eigenmodes('ClassifyContext', struct('isEigenView',false,'hasSourceModes',false));
assert(strcmp(c.kind,'none') && ~c.selectEnabled && ~c.activeEnabled && ~c.forceActive, 'none -> all off');

% eigenmode view -> selection on, Active off (N/A), force isActive=0
c = panel_eigenmodes('ClassifyContext', struct('isEigenView',true,'hasSourceModes',false));
assert(strcmp(c.kind,'view') && c.selectEnabled && ~c.activeEnabled, 'view -> select on, Active off');

% source map with modes -> selection + Active on
c = panel_eigenmodes('ClassifyContext', struct('isEigenView',false,'hasSourceModes',true));
assert(strcmp(c.kind,'source') && c.selectEnabled && c.activeEnabled, 'source -> select + Active on');

fprintf('ALL TESTS PASSED: test_eigenmode_lever_lifecycle\n');
end
```

- [ ] **Step 2: Run it — verify it FAILS** (`ClassifyContext` not implemented).

- [ ] **Step 3: Add `ClassifyContext` and rewire `UpdatePanel`** in `panel_eigenmodes.m`.

Add the pure classifier:
```matlab
%% ===== PURE: panel context from front-figure facts =====
function c = ClassifyContext(facts) %#ok<DEFNU>
    c = struct('kind','none','selectEnabled',false,'activeEnabled',false,'forceActive',false);
    if facts.isEigenView
        c.kind='view';   c.selectEnabled=true;  c.activeEnabled=false;
    elseif facts.hasSourceModes
        c.kind='source'; c.selectEnabled=true;  c.activeEnabled=true;
    end
end
```

Rewrite `UpdatePanel(hFig)` to: resolve the figure (default to the current 3D figure when omitted/empty), gather facts (`isEigenView` from appdata, `hasSourceModes`/`SurfaceFile` from `GetFigureSurfaceWithModes`), classify, then enable/disable the selection controls and the Active checkbox accordingly; when not a source context force `isActive=0`; when `none`, set the readout to "no eigenmode view". For an eligible view, reset state on surface/`K_paired` change and `RefreshControls`. Concretely:

```matlab
function UpdatePanel(hFig) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'EigenModes');
    if isempty(ctrl), return; end
    if (nargin < 1) || isempty(hFig) || ~ishandle(hFig)
        hFig = bst_figures('GetCurrentFigure', '3D');
    end
    isEigenView = ~isempty(hFig) && ishandle(hFig) && ~isempty(getappdata(hFig, 'EigenView'));
    SurfaceFile = '';
    if isEigenView
        ev = getappdata(hFig, 'EigenView'); SurfaceFile = ev.SurfaceFile;
    elseif ~isempty(hFig) && ishandle(hFig)
        SurfaceFile = GetFigureSurfaceWithModes(hFig);
    end
    facts = struct('isEigenView', isEigenView, 'hasSourceModes', ~isEigenView && ~isempty(SurfaceFile));
    c = ClassifyContext(facts);
    % Enable/disable the selection controls and the Active toggle
    SetSelectEnabled(ctrl, c.selectEnabled);
    ctrl.jCheckActive.setEnabled(c.activeEnabled);
    ctrl.jCheckActive.setVisible(c.activeEnabled);
    if ~strcmp(c.kind, 'source')
        SetActive(0);                          % no stale filtering off-source
    end
    if strcmp(c.kind, 'none')
        ctrl.jLabelReadout.setText('no eigenmode view');
        return;
    end
    % Populate from the eligible surface
    [Eig, ~] = in_tess_eigenmodes(SurfaceFile);
    K = double(max(Eig.CompRank));
    st = GetState();
    if ~file_compare(st.SurfaceFile, SurfaceFile) || (st.nModes ~= K)
        ResetState(SurfaceFile, K);
        if isEigenView
            SetWindowShape('single'); SetCurrentMode(1);
        end
    end
    ctrl.jSliderLo.setMaximum(K); ctrl.jSliderHi.setMaximum(K);
    RefreshControls();
end
```

Add the helper `SetSelectEnabled` (enables only the selection controls — sliders, radios, band label — leaving `jCheckActive`/`jPanelTop` managed separately):
```matlab
function SetSelectEnabled(ctrl, isOn)
    sel = {'jSliderLo','jSliderHi','jLabelBand','jRadioSingle','jRadioBox','jRadioTaper','jRadioGain'};
    for i = 1:numel(sel)
        if isfield(ctrl, sel{i}) && isa(ctrl.(sel{i}), 'javax.swing.JComponent')
            ctrl.(sel{i}).setEnabled(logical(isOn));
        end
    end
end
```
(Keep the old `SetPanelEnabled` if still referenced, or remove it if now unused.)

- [ ] **Step 4: Drive `UpdatePanel` on figure close / current-3D clearing** in `toolbox/core/bst_figures.m`. Find where figures are deleted/closed and the current figure is recomputed (search for `panel_surface('UpdatePanel')` calls that run on close / current-figure change). Next to each such call that fires on close or current-3D change, add (gated, so it is cheap when the tab is absent):
```matlab
        if gui_brainstorm('isTabVisible', 'EigenModes')
            panel_eigenmodes('UpdatePanel');
        end
```
At minimum ensure it runs when a 3D figure closes and when the current 3D figure changes to a non-eligible figure, so the panel goes inert. (The Type3D focus hook from the lever already covers focus changes to eligible figures.)

- [ ] **Step 5: Run the lifecycle test + full lever suite**:
```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_lifecycle.m')
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_weights.m')
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_state.m')
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_paired.m')
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_integration.m')
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_panel.m')
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_viewer_synth.m')
```
All `ALL TESTS PASSED`. `checkcode` both modified files — no new errors.

- [ ] **Step 6: Commit**
```bash
git add toolbox/gui/panel_eigenmodes.m toolbox/core/bst_figures.m dev/tests/test_eigenmode_lever_lifecycle.m
git commit -m "Eigenmode lever: inert when no eligible view; gate Active to source context"
```

---

## Task 4: End-to-end live test

**Files:** Create `dev/tests/test_eigenmode_viewer_e2e.m`.

- [ ] **Step 1: Write the e2e test**:

```matlab
function test_eigenmode_viewer_e2e(SurfaceFile)
% Live: open View Eigenmodes, step modes via the lever, superpose a band.
% USAGE: test_eigenmode_viewer_e2e(SurfaceFile)   % cortex with eigenmodes
if ~brainstorm('status'); brainstorm nogui; end
[~, isComputed] = in_tess_eigenmodes(SurfaceFile);
if ~isComputed
    process_eigenmodes('Compute', SurfaceFile, 200, 'barycentric', true, false, false);
end
hFig = view_eigenmodes(SurfaceFile);
assert(~isempty(hFig), 'viewer failed to open');
ev = getappdata(hFig, 'EigenView');
assert(~isempty(ev) && isfield(ev,'PairedGrid'), 'figure tagged EigenView with PairedGrid');

% Mode 1 (single) is shown
TessInfo = getappdata(hFig, 'Surface');
c1 = TessInfo(1).Data;
assert(max(abs(c1 - ev.PairedGrid(:,1))) < 1e-6, 'mode 1 shown initially');

% Step to mode 2 via the lever (as the arrow keys do)
panel_eigenmodes('SetCurrentMode', 2);
TessInfo = getappdata(hFig, 'Surface');
c2 = TessInfo(1).Data;
assert(max(abs(c2 - ev.PairedGrid(:,2))) < 1e-6, 'stepping shows mode 2');

% Superpose a band [1,5] (box) -> sum of paired columns 1..5
panel_eigenmodes('SetWindowShape', 'box');
panel_eigenmodes('SetBand', 1, 5);
TessInfo = getappdata(hFig, 'Surface');
cb = TessInfo(1).Data;
assert(max(abs(cb - sum(ev.PairedGrid(:,1:5),2))) < 1e-6, 'band -> superposition');

close(hFig);
fprintf('ALL TESTS PASSED: test_eigenmode_viewer_e2e\n');
end
```

- [ ] **Step 2: Run it on real data (best effort).** Probe the protocol for a cortex with eigenmodes (`bst_get('ProtocolSubjects')`; a `tess_cortex_*` with eigenmodes, or compute them). Run `test_eigenmode_viewer_e2e('<SurfaceFile>')`. If GUI/data unavailable, document it and rely on the headless suite as the gate (run all eight lever/viewer tests).

- [ ] **Step 3: `checkcode` the e2e file; commit**
```bash
git add dev/tests/test_eigenmode_viewer_e2e.m
git commit -m "Eigenmode viewer: end-to-end live test (step + superpose via the lever)"
```

## Done criteria
- A low paired-band keeps both hemispheres (symmetry) — `test_eigenmode_lever_paired` green.
- `view_eigenmodes` shows `Φ·w` (single mode or superposition), driven by the panel and arrows; no `panel_time` hijack.
- The panel is inert when no eligible view is in front; the Active toggle is hidden/disabled in viewer context.
- All headless tests pass; e2e passes live (or is documented data-gated).
