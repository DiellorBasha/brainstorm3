# panel_eigenfilter_options Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a standalone Brainstorm GUI options panel for the eigen-domain spatial filter (`bst_eigen` `Method='filter'`), with a live spectral-response `h(λ)` preview, built like `panel_timefreq_options`.

**Architecture:** Rename the convention-violating kernel widget `bst_eigfilter_panel.m` → `panel_eigenfilter_design.m` and give it a reusable `DrawResponse` plot verb. Build a new `panel_eigenfilter_options.m` (`CreatePanel(EigenFile)` / `GetPanelContents`) that embeds the design widget and toggles a live `h(λ)` figure. No changes to the `bst_eigfilter_*` kernel library.

**Tech Stack:** MATLAB, Brainstorm GUI toolkit (`gui_river`, `gui_component`, `bst_mutex`, `BstPanel`, `macro_method` dispatch), Java Swing. Validation via the MATLAB MCP (`evaluate_matlab_code`) + `checkcode`.

## Global Constraints

- **No library changes.** Do NOT touch `toolbox/eigen/eigfilter/` or any `bst_eigfilter_*` function. Only the three GUI files in this plan change.
- **Brainstorm panel contract:** every panel file begins with `eval(macro_method);` dispatch; public verbs are `CreatePanel` / `GetPanelContents`; controls are returned via `BstPanel(panelName, jScroll, ctrl)` and re-fetched with `bst_get('PanelControls', panelName)`.
- **License header:** every new/renamed `.m` file keeps the standard Brainstorm GPL header block and an `% Authors: Diellor Basha, 2026` line (copy verbatim from `panel_eigenmodes_compute.m`).
- **Naming:** GUI panels use the `panel_*` prefix (never `bst_*…_panel`).
- **No `clear`** in any MATLAB validation snippet (it wipes `GlobalData` and hangs the live Brainstorm session); use `rehash` if needed — edited `.m` files auto-reload.
- **Commits:** work on a feature branch (Task 0), never directly on `development`. Commit per task.

---

### Task 0: Create the feature branch

**Files:** none (git only)

- [ ] **Step 1: Branch off development**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git checkout -b feature/panel-eigenfilter-options
```

Expected: `Switched to a new branch 'feature/panel-eigenfilter-options'`

---

### Task 1: Rename the kernel widget and update its caller

Rename `bst_eigfilter_panel.m` → `panel_eigenfilter_design.m` (function name + header + API doc examples), and repoint its only external caller, `panel_helmholtz.m`. Pure mechanical change — no behavior difference.

**Files:**
- Rename: `toolbox/gui/bst_eigfilter_panel.m` → `toolbox/gui/panel_eigenfilter_design.m`
- Modify: `toolbox/gui/panel_eigenfilter_design.m` (function name, header, doc examples)
- Modify: `toolbox/gui/panel_helmholtz.m` (6 call sites)

**Interfaces:**
- Produces: `panel_eigenfilter_design('Kernels')` → `[keys, displays]`; `('CurrentKernel', jKernel, keys)` → `key`; `('BuildSliders', jParams, kernelKey, Lambda, onSettle)`; `('ParamNames', jParams)` → cellstr; `('ReadParams', jParams, Lambda)` → params struct. (All identical to the old `bst_eigfilter_panel` verbs.)

- [ ] **Step 1: Rename the file (preserve history)**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git mv toolbox/gui/bst_eigfilter_panel.m toolbox/gui/panel_eigenfilter_design.m
```

- [ ] **Step 2: Update the function name + header in the renamed file**

In `toolbox/gui/panel_eigenfilter_design.m`, change the first line and the doc header. Old:

```matlab
function varargout = bst_eigfilter_panel(varargin)
% BST_EIGFILTER_PANEL: Shared "Filter kernel" UI section (kernel dropdown + mode-index
% scale sliders), reused by panel_helmholtz.
%
% API (dispatched via macro_method):
%   [keys, displays] = bst_eigfilter_panel('Kernels')
%   key   = bst_eigfilter_panel('CurrentKernel', jKernel, keys)
%   bst_eigfilter_panel('BuildSliders', jParams, kernelKey, Lambda, onSettle)
%   names = bst_eigfilter_panel('ParamNames', jParams)
%   params = bst_eigfilter_panel('ReadParams', jParams, Lambda)
```

New:

```matlab
function varargout = panel_eigenfilter_design(varargin)
% PANEL_EIGENFILTER_DESIGN: Shared eigenfilter "design" UI section (kernel dropdown +
% mode-index scale sliders), reused by panel_helmholtz and panel_eigenfilter_options.
% Also renders the filter's spectral response h(lambda) (DrawResponse).
%
% API (dispatched via macro_method):
%   [keys, displays] = panel_eigenfilter_design('Kernels')
%   key   = panel_eigenfilter_design('CurrentKernel', jKernel, keys)
%   panel_eigenfilter_design('BuildSliders', jParams, kernelKey, Lambda, onSettle)
%   names = panel_eigenfilter_design('ParamNames', jParams)
%   params = panel_eigenfilter_design('ReadParams', jParams, Lambda)
%   panel_eigenfilter_design('DrawResponse', hAxes, kernelName, params, Lambda)
```

Leave every other function in the file (`Kernels`, `CurrentKernel`, `BuildSliders`, `ParamNames`, `ReadParams`, and all `i_*` internals) unchanged. The body uses `eval(macro_method)` dispatch and has no self-name calls, so nothing else needs editing.

- [ ] **Step 3: Repoint the caller `panel_helmholtz.m`**

Replace all 6 occurrences of `bst_eigfilter_panel` with `panel_eigenfilter_design`:

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
sed -i '' 's/bst_eigfilter_panel/panel_eigenfilter_design/g' toolbox/gui/panel_helmholtz.m
```

- [ ] **Step 4: Verify no stale references remain anywhere**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
grep -rn "bst_eigfilter_panel" toolbox/ dev/
```

Expected: no output (zero matches).

- [ ] **Step 5: Lint both files**

Via the MATLAB MCP `evaluate_matlab_code`:

```matlab
disp(checkcode('toolbox/gui/panel_eigenfilter_design.m'))
disp(checkcode('toolbox/gui/panel_helmholtz.m'))
```

Expected: no real errors (Brainstorm-idiom warnings such as `#ok` suppressions are acceptable; there must be no "undefined function/variable" or syntax errors).

- [ ] **Step 6: Commit**

```bash
git add toolbox/gui/panel_eigenfilter_design.m toolbox/gui/panel_helmholtz.m
git commit -m "refactor(gui): rename bst_eigfilter_panel -> panel_eigenfilter_design (panel_ convention)"
```

---

### Task 2: Add the `DrawResponse` verb to `panel_eigenfilter_design`

Add a stateless plot verb that evaluates the kernel on `Lambda` and plots gain vs mode index. Testable without a running Brainstorm (needs only `bst_eigfilter_kernel`/`bst_eigfilter_evaluate` on the path).

**Files:**
- Modify: `toolbox/gui/panel_eigenfilter_design.m` (append one function)
- Test (scratch): `dev/tests/test_drawresponse.m`

**Interfaces:**
- Consumes: `bst_eigfilter_kernel(name, params)` → handle (or cell bank); `bst_eigfilter_evaluate(g, Lambda)` → `[K x 1]` gain.
- Produces: `panel_eigenfilter_design('DrawResponse', hAxes, kernelName, params, Lambda)` — draws one line into `hAxes`; returns nothing.

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_drawresponse.m`:

```matlab
function test_drawresponse()
% Headless check that DrawResponse plots a gain curve into an axes.
Lambda = linspace(0.01, 5, 50)';
hFig = figure('Visible','off');
hAx  = axes('Parent', hFig);
panel_eigenfilter_design('DrawResponse', hAx, 'heat', struct('t',0.1), Lambda);
hLine = findobj(hAx, 'Type', 'line');
assert(~isempty(hLine), 'DrawResponse drew no line.');
yd = get(hLine(1), 'YData');
assert(numel(yd) == numel(Lambda), 'Gain length must equal numel(Lambda).');
assert(all(yd <= 1 + 1e-9) && all(yd >= 0), 'Heat gain must be in [0,1].');
close(hFig);
disp('test_drawresponse PASSED');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Via the MATLAB MCP `run_matlab_file` on `dev/tests/test_drawresponse.m`.
Expected: error — unknown command/verb `DrawResponse` (dispatch falls through in `macro_method`).

- [ ] **Step 3: Implement `DrawResponse`**

Append this function to `toolbox/gui/panel_eigenfilter_design.m` (after the existing verb functions, before/among the `i_*` internals — any top-level position works):

```matlab
function DrawResponse(hAxes, kernelName, params, Lambda) %#ok<DEFNU>
% Plot the filter's spectral response h(lambda) (gain vs mode index) into hAxes.
    Lambda = double(Lambda(:));
    if isempty(Lambda)
        cla(hAxes);
        return;
    end
    g = bst_eigfilter_kernel(kernelName, params);
    if iscell(g)            % sliders yield single-scale params; guard the bank case
        g = g{1};
    end
    h = bst_eigfilter_evaluate(g, Lambda);
    k = (1:numel(Lambda))';
    plot(hAxes, k, h, 'LineWidth', 2);
    set(hAxes, 'XGrid', 'on', 'YGrid', 'on');
    xlabel(hAxes, 'Mode index k');
    ylabel(hAxes, 'Gain h(\lambda)');
    title(hAxes, sprintf('Spectral response: %s', kernelName), 'Interpreter', 'none');
    xlim(hAxes, [1, max(2, numel(Lambda))]);
end
```

- [ ] **Step 4: Run the test to verify it passes**

Via the MATLAB MCP `run_matlab_file` on `dev/tests/test_drawresponse.m`.
Expected: `test_drawresponse PASSED`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_eigenfilter_design.m dev/tests/test_drawresponse.m
git commit -m "feat(gui): add DrawResponse spectral-response plot verb to panel_eigenfilter_design"
```

---

### Task 3: Create `panel_eigenfilter_options` core (CreatePanel + GetPanelContents)

The new options panel without the live Display yet: eigen-basis info, embedded design widget, comment field, OK/Cancel, and `GetPanelContents` returning a `bst_eigen`-ready OPTIONS struct.

**Files:**
- Create: `toolbox/gui/panel_eigenfilter_options.m`
- Test (scratch): `dev/tests/test_panel_eigenfilter_options.m`

**Interfaces:**
- Consumes: `in_bst_eigen(EigenFile)` → `EigenMat` with `.Variant`, `.Lambda{h}`; `panel_eigenfilter_design` verbs from Task 1.
- Produces: `panel_eigenfilter_options('CreatePanel', EigenFile)` → `[bstPanel, panelName]`; `panel_eigenfilter_options('GetPanelContents')` → struct with fields `Method`(='filter'), `EigenFile`, `KernelName`, `KernelParams`, `Comment`.

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_panel_eigenfilter_options.m`:

```matlab
function test_panel_eigenfilter_options(EigenFile)
% Build the panel from a real eigen_ node, read it back, assert OPTIONS fields.
% Pass an eigen_ file relative to the protocol (e.g. from bst_get('AnyFile', ...)).
if nargin < 1 || isempty(EigenFile)
    error('Pass the path to an eigen_ node file.');
end
[bstPanel, panelName] = panel_eigenfilter_options('CreatePanel', EigenFile);
assert(~isempty(bstPanel), 'CreatePanel returned empty.');
gui_show(bstPanel, 'JavaWindow', panelName, 0, 0, 0);
drawnow;
s = panel_eigenfilter_options('GetPanelContents');
gui_hide(panelName);
assert(strcmp(s.Method, 'filter'), 'Method must be ''filter''.');
assert(ischar(s.EigenFile) && ~isempty(s.EigenFile), 'EigenFile must be set.');
assert(ischar(s.KernelName) && ~isempty(s.KernelName), 'KernelName must be set.');
assert(isstruct(s.KernelParams), 'KernelParams must be a struct.');
assert(ischar(s.Comment), 'Comment must be a char.');
disp('test_panel_eigenfilter_options PASSED');
end
```

- [ ] **Step 2: Run the test to verify it fails**

First locate a real eigen_ node, then run, via the MATLAB MCP `evaluate_matlab_code`:

```matlab
% find any eigen_ node in the current protocol
ProtocolInfo = bst_get('ProtocolInfo');
d = dir(fullfile(ProtocolInfo.STUDIES, '**', 'eigen_*.mat'));
EigenFile = file_short(fullfile(d(1).folder, d(1).name));
test_panel_eigenfilter_options(EigenFile)
```

Expected: error — `Undefined function 'panel_eigenfilter_options'`.

- [ ] **Step 3: Implement the panel file**

Create `toolbox/gui/panel_eigenfilter_options.m` (copy the GPL header + `Authors` line from `panel_eigenmodes_compute.m`):

```matlab
function varargout = panel_eigenfilter_options(varargin)
% PANEL_EIGENFILTER_OPTIONS: Options for the eigen-domain spatial filter (bst_eigen 'filter').
%
% USAGE:  [bstPanel, panelName] = panel_eigenfilter_options('CreatePanel', EigenFile)
%                            s  = panel_eigenfilter_options('GetPanelContents')
%
% The spatial analogue of panel_timefreq_options. EigenFile is an eigen_ node; the panel
% reads its eigenvalues (Lambda) to drive the kernel sliders, embeds the shared
% panel_eigenfilter_design widget, and (Display toggle) previews the spectral response.
% A future process_eigenfilter attaches via an 'editpref' option; CreatePanel is
% arg-type aware so a (sProcess, sInputs) branch can be added without a rewrite.

% <<COPY THE STANDARD BRAINSTORM GPL HEADER BLOCK HERE, verbatim from panel_eigenmodes_compute.m>>
%
% Authors: Diellor Basha, 2026

eval(macro_method);
end


%% ===== CREATE PANEL =====
function [bstPanelNew, panelName] = CreatePanel(EigenFile) %#ok<DEFNU>
    panelName = 'EigenfilterOptions';
    import java.awt.*;
    import javax.swing.*;

    % Arg-type aware: char/file => EigenFile path (now); struct => future (sProcess,sInputs)
    if isstruct(EigenFile)
        bst_error('panel_eigenfilter_options: the (sProcess, sInputs) path is not implemented yet; pass an eigen_ file.', 'Eigenfilter options', 0);
        bstPanelNew = []; panelName = []; return;
    end
    if isempty(EigenFile) || ~ischar(EigenFile)
        bst_error('panel_eigenfilter_options: a valid eigen_ file is required.', 'Eigenfilter options', 0);
        bstPanelNew = []; panelName = []; return;
    end

    % Load the eigenbasis (the spatial axis)
    EigenMat = in_bst_eigen(EigenFile);
    Lambda = [];
    for h = 1:numel(EigenMat.Lambda)
        if ~isempty(EigenMat.Lambda{h})
            Lambda = EigenMat.Lambda{h}(:);
            break;
        end
    end
    if isempty(Lambda)
        bst_error('panel_eigenfilter_options: the eigen_ node has no eigenvalues.', 'Eigenfilter options', 0);
        bstPanelNew = []; panelName = []; return;
    end
    K     = numel(Lambda);
    nHemi = sum(~cellfun(@isempty, EigenMat.Lambda));

    % Display figure handle (shared across nested callbacks)
    hFigResp = [];

    % ===== MAIN PANEL =====
    jPanelNew = gui_river([5,5], [10,15,12,10]);

    % ===== EIGEN BASIS INFO =====
    jPanelInfo = gui_river([2,2], [0,10,10,10], 'Eigen basis');
        gui_component('label', jPanelInfo, '',   sprintf('Variant:  %s', EigenMat.Variant), [], [], [], []);
        gui_component('label', jPanelInfo, 'br', sprintf('Modes:  %d      Hemispheres:  %d', K, nHemi), [], [], [], []);
    jPanelNew.add('br hfill', jPanelInfo);

    % ===== FILTER DESIGN (reuse the shared design widget) =====
    jPanelDes = gui_river([2,2], [0,10,12,10], 'Filter design');
        [keys, displays] = panel_eigenfilter_design('Kernels');
        gui_component('label', jPanelDes, '', 'Kernel: ', [], [], [], []);
        jKernel = gui_component('combobox', jPanelDes, 'tab hfill', [], {displays}, [], [], []);
        iHeat = find(strcmp(keys, 'heat'), 1);
        if ~isempty(iHeat); jKernel.setSelectedIndex(iHeat-1); end
        jParams = gui_river([2,2], [0,2,0,2]);
        jPanelDes.add('br hfill', jParams);
        panel_eigenfilter_design('BuildSliders', jParams, ...
            panel_eigenfilter_design('CurrentKernel', jKernel, keys), Lambda, @() UpdateResponse());
        java_setcb(jKernel, 'ActionPerformedCallback', @(hh,ee) OnKernel());
    jPanelNew.add('br hfill', jPanelDes);

    % ===== DISPLAY TOGGLE =====
    jPanelDisp = gui_river([2,2], [0,10,8,10], 'Display');
        jToggleDisp = gui_component('toggle', jPanelDisp, '', 'Show spectral response', [], [], @(hh,ee) ToggleDisplay());
    jPanelNew.add('br hfill', jPanelDisp);

    % ===== OUTPUT COMMENT =====
    jPanelCom = gui_river([2,2], [0,10,8,10], 'Output');
        gui_component('label', jPanelCom, '', 'Comment: ', [], [], [], []);
        jTextComment = gui_component('text', jPanelCom, 'tab hfill', '', [], [], [], []);
    jPanelNew.add('br hfill', jPanelCom);

    % ===== OK / CANCEL =====
    gui_component('button', jPanelNew, 'br right', 'Cancel', [], [], @ButtonCancel_Callback, []);
    gui_component('button', jPanelNew, [],         'OK',     [], [], @ButtonOk_Callback, []);

    % ===== ASSEMBLE =====
    jPanelScroll = javax.swing.JScrollPane(jPanelNew);
    bst_mutex('create', panelName);
    ctrl = struct('jKernel',      jKernel, ...
                  'KernelKeys',   {keys}, ...
                  'jParams',      jParams, ...
                  'jToggleDisp',  jToggleDisp, ...
                  'jTextComment', jTextComment, ...
                  'Lambda',       Lambda, ...
                  'EigenFile',    EigenFile);
    bstPanelNew = BstPanel(panelName, jPanelScroll, ctrl);

    %% ===== NESTED CALLBACKS =====
    function OnKernel(varargin)
        key = panel_eigenfilter_design('CurrentKernel', jKernel, keys);
        panel_eigenfilter_design('BuildSliders', jParams, key, Lambda, @() UpdateResponse());
        UpdateResponse();
    end

    function ToggleDisplay(varargin)
        if jToggleDisp.isSelected()
            if isempty(hFigResp) || ~ishandle(hFigResp)
                hFigResp = figure('MenuBar','none', 'Toolbar','none', 'NumberTitle','off', ...
                                  'Name','Eigenfilter spectral response', 'Pointer','arrow');
            end
            UpdateResponse();
        else
            if ~isempty(hFigResp) && ishandle(hFigResp); close(hFigResp); end
            hFigResp = [];
        end
    end

    function UpdateResponse(varargin)
        if isempty(hFigResp) || ~ishandle(hFigResp); return; end
        key    = panel_eigenfilter_design('CurrentKernel', jKernel, keys);
        params = panel_eigenfilter_design('ReadParams', jParams, Lambda);
        hAxes  = findobj(hFigResp, 'Type', 'axes');
        if isempty(hAxes)
            hAxes = axes('Parent', hFigResp);
        else
            hAxes = hAxes(1);
        end
        panel_eigenfilter_design('DrawResponse', hAxes, key, params, Lambda);
    end

    function ButtonCancel_Callback(varargin)
        if ~isempty(hFigResp) && ishandle(hFigResp); close(hFigResp); end
        gui_hide(panelName);
    end

    function ButtonOk_Callback(varargin)
        if ~isempty(hFigResp) && ishandle(hFigResp); close(hFigResp); end
        bst_mutex('release', panelName);
    end
end


%% ===== GET PANEL CONTENTS =====
function s = GetPanelContents() %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'EigenfilterOptions');
    key    = panel_eigenfilter_design('CurrentKernel', ctrl.jKernel, ctrl.KernelKeys);
    params = panel_eigenfilter_design('ReadParams', ctrl.jParams, ctrl.Lambda);
    s = struct();
    s.Method       = 'filter';
    s.EigenFile    = ctrl.EigenFile;
    s.KernelName   = key;
    s.KernelParams = params;
    s.Comment      = char(ctrl.jTextComment.getText());
end
```

NOTE for the implementer: replace the `<<COPY THE STANDARD BRAINSTORM GPL HEADER BLOCK HERE>>` line with the exact 17-line GPL block found between the `% @===` markers in `toolbox/gui/panel_eigenmodes_compute.m`.

- [ ] **Step 4: Run the test to verify it passes**

Re-run Step 2's snippet (locate eigen_ node + `test_panel_eigenfilter_options(EigenFile)`).
Expected: `test_panel_eigenfilter_options PASSED` (a panel briefly shows and hides).

- [ ] **Step 5: Lint**

```matlab
disp(checkcode('toolbox/gui/panel_eigenfilter_options.m'))
```

Expected: no real errors.

- [ ] **Step 6: Commit**

```bash
git add toolbox/gui/panel_eigenfilter_options.m dev/tests/test_panel_eigenfilter_options.m
git commit -m "feat(gui): add panel_eigenfilter_options (CreatePanel/GetPanelContents)"
```

---

### Task 4: End-to-end validation — drive `bst_eigen` from the panel

Confirm the OPTIONS produced by the panel actually run the filter and write a `results_eigenfilter` file. No new product code — this is the integration gate (and the live-Display visual check).

**Files:**
- Test (scratch): `dev/tests/test_eigenfilter_end_to_end.m`

**Interfaces:**
- Consumes: `panel_eigenfilter_options('GetPanelContents')` OPTIONS; `bst_eigen(Data, OPTIONS)`.

- [ ] **Step 1: Write the end-to-end test**

Create `dev/tests/test_eigenfilter_end_to_end.m`:

```matlab
function OutputFiles = test_eigenfilter_end_to_end(ResultsFile, EigenFile)
% Build OPTIONS from the panel, run the eigen filter, assert a results file is returned.
% ResultsFile : a source map (results_*.mat) on the same surface as EigenFile.
% EigenFile   : an eigen_ node.
[bstPanel, panelName] = panel_eigenfilter_options('CreatePanel', EigenFile);
gui_show(bstPanel, 'JavaWindow', panelName, 0, 0, 0);
drawnow;
OPTIONS = panel_eigenfilter_options('GetPanelContents');
gui_hide(panelName);
OPTIONS.iTargetStudy = 'NoSave';   % return contents, do not write to DB
[OutputFiles, Messages, isError] = bst_eigen(ResultsFile, OPTIONS);
assert(~isError, 'bst_eigen reported an error: %s', Messages);
assert(~isempty(OutputFiles), 'bst_eigen produced no output.');
FileMat = OutputFiles{1};
assert(isfield(FileMat, 'ImageGridAmp') && ~isempty(FileMat.ImageGridAmp), ...
    'Filtered source map is empty.');
disp('test_eigenfilter_end_to_end PASSED');
end
```

- [ ] **Step 2: Run it against real files**

Via the MATLAB MCP `evaluate_matlab_code` (pick a `results_*` file and matching `eigen_` node from the loaded protocol):

```matlab
ProtocolInfo = bst_get('ProtocolInfo');
dE = dir(fullfile(ProtocolInfo.STUDIES, '**', 'eigen_*.mat'));
EigenFile = file_short(fullfile(dE(1).folder, dE(1).name));
% choose a results file on the same surface (Dirac/unconstrained for 3-vector variants)
ResultsFile = '<fill in a results_*.mat from the same subject/surface>';
test_eigenfilter_end_to_end(ResultsFile, EigenFile)
```

Expected: `test_eigenfilter_end_to_end PASSED`.

- [ ] **Step 3: Visual check of the live Display**

Via the MATLAB MCP `evaluate_matlab_code` — open the panel, toggle Display, move a slider, confirm the `h(λ)` curve redraws:

```matlab
[bstPanel, panelName] = panel_eigenfilter_options('CreatePanel', EigenFile);
gui_show(bstPanel, 'JavaWindow', panelName, 0, 0, 0);
% In the panel: click "Show spectral response", switch kernels, drag sliders.
```

Expected: a figure titled "Eigenfilter spectral response" shows a gain curve that updates when the kernel/sliders change. Close the panel when done.

- [ ] **Step 4: Commit**

```bash
git add dev/tests/test_eigenfilter_end_to_end.m
git commit -m "test(gui): end-to-end panel_eigenfilter_options -> bst_eigen filter"
```

---

## Self-Review

- **Spec coverage:** Rename + caller update (Task 1) ✓; `DrawResponse` (Task 2) ✓; `CreatePanel(EigenFile)` arg-type-aware + `GetPanelContents` OPTIONS (Task 3) ✓; Display toggle + live `h(λ)` (Task 3 UI + Task 4 visual check) ✓; error handling (CreatePanel guards, DrawResponse empty-Lambda guard) ✓; validation via MCP + checkcode ✓. No-library-change constraint respected (only 3 GUI files) ✓.
- **Placeholder scan:** the only intentional fill-in is the GPL header copy (Task 3 Step 3 note) and the e2e `ResultsFile` choice (environment-specific) — both are explicit instructions, not vague TODOs.
- **Type consistency:** `panel_eigenfilter_design` verb names match across Tasks 1–4; `ctrl` field names set in CreatePanel match those read in `GetPanelContents` (`jKernel`, `KernelKeys`, `jParams`, `jTextComment`, `Lambda`, `EigenFile`); OPTIONS fields (`Method`/`EigenFile`/`KernelName`/`KernelParams`/`Comment`) match `bst_eigen`'s `Def_OPTIONS`.
```
