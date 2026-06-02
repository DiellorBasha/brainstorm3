# Eigenmode Time Series Viewer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an eigenmode time series viewer that plots one trace per eigenmode coefficient `θₖ(t)` in the standard butterfly/column layout, tracking the `EigenModes` panel band live.

**Architecture:** A new `view_eigenmodes_timeseries.m` resolves the study's head model + surface eigenmodes, runs `bst_eigenmodes_transform` once to cache the full coefficient matrix `Theta [K_raw × nTime]`, selects the panel-selected band's rows (each paired rank → its left/right raw columns = two traces), and hands them to the existing `view_timeseries_matrix`. A new dispatch branch in `bst_figures('FireModesChanged')` re-selects cached rows on band change; a new menu item in the 3D figure popup launches it.

**Tech Stack:** MATLAB, Brainstorm GUI conventions (`eval(macro_method)` string dispatch, `view_timeseries_matrix`, `in_tess_eigenmodes`, `bst_eigenmodes_transform`).

---

## File Structure

| File | Responsibility | Action |
|------|----------------|--------|
| `toolbox/gui/view_eigenmodes_timeseries.m` | Entry point, pure band→trace selection, transform/cache, live refresh callback | **Create** |
| `toolbox/core/bst_figures.m` (`FireModesChanged`, ~line 1010) | Dispatch band changes to the time series client | **Modify** |
| `toolbox/gui/figure_3d.m` (`DisplayFigurePopup`, ~line 1655) | "Eigenmode time series" launch menu item + guard | **Modify** |
| `dev/tests/test_view_eigenmodes_timeseries_pure.m` | Unit test for the pure band→trace selector | **Create** |

**Key reference signatures (already in the codebase — do not change):**
- `[Kernel, Info] = bst_eigenmodes_transform(Gain, Phi)` — `Gain` is `[nch × nVert]` constrained, restricted to channels; `Phi` is `[nVert × K]`; `Theta = Kernel * Data(iCh,:)`.
- `[Eigenmodes, isComputed] = in_tess_eigenmodes(SurfaceFile)` — returns struct with `.Vectors [nVert×nModes]`, `.Values`, `.Component`, `.CompRank`, `.nComponents`, `.nModes`.
- `[hFig, iDS, iFig] = view_timeseries_matrix(BaseFiles, F, TimeVector, Modality, AxesLabels, LinesLabels, LinesColor, hFig, Std, DisplayUnits)` — `F` is a **cell** `{[nRows×nTime]}`; `LinesLabels` is a **cell** `{1×nRows}` of strings; `LinesColor` is a **cell** `{1×nRows}` of `[1×3]` RGB. Passing an existing `hFig` redraws in place.
- `panel_eigenmodes('GetState')` → struct with `.SurfaceFile`, `.nModes` (= K_paired), `.Band` `[kLo kHi]`. `panel_eigenmodes('ResetState', SurfaceFile, K)` initializes a default band. `panel_eigenmodes('RefreshControls')` reflects state into the panel controls.

---

## Task 1: Pure band→trace selector + color helper (with dispatch skeleton)

**Files:**
- Create: `toolbox/gui/view_eigenmodes_timeseries.m`
- Test: `dev/tests/test_view_eigenmodes_timeseries_pure.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_view_eigenmodes_timeseries_pure.m`:

```matlab
function test_view_eigenmodes_timeseries_pure
% Verify the band->trace selector: each paired rank in [kLo,kHi] yields its
% left column then its right column as two traces; single-component data
% yields one unlabelled-hemisphere trace per rank.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% --- Paired case: 2 components, ranks 1..2 each ---
% Raw columns:  1=L/r1, 2=L/r2, 3=R/r1, 4=R/r2
Component = [1;1;2;2];
CompRank  = [1;2;1;2];

[iRows, Labels, Hemi] = view_eigenmodes_timeseries('GetBandTraces', Component, CompRank, 1, 2);
assert(isequal(iRows(:)', [1 3 2 4]), 'Order must be L then R, per rank ascending.');
assert(isequal(Labels, {'Mode 1 L','Mode 1 R','Mode 2 L','Mode 2 R'}), 'Paired labels wrong.');
assert(isequal(Hemi(:)', [1 2 1 2]), 'Hemisphere ids wrong.');

% Sub-band selects a single rank -> 2 traces
[iRows2, Labels2] = view_eigenmodes_timeseries('GetBandTraces', Component, CompRank, 2, 2);
assert(isequal(iRows2(:)', [2 4]), 'Single-rank band must give that rank''s L,R columns.');
assert(isequal(Labels2, {'Mode 2 L','Mode 2 R'}), 'Single-rank labels wrong.');

% --- Single-component case: labels drop the hemisphere suffix ---
Comp1 = [1;1;1];
Rank1 = [1;2;3];
[iRows3, Labels3, Hemi3] = view_eigenmodes_timeseries('GetBandTraces', Comp1, Rank1, 1, 2);
assert(isequal(iRows3(:)', [1 2]), 'Single-component rows wrong.');
assert(isequal(Labels3, {'Mode 1','Mode 2'}), 'Single-component labels must omit L/R.');
assert(isequal(Hemi3(:)', [0 0]), 'Single-component hemisphere id must be 0.');

% --- Asymmetric rank (rank 2 only on left) ---
CompA = [1;1;2];
RankA = [1;2;1];
[iRowsA, LabelsA] = view_eigenmodes_timeseries('GetBandTraces', CompA, RankA, 1, 2);
assert(isequal(iRowsA(:)', [1 3 2]), 'Asymmetric ordering wrong.');
assert(isequal(LabelsA, {'Mode 1 L','Mode 1 R','Mode 2 L'}), 'Asymmetric labels wrong.');

% --- Color helper: hemisphere -> distinct RGB cell array ---
colors = view_eigenmodes_timeseries('HemiColors', [1 2 0]);
assert(iscell(colors) && numel(colors) == 3, 'HemiColors must return a 1x3 cell.');
assert(~isequal(colors{1}, colors{2}), 'Left and right must differ in color.');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run (MATLAB MCP `run_matlab_test_file`, or from MATLAB):
```matlab
run('dev/tests/test_view_eigenmodes_timeseries_pure.m')
```
Expected: FAIL — `Undefined function 'view_eigenmodes_timeseries'` (file not created yet).

- [ ] **Step 3: Create the file with dispatch + the pure helpers**

Create `toolbox/gui/view_eigenmodes_timeseries.m`:

```matlab
function varargout = view_eigenmodes_timeseries(varargin)
% VIEW_EIGENMODES_TIMESERIES: Plot eigenmode coefficients theta_k(t) over time.
%
% USAGE:  hFig = view_eigenmodes_timeseries(DataFile)
%         [iRows, Labels, Hemi] = view_eigenmodes_timeseries('GetBandTraces', Component, CompRank, kLo, kHi)
%         colors = view_eigenmodes_timeseries('HemiColors', Hemi)
%         view_eigenmodes_timeseries('ModesChangedCallback', hFig)
%
% One trace per eigenmode coefficient (sensor->mode transform), for the paired
% ranks in the EigenModes panel band. Each paired rank yields a left and a right
% trace. The figure tracks the panel band live via bst_figures('FireModesChanged').
%
% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Diellor Basha, 2026

if (nargin >= 1) && ischar(varargin{1}) && any(strcmp(varargin{1}, {'GetBandTraces','HemiColors','ModesChangedCallback'}))
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
[varargout{1:nargout}] = ViewFigure(varargin{:});
end


%% ===== PURE: band (paired-rank) -> raw-column traces =====
% For each paired rank k in kLo:kHi, emit its left column(s) then right column(s).
% Labels carry an L/R suffix only when the data has two components.
function [iRows, Labels, Hemi] = GetBandTraces(Component, CompRank, kLo, kHi)
    Component = Component(:);
    CompRank  = CompRank(:);
    isPaired  = any(Component == 2);
    iRows = [];
    Labels = {};
    Hemi = [];
    for k = kLo:kHi
        % Left (component 1) then right (component 2)
        for c = find(CompRank == k & Component == 1)'
            iRows(end+1) = c; %#ok<AGROW>
            Hemi(end+1)  = 1; %#ok<AGROW>
            if isPaired
                Labels{end+1} = sprintf('Mode %d L', k); %#ok<AGROW>
            else
                Labels{end+1} = sprintf('Mode %d', k);   %#ok<AGROW>
                Hemi(end)     = 0;
            end
        end
        for c = find(CompRank == k & Component == 2)'
            iRows(end+1) = c; %#ok<AGROW>
            Hemi(end+1)  = 2; %#ok<AGROW>
            Labels{end+1} = sprintf('Mode %d R', k); %#ok<AGROW>
        end
    end
end


%% ===== PURE: hemisphere id -> RGB cell array (left warm, right cool) =====
function colors = HemiColors(Hemi)
    Hemi = Hemi(:)';
    cL = [0.85 0.33 0.10];   % left  = warm orange
    cR = [0.00 0.45 0.74];   % right = cool blue
    c0 = [0.20 0.20 0.20];   % single-component = neutral
    colors = cell(1, numel(Hemi));
    for i = 1:numel(Hemi)
        switch Hemi(i)
            case 1, colors{i} = cL;
            case 2, colors{i} = cR;
            otherwise, colors{i} = c0;
        end
    end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run:
```matlab
run('dev/tests/test_view_eigenmodes_timeseries_pure.m')
```
Expected: `ALL TESTS PASSED`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/view_eigenmodes_timeseries.m dev/tests/test_view_eigenmodes_timeseries_pure.m
git commit -m "Eigenmode time series: pure band->trace selector + color helper"
```

---

## Task 2: Entry point — transform, cache, first plot

**Files:**
- Modify: `toolbox/gui/view_eigenmodes_timeseries.m` (add `ViewFigure`)

This task has no automated unit test (it does DB/GUI I/O); it is validated by the smoke test in Task 5. Implement it directly, then sanity-check that the file still loads (`which view_eigenmodes_timeseries`).

- [ ] **Step 1: Add the `ViewFigure` function**

Append to `toolbox/gui/view_eigenmodes_timeseries.m` (before the final `GetBandTraces`/`HemiColors`, or after — MATLAB local-function order is free):

```matlab
%% ===== GUI: build the eigenmode coefficient time series figure =====
function hFig = ViewFigure(DataFile)
    hFig = [];
    if isempty(DataFile)
        bst_error('No data file provided.', 'Eigenmode time series', 0);
        return;
    end
    % ----- Study + head model + surface -----
    [sStudy, iStudy] = bst_get('AnyFile', DataFile); %#ok<ASGLU>
    if isempty(sStudy) || ~isfield(sStudy, 'iHeadModel') || isempty(sStudy.iHeadModel) || (sStudy.iHeadModel < 1)
        bst_error('No head model available for this study.', 'Eigenmode time series', 0);
        return;
    end
    HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;
    HeadModelMat  = in_bst_headmodel(HeadModelFile, 0, 'HeadModelType', 'SurfaceFile');
    if ~strcmpi(HeadModelMat.HeadModelType, 'surface')
        bst_error('Eigenmode transform requires a surface head model.', 'Eigenmode time series', 0);
        return;
    end
    SurfaceFile = HeadModelMat.SurfaceFile;

    % ----- Eigenmodes -----
    [Eig, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if ~isComputed
        bst_error(['No eigenmodes on this surface.' 10 'Run "Compute eigenmodes" first.'], 'Eigenmode time series', 0);
        return;
    end

    % ----- Constrained gain (fixed orientation: [nch x nVert]) -----
    HM   = in_bst_headmodel(HeadModelFile, 1);
    Gain = double(HM.Gain);
    if size(Gain, 2) ~= size(Eig.Vectors, 1)
        bst_error(sprintf(['Head model has %d vertices but eigenmodes have %d.' 10 ...
            'Recompute the head model.'], size(Gain,2), size(Eig.Vectors,1)), 'Eigenmode time series', 0);
        return;
    end

    % ----- Channels + recordings -----
    ChannelFile = bst_get('ChannelFileForStudy', sStudy.FileName);
    if isempty(ChannelFile)
        bst_error('No channel file found.', 'Eigenmode time series', 0);
        return;
    end
    ChannelMat = in_bst_channel(ChannelFile);
    DataMat    = in_bst_data(DataFile);
    if isstruct(DataMat.F)
        bst_error('Eigenmode time series requires imported (non-raw) recordings.', 'Eigenmode time series', 0);
        return;
    end
    if isfield(DataMat, 'ChannelFlag') && ~isempty(DataMat.ChannelFlag)
        ChannelFlag = DataMat.ChannelFlag;
    else
        ChannelFlag = ones(length(ChannelMat.Channel), 1);
    end
    iCh = good_channel(ChannelMat.Channel, ChannelFlag, 'MEG');
    if isempty(iCh)
        iCh = good_channel(ChannelMat.Channel, ChannelFlag, 'EEG');
    end
    if isempty(iCh)
        bst_error('No good MEG or EEG channels found.', 'Eigenmode time series', 0);
        return;
    end

    % ----- Transform: full coefficient matrix Theta [K_raw x nTime] -----
    nCh   = numel(iCh);
    K_raw = min(nCh, double(Eig.nModes));
    Phi   = double(Eig.Vectors(:, 1:K_raw));
    [Kernel, ~] = bst_eigenmodes_transform(Gain(iCh, :), Phi);   % [K_raw x nCh]
    Theta = Kernel * double(DataMat.F(iCh, :));                  % [K_raw x nTime]

    % ----- Cache everything the live refresh needs -----
    cache = struct( ...
        'SurfaceFile', SurfaceFile, ...
        'DataFile',    DataFile, ...
        'Theta',       Theta, ...
        'Component',   Eig.Component(1:K_raw), ...
        'CompRank',    Eig.CompRank(1:K_raw), ...
        'TimeVector',  DataMat.Time);

    % ----- Ensure the lever is initialized for this surface (paired ranks) -----
    Kp = double(max(cache.CompRank));
    st = panel_eigenmodes('GetState');
    if ~file_compare(st.SurfaceFile, SurfaceFile) || (st.nModes ~= Kp)
        panel_eigenmodes('ResetState', SurfaceFile, Kp);
    end
    band = panel_eigenmodes('GetState');
    band = band.Band;

    % ----- First plot -----
    [iRows, Labels, Hemi] = GetBandTraces(cache.Component, cache.CompRank, band(1), band(2));
    if isempty(iRows)
        bst_error('No eigenmodes in the selected band.', 'Eigenmode time series', 0);
        return;
    end
    F      = cache.Theta(iRows, :);
    colors = HemiColors(Hemi);
    hFig = view_timeseries_matrix(DataFile, {F}, cache.TimeVector, '', {'Eigenmode coefficients'}, Labels, colors, []);
    if isempty(hFig)
        return;
    end
    set(hFig, 'Name', ['Eigenmode time series: ' SurfaceFile]);
    setappdata(hFig, 'EigenTimeSeries', cache);

    % ----- Show + sync the panel -----
    gui_brainstorm('ShowToolTab', 'EigenModes');
    try
        panel_eigenmodes('RefreshControls');
    catch
        % Non-fatal: the panel still works, controls just won't pre-sync.
    end
end
```

- [ ] **Step 2: Verify the file parses and loads**

Run:
```matlab
clear functions; which('view_eigenmodes_timeseries')
```
Expected: prints the path to `toolbox/gui/view_eigenmodes_timeseries.m` with no parse error.

- [ ] **Step 3: Re-run the pure test (no regression)**

Run:
```matlab
run('dev/tests/test_view_eigenmodes_timeseries_pure.m')
```
Expected: `ALL TESTS PASSED` (adding `ViewFigure` must not break the dispatch).

- [ ] **Step 4: Commit**

```bash
git add toolbox/gui/view_eigenmodes_timeseries.m
git commit -m "Eigenmode time series: entry point (transform, cache, first plot)"
```

---

## Task 3: Live tracking — refresh callback + FireModesChanged dispatch

**Files:**
- Modify: `toolbox/gui/view_eigenmodes_timeseries.m` (add `ModesChangedCallback`)
- Modify: `toolbox/core/bst_figures.m` (`FireModesChanged`, per-figure loop ~line 1014)

- [ ] **Step 1: Add `ModesChangedCallback` to the viewer**

Append to `toolbox/gui/view_eigenmodes_timeseries.m`:

```matlab
%% ===== Re-select the band's rows and redraw on a lever change =====
function ModesChangedCallback(hFig) %#ok<DEFNU>
    cache = getappdata(hFig, 'EigenTimeSeries');
    if isempty(cache)
        return;
    end
    st   = panel_eigenmodes('GetState');
    band = st.Band;
    [iRows, Labels, Hemi] = GetBandTraces(cache.Component, cache.CompRank, band(1), band(2));
    if isempty(iRows)
        return;
    end
    F      = cache.Theta(iRows, :);
    colors = HemiColors(Hemi);
    % Redraw into the same figure (view_timeseries_matrix replots when hFig is given)
    view_timeseries_matrix(cache.DataFile, {F}, cache.TimeVector, '', {'Eigenmode coefficients'}, Labels, colors, hFig);
    % Re-assert our tag (redraw may overwrite figure appdata)
    setappdata(hFig, 'EigenTimeSeries', cache);
end
```

- [ ] **Step 2: Read the current dispatch loop to place the new branch**

Run:
```bash
sed -n '1010,1035p' toolbox/core/bst_figures.m
```
Expected: shows the per-figure loop with the line
`if strcmpi(get(sFig.hFigure, 'Visible'), 'off') || ~strcmpi(sFig.Id.Type, '3DViz')` followed by the `'Surface'` match and the `EigenView`/`UpdateSurfaceData` branch.

- [ ] **Step 3: Insert the time series client branch**

In `toolbox/core/bst_figures.m`, inside `FireModesChanged`'s per-figure loop, replace the combined skip line:

```matlab
            sFig = GlobalData.DataSet(iDS).Figure(iFig);
            if strcmpi(get(sFig.hFigure, 'Visible'), 'off') || ~strcmpi(sFig.Id.Type, '3DViz')
                continue;
            end
```

with this (split the visibility guard from the 3D-only guard, and add the time series client first):

```matlab
            sFig = GlobalData.DataSet(iDS).Figure(iFig);
            if strcmpi(get(sFig.hFigure, 'Visible'), 'off')
                continue;
            end
            % Eigenmode time series client (not a 3D figure): re-select band rows.
            etsInfo = getappdata(sFig.hFigure, 'EigenTimeSeries');
            if ~isempty(etsInfo)
                if file_compare(etsInfo.SurfaceFile, SurfaceFile)
                    view_eigenmodes_timeseries('ModesChangedCallback', sFig.hFigure);
                end
                continue;
            end
            if ~strcmpi(sFig.Id.Type, '3DViz')
                continue;
            end
```

- [ ] **Step 4: Verify the dispatch edit parses**

Run:
```matlab
clear functions; bst_figures('FireModesChanged');
```
Expected: returns with no error (it's a no-op when no eigenmode state/figures exist — the early `return` guards in `FireModesChanged` handle the empty case).

- [ ] **Step 5: Re-run the pure test (no regression)**

Run:
```matlab
run('dev/tests/test_view_eigenmodes_timeseries_pure.m')
```
Expected: `ALL TESTS PASSED`.

- [ ] **Step 6: Commit**

```bash
git add toolbox/gui/view_eigenmodes_timeseries.m toolbox/core/bst_figures.m
git commit -m "Eigenmode time series: live band tracking via FireModesChanged"
```

---

## Task 4: Launch menu item in the 3D figure popup

**Files:**
- Modify: `toolbox/gui/figure_3d.m` (`DisplayFigurePopup`, after the "View sources" item ~line 1665)

- [ ] **Step 1: Read the insertion point**

Run:
```bash
sed -n '1648,1690p' toolbox/gui/figure_3d.m
```
Expected: shows the `if ~isempty(DataFile) && ~ismember(Modality, ...)` block adding `Recordings`, `Topography`, `View sources`, ending before `% ==== MENU: 2DLAYOUT ====`.

- [ ] **Step 2: Add the menu item**

In `toolbox/gui/figure_3d.m`, immediately after the "View SOURCES" `gui_component` block (the one calling `bst_figures('ViewResults',hFig)`) and still inside the `if ~isempty(DataFile) && ~ismember(Modality, {...})` block, insert:

```matlab
        % === Eigenmode time series ===
        % Enabled when the study has a surface head model with computed eigenmodes.
        if ~isempty(sStudy) && isfield(sStudy, 'iHeadModel') && ~isempty(sStudy.iHeadModel) && (sStudy.iHeadModel >= 1)
            HmFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;
            HmMat  = in_bst_headmodel(HmFile, 0, 'HeadModelType', 'SurfaceFile');
            if strcmpi(HmMat.HeadModelType, 'surface')
                [~, isEig] = in_tess_eigenmodes(HmMat.SurfaceFile);
                if isEig
                    gui_component('MenuItem', jPopup, [], 'Eigenmode time series', IconLoader.ICON_TS_DISPLAY, [], @(h,ev)bst_call(@view_eigenmodes_timeseries, DataFile));
                end
            end
        end
```

Note: `sStudy` is already defined in this block (`sStudy = bst_get('AnyFile', DataFile);` a few lines above). `IconLoader` is imported at the top of `DisplayFigurePopup`.

- [ ] **Step 3: Verify it parses**

Run:
```matlab
clear functions; which('figure_3d')
```
Expected: prints the path with no parse error.

- [ ] **Step 4: Commit**

```bash
git add toolbox/gui/figure_3d.m
git commit -m "Eigenmode time series: launch item in 3D figure popup"
```

---

## Task 5: End-to-end smoke validation

**Files:**
- Create: `dev/tests/test_eigenmode_timeseries_e2e.m`

This is a guarded smoke harness (skips cleanly when no suitable protocol/data is loaded), mirroring `dev/tests/test_eigenmode_viewer_e2e.m`. It is run manually against a real protocol.

- [ ] **Step 1: Write the smoke harness**

Create `dev/tests/test_eigenmode_timeseries_e2e.m`:

```matlab
function test_eigenmode_timeseries_e2e
% Smoke test: open the eigenmode time series viewer on the currently selected
% data file, confirm a figure with the expected number of traces for the
% current band, then move the band and confirm the trace count changes.
% Requires: a loaded protocol with an imported data file whose study has a
% surface head model + computed eigenmodes. Skips cleanly otherwise.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% --- Resolve a usable data file from the current protocol ---
sStudy = bst_get('StudyWithCondition', []); %#ok<NASGU>
[sStudies, iStudies] = bst_get('ProtocolStudies'); %#ok<ASGLU>
DataFile = '';
for iS = 1:numel(iStudies)
    s = bst_get('Study', iStudies(iS));
    if ~isempty(s) && isfield(s,'iHeadModel') && ~isempty(s.iHeadModel) && (s.iHeadModel >= 1) ...
            && isfield(s,'Data') && ~isempty(s.Data)
        hm = in_bst_headmodel(s.HeadModel(s.iHeadModel).FileName, 0, 'HeadModelType', 'SurfaceFile');
        if strcmpi(hm.HeadModelType, 'surface')
            [~, isEig] = in_tess_eigenmodes(hm.SurfaceFile);
            if isEig
                DataFile = s.Data(1).FileName;
                break;
            end
        end
    end
end
if isempty(DataFile)
    disp('SKIP: no study with surface head model + eigenmodes + imported data.');
    return;
end

% --- Open the viewer ---
hFig = view_eigenmodes_timeseries(DataFile);
assert(~isempty(hFig) && ishandle(hFig), 'Viewer figure was not created.');
cache = getappdata(hFig, 'EigenTimeSeries');
assert(~isempty(cache) && isfield(cache, 'Theta'), 'Figure missing EigenTimeSeries cache.');

% --- Expected trace count for the current band ---
st = panel_eigenmodes('GetState');
[iRows0] = view_eigenmodes_timeseries('GetBandTraces', cache.Component, cache.CompRank, st.Band(1), st.Band(2));
nLines = numel(findobj(hFig, 'Type', 'line', 'Tag', 'DataLine'));
assert(nLines >= numel(iRows0) && nLines > 0, sprintf('Expected ~%d traces, found %d.', numel(iRows0), nLines));

% --- Narrow the band to a single rank, confirm tracking redraws fewer traces ---
panel_eigenmodes('SetWindowShape', 'single');
panel_eigenmodes('SetCurrentMode', st.Band(1));
[iRows1] = view_eigenmodes_timeseries('GetBandTraces', cache.Component, cache.CompRank, ...
    panel_eigenmodes('GetState').Band(1), panel_eigenmodes('GetState').Band(2));
assert(numel(iRows1) <= numel(iRows0), 'Single-mode band should not increase trace count.');

disp('ALL TESTS PASSED');
close(hFig);
end
```

- [ ] **Step 2: Run the smoke harness**

Open Brainstorm with a protocol that has eigenmodes + a head model + imported recordings, then run:
```matlab
run('dev/tests/test_eigenmode_timeseries_e2e.m')
```
Expected: `ALL TESTS PASSED` (or `SKIP: ...` if no suitable data — in which case validate manually: right-click a 3D source figure → "Eigenmode time series", confirm traces appear, toggle Butterfly/Column from the figure menu, and move the EigenModes slider to see traces update).

- [ ] **Step 3: Commit**

```bash
git add dev/tests/test_eigenmode_timeseries_e2e.m
git commit -m "Eigenmode time series: end-to-end smoke harness"
```

---

## Self-Review Notes

- **Spec coverage:** Component 1 (entry point) → Task 2; Component 2 (band→traces) → Task 1; Component 3 (display via `view_timeseries_matrix`) → Tasks 2/3; Component 4 (live tracking) → Task 3; Component 5 (launch menu) → Task 4. Edge cases (no eigenmodes, no head model, raw data, single component, empty band, channel selection) → guarded in Task 2 / tested in Task 1. Testing section → Tasks 1 & 5.
- **Raw-data scope:** the spec assumed imported recordings; Task 2 explicitly errors on raw (`isstruct(DataMat.F)`), which is called out here rather than left implicit.
- **Color format:** `LinesColor` is a cell array of `[1×3]` RGB (confirmed against `view_timeseries_matrix` header and `view_scouts` usage), not an `[n×3]` matrix.
