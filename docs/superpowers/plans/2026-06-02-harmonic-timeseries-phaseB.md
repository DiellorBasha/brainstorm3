# Harmonic Eigenmode Inverse — Phase B (Consistent time series)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Re-point the eigenmode time series viewer to read a Harmonic results node's stored operator (`EigenKernel = M̃`), launch it from the node's database-tree menu, and remove the old figure-popup launch — so the traces are whitened, unregularized, and exactly consistent with the displayed Harmonic cortex map (`Φ·θ = ImagingKernel·D`).

**Architecture:** `view_eigenmodes_timeseries(ResultsFile)` loads the Harmonic node (`Function='eigenmode_harmonic'`) via `in_bst_results`, reads `M̃`/`GoodChannel`/`DataFile`/`SurfaceFile`, resolves the recordings dataset, reads the current window `D`, and caches `Theta = M̃·D`. All downstream band/page/reload tracking (built in the time-series feature) is unchanged. A tree-menu item launches it; the `figure_3d` popup item is removed.

**Tech Stack:** MATLAB; Brainstorm `in_bst_results`, `bst_memory` (LoadDataFile/GetRecordingsValues/GetTimeVector), `tree_callbacks` popup, `view_timeseries_matrix`, `in_tess_eigenmodes`.

**Depends on:** Phase A (merged) — Harmonic nodes carry `EigenKernel` (`[K×nGoodCh]`), `ImagingKernel = Φ·M̃` (`[nVert×nGoodCh]`), `Function='eigenmode_harmonic'`, `GoodChannel` (index vector), `DataFile`, `SurfaceFile`.

---

## Reference facts (verified)

- `view_eigenmodes_timeseries.m` dispatch (line ~33) routes string-first calls to `GetBandTraces/HemiColors/ModesChangedCallback/CurrentTimeChangedCallback/ReloadCallback`; otherwise `ViewFigure(varargin{:})`. `ViewFigure` currently takes a `DataFile` and computes its own unwhitened transform — THIS plan replaces its body to take a `ResultsFile`.
- The viewer's cache fields consumed downstream: `SurfaceFile, DataFile, Kernel, GoodChannel, Component, CompRank, Theta, TimeVector, WindowTime, Band`. Helpers `BandData`, `RefreshTraces`, `PlotBand`, `SyncWindow`, `ModesChangedCallback`, `CurrentTimeChangedCallback`, `ReloadCallback` all read `cache.Kernel`/`cache.GoodChannel`/`cache.Theta`/etc. — they stay UNCHANGED; only `ViewFigure` changes its data source.
- `Results = in_bst_results(ResultsFile, LoadFull, FieldsList...)` — with `LoadFull=0` it does NOT multiply kernel×data; pass specific field names to read just those.
- `bst_memory('GetRecordingsValues', iDS, iChannel, 'UserTimeWindow', 0)` reads the current displayed window for the good channels, unscaled (physical units, matching the kernel).
- `tree_callbacks.m` results-node popup: `case {'results','link'}` at line ~1801. In scope: `filenameRelative` (the results file), `sStudy`, `iResult`, `sSubject`. A "Cortical activations" submenu `jMenuActivations` is built (~line 1823); the "View eigenspectrum" item is added to it (~line 1836).
- `figure_3d.m` lines ~1667–1685: the `% === Eigenmode time series ===` block (the launch added in the time-series feature) — to be removed.

---

## File structure

| File | Responsibility | Action |
|------|----------------|--------|
| `toolbox/gui/view_eigenmodes_timeseries.m` | `ViewFigure(ResultsFile)` reads the Harmonic node | Modify (rewrite `ViewFigure` body) |
| `toolbox/tree/tree_callbacks.m` | "Eigenmode time series" on a Harmonic results node | Modify |
| `toolbox/gui/figure_3d.m` | remove the old popup launch | Modify |
| `dev/tests/test_eigenmode_timeseries_e2e.m` | open via a Harmonic node + existing tracking asserts | Modify |

---

## Task 1: Re-point `ViewFigure` to read the Harmonic node

**Files:**
- Modify: `toolbox/gui/view_eigenmodes_timeseries.m`

- [ ] **Step 1: Replace the entire `ViewFigure` function** (from `function hFig = ViewFigure(DataFile)` through its closing `end`) with:

```matlab
%% ===== GUI: build the eigenmode coefficient time series from a Harmonic node =====
function hFig = ViewFigure(ResultsFile)
    global GlobalData;
    hFig = [];
    if isempty(ResultsFile) || ~ischar(ResultsFile)
        bst_error('Open this from an "Eigenmode HARMONIC" results node.', 'Eigenmode time series', 0);
        return;
    end
    % ----- Load the Harmonic node (kernel only: LoadFull=0 does NOT multiply by data) -----
    ResMat = in_bst_results(ResultsFile, 0, 'Function', 'EigenKernel', 'GoodChannel', 'DataFile', 'SurfaceFile');
    if isempty(ResMat) || ~isfield(ResMat, 'Function') || ~strcmpi(ResMat.Function, 'eigenmode_harmonic') ...
            || ~isfield(ResMat, 'EigenKernel') || isempty(ResMat.EigenKernel)
        bst_error(['This is not an Eigenmode HARMONIC results node.' 10 ...
                   'Compute sources with method "Harmonic (eigenmodes)" first.'], 'Eigenmode time series', 0);
        return;
    end
    Mtilde      = double(ResMat.EigenKernel);     % [K x nGoodCh]
    GoodChannel = ResMat.GoodChannel;             % index vector into the channel file
    DataFile    = ResMat.DataFile;
    SurfaceFile = ResMat.SurfaceFile;
    if isempty(DataFile)
        bst_error('This Harmonic node has no associated recordings to project.', 'Eigenmode time series', 0);
        return;
    end

    % ----- Eigenmodes for paired-rank trace mapping (first K columns) -----
    [Eig, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if ~isComputed
        bst_error('No eigenmodes on this surface.', 'Eigenmode time series', 0);
        return;
    end
    K = size(Mtilde, 1);
    if size(Eig.Vectors, 2) < K
        bst_error('Eigenmode count is smaller than the Harmonic kernel rank.', 'Eigenmode time series', 0);
        return;
    end

    % ----- Load the recordings dataset (raw or imported); read the CURRENT window -----
    iDS = bst_memory('LoadDataFile', DataFile);
    if isempty(iDS)
        bst_error('Could not load the recordings.', 'Eigenmode time series', 0);
        return;
    end
    F = bst_memory('GetRecordingsValues', iDS, GoodChannel, 'UserTimeWindow', 0);   % [nGoodCh x nTime]
    [TimeVector, ~] = bst_memory('GetTimeVector', iDS, [], 'UserTimeWindow');
    Theta = Mtilde * F;                                                             % [K x nTime]

    % ----- Cache (M̃ + GoodChannel let SyncWindow recompute on a window change) -----
    cache = struct( ...
        'SurfaceFile', SurfaceFile, ...
        'DataFile',    DataFile, ...
        'Kernel',      Mtilde, ...
        'GoodChannel', GoodChannel, ...
        'Component',   Eig.Component(1:K), ...
        'CompRank',    Eig.CompRank(1:K), ...
        'Theta',       Theta, ...
        'TimeVector',  TimeVector, ...
        'WindowTime',  GlobalData.DataSet(iDS).Measures.Time);  % [t0 t1]; window-change detection

    % ----- Ensure the lever is initialized for this surface (paired ranks) -----
    Kp = double(max(cache.CompRank));
    st = panel_eigenmodes('GetState');
    if ~file_compare(st.SurfaceFile, SurfaceFile) || (st.nModes ~= Kp)
        panel_eigenmodes('ResetState', SurfaceFile, Kp);
    end
    band = panel_eigenmodes('GetState');
    band = band.Band;
    cache.Band = band;   % last-rendered band; lets a reload replot without the panel

    % ----- First plot -----
    [Fband, Labels, colors] = BandData(cache, band);
    if isempty(Fband)
        bst_error('No eigenmodes in the selected band.', 'Eigenmode time series', 0);
        return;
    end
    hFig = view_timeseries_matrix(DataFile, {Fband}, cache.TimeVector, '', {'Eigenmode coefficients'}, Labels, colors, []);
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

Do NOT change the dispatch line or any other subfunction (`BandData`, `RefreshTraces`, `PlotBand`, `SyncWindow`, `ModesChangedCallback`, `CurrentTimeChangedCallback`, `ReloadCallback`, `GetBandTraces`, `HemiColors`). They already operate on `cache.Kernel`/`cache.GoodChannel`/`cache.Theta`/`cache.Component`/`cache.CompRank`/`cache.Band`, which this new `ViewFigure` still populates with the same shapes (`Kernel` is now `M̃` `[K×nGoodCh]`; `Theta` is `[K×nTime]`; `GoodChannel` an index vector; `Component`/`CompRank` length `K`).

- [ ] **Step 2: Verify it parses + pure test unaffected.** (MATLAB MCP)
  - `check_matlab_code` on `toolbox/gui/view_eigenmodes_timeseries.m` → no new issues.
  - `clear functions; cd('/Users/diellorbasha/workspace/research/code/brainstorm3'); if ~brainstorm('status'); brainstorm nogui; end; which('view_eigenmodes_timeseries')` → path, no parse error.
  - `run('dev/tests/test_view_eigenmodes_timeseries_pure.m')` → `ALL TESTS PASSED` (the pure helpers are untouched).

- [ ] **Step 3: Commit**
```bash
git add toolbox/gui/view_eigenmodes_timeseries.m
git commit -m "Harmonic time series: viewer reads the Harmonic node's EigenKernel (consistent, whitened)"
```

---

## Task 2: Launch from the Harmonic results node (DB tree)

**Files:**
- Modify: `toolbox/tree/tree_callbacks.m`

- [ ] **Step 1: Add the gated menu item.** In the `case {'results', 'link'}` block (~line 1801), find the "View eigenspectrum" item inside `jMenuActivations`:
```matlab
                    % === VIEW EIGENSPECTRUM ===
                    if ismember(sStudy.Result(iResult).HeadModelType, {'surface'}) && ~isempty(sSubject) && ~isempty(sSubject.iCortex)
                        gui_component('MenuItem', jMenuActivations, [], 'View eigenspectrum', IconLoader.ICON_TIMEFREQ, [], @(h,ev)view_eigenmode_spectrum(filenameRelative));
                    end
```
Immediately AFTER that `end`, insert:
```matlab
                    % === EIGENMODE TIME SERIES (Harmonic nodes only) ===
                    if (length(bstNodes) == 1) && ismember(sStudy.Result(iResult).HeadModelType, {'surface'})
                        isHarm = false;
                        try
                            ResHdr = in_bst_results(filenameRelative, 0, 'Function');
                            isHarm = ~isempty(ResHdr) && isfield(ResHdr, 'Function') && strcmpi(ResHdr.Function, 'eigenmode_harmonic');
                        catch
                            isHarm = false;
                        end
                        if isHarm
                            gui_component('MenuItem', jMenuActivations, [], 'Eigenmode time series', IconLoader.ICON_TS_DISPLAY, [], @(h,ev)bst_call(@view_eigenmodes_timeseries, filenameRelative));
                        end
                    end
```
This loads only the `Function` field (cheap), guards against unreadable files, and only shows the item for `eigenmode_harmonic` nodes. `filenameRelative`, `sStudy`, `iResult`, `bstNodes`, `IconLoader` are all in scope here.

- [ ] **Step 2: Verify it parses.** (MATLAB MCP)
  - `check_matlab_code` on `toolbox/tree/tree_callbacks.m` → no NEW issues from the added lines (the file has pre-existing warnings).
  - `clear functions; which('tree_callbacks')` → path, no parse error.

- [ ] **Step 3: Commit**
```bash
git add toolbox/tree/tree_callbacks.m
git commit -m "Harmonic time series: launch from the Harmonic results node (DB tree menu)"
```

---

## Task 3: Remove the old figure-popup launch

**Files:**
- Modify: `toolbox/gui/figure_3d.m`

- [ ] **Step 1: Remove the block.** In `DisplayFigurePopup`, delete the entire `% === Eigenmode time series ===` block (it currently sits right after the "View sources" item and before `% === VIEW SPECTRUM ===`). Remove exactly this block:
```matlab
        % === Eigenmode time series ===
        % Enabled when the study has a surface head model with computed eigenmodes.
        % Wrapped in try/catch: a missing/corrupt head model or surface file must
        % suppress this item silently, never break the whole popup.
        if ~isempty(sStudy) && isfield(sStudy, 'iHeadModel') && ~isempty(sStudy.iHeadModel) ...
                && (sStudy.iHeadModel >= 1) && (length(sStudy.HeadModel) >= sStudy.iHeadModel)
            try
                HmFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;
                HmMat  = in_bst_headmodel(HmFile, 0, 'HeadModelType', 'SurfaceFile');
                if strcmpi(HmMat.HeadModelType, 'surface')
                    [~, isEig] = in_tess_eigenmodes(HmMat.SurfaceFile);
                    if isEig
                        gui_component('MenuItem', jPopup, [], 'Eigenmode time series', IconLoader.ICON_TS_DISPLAY, [], @(h,ev)bst_call(@view_eigenmodes_timeseries, DataFile));
                    end
                end
            catch
                % Non-fatal: skip the item if the head model / surface can't be read.
            end
        end
```
Leave the surrounding "View sources" item and "VIEW SPECTRUM" block untouched.

- [ ] **Step 2: Verify it parses.** (MATLAB MCP)
  - `check_matlab_code` on `toolbox/gui/figure_3d.m` → no new issues; `grep` confirms `view_eigenmodes_timeseries` no longer appears in `figure_3d.m`.
  - `clear functions; which('figure_3d')` → path, no parse error.

- [ ] **Step 3: Commit**
```bash
git add toolbox/gui/figure_3d.m
git commit -m "Harmonic time series: remove the superseded 3D-figure-popup launch"
```

---

## Task 4: Update the e2e to open via a Harmonic node

**Files:**
- Modify: `dev/tests/test_eigenmode_timeseries_e2e.m`

The existing e2e opened the viewer with `view_eigenmodes_timeseries(DataFile)` — that signature is gone. Update it to compute a temporary Harmonic node, open the viewer from it, run the existing tracking assertions, then delete the temp node.

- [ ] **Step 1: Replace the data-resolution + open block.** The current test finds a `DataFile` and calls `hFig = view_eigenmodes_timeseries(DataFile)`. Replace the discovery + open section (everything from the `sProtocol = bst_get('ProtocolStudies')` discovery through the `hFig = view_eigenmodes_timeseries(...)` line, but BEFORE the `cache = getappdata(...)` assertions) with this — it additionally requires a noise cov and builds+saves a temporary Harmonic node:

```matlab
sProtocol = bst_get('ProtocolStudies');
if isempty(sProtocol) || ~isfield(sProtocol, 'Study') || isempty(sProtocol.Study)
    disp('SKIP: no protocol loaded.');
    return;
end
% Find a study with surface head model + eigenmodes + noise cov + a data file
found = [];
for iS = 1:numel(sProtocol.Study)
    s = sProtocol.Study(iS);
    if ~isfield(s,'iHeadModel') || isempty(s.iHeadModel) || (s.iHeadModel < 1) ...
            || (length(s.HeadModel) < s.iHeadModel) || ~isfield(s,'Data') || isempty(s.Data) ...
            || ~isfield(s,'NoiseCov') || isempty(s.NoiseCov) || isempty(s.NoiseCov(1).FileName)
        continue;
    end
    try
        hm = in_bst_headmodel(s.HeadModel(s.iHeadModel).FileName, 0, 'HeadModelType', 'SurfaceFile');
        if strcmpi(hm.HeadModelType, 'surface')
            [~, isEig] = in_tess_eigenmodes(hm.SurfaceFile);
            if isEig
                found = struct('hmFile', s.HeadModel(s.iHeadModel).FileName, ...
                    'surf', hm.SurfaceFile, 'ncFile', s.NoiseCov(1).FileName, ...
                    'dataFile', s.Data(1).FileName);
                break;
            end
        end
    catch
    end
end
if isempty(found)
    disp('SKIP: no study with surface head model + eigenmodes + noise cov + data.');
    return;
end
[sStudyT, iStudyT] = bst_get('AnyFile', found.dataFile);

% Build a temporary Harmonic node via the engine (mirrors what the GUI saves)
HM = in_bst_headmodel(found.hmFile, 1);                 % constrained [nch x nVert]
goodMask = all(isfinite(double(HM.Gain)), 2);
[InvE, errE] = bst_inverse_eigenmodes(found.hmFile, found.surf, found.ncFile, ...
    'Method', 'harmonic', 'nModes', 0, 'GoodChannel', goodMask);
assert(isempty(errE), ['harmonic inverse failed: ' errE]);
Kh   = InvE.nModes;
[EigT, ~] = in_tess_eigenmodes(found.surf);
PhiT = double(EigT.Vectors(:, 1:Kh));
ResMat = db_template('resultsmat');
ResMat.ImagingKernel = PhiT * InvE.ImagingKernel;       % [nVert x nGoodCh]
ResMat.ImageGridAmp  = [];
ResMat.nComponents   = 1;
ResMat.Function      = 'eigenmode_harmonic';
ResMat.EigenKernel   = InvE.ImagingKernel;              % M̃ [K x nGoodCh]
ResMat.GoodChannel   = InvE.GoodChannel;
ResMat.Whitener      = InvE.Whitener;
ResMat.SurfaceFile   = found.surf;
ResMat.HeadModelFile = found.hmFile;
ResMat.HeadModelType = 'surface';
ResMat.DataFile      = found.dataFile;
ResMat.Comment       = 'TEST Eigenmode HARMONIC';
ResMat.Time          = [];
StudyDir   = bst_fileparts(file_fullpath(sStudyT.FileName));
NodeFile   = bst_process('GetNewFilename', StudyDir, 'results_TESTharmonic');
bst_save(NodeFile, ResMat, 'v6');
db_add_data(iStudyT, NodeFile, ResMat);
nodeRel = file_short(NodeFile);
% Always clean up the temp node + figure, however the test exits
cleanupNode = onCleanup(@() cleanupHarmonic(NodeFile, iStudyT));

% Open the viewer FROM the Harmonic node
hFig = view_eigenmodes_timeseries(nodeRel);
```

- [ ] **Step 2: Add a `θ = M̃·D` consistency assertion** right after the existing `cache` integrity asserts (after the line that asserts `size(cache.Theta,2) == numel(cache.TimeVector)`), insert:

```matlab
% Consistency: the cached Theta must equal M̃ * D for the current window
iDSt = bst_memory('LoadDataFile', found.dataFile);
Dwin = bst_memory('GetRecordingsValues', iDSt, cache.GoodChannel, 'UserTimeWindow', 0);
assert(max(abs(cache.Theta(:) - reshape(InvE.ImagingKernel * Dwin, [], 1))) < 1e-6 * (1 + max(abs(cache.Theta(:)))), ...
    'Cached Theta must equal EigenKernel * D(window).');
```

- [ ] **Step 3: Add the cleanup helper** at the END of the test file (after the existing `restorePanel` local function), add:

```matlab
% Delete the temporary Harmonic node and reload the study (DB stays consistent).
function cleanupHarmonic(NodeFile, iStudy)
    try
        if exist(file_fullpath(NodeFile), 'file')
            file_delete(file_fullpath(NodeFile), 1);
            db_reload_studies(iStudy);
        end
    catch
        % Best-effort cleanup.
    end
end
```

(The existing `onCleanup(@() restorePanel(hFig, st0))` and band/page/reload assertions remain; they now run against the node-launched figure. Keep them all.)

- [ ] **Step 4: Run the e2e** → expect `ALL TESTS PASSED` (if suitable data) or a `SKIP:` line. Must not error, and must clean up the temp node.
```matlab
run('dev/tests/test_eigenmode_timeseries_e2e.m')
```

- [ ] **Step 5: Commit**
```bash
git add dev/tests/test_eigenmode_timeseries_e2e.m
git commit -m "Harmonic time series: e2e opens via a Harmonic node + asserts theta = M̃·D"
```

---

## Self-Review Notes

- **Spec coverage:** unit 5 (viewer reads node) → Task 1; unit 6 (tree launch + remove figure popup) → Tasks 2–3. The consistency contract `Φ·θ = ImagingKernel·D` is exercised by Task 4's `θ = M̃·D` assertion (and `ImagingKernel = Φ·M̃` from Phase A).
- **Unchanged tracking machinery:** `SyncWindow`/`Reload`/`ModesChanged`/`CurrentTimeChanged` keep working because `cache.Kernel` is still the operator applied to the window's recordings (now `M̃` instead of the old transform) and all other cache fields keep their shapes.
- **GoodChannel:** the node stores an index vector; `GetRecordingsValues(iDS, GoodChannel, …)` reads exactly those channels, and `M̃` columns correspond to them (built on `Gain(GoodChannel,:)` whitened) — so `M̃·D` is dimensionally and physically consistent.
- **Cleanup discipline:** the temp node is removed with `file_delete` + `db_reload_studies` (never a raw `.mat` delete), per the project's DB-deletion rule.
