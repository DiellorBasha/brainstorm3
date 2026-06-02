# Eigenmode Time Series — Raw / Lazy-Load Enhancement Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Make the eigenmode time series viewer work on raw/continuous recordings by reading the figure's **current displayed time window** from the loaded dataset (the same lazy path "Display on cortex" uses), instead of `in_bst_data` on the raw link — and auto-refresh when the displayed window changes (raw page scroll).

**Architecture:** `ViewFigure` resolves the recordings dataset `iDS` from the data file (`bst_memory('LoadDataFile')`), builds the transform kernel once (`Kernel = pinv(Gain(iCh,:)·Phi)`), then reads the current window via `bst_memory('GetRecordingsValues', iDS, iCh, 'UserTimeWindow', 0)` (unscaled — matches lead-field units; lazily loads the raw page) and `Theta = Kernel·F`. The cache holds `Kernel`, `GoodChannel`, the window `Theta`/`TimeVector`, and the window bounds. A new `FireCurrentTimeChanged` branch calls a `CurrentTimeChangedCallback` that re-reads + recomputes `Theta` (kernel is cached) only when the window bounds change.

**Tech Stack:** MATLAB; `bst_memory` (GetRecordingsValues / GetTimeVector / LoadDataFile), `bst_figures` (FireCurrentTimeChanged), existing `bst_eigenmodes_transform`.

---

## Background: verified APIs

- `[F] = bst_memory('GetRecordingsValues', iDS, iChannel, iTime, isGradMagScale)` — returns `Measures.F(iChannel, iTimeWindow)`; lazily calls `LoadRecordingsMatrix(iDS)` if `Measures.F` is empty (raw → reads the current page). `iTime='UserTimeWindow'` → the displayed window; `isGradMagScale=0` → no display scaling (physical units, consistent with the gain).
- `[TimeVector, iTime] = bst_memory('GetTimeVector', iDS, [], 'UserTimeWindow')` — the window's time vector.
- `iDS = bst_memory('LoadDataFile', DataFile)` — dataset index for a data file (returns the already-loaded one if present). Used by `view_timeseries_matrix` internally, so the TS figure lands in this same `iDS`.
- `GlobalData.DataSet(iDS).Measures.Time` — `[t0 t1]` bounds of the loaded window (for raw, the current page); used for change detection.
- `GlobalData.DataSet(iDS).Channel`, `.Measures.ChannelFlag` — channel array + current good/bad flags.
- `FireCurrentTimeChanged` (bst_figures.m:944) — per-figure `switch` on `Id.Type`; `ResultsTimeSeries` already routes to `figure_timeseries('CurrentTimeChangedCallback', iDS, iFig)` for the cursor. We ADD our re-read after the switch (not replacing it).

---

## Task 1: Lazy current-window read in `ViewFigure` (+ refactor cache/helpers)

**Files:**
- Modify: `toolbox/gui/view_eigenmodes_timeseries.m`

Replace the data-loading and plotting internals of `ViewFigure`, refactor the band→plot logic into shared helpers, and update `ModesChangedCallback`. The pure helpers `GetBandTraces`/`HemiColors` are unchanged.

- [ ] **Step 1: Update the dispatch list** to add `CurrentTimeChangedCallback` and a shared `RefreshTraces` is internal (not dispatched). Change line 33's name list from
  `{'GetBandTraces','HemiColors','ModesChangedCallback'}` to
  `{'GetBandTraces','HemiColors','ModesChangedCallback','CurrentTimeChangedCallback'}`.

- [ ] **Step 2: Replace `ViewFigure`** (the entire current subfunction) with this version:

```matlab
%% ===== GUI: build the eigenmode coefficient time series figure =====
function hFig = ViewFigure(DataFile)
    global GlobalData;
    hFig = [];
    if isempty(DataFile) || ~ischar(DataFile)
        bst_error('No data file provided.', 'Eigenmode time series', 0);
        return;
    end
    % ----- Study + head model + surface -----
    [sStudy, ~] = bst_get('AnyFile', DataFile);
    if isempty(sStudy) || ~isfield(sStudy, 'iHeadModel') || isempty(sStudy.iHeadModel) ...
            || (sStudy.iHeadModel < 1) || (length(sStudy.HeadModel) < sStudy.iHeadModel)
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

    % ----- Load the recordings dataset (raw or imported); read the CURRENT window -----
    % This is the lazy path used by "Display on cortex": GetRecordingsValues loads
    % the current page on demand for raw files and returns the displayed window.
    iDS = bst_memory('LoadDataFile', DataFile);
    if isempty(iDS)
        bst_error('Could not load the recordings.', 'Eigenmode time series', 0);
        return;
    end
    Channels    = GlobalData.DataSet(iDS).Channel;
    ChannelFlag = GlobalData.DataSet(iDS).Measures.ChannelFlag;
    if isempty(Channels)
        bst_error('No channels found for this recording.', 'Eigenmode time series', 0);
        return;
    end
    if isempty(ChannelFlag)
        ChannelFlag = ones(length(Channels), 1);
    end
    iCh = good_channel(Channels, ChannelFlag, 'MEG');
    if isempty(iCh)
        iCh = good_channel(Channels, ChannelFlag, 'EEG');
    end
    if isempty(iCh)
        bst_error('No good MEG or EEG channels found.', 'Eigenmode time series', 0);
        return;
    end

    % ----- Transform kernel over ALL raw modes (complete pairing; rank-safe pinv) -----
    K_raw = double(Eig.nModes);
    Phi   = double(Eig.Vectors(:, 1:K_raw));
    [Kernel, ~] = bst_eigenmodes_transform(Gain(iCh, :), Phi);   % [K_raw x nCh]

    % ----- Current-window recordings -> coefficients (unscaled: isGradMagScale=0) -----
    F = bst_memory('GetRecordingsValues', iDS, iCh, 'UserTimeWindow', 0);   % [nCh x nTime]
    [TimeVector, ~] = bst_memory('GetTimeVector', iDS, [], 'UserTimeWindow');
    Theta = Kernel * F;                                                     % [K_raw x nTime]

    % ----- Cache (kernel + channels let us recompute on window change) -----
    cache = struct( ...
        'SurfaceFile', SurfaceFile, ...
        'DataFile',    DataFile, ...
        'Kernel',      Kernel, ...
        'GoodChannel', iCh, ...
        'Component',   Eig.Component(1:K_raw), ...
        'CompRank',    Eig.CompRank(1:K_raw), ...
        'Theta',       Theta, ...
        'TimeVector',  TimeVector, ...
        'WindowTime',  GlobalData.DataSet(iDS).Measures.Time);

    % ----- Ensure the lever is initialized for this surface (paired ranks) -----
    Kp = double(max(cache.CompRank));
    st = panel_eigenmodes('GetState');
    if ~file_compare(st.SurfaceFile, SurfaceFile) || (st.nModes ~= Kp)
        panel_eigenmodes('ResetState', SurfaceFile, Kp);
    end
    band = panel_eigenmodes('GetState');
    band = band.Band;

    % ----- First plot -----
    [F0, Labels, colors] = BandData(cache, band);
    if isempty(F0)
        bst_error('No eigenmodes in the selected band.', 'Eigenmode time series', 0);
        return;
    end
    hFig = view_timeseries_matrix(DataFile, {F0}, cache.TimeVector, '', {'Eigenmode coefficients'}, Labels, colors, []);
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


%% ===== PURE-ish: band -> {F, Labels, colors} from a cache's current Theta =====
function [F, Labels, colors] = BandData(cache, band)
    [iRows, Labels, Hemi] = GetBandTraces(cache.Component, cache.CompRank, band(1), band(2));
    if isempty(iRows)
        F = []; colors = {};
        return;
    end
    F      = cache.Theta(iRows, :);
    colors = HemiColors(Hemi);
end


%% ===== Re-slice the band and redraw into the existing figure =====
function RefreshTraces(hFig)
    cache = getappdata(hFig, 'EigenTimeSeries');
    if isempty(cache)
        return;
    end
    st = panel_eigenmodes('GetState');
    if ~file_compare(st.SurfaceFile, cache.SurfaceFile)
        return;   % panel currently driving a different surface
    end
    [F, Labels, colors] = BandData(cache, st.Band);
    if isempty(F)
        return;
    end
    view_timeseries_matrix(cache.DataFile, {F}, cache.TimeVector, '', {'Eigenmode coefficients'}, Labels, colors, hFig);
    % Restore our title (view_timeseries_matrix calls UpdateFigureName, which
    % overwrites it). The EigenTimeSeries appdata survives the redraw; re-assert
    % it defensively in case view_timeseries_matrix's behaviour ever changes.
    setappdata(hFig, 'EigenTimeSeries', cache);
    set(hFig, 'Name', ['Eigenmode time series: ' cache.SurfaceFile]);
end
```

- [ ] **Step 3: Replace `ModesChangedCallback`** with a thin wrapper:

```matlab
%% ===== Lever changed: re-slice the band (no data re-read) =====
function ModesChangedCallback(hFig) %#ok<DEFNU>
    RefreshTraces(hFig);
end
```

- [ ] **Step 4: Verify** (MATLAB MCP):
  - `check_matlab_code` on the file → no new issues.
  - `clear functions; cd('/Users/diellorbasha/workspace/research/code/brainstorm3'); if ~brainstorm('status'); brainstorm nogui; end; which('view_eigenmodes_timeseries')` → path, no parse error.
  - `run('dev/tests/test_view_eigenmodes_timeseries_pure.m')` → `ALL TESTS PASSED` (pure helpers unaffected).

- [ ] **Step 5: Commit**
```bash
git add toolbox/gui/view_eigenmodes_timeseries.m
git commit -m "Eigenmode time series: lazy current-window read (raw + imported)"
```

---

## Task 2: Auto-follow the displayed window (`CurrentTimeChangedCallback` + dispatch)

**Files:**
- Modify: `toolbox/gui/view_eigenmodes_timeseries.m` (add `CurrentTimeChangedCallback`)
- Modify: `toolbox/core/bst_figures.m` (`FireCurrentTimeChanged`)

- [ ] **Step 1: Add `CurrentTimeChangedCallback`** to `view_eigenmodes_timeseries.m`:

```matlab
%% ===== Displayed window changed (e.g. raw page scroll): re-read + recompute =====
function CurrentTimeChangedCallback(hFig, iDS) %#ok<DEFNU>
    global GlobalData;
    cache = getappdata(hFig, 'EigenTimeSeries');
    if isempty(cache)
        return;
    end
    if (nargin < 2) || isempty(iDS) || (iDS < 1) || (iDS > numel(GlobalData.DataSet))
        return;
    end
    % Only re-read when the displayed window bounds actually changed (a raw page
    % scroll). Moving the time cursor within the same page leaves Measures.Time
    % unchanged, so this is a cheap no-op on every cursor move.
    Win = GlobalData.DataSet(iDS).Measures.Time;
    if isequal(Win, cache.WindowTime)
        return;
    end
    F = bst_memory('GetRecordingsValues', iDS, cache.GoodChannel, 'UserTimeWindow', 0);
    [TimeVector, ~] = bst_memory('GetTimeVector', iDS, [], 'UserTimeWindow');
    cache.Theta      = cache.Kernel * F;
    cache.TimeVector = TimeVector;
    cache.WindowTime = Win;
    setappdata(hFig, 'EigenTimeSeries', cache);
    RefreshTraces(hFig);
end
```

- [ ] **Step 2: Read the dispatch loop** in `bst_figures.m`:
```bash
sed -n '949,990p' toolbox/core/bst_figures.m
```
Confirm the per-figure `switch (sFig.Id.Type)` ends with `end` before the closing `end` of the `iFig` loop.

- [ ] **Step 3: Insert the re-read after the switch.** In `FireCurrentTimeChanged`, immediately AFTER the `switch (sFig.Id.Type) ... end` block and BEFORE the loop's closing `end`, add:

```matlab
            % Eigenmode time series: re-read recordings when the displayed window
            % changes (raw page scroll). Runs in addition to the cursor update above.
            if ~isempty(getappdata(sFig.hFigure, 'EigenTimeSeries'))
                view_eigenmodes_timeseries('CurrentTimeChangedCallback', sFig.hFigure, iDS);
            end
```

(Match the existing indentation inside the `iFig` loop.)

- [ ] **Step 4: Verify** (MATLAB MCP):
  - `check_matlab_code` on both files → no new issues.
  - `clear functions; cd('/Users/diellorbasha/workspace/research/code/brainstorm3'); if ~brainstorm('status'); brainstorm nogui; end; bst_figures('FireCurrentTimeChanged'); disp('FireCurrentTimeChanged OK')` → prints OK (no-op on empty state).
  - `run('dev/tests/test_view_eigenmodes_timeseries_pure.m')` → `ALL TESTS PASSED`.

- [ ] **Step 5: Commit**
```bash
git add toolbox/gui/view_eigenmodes_timeseries.m toolbox/core/bst_figures.m
git commit -m "Eigenmode time series: auto-follow displayed window on FireCurrentTimeChanged"
```

---

## Task 3: Tests + docs

**Files:**
- Modify: `dev/tests/test_eigenmode_timeseries_e2e.m`
- Modify: `docs/superpowers/specs/2026-06-02-eigenmode-timeseries-design.md`

- [ ] **Step 1: Relax the e2e harness to allow raw data.** In `test_eigenmode_timeseries_e2e.m`, remove the inner `@raw`-skipping loop (raw is now supported) and pick the first data file:
  - Replace the inner `for iD = 1:numel(s.Data) ... end` block (and the `if ~isempty(DataFile), break; end`) with:
    ```matlab
                DataFile = s.Data(1).FileName;
                break;
    ```
  - Update the comment above it accordingly (raw is now supported via the lazy current-window read).
  Keep everything else (cleanup/onCleanup, band-collapse asserts) intact.

- [ ] **Step 2: Run the e2e** → expect `ALL TESTS PASSED` or a `SKIP:` line, never an error. Paste the result.

- [ ] **Step 3: Update the design spec.** In `2026-06-02-eigenmode-timeseries-design.md`:
  - In the **Data** section, change the source from "`in_bst_data` on the imported file" to "the figure's current displayed time window via `bst_memory('GetRecordingsValues', …)`, which lazily loads the raw page — works for both imported and raw recordings."
  - In **Out of scope**, REMOVE the line excluding raw recordings (now supported).
  - Add a short **Live tracking** note that the viewer follows the displayed window via `FireCurrentTimeChanged` (re-reads only on window-bounds change), in addition to the `FireModesChanged` band tracking.

- [ ] **Step 4: Commit**
```bash
git add dev/tests/test_eigenmode_timeseries_e2e.m docs/superpowers/specs/2026-06-02-eigenmode-timeseries-design.md
git commit -m "Eigenmode time series: e2e allows raw + spec reflects lazy window read"
```

---

## Self-Review Notes

- **Signature unchanged:** `view_eigenmodes_timeseries(DataFile)` still takes the data file; the menu in `figure_3d.m` is untouched. The viewer now resolves the recordings dataset itself and reads the current window — so raw works without changing the launch site.
- **Units:** `isGradMagScale=0` keeps recordings in physical units consistent with the lead field (the inverse/display path multiplies the un-rescaled `Measures.F` by the kernel; we mirror that).
- **Cheap cursor moves:** `CurrentTimeChangedCallback` returns immediately unless `Measures.Time` (window bounds) changed, so ordinary cursor motion within a page costs one struct compare.
- **Kernel cached:** window changes recompute only `Theta = Kernel·F`, not the SVD.
