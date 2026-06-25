# Dynamics Bidirectional Time/Frequency Focus — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give the Dynamics panel an ephemeral "focus" state for the time and frequency axes, bound bidirectionally to the native selection box on the recording time-series figure and the PSD spectrum figure.

**Architecture:** The panel drives the native selection on each figure (panel→view) and listens for user edits via a `NotifySelection` dispatcher the figures call from their mouse-up handlers (view→panel). A re-entrancy guard (`i_driving`) prevents feedback loops. A PSD lifecycle helper guarantees a spectrum figure exists when a band is selected.

**Tech Stack:** MATLAB, Brainstorm GUI (`panel_bst_dynamics.m`, `figure_timeseries.m`, `figure_spectrum.m`), Brainstorm process system (`process_psd`), MATLAB MCP for verification.

## Global Constraints

- No unit-test framework exists in Brainstorm. Verify each task with: (1) static lint via the MATLAB MCP `check_matlab_code` on every edited file (zero NEW errors), (2) assertable MATLAB snippets for pure helpers, (3) a scripted/manual Brainstorm session for GUI behavior. The final task is a manual end-to-end walkthrough.
- Validation recording: **`S01_AEF_01_notch`**, alpha band **7–13 Hz**, right parieto-occipital ~10.55 Hz burst at ~22.6 s.
- PSD auto-compute: `process_psd`, averaged, `Measure='magnitude'`, MEG sensors, spectrum X-axis fixed to **0–60 Hz**. Triggered ONLY on band-preset selection when no PSD figure/file exists.
- Time selection is **single-window**; multi-window navigation is served by the Atoms tree, not by drawing N boxes.
- Figure-file hooks MUST be inert for non-Dynamics users: gate every hook on a cheap `~isempty(getappdata(0,'DynamicsTarget'))` check before calling into the panel, and wrap the call in `try/catch`.
- The panel programmatically drives selections via `figure_timeseries('SetTimeSelectionManual', ...)` and `figure_spectrum('SetFreqSelection', ...)`. These do NOT route through mouse-up, so they cannot echo; the `i_driving` guard is belt-and-suspenders for any redraw-triggered path.
- Follow existing panel idioms: `[ctrl,st]=i_cs()`, `setappdata(0,'DynamicsTarget',st)`, `i_field(st,name,default)`, `bst_get('PanelControls','Dynamics')`.

---

## File Structure

- **Modify** `toolbox/gui/panel_bst_dynamics.m` — all panel-side logic: guard, dispatcher, ownership helpers, PSD lifecycle, frequency drive/sync, time drive/sync, tree staged-window nodes. (Single file, matches existing convention of keeping the whole panel together.)
- **Modify** `toolbox/gui/figure_spectrum.m` — one hook in `FigureMouseUpCallback` (freq sync-back).
- **Modify** `toolbox/gui/figure_timeseries.m` — one hook in `FigureMouseUpCallback` (time sync-back).

---

## Task 1: Focus infrastructure — guard, dispatcher, ownership, band match

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (add helpers near `i_cs` at line 544; extend `SetTarget` struct at line 881)

**Interfaces:**
- Produces:
  - `i_driving(tf)` set / `i_driving()` get → logical re-entrancy flag.
  - `panel_bst_dynamics('NotifySelection', hFig, axis, range)` — `axis ∈ {'time','freq'}`, `range=[lo hi]`.
  - `i_band_match(lo,hi)` → preset name (`'alpha'`…) or `'custom'`.
  - `i_owns_rec(st,hFig)`, `i_owns_spec(st,hFig)`, `i_rec_figure(st)` → logical / handle.
  - `st.hSpec`, `st.focusTime`, `st.detSel` fields initialized to `[]`.

- [ ] **Step 1: Initialize new state fields in `SetTarget`**

In `toolbox/gui/panel_bst_dynamics.m`, line 881-883, add the three fields to the struct:

```matlab
    setappdata(0, 'DynamicsTarget', struct('hFig',hFig, 'T',T, 'file',file, 'curGroup',0, ...
        'nodeList',{ {} }, 'nodeInfo',[], 'occMap',[], 'Lambda',[], 'showPhase',[1 1 1 1], ...
        'hSpec',[], 'focusTime',[], 'detSel',[], ...
        'nav', bst_dynamics('NewGroup', 'cursor')));
```

- [ ] **Step 2: Add the focus infrastructure helpers**

In `toolbox/gui/panel_bst_dynamics.m`, immediately after `i_cs` (after line 548), insert:

```matlab
%% ===== BIDIRECTIONAL FOCUS: re-entrancy guard =====
% True while the panel is DRIVING a figure selection (panel->view), so a redraw-triggered
% hook cannot echo back (view->panel) and create a feedback loop.
function tf = i_driving(varargin)
    persistent FLAG;
    if isempty(FLAG), FLAG = false; end
    if (nargin >= 1), FLAG = logical(varargin{1}); end
    tf = FLAG;
end

%% ===== NOTIFY SELECTION (view -> panel) =====
% Called by figure_timeseries (axis='time') / figure_spectrum (axis='freq') on mouse-up
% when the user edited a native selection box. range=[lo hi] in seconds (time) or Hz (freq).
% No-op unless a Dynamics session owns the notifying figure and we are not mid-drive.
function NotifySelection(hFig, axis, range) %#ok<DEFNU>
    if i_driving(), return; end
    st = getappdata(0, 'DynamicsTarget');
    if isempty(st), return; end
    if isempty(range) || (numel(range) < 2) || any(~isfinite(range)), return; end
    range = sort(double(range(:)'));
    switch axis
        case 'freq'
            if ~i_owns_spec(st, hFig), return; end
            i_sync_freq(st, range);
        case 'time'
            if ~i_owns_rec(st, hFig), return; end
            i_sync_time(st, range);
    end
end

%% ===== FOCUS OWNERSHIP / FIGURE LOOKUP =====
% The recording time-series figure for this Dynamics session (matches st.T.DataFile);
% falls back to any DataTimeSeries figure (the time selection is linked across them anyway).
function hFig = i_rec_figure(st)
    global GlobalData; %#ok<TLEV>
    hFig = [];
    if isempty(st) || ~isfield(st,'T') || isempty(st.T) || isempty(st.T.DataFile), return; end
    hAll = bst_figures('GetFiguresByType', {'DataTimeSeries'});
    for h = hAll(:)'
        [~,~,iDS] = bst_figures('GetFigure', h);
        if ~isempty(iDS) && ~isempty(GlobalData.DataSet(iDS).DataFile) && file_compare(GlobalData.DataSet(iDS).DataFile, st.T.DataFile)
            hFig = h;  return;
        end
    end
    if ~isempty(hAll), hFig = hAll(1); end
end
function tf = i_owns_rec(st, hFig)
    global GlobalData; %#ok<TLEV>
    tf = false;
    if isempty(hFig) || ~ishandle(hFig) || isempty(st.T.DataFile), return; end
    [~,~,iDS] = bst_figures('GetFigure', hFig);
    if isempty(iDS) || isempty(GlobalData.DataSet(iDS).DataFile), return; end
    tf = file_compare(GlobalData.DataSet(iDS).DataFile, st.T.DataFile);
end
function tf = i_owns_spec(st, hFig)
    tf = isfield(st,'hSpec') && ~isempty(st.hSpec) && ishandle(st.hSpec) && isequal(double(st.hSpec), double(hFig));
end

%% ===== BAND MATCH: a [lo hi] range -> preset name or 'custom' =====
function nm = i_band_match(lo, hi)
    nm = 'custom';
    b = i_bands();
    for k = 1:size(b,1)
        if (abs(b{k,2}(1) - lo) < 0.51) && (abs(b{k,2}(2) - hi) < 0.51), nm = b{k,1};  return; end
    end
end
```

- [ ] **Step 3: Lint the edited file (MATLAB MCP)**

Run: `check_matlab_code` on `toolbox/gui/panel_bst_dynamics.m`
Expected: no NEW errors (pre-existing Brainstorm-idiom warnings such as `#ok` suppressions are acceptable). `i_sync_freq`/`i_sync_time` are referenced but defined in later tasks — a "function not found" is NOT reported by checkcode (it is static, single-file), so expect none.

- [ ] **Step 4: Assert the pure helpers (MATLAB MCP)**

Run (`evaluate_matlab_code`):

```matlab
assert(strcmp(panel_bst_dynamics('i_band_match', 8, 13), 'alpha'), 'alpha match');
assert(strcmp(panel_bst_dynamics('i_band_match', 8.4, 12.7), 'alpha'), 'alpha tolerance');
assert(strcmp(panel_bst_dynamics('i_band_match', 9, 11), 'custom'), 'custom');
panel_bst_dynamics('i_driving', true);  assert(panel_bst_dynamics('i_driving'), 'guard on');
panel_bst_dynamics('i_driving', false); assert(~panel_bst_dynamics('i_driving'), 'guard off');
disp('TASK1_OK');
```

Expected: prints `TASK1_OK` with no assertion error. (Local subfunctions are reachable via `eval(macro_method)` dispatch.)

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m
git commit -m "feat(dynamics): focus infra - guard, NotifySelection dispatcher, ownership, band match

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: PSD lifecycle helper

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (add helpers after the Task 1 block)

**Interfaces:**
- Consumes: `st.T.DataFile`.
- Produces: `i_ensure_psd(st)` → spectrum figure handle (`[]` on failure); `i_find_psd_file(DataFile)` → timefreq file path or `''`; `i_compute_psd(DataFile)` → timefreq file path or `''`.

- [ ] **Step 1: Add the PSD lifecycle helpers**

In `toolbox/gui/panel_bst_dynamics.m`, after the Task 1 helpers, insert:

```matlab
%% ===== PSD LIFECYCLE: find-or-open-or-compute the spectrum figure for this recording =====
function hSpec = i_ensure_psd(st)
    hSpec = [];
    DataFile = st.T.DataFile;
    if isempty(DataFile), return; end
    % 1) already-open Spectrum figure for this recording?
    hAll = bst_figures('GetFiguresByType', 'Spectrum');
    for h = hAll(:)'
        TfInfo = getappdata(h, 'Timefreq');
        if isempty(TfInfo) || ~isfield(TfInfo,'FileName') || isempty(TfInfo.FileName), continue; end
        [~,~,~,~,sTf] = bst_get('TimefreqFile', TfInfo.FileName);
        if ~isempty(sTf) && ~isempty(sTf.DataFile) && file_compare(sTf.DataFile, DataFile)
            hSpec = h;  i_fix_spec_xlim(hSpec);  return;
        end
    end
    % 2) precomputed PSD timefreq file for this recording?
    TfFile = i_find_psd_file(DataFile);
    % 3) else compute one
    if isempty(TfFile), TfFile = i_compute_psd(DataFile); end
    if isempty(TfFile), return; end
    hSpec = view_spectrum(TfFile, 'Spectrum');
    i_fix_spec_xlim(hSpec);
end

% Fix the spectrum X-axis to 0-60 Hz (focus convention).
function i_fix_spec_xlim(hSpec)
    if isempty(hSpec) || ~ishandle(hSpec), return; end
    hAxes = findobj(hSpec, '-depth', 1, 'Tag', 'AxesGraph');
    if ~isempty(hAxes), try, set(hAxes, 'XLim', [0 60]); catch, end; end %#ok<CTCH>
end

% Find an existing PSD timefreq file associated with the recording (Comment contains "PSD").
function TfFile = i_find_psd_file(DataFile)
    TfFile = '';
    sStudy = bst_get('AnyFile', DataFile);
    if isempty(sStudy) || ~isfield(sStudy,'Timefreq') || isempty(sStudy.Timefreq), return; end
    for i = 1:numel(sStudy.Timefreq)
        sT = sStudy.Timefreq(i);
        if ~isempty(sT.DataFile) && file_compare(sT.DataFile, DataFile) && ~isempty(regexpi(sT.Comment, 'PSD', 'once'))
            TfFile = sT.FileName;  return;
        end
    end
end

% Compute an averaged magnitude PSD over the recording (MEG), return its timefreq file path.
function TfFile = i_compute_psd(DataFile)
    TfFile = '';
    sIn = bst_process('GetInputStruct', DataFile);
    if isempty(sIn), return; end
    sOut = bst_process('CallProcess', 'process_psd', sIn, [], ...
        'timewindow',  [], ...
        'win_length',  4, ...
        'win_overlap', 50, ...
        'units',       'physical', ...
        'sensortypes', 'MEG', ...
        'win_std',     0, ...
        'edit', struct('Comment','PSD: Dynamics focus', 'TimeBands',[], 'Freqs',[], ...
                       'ClusterFuncTime','none', 'Measure','magnitude', 'Output','all', 'SaveKernel',0));
    if ~isempty(sOut) && isfield(sOut,'FileName') && ~isempty(sOut(1).FileName), TfFile = sOut(1).FileName; end
end
```

- [ ] **Step 2: Lint (MATLAB MCP)**

Run: `check_matlab_code` on `toolbox/gui/panel_bst_dynamics.m`
Expected: no NEW errors.

- [ ] **Step 3: Verify in a live session (MATLAB MCP)**

Pre-req: Brainstorm running with `S01_AEF_01_notch` loaded and its recording time-series open (the Dynamics session targets it). Run (`evaluate_matlab_code`):

```matlab
st = getappdata(0, 'DynamicsTarget');
assert(~isempty(st) && ~isempty(st.T.DataFile), 'Dynamics session must be active');
hSpec = panel_bst_dynamics('i_ensure_psd', st);
assert(~isempty(hSpec) && ishandle(hSpec), 'ensure_psd returned a figure');
TfInfo = getappdata(hSpec, 'Timefreq');
assert(strcmpi(TfInfo.DisplayMode,'Spectrum'), 'spectrum mode');
hAxes = findobj(hSpec,'-depth',1,'Tag','AxesGraph');
assert(isequal(get(hAxes,'XLim'), [0 60]), 'xlim 0-60');
disp('TASK2_OK');
```

Expected: a spectrum figure opens (or is reused), prints `TASK2_OK`.

- [ ] **Step 4: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m
git commit -m "feat(dynamics): PSD lifecycle helper (find-or-open-or-compute spectrum)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Frequency DRIVE (panel → spectrum)

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (`i_drive`, freq case, lines 239-247; add overlay helpers)

**Interfaces:**
- Consumes: `i_ensure_psd`, `i_driving`, `figure_spectrum('SetFreqSelection', hSpec, [lo hi])`, `figure_spectrum('DrawSelection', hSpec)`.
- Produces: `i_freq_overlay(st,lo,hi)` → st; `i_freq_overlay_clear(st)` → st.

- [ ] **Step 1: Rewrite the `freq` case of `i_drive`**

In `toolbox/gui/panel_bst_dynamics.m`, replace the `case 'freq'` block (lines 239-247) with:

```matlab
        case 'freq'
            if isfinite(loc.center) && (loc.extent>0)
                lo = loc.center - loc.extent;  hi = loc.center + loc.extent;
                panel_filter('SetFilters', 1, hi, 1, lo, 0, [], 0, 1);
                st.curBand = [lo hi];  st.curBandName = i_freq_name(ctrl);
                st = i_freq_overlay(st, lo, hi);                 % ensure PSD + drive its freq selection
            else
                panel_filter('SetFilters', 0, [], 0, [], 0, [], 0, 0);
                st.curBand = [];  st.curBandName = '';
                st = i_freq_overlay_clear(st);                   % clear the spectrum band strip
            end
```

(The `setappdata(0,'DynamicsTarget',st)` at the end of `i_drive`, line 259, persists the updated `st.hSpec`.)

- [ ] **Step 2: Add the overlay helpers**

After the PSD lifecycle helpers (Task 2), insert:

```matlab
%% ===== FREQ OVERLAY: drive the spectrum band strip (panel -> view) =====
function st = i_freq_overlay(st, lo, hi)
    hSpec = i_ensure_psd(st);
    st.hSpec = hSpec;
    if isempty(hSpec) || ~ishandle(hSpec), return; end
    i_driving(true);
    try, figure_spectrum('SetFreqSelection', hSpec, [lo hi]); catch, end %#ok<CTCH>
    i_driving(false);
end
function st = i_freq_overlay_clear(st)
    if isfield(st,'hSpec') && ~isempty(st.hSpec) && ishandle(st.hSpec)
        i_driving(true);
        try
            setappdata(st.hSpec, 'GraphSelection', []);          % [] clears WITHOUT prompting (SetFreqSelection([]) would prompt)
            figure_spectrum('DrawSelection', st.hSpec);
        catch
        end
        i_driving(false);
    end
end
```

- [ ] **Step 3: Lint (MATLAB MCP)**

Run: `check_matlab_code` on `toolbox/gui/panel_bst_dynamics.m`
Expected: no NEW errors.

- [ ] **Step 4: Verify drive in a live session (MATLAB MCP)**

Pre-req: Dynamics session active on `S01_AEF_01_notch`. Run:

```matlab
ctrl = bst_get('PanelControls','Dynamics');
ctrl.jFreqBand.setSelectedItem('alpha');  drawnow;     % fires OnFreqPreset -> i_drive('freq')
panel_bst_dynamics('OnFreqPreset');  drawnow;          % deterministic re-trigger
st = getappdata(0,'DynamicsTarget');
assert(~isempty(st.hSpec) && ishandle(st.hSpec), 'spectrum opened');
gs = getappdata(st.hSpec, 'GraphSelection');
assert(~isempty(gs) && abs(min(gs)-8)<0.6 && abs(max(gs)-13)<0.6, sprintf('band strip 8-13, got [%g %g]', min(gs), max(gs)));
% clear path:
ctrl.jFreqBand.setSelectedItem('none');  panel_bst_dynamics('OnFreqPreset');  drawnow;
gs2 = getappdata(st.hSpec, 'GraphSelection');
assert(isempty(gs2), 'band strip cleared on none');
disp('TASK3_OK');
```

Expected: spectrum shows an 8–13 Hz strip, then clears; prints `TASK3_OK`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m
git commit -m "feat(dynamics): frequency drive - band preset paints the PSD band strip

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Frequency sync-back (spectrum → panel)

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (add `i_sync_freq`)
- Modify: `toolbox/gui/figure_spectrum.m` (`FigureMouseUpCallback`, after line 404)

**Interfaces:**
- Consumes: `NotifySelection` (Task 1), `i_band_match`, `i_driving`.
- Produces: `i_sync_freq(st, range)`.

- [ ] **Step 1: Add `i_sync_freq` to the panel**

In `toolbox/gui/panel_bst_dynamics.m`, after the freq overlay helpers, insert:

```matlab
%% ===== FREQ SYNC-BACK: a user edit on the spectrum updates the panel (view -> panel) =====
function i_sync_freq(st, range)
    ctrl = bst_get('PanelControls', 'Dynamics');  if isempty(ctrl), return; end
    lo = range(1);  hi = range(2);
    if (hi - lo) < 1e-6, return; end
    % reflect the dragged band in the panel fields
    ctrl.jFreqC.setText(num2str((lo+hi)/2));
    ctrl.jFreqW.setText(num2str((hi-lo)/2));
    % apply the matching time-series bandpass directly (combo cascade may not fire if value is unchanged)
    panel_filter('SetFilters', 1, hi, 1, lo, 0, [], 0, 1);
    nm = i_band_match(lo, hi);
    st.curBand = [lo hi];
    if strcmpi(nm,'custom'), st.curBandName = ''; else, st.curBandName = nm; end
    setappdata(0, 'DynamicsTarget', st);
    % reflect the band name in the combobox (display); any cascade is idempotent and overlay-guarded
    try, ctrl.jFreqBand.setSelectedItem(nm); catch, end %#ok<CTCH>
end
```

- [ ] **Step 2: Add the mouse-up hook in `figure_spectrum`**

In `toolbox/gui/figure_spectrum.m`, in `FigureMouseUpCallback`, insert immediately after line 404 (the `end` closing the `if ~hasMoved && ~isempty(MouseStatus)` block) and before line 406 (`% Reset MouseMove callbacks`):

```matlab
    % Bidirectional Dynamics focus: report a completed user freq-selection to the Dynamics panel.
    if hasMoved && ~isempty(getappdata(0, 'DynamicsTarget'))
        TfInfo = getappdata(hFig, 'Timefreq');
        GraphSelection = getappdata(hFig, 'GraphSelection');
        if ~isempty(TfInfo) && strcmpi(TfInfo.DisplayMode,'Spectrum') && ~isempty(GraphSelection) && all(isfinite(GraphSelection))
            try, panel_bst_dynamics('NotifySelection', hFig, 'freq', GraphSelection); catch, end %#ok<CTCH>
        end
    end
```

- [ ] **Step 3: Lint both files (MATLAB MCP)**

Run: `check_matlab_code` on `toolbox/gui/panel_bst_dynamics.m` and `toolbox/gui/figure_spectrum.m`
Expected: no NEW errors.

- [ ] **Step 4: Verify sync-back + no-loop in a live session (MATLAB MCP)**

Pre-req: Dynamics session active; alpha selected (Task 3 leaves `st.hSpec` open). Run (simulates a user freq-selection by setting the appdata + calling the hook target directly):

```matlab
st = getappdata(0,'DynamicsTarget');  hSpec = st.hSpec;
% simulate a user-dragged 18-24 Hz selection landing in mouse-up:
setappdata(hSpec, 'GraphSelection', [18 24]);
panel_bst_dynamics('NotifySelection', hSpec, 'freq', [18 24]);  drawnow;
ctrl = bst_get('PanelControls','Dynamics');
c = str2double(char(ctrl.jFreqC.getText()));  w = str2double(char(ctrl.jFreqW.getText()));
assert(abs(c-21)<1e-6 && abs(w-3)<1e-6, sprintf('fields center/half got %g/%g', c, w));
assert(strcmpi(char(ctrl.jFreqBand.getSelectedItem()),'custom'), 'combo custom');
% no-loop: with the guard set, a notify must be ignored
panel_bst_dynamics('i_driving', true);
setappdata(hSpec, 'GraphSelection', [1 2]);
panel_bst_dynamics('NotifySelection', hSpec, 'freq', [1 2]);
panel_bst_dynamics('i_driving', false);
c2 = str2double(char(ctrl.jFreqC.getText()));
assert(abs(c2-21)<1e-6, 'guarded notify ignored');
disp('TASK4_OK');
```

Expected: panel fields read center 21 / half-width 3, combo `custom`, guarded notify ignored; prints `TASK4_OK`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m toolbox/gui/figure_spectrum.m
git commit -m "feat(dynamics): frequency sync-back - dragging the PSD band updates the panel

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Time DRIVE (panel → time-series), saved windows + detect

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (add `i_focus_time`; `OnDetect` after line 410; `TreeSel_Callback` window/atom branches lines 960-971)

**Interfaces:**
- Consumes: `i_rec_figure`, `i_driving`, `figure_timeseries('SetTimeSelectionManual', hRec, [t1 t2])`.
- Produces: `i_focus_time(st, win)` where `win=[tStart tEnd]`.

- [ ] **Step 1: Add `i_focus_time`**

In `toolbox/gui/panel_bst_dynamics.m`, after `i_sync_freq`, insert:

```matlab
%% ===== TIME FOCUS: drive the recording's Time Selection box (panel -> view) =====
function i_focus_time(st, win)
    if isempty(win) || (numel(win) < 2) || any(~isfinite(win(1:2))), return; end
    hRec = i_rec_figure(st);
    if isempty(hRec) || ~ishandle(hRec), return; end
    i_driving(true);
    try, figure_timeseries('SetTimeSelectionManual', hRec, [win(1) win(2)]); catch, end %#ok<CTCH>
    i_driving(false);
end
```

- [ ] **Step 2: Drive the box to the first window on Detect**

In `OnDetect`, after `i_apply(st);` (line 410) and before the `bst_progress('text', ...)` (line 411), insert:

```matlab
    [~, st] = i_cs();                                            % i_apply rewrote the target
    i_focus_time(st, [evt(1,1), evt(2,1)]);                      % focus the FIRST detected window
```

- [ ] **Step 3: Drive the box when a saved window / atom is selected in the tree**

In `TreeSel_Callback`, the `'window'` branch (line 960-963) currently ends with `i_jump(...)`. Replace those lines:

```matlab
        elseif strcmp(info.kind, 'window')
            [rows, occMap] = i_window_atoms(st.T, info.g, info.w, i_field(st,'showPhase',[1 1 1 1]));
            for k = 1:numel(rows), model.addElement(rows{k}); end
            i_jump(st.T.Groups(info.g).times(1, info.w));   % selecting a window jumps to its onset
            i_focus_time(st, st.T.Groups(info.g).times(:, info.w)');   % and focuses the [onset offset] box
            st.detSel = [];                                  % saved window -> no staged-edit target
```

- [ ] **Step 4: Lint (MATLAB MCP)**

Run: `check_matlab_code` on `toolbox/gui/panel_bst_dynamics.m`
Expected: no NEW errors.

- [ ] **Step 5: Verify in a live session (MATLAB MCP)**

Pre-req: Dynamics session active on `S01_AEF_01_notch`, alpha band selected. Run:

```matlab
panel_bst_dynamics('OnDetect');  drawnow;
st = getappdata(0,'DynamicsTarget');
hRec = panel_bst_dynamics('i_rec_figure', st);
gs = getappdata(hRec, 'GraphSelection');
assert(~isempty(gs) && numel(gs)==2 && all(isfinite(gs)), 'time box set on detect');
assert(gs(2) > gs(1), 'box has positive extent');
disp(sprintf('TASK5_OK box=[%.3f %.3f]', gs(1), gs(2)));
```

Expected: after Detect the recording shows a highlighted window box; prints `TASK5_OK` with the box bounds.

- [ ] **Step 6: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m
git commit -m "feat(dynamics): time drive - Detect and saved-window selection set the Time Selection box

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: Staged-detection per-window navigation

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (`BuildTree` detection-mirror block lines 916-930; `TreeSel_Callback` detevt branch lines 972-976)

**Interfaces:**
- Consumes: `i_focus_time` (Task 5), `panel_record('GetEvents', [], 1)`.
- Produces: tree node kinds `'detwinroot'` / `'detwin'`; `st.detSel = [ie w]` set on staged-window selection.

- [ ] **Step 1: Emit per-window leaves for the staged band-window event**

In `BuildTree`, replace the detection-event loop (lines 923-928):

```matlab
            for ie = find(isDet)
                e = evs(ie);
                isWin = ~isempty(regexp(e.label, '\([0-9.]+-[0-9.]+ Hz\)$', 'once')) && (size(e.times,1) == 2);
                if isWin
                    winNode = DefaultMutableTreeNode(sprintf('%s  (%d)', e.label, size(e.times,2)));
                    detNode.add(winNode);
                    nodeList{end+1} = winNode;  nodeInfo(end+1) = struct('kind','detwinroot','g',ie,'w',0); %#ok<AGROW>
                    for w = 1:size(e.times,2)
                        leaf = DefaultMutableTreeNode(sprintf(' %.3f - %.3f s', e.times(1,w), e.times(2,w)));
                        winNode.add(leaf);
                        nodeList{end+1} = leaf;  nodeInfo(end+1) = struct('kind','detwin','g',ie,'w',w); %#ok<AGROW>
                    end
                else
                    leaf = DefaultMutableTreeNode(sprintf('%s  (%d)', e.label, size(e.times,2)));
                    detNode.add(leaf);
                    nodeList{end+1} = leaf;  nodeInfo(end+1) = struct('kind','detevt','g',ie,'w',0); %#ok<AGROW>
                end
            end
```

- [ ] **Step 2: Handle the new node kinds in `TreeSel_Callback`**

In `TreeSel_Callback`, replace the `'detevt'` branch (lines 972-976):

```matlab
        elseif strcmp(info.kind, 'detevt')
            evs = panel_record('GetEvents', [], 1);
            if (info.g <= numel(evs)) && ~isempty(evs(info.g).times)
                i_jump(evs(info.g).times(1,1));
            end
            st.detSel = [];
        elseif strcmp(info.kind, 'detwinroot')
            evs = panel_record('GetEvents', [], 1);
            if (info.g <= numel(evs)) && ~isempty(evs(info.g).times)
                i_focus_time(st, evs(info.g).times(:,1)');
                st.detSel = [info.g, 1];
            end
        elseif strcmp(info.kind, 'detwin')
            evs = panel_record('GetEvents', [], 1);
            if (info.g <= numel(evs)) && (info.w <= size(evs(info.g).times,2))
                win = evs(info.g).times(:, info.w)';
                i_jump(win(1));
                i_focus_time(st, win);
                st.detSel = [info.g, info.w];   % remember for time-drag sync-back (Task 7)
            end
```

Also set `st.detSel = []` in the other selection branches so a stale staged-window target never lingers. In the same callback, in the `'stack'` branch (line 956) and `'atom'` branch (line 964) add `st.detSel = [];` as the first statement of each branch.

- [ ] **Step 3: Lint (MATLAB MCP)**

Run: `check_matlab_code` on `toolbox/gui/panel_bst_dynamics.m`
Expected: no NEW errors.

- [ ] **Step 4: Verify staged navigation in a live session (MATLAB MCP)**

Pre-req: Detect has been run (Task 5), so a staged band-window event with ≥2 windows exists. Run:

```matlab
evs = panel_record('GetEvents', [], 1);
iWin = find(~cellfun(@isempty, regexp({evs.label}, '\([0-9.]+-[0-9.]+ Hz\)$', 'once')), 1);
assert(~isempty(iWin) && size(evs(iWin).times,2) >= 2, 'need >=2 staged windows');
st = getappdata(0,'DynamicsTarget');
% drive the SECOND staged window via the selection handler path:
st.detSel = [];  setappdata(0,'DynamicsTarget',st);
panel_bst_dynamics('i_focus_time', st, evs(iWin).times(:,2)');
hRec = panel_bst_dynamics('i_rec_figure', st);
gs = getappdata(hRec, 'GraphSelection');
assert(abs(gs(1)-evs(iWin).times(1,2))<0.05 && abs(gs(2)-evs(iWin).times(2,2))<0.05, 'box on 2nd window');
disp('TASK6_OK');
```

Expected: the Detection node now expands into per-window leaves; driving the 2nd window moves the box; prints `TASK6_OK`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m
git commit -m "feat(dynamics): staged-detection per-window navigation in the Atoms tree

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: Time sync-back (time-series → panel)

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (add `i_sync_time`)
- Modify: `toolbox/gui/figure_timeseries.m` (`FigureMouseUpCallback`, `elseif hasMoved` branch, after line 667)

**Interfaces:**
- Consumes: `NotifySelection` (Task 1), `st.detSel` (Task 6), `panel_record('SetEvents'/'ReplotEvents'/'GetCurrentDataset')`.
- Produces: `i_sync_time(st, range)`.

- [ ] **Step 1: Add `i_sync_time` to the panel**

In `toolbox/gui/panel_bst_dynamics.m`, after `i_focus_time`, insert:

```matlab
%% ===== TIME SYNC-BACK: a user edit of the Time Selection updates the panel (view -> panel) =====
% Records the active focus window. If a STAGED detection window is selected (st.detSel),
% rewrites that window's [onset offset] in the staged event (pre-save adjust). Never mutates
% already-saved atoms.
function i_sync_time(st, range)
    st.focusTime = range;
    setappdata(0, 'DynamicsTarget', st);
    if ~isfield(st,'detSel') || isempty(st.detSel), return; end
    ie = st.detSel(1);  w = st.detSel(2);
    evs = panel_record('GetEvents', [], 1);
    if (ie > numel(evs)) || (w > size(evs(ie).times,2)) || (size(evs(ie).times,1) ~= 2), return; end
    global GlobalData; %#ok<TLEV>
    iDS = panel_record('GetCurrentDataset');
    wasMod = ~isempty(iDS) && ~isempty(GlobalData.DataSet(iDS).Measures.sFile) && GlobalData.DataSet(iDS).Measures.isModified;
    sEvent = evs(ie);
    sEvent.times(:, w) = sort(range(:));
    i_driving(true);
    try
        panel_record('SetEvents', sEvent, ie);
        panel_record('ReplotEvents');
    catch
    end
    i_driving(false);
    if ~isempty(iDS) && ~wasMod, GlobalData.DataSet(iDS).Measures.isModified = 0; end   % staged edit must not dirty the recording
end
```

- [ ] **Step 2: Add the mouse-up hook in `figure_timeseries`**

In `toolbox/gui/figure_timeseries.m`, `FigureMouseUpCallback`, replace the `elseif hasMoved` branch (lines 663-668):

```matlab
    % If time selection was defined: check if its length is non-zero
    elseif hasMoved
        GraphSelection = getappdata(hFig, 'GraphSelection');
        if (length(GraphSelection) == 2) && (GraphSelection(1) == GraphSelection(2))
            SetTimeSelectionLinked(hFig, []);
        elseif (length(GraphSelection) == 2) && all(isfinite(GraphSelection)) && ~isempty(getappdata(0, 'DynamicsTarget'))
            % Bidirectional Dynamics focus: report a completed user time-selection to the panel.
            try, panel_bst_dynamics('NotifySelection', hFig, 'time', GraphSelection); catch, end %#ok<CTCH>
        end
    end
```

- [ ] **Step 3: Lint both files (MATLAB MCP)**

Run: `check_matlab_code` on `toolbox/gui/panel_bst_dynamics.m` and `toolbox/gui/figure_timeseries.m`
Expected: no NEW errors.

- [ ] **Step 4: Verify sync-back + staged mutation in a live session (MATLAB MCP)**

Pre-req: Detect run; select a staged window so `st.detSel` is set. Run:

```matlab
evs0 = panel_record('GetEvents', [], 1);
iWin = find(~cellfun(@isempty, regexp({evs0.label}, '\([0-9.]+-[0-9.]+ Hz\)$', 'once')), 1);
st = getappdata(0,'DynamicsTarget');  st.detSel = [iWin, 1];  setappdata(0,'DynamicsTarget',st);
hRec = panel_bst_dynamics('i_rec_figure', st);
newWin = evs0(iWin).times(:,1)' + 0.010;             % nudge the 1st window by +10 ms
% simulate the completed drag landing in mouse-up:
setappdata(hRec, 'GraphSelection', newWin);
panel_bst_dynamics('NotifySelection', hRec, 'time', newWin);  drawnow;
st2 = getappdata(0,'DynamicsTarget');
assert(isequal(round(st2.focusTime*1000), round(newWin*1000)), 'focusTime recorded');
evs1 = panel_record('GetEvents', [], 1);
assert(max(abs(evs1(iWin).times(:,1) - sort(newWin(:)))) < 1e-6, 'staged window 1 rewritten');
% no-loop guard:
panel_bst_dynamics('i_driving', true);
panel_bst_dynamics('NotifySelection', hRec, 'time', [0 0.001]);
panel_bst_dynamics('i_driving', false);
st3 = getappdata(0,'DynamicsTarget');
assert(isequal(st3.focusTime, st2.focusTime), 'guarded notify ignored');
disp('TASK7_OK');
```

Expected: `focusTime` updates, the staged window's onset/offset shift, guarded notify is ignored; prints `TASK7_OK`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m toolbox/gui/figure_timeseries.m
git commit -m "feat(dynamics): time sync-back - dragging the Time Selection adjusts the staged window

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 8: End-to-end manual validation

**Files:** none (validation + notes only)

- [ ] **Step 1: Walk the full focus loop on the alpha example block**

In the live Brainstorm session (`S01_AEF_01_notch`, recording open, Dynamics panel docked):

1. Pick **alpha** in the Frequency combo. Expect: time-series bandpass applies; a PSD spectrum opens (auto-computed if absent) with X-axis 0–60 Hz and an 8–13 Hz vertical band strip.
2. Drag the band strip on the spectrum to ~18–24 Hz and release. Expect: Frequency fields read center 21 / half-width 3, combo flips to `custom`, the time-series bandpass follows. No flicker/loop.
3. Re-pick **alpha** (re-assert focus). Expect: strip and filter return to 8–13 Hz.
4. Click **Detect**. Expect: the recording shows a Time Selection box on the FIRST alpha window; the Atoms tree shows "Detection (events) [unsaved]" expanding to per-window leaves.
5. Select successive window leaves. Expect: the single box moves to each window's `[onset offset]`; clicking elsewhere with no drag clears the box (defocus — by design).
6. Select a window leaf, then drag the Time Selection box edges and release. Expect: that staged window's bounds update (the leaf label refreshes on the next tree rebuild) without dirtying the recording.

- [ ] **Step 2: Confirm non-Dynamics inertness**

Close the Dynamics session (the `x` button). Open a normal recording + a normal PSD spectrum (no Dynamics). Drag a time selection and a frequency selection. Expect: no errors, identical behavior to before this work (the hooks early-out on the empty `DynamicsTarget` check).

- [ ] **Step 3: Record the validation result**

Append a short "Validation" note (date, pass/fail per step, any deviations) to `dev/2026-06-25-dynamics-focus-state-design.md`, then commit:

```bash
git add dev/2026-06-25-dynamics-focus-state-design.md
git commit -m "docs(dynamics): record focus-state end-to-end validation on the alpha block

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-Review Notes

- **Spec coverage:** §2 native-selection binding → Tasks 3/5 (drive) + 4/7 (sync). §3 freq drive (filter + ensure-PSD + strip) → Tasks 2/3; freq sync → Task 4; time drive (first window + list nav) → Tasks 5/6; time sync → Task 7. §4 components: `NotifySelection` → Task 1; sync hooks → Tasks 4/7; re-entrancy guard → Task 1; PSD lifecycle → Task 2. §5 persistence (no persistent patch; re-assert on selection) → native semantics used throughout; clearing = defocus verified in Task 8 step 5. §6 validation → Task 8.
- **Refinement vs spec:** hooks live in each figure's `FigureMouseUpCallback` (fire once on completion), not in `SetTimeSelectionLinked`/`DrawSelection` (which fire per motion tick). Same outcome, no lag.
- **Ambiguity resolved:** "currently-selected staged atom group window" = the staged band-window detection event's window selected in the tree (`st.detSel`); saved windows are drive-only (never mutated).
- **Type consistency:** `i_focus_time(st,win)` win=`[tStart tEnd]`; `NotifySelection(hFig,axis,range)` range=`[lo hi]`; `st.detSel=[ie w]`; `st.hSpec` handle; `st.focusTime=[lo hi]` — consistent across Tasks 1/5/6/7.
