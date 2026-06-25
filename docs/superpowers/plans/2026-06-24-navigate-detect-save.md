# Navigate / Detect / Save Contract Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the atom panel's persistence boundary explicit — Detect stages a time skeleton as Brainstorm events (rendered on the time series + Atoms tree, navigable), and only the Save actions write the atom table, recording time and frequency as numeric tensor indices.

**Architecture:** `OnDetect` stops writing `st.T`; it installs the refphase result as a Brainstorm event group on the recording (via the `panel_record` events API) — the Detect staging. A new `OnSaveDetection` converts those events into atom groups in `st.T`, stamping each with a numeric freq Localization via `bst_atom`. `OnSaveCursor` commits the live cursor (`st.nav`) as one atom; `OnLoadAtom` round-trips a saved atom back into the navigator blocks. `OnRecord`/`OnCaptureRegion` are unchanged writers.

**Tech Stack:** MATLAB, Brainstorm GUI; `panel_record` events API; `process_evt_refphase` (Compute); `bst_atom` / `bst_dynamics`; tested headless under `brainstorm nogui`.

## Global Constraints

- Phase 4 of `docs/superpowers/specs/2026-06-24-atom-tensor-architecture-analysis.md`; spec `docs/superpowers/specs/2026-06-24-navigate-detect-save-design.md`.
- **Three buffers:** `st.nav` (transient cursor; Navigate) · detection **events** on the recording (Detect staging, not atoms) · `st.T` (saved atom table; only persisted buffer). **Only** the writers — Record, Capture, **Save cursor**, **Save detection** — write `st.T`. Detect writes only events.
- **Detect → events:** `OnDetect` installs the band-window + 4 phase markers as a Brainstorm event group on the recording (auto-renders on the time series; mirrored in the Atoms tree; navigable). No `AddGroup(st.T,…)`, no auto-save in `OnDetect`.
- **Frequency is numeric on the atom:** Save stamps each atom's freq Localization `(center, extent)` via `bst_atom('Set', G, 'freq', …)`, sourced from **the band the detection was run at** (read from the detection event group's stored band — not the possibly-drifted live navigator band), plus `bandName` for the label.
- **Save cursor** = one atom from `st.nav` + the measured descriptor (distinct from Record's extrema). **Load atom** = explicit action (Atoms menu / double-click), not auto-on-select.
- Source detection stays **manual** this phase (navigate to a detection event → Region tool / Capture / Record). The Φ/Ψ-extrema+saddles source detector is a separate future phase.
- Event-group naming (matches today's labels): extended `"<band> (lo-hi Hz)"`, simple `"<band>_peak/_trough/_rising/_falling"`. Phase numeric value via `process_evt_refphase('PhaseValue', name)`.
- Do not start implementation on `development`; this work is on branch `feature/navigate-detect-save` (spec committed there).
- Tests run headless under **`brainstorm nogui` (GuiLevel 0)** — NOT `server` (-1). Do not `clear`/restart Brainstorm or close all figures; `rehash` and re-run. The unconstrained-kernel fixture backs the figure-gated tests.

---

### Task 1: Detect → detection events (stop writing st.T) + tree mirror

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (`OnDetect` rewrite; `i_detect_events` helper; `BuildTree` mirrors the detection event group)
- Test: `dev/test_detect_save.m` (new — T1: detect installs events, writes nothing to st.T)

**Interfaces:**
- Consumes: `process_evt_refphase('Compute'/'PhaseValue', …)`; `panel_record('SetEvents'/'GetEvents'/'UpdateEventsList')`; `i_load_meg` / `i_bands` (existing locals); `st.curBand`/`st.curBandName` (set by the navigator freq block).
- Produces:
  - `OnDetect` installs a detection event group (extended `<band> (lo-hi Hz)` + 4 phase simple events) on the loaded recording; writes nothing to `st.T`.
  - `i_detect_event_group(bandName)` → the event labels for a band: `{winLabel, peakL, troughL, risingL, fallingL}` (used by Task 2/Task 1 tree).
  - `BuildTree` shows a top node `Detection (events)` mirroring the detection event group (read from the recording), navigable (selecting a marker jumps the time cursor).

- [ ] **Step 1: Write the failing test**

Create `dev/test_detect_save.m`:

```matlab
function test_detect_save()
% TEST_DETECT_SAVE: Navigate/Detect/Save contract — Detect stages events, only Save writes st.T.
%
% USAGE:  test_detect_save   % Brainstorm running in nogui (GuiLevel 0)
%
% Authors: Diellor Basha, 2026

    PF = {'FAIL','PASS'};  pass = true;

    [linkFile, relData] = i_find_kernel_ds();
    if isempty(linkFile)
        fprintf('SKIPPED (no unconstrained kernel link)\n');
        fprintf('\n==== SUITE: %s ====\n', PF{pass+1});  return;
    end
    sStudy = bst_get('DataFile', relData);
    R = '';  for j=1:numel(sStudy.Result), if ~isempty(regexp(sStudy.Result(j).Comment,'MN: MEG\(Unconstr\)','once')) && ~isempty(regexp(sStudy.Result(j).FileName,'KERNEL','once')), R = sStudy.Result(j).FileName; break; end; end
    hFig = view_dynamics('FromResult', R);  drawnow;
    ctrl = bst_get('PanelControls', 'Dynamics');

    % select the alpha band in the navigator (sets st.curBand)
    ctrl.jFreqBand.setSelectedItem('alpha');  panel_bst_dynamics('OnFreqPreset');  drawnow;

    % ---------- T1: Detect installs events, writes NOTHING to st.T ----------
    st0 = getappdata(0,'DynamicsTarget');  nG0 = st0.T.nGroups;
    panel_bst_dynamics('OnDetect');  drawnow;
    evs = panel_record('GetEvents', [], 1);                       % all events on the recording
    haveWin  = ~isempty(evs) && any(~cellfun(@isempty, regexp({evs.label}, '^alpha \(', 'once')));
    havePeak = ~isempty(evs) && any(strcmpi({evs.label}, 'alpha_peak'));
    st1 = getappdata(0,'DynamicsTarget');
    ok1 = haveWin && havePeak && (st1.T.nGroups == nG0);          % events present, st.T unchanged
    fprintf('T1 detect->events: win=%d peak=%d stT_unchanged=%d => %s\n', haveWin, havePeak, (st1.T.nGroups==nG0), PF{ok1+1});
    pass = pass && ok1;

    if ishandle(hFig), close(hFig); end
    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end

function [linkFile, relData] = i_find_kernel_ds()
    linkFile = '';
    relData = 'Subject01/S01_AEF_20131218_01_notch/data_block001_02.mat';
    [sStudy, ~] = bst_get('DataFile', relData);
    if isempty(sStudy), return; end
    comments = {sStudy.Result.Comment};  fnames = {sStudy.Result.FileName};
    isMN = ~cellfun(@isempty, regexp(comments, 'MN: MEG\(Unconstr\)', 'once')) & ...
           ~cellfun(@isempty, regexp(fnames,   'KERNEL', 'once'));
    for j = find(isMN)
        try
            r = in_bst_results(fnames{j}, 0, 'nComponents','ImagingKernel');
            if (r.nComponents==3) && ~isempty(r.ImagingKernel), linkFile = ['link|' fnames{j} '|' relData];  return; end
        catch
        end
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run (MATLAB, `brainstorm nogui`):
```matlab
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3/dev'); rehash; test_detect_save
```
Expected: FAIL — `OnDetect` still writes `st.T` (so `stT_unchanged=0`) and installs no events (`win=0`).

- [ ] **Step 3: Rewrite `OnDetect` to install events instead of writing st.T**

In `toolbox/gui/panel_bst_dynamics.m`, replace the body of `OnDetect` **after** the `[evt, markers] = process_evt_refphase('Compute', …)` call and the `if isempty(evt) … return; end` guard — i.e. replace the block from `% rebuild the band's groups` through the final `bst_progress('text', …)` — with the event-install path:

```matlab
    % Stage the detection as a Brainstorm EVENT group on the recording (NOT atoms).
    % The events auto-render on the time series + mirror in the tree; Save converts them to atoms.
    i_detect_events(bandName, band, evt, markers);
    i_apply(st);                                                  % refresh tree (mirrors the detection events)
    bst_progress('text', sprintf('Detected %d %s windows (events; not yet saved)', size(evt,2), bandName));
```

Add the `i_detect_events` helper (place after `i_remove_band`). It builds the 5 event structs and installs them via the `panel_record` events API, replacing any prior same-label detection events:

```matlab
% Install the refphase detection as a Brainstorm event group on the loaded recording
% (extended band-window + 4 phase-marker simple events). Replaces prior same-label events.
function i_detect_events(bandName, band, evt, markers)
    winLabel = sprintf('%s (%g-%g Hz)', bandName, band(1), band(2));
    defs = { winLabel,                    evt,             [0.5 0.5 0.5]; ...
             [bandName '_peak'],          markers.peak,    [0.90 0.10 0.10]; ...
             [bandName '_trough'],        markers.trough,  [0.10 0.25 0.90]; ...
             [bandName '_rising'],        markers.rising,  [0.10 0.70 0.20]; ...
             [bandName '_falling'],       markers.falling, [0.95 0.55 0.10] };
    for k = 1:size(defs,1)
        label = defs{k,1};  times = defs{k,2};  color = defs{k,3};
        if isempty(times), continue; end
        sEvent = db_template('event');
        sEvent.label  = label;
        sEvent.color  = color;
        sEvent.times  = times;
        sEvent.epochs = ones(1, size(times,2));
        % find-or-replace by label (replace prior detection for this band on re-detect)
        all = panel_record('GetEvents', [], 1);
        iEvt = [];
        if ~isempty(all), iEvt = find(strcmpi({all.label}, label), 1); end
        if isempty(iEvt), iEvt = numel(all) + 1; end
        panel_record('SetEvents', sEvent, iEvt);
    end
    panel_record('UpdateEventsList');
end
```

⚠️ **Verify-while-implementing (events render on the time series):** `SetEvents`+`UpdateEventsList` updates the in-memory events + the Record-panel list. Confirm the **time-series figure draws the event markers**. If it does not refresh automatically, append a figure reload after `UpdateEventsList` in `i_detect_events`: `bst_figures('ReloadFigures', [], 0);` (the cheap per-figure reload used by `panel_record`). Use the minimal call that makes the markers appear; note what you used in the report.

In `BuildTree`, add a `Detection (events)` mirror node. After the existing top-level group loop builds `root`, append:

```matlab
    % mirror the detection event group (the unsaved time skeleton) as a distinct node
    evs = panel_record('GetEvents', [], 1);
    if ~isempty(evs)
        isDet = ~cellfun(@isempty, regexp({evs.label}, '_(peak|trough|rising|falling)$|\([0-9.]+-[0-9.]+ Hz\)$', 'once'));
        if any(isDet)
            detNode = DefaultMutableTreeNode('Detection (events)  [unsaved]');
            root.add(detNode);
            nodeList{end+1} = detNode;  nodeInfo(end+1) = struct('kind','detroot','g',0,'w',0); %#ok<AGROW>
            for ie = find(isDet)
                e = evs(ie);
                leaf = DefaultMutableTreeNode(sprintf('%s  (%d)', e.label, size(e.times,2)));
                detNode.add(leaf);
                nodeList{end+1} = leaf;  nodeInfo(end+1) = struct('kind','detevt','g',ie,'w',0); %#ok<AGROW>
            end
        end
    end
```

And in `TreeSel_Callback`, handle the `detevt` kind — selecting a detection-event leaf jumps the time cursor to its first time (read-only, no occurrence list):

```matlab
        elseif strcmp(info.kind, 'detevt')
            evs = panel_record('GetEvents', [], 1);
            if (info.g <= numel(evs)) && ~isempty(evs(info.g).times)
                i_jump(evs(info.g).times(1,1));
            end
```

(`detroot` selection is a no-op — it just groups the leaves.)

- [ ] **Step 4: Run test to verify it passes**

Run:
```matlab
rehash; test_detect_save
```
Expected: `T1 detect->events: win=1 peak=1 stT_unchanged=1 => PASS`, `==== SUITE: PASS ====`. Also eyeball that the `alpha` events render on the open time-series figure.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m dev/test_detect_save.m
git commit -m "feat(dynamics): Detect stages refphase as Brainstorm events (no st.T write) + tree mirror"
```

---

### Task 2: Save detection + Clear detection (events → atoms)

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (`OnSaveDetection`, `OnClearDetection`; a "Save detection"/"Clear preview" button)
- Modify: `dev/test_detect_save.m` (T2: save promotes; T3: clear removes); `dev/test_dynamics_atoms.m` (split T5)

**Interfaces:**
- Consumes: `panel_record('GetEvents'/'SetEvents'/'UpdateEventsList')`; `bst_dynamics('NewGroup'/'AddGroup')`; `bst_atom('Set'/'NewLoc')`; `i_remove_band` (existing); `process_evt_refphase('PhaseValue')`.
- Produces:
  - `OnSaveDetection` — converts the detection event group into atom groups in `st.T` (band-window extended + 4 phase children, times from events), each stamped with a numeric freq Localization; `i_apply`; auto-save.
  - `OnClearDetection` — removes the detection event group from the recording (discard, no commit).

- [ ] **Step 1: Write the failing test**

In `dev/test_detect_save.m`, before the `if ishandle(hFig)` cleanup, add T2 and T3:

```matlab
    % ---------- T2: Save detection promotes events -> atoms (numeric freq) ----------
    panel_bst_dynamics('OnSaveDetection');  drawnow;
    st = getappdata(0,'DynamicsTarget');  Td = st.T;  parents = {Td.Groups.parent};
    gW = find(cellfun(@isempty,parents) & arrayfun(@(k) strcmp(Td.Groups(k).bandName,'alpha') && size(Td.Groups(k).times,1)==2, 1:Td.nGroups), 1);
    chN = [];  if ~isempty(gW), chN = find(strcmpi(parents, Td.Groups(gW).label)); end
    lf = bst_atom('Get', Td.Groups(gW), 'freq');                  % numeric freq stamped
    ok2 = ~isempty(gW) && (numel(chN)==4) && (abs(lf.center-10.5)<1e-6) && (abs(lf.extent-2.5)<1e-6);
    fprintf('T2 save detection: window=%d children=%d freqC=%g freqW=%g => %s\n', ~isempty(gW), numel(chN), lf.center, lf.extent, PF{ok2+1});
    pass = pass && ok2;

    % ---------- T3: Clear removes the detection events ----------
    panel_bst_dynamics('OnDetect');  drawnow;                     % re-stage events
    panel_bst_dynamics('OnClearDetection');  drawnow;
    evs2 = panel_record('GetEvents', [], 1);
    stillDet = ~isempty(evs2) && any(strcmpi({evs2.label}, 'alpha_peak'));
    ok3 = ~stillDet;
    fprintf('T3 clear detection: detectionEventsGone=%d => %s\n', ~stillDet, PF{ok3+1});
    pass = pass && ok3;
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```matlab
rehash; test_detect_save
```
Expected: FAIL at T2 — `OnSaveDetection` is an unknown verb (the events were staged in T1 but nothing promotes them).

- [ ] **Step 3: Implement `OnSaveDetection` + `OnClearDetection`**

Add to `toolbox/gui/panel_bst_dynamics.m`:

```matlab
%% ===== SAVE DETECTION: convert the detection event group -> atoms in st.T =====
function OnSaveDetection() %#ok<DEFNU>
    [~, st] = i_cs();  if isempty(st), return; end
    evs = panel_record('GetEvents', [], 1);
    if isempty(evs), java_dialog('msgbox', 'No detection to save. Run Detect first.', 'Save detection');  return; end
    % find the band-window event ("<band> (lo-hi Hz)") + its 4 phase events
    iWin = find(~cellfun(@isempty, regexp({evs.label}, '^.+ \([0-9.]+-[0-9.]+ Hz\)$', 'once')), 1);
    if isempty(iWin), java_dialog('msgbox', 'No detection window event found.', 'Save detection');  return; end
    winLabel = evs(iWin).label;
    tok = regexp(winLabel, '^(.+) \(([0-9.]+)-([0-9.]+) Hz\)$', 'tokens', 'once');
    bandName = tok{1};  band = [str2double(tok{2}), str2double(tok{3})];
    DataFile = st.T.DataFile;
    if isempty(DataFile)
        for g=1:numel(st.T.Groups), if ~isempty(st.T.Groups(g).DataFile), DataFile = st.T.Groups(g).DataFile; break; end; end
    end
    % numeric freq Localization from the band the detection was RUN at (from the event label)
    fLoc = bst_atom('NewLoc', 'freq');  fLoc.center = mean(band);  fLoc.extent = (band(2)-band(1))/2;  fLoc.label = bandName;
    % rebuild this band's saved groups
    st.T = i_remove_band(st.T, bandName);
    W = bst_dynamics('NewGroup', winLabel);
    W.type='extended';  W.times=evs(iWin).times;  W.epochs=ones(1,size(evs(iWin).times,2));  W.bandName=bandName;
    W.color=[0.5 0.5 0.5];  W.SurfaceFile=st.T.SurfaceFile;  W.DataFile=DataFile;
    W = bst_atom('Set', W, 'freq', 1, fLoc);                      % numeric freq tensor index
    st.T = bst_dynamics('AddGroup', st.T, W);
    phaseColors = struct('peak',[0.90 0.10 0.10],'trough',[0.10 0.25 0.90],'rising',[0.10 0.70 0.20],'falling',[0.95 0.55 0.10]);
    for ph = {'peak','trough','rising','falling'}
        lbl = [bandName '_' ph{1}];  ie = find(strcmpi({evs.label}, lbl), 1);
        if isempty(ie) || isempty(evs(ie).times), continue; end
        P = bst_dynamics('NewGroup', lbl);
        P.type='simple';  P.parent=winLabel;  P.phase=process_evt_refphase('PhaseValue', ph{1});
        P.times=evs(ie).times(1,:);  P.epochs=ones(1,size(evs(ie).times,2));  P.bandName=bandName;
        P.color=phaseColors.(ph{1});  P.SurfaceFile=st.T.SurfaceFile;  P.DataFile=DataFile;
        P = bst_atom('Set', P, 'freq', 1, fLoc);                  % numeric freq on each child
        st.T = bst_dynamics('AddGroup', st.T, P);
    end
    if isempty(st.T.DataFile),    st.T.DataFile    = DataFile;      end
    if isempty(st.T.SurfaceFile), st.T.SurfaceFile = W.SurfaceFile; end
    i_apply(st);
    if ~isempty(st.file), try, bst_dynamics('Save', st.file, st.T); catch, end; end %#ok<CTCH>
    bst_progress('text', sprintf('Saved %s detection as atoms', bandName));
end

%% ===== CLEAR DETECTION: remove the detection event group (discard) =====
function OnClearDetection() %#ok<DEFNU>
    evs = panel_record('GetEvents', [], 1);
    if isempty(evs), return; end
    isDet = ~cellfun(@isempty, regexp({evs.label}, '_(peak|trough|rising|falling)$|\([0-9.]+-[0-9.]+ Hz\)$', 'once'));
    keep = evs(~isDet);
    if isempty(keep), keep = repmat(db_template('event'), 1, 0); end   % empty event array, not []
    panel_record('SetEvents', keep);                             % replace the whole event set with the non-detection ones
    panel_record('UpdateEventsList');
    [~, st] = i_cs();  if ~isempty(st), i_apply(st); end
end
```

⚠️ **Verify-while-implementing:** `panel_record('SetEvents', keep)` (single-arg) replaces the whole events structure — confirm this clears the detection events and refreshes the time series (add `bst_figures('ReloadFigures', [], 0)` if needed, matching whatever Task 1 used).

- [ ] **Step 4: Add the Save-detection / Clear buttons**

In `CreatePanel`, the Actions row (added in Phase 3) currently has Detect / Record / Capture. Add the Save-detection + Clear buttons after Detect:

```matlab
    gui_component('button', jAct, 'br hfill', 'Save detection', [], 'Promote the staged detection events into saved atoms (with numeric frequency)', @(h,e)bst_call(@OnSaveDetection));
    gui_component('button', jAct, 'br hfill', 'Clear preview',  [], 'Discard the staged detection events without saving', @(h,e)bst_call(@OnClearDetection));
```

- [ ] **Step 5: Split T5 in the atom suite**

In `dev/test_dynamics_atoms.m`, T5 currently calls `OnDetect` and asserts `st.T` gains the band-window + children. Under the new contract `OnDetect` stages events and `OnSaveDetection` writes `st.T`. Replace the T5 block's detect-then-assert with detect → assert events + st.T band empty, then save → assert the groups:

```matlab
    % T5: Detect stages events (no st.T write); Save detection promotes them (numeric freq)
    st0 = getappdata(0,'DynamicsTarget');  haveBand = ~isempty(st0.curBand);   % alpha, from T4
    panel_bst_dynamics('OnDetect');  drawnow;
    evsd = panel_record('GetEvents', [], 1);
    stagedOK = ~isempty(evsd) && any(strcmpi({evsd.label},'alpha_peak'));
    stMid = getappdata(0,'DynamicsTarget');
    bandEmptyMid = isempty(find(cellfun(@isempty,{stMid.T.Groups.parent}) & arrayfun(@(k) strcmp(stMid.T.Groups(k).bandName,'alpha') && size(stMid.T.Groups(k).times,1)==2, 1:stMid.T.nGroups), 1));
    panel_bst_dynamics('OnSaveDetection');  drawnow;
    st = getappdata(0,'DynamicsTarget');  Td = st.T;  parents = {Td.Groups.parent};
    gW = find(cellfun(@isempty,parents) & arrayfun(@(k) strcmp(Td.Groups(k).bandName,'alpha') && size(Td.Groups(k).times,1)==2, 1:Td.nGroups), 1);
    chN = [];  if ~isempty(gW), chN = find(strcmpi(parents, Td.Groups(gW).label)); end
    lf5 = []; if ~isempty(gW), lf5 = bst_atom('Get', Td.Groups(gW), 'freq'); end
    freqNum = ~isempty(lf5) && isfinite(lf5.center) && isfinite(lf5.extent);
    streamKept = any(arrayfun(@(k) strcmp(Td.Groups(k).Function,'stream'), 1:Td.nGroups));   % T4 record group survives
    ok5 = haveBand && stagedOK && bandEmptyMid && ~isempty(gW) && (numel(chN)==4) && freqNum && streamKept;
    fprintf('T5 detect/save: staged=%d midEmpty=%d window=%d children=%d freqNum=%d recordKept=%d => %s\n', stagedOK, bandEmptyMid, ~isempty(gW), numel(chN), freqNum, streamKept, PF{ok5+1});
    pass = pass && ok5;
```

- [ ] **Step 6: Run the tests**

Run:
```matlab
rehash; test_detect_save
test_dynamics_atoms
```
Expected: `test_detect_save` `==== SUITE: PASS ====` (T1–T3); `test_dynamics_atoms` `==== SUITE: PASS ====` (T1–T8, T5 split). T6–T8 (region/filter/capture) unaffected.

- [ ] **Step 7: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m dev/test_detect_save.m dev/test_dynamics_atoms.m
git commit -m "feat(dynamics): Save detection (events->atoms, numeric freq) + Clear; split atom-suite T5"
```

---

### Task 3: Save cursor as one atom

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (`OnSaveCursor`; a "Save cursor" button)
- Modify: `dev/test_detect_save.m` (T4)

**Interfaces:**
- Consumes: `st.nav` (the navigator cursor group); `bst_atom('Get', st.nav, axis)`; `bst_dynamics('NewGroup'/'AddGroup')`; `i_find_group`/`i_disp_band`/`i_op_color` (existing); `HelmholtzState` for the descriptor.
- Produces: `OnSaveCursor` — appends exactly one occurrence (the 4-D cursor + measured descriptor) into the `(bandName, Function)` group in `st.T`.

- [ ] **Step 1: Write the failing test**

In `dev/test_detect_save.m`, before cleanup, add T4:

```matlab
    % ---------- T4: Save cursor commits ONE atom from st.nav ----------
    ctrl.jFreqC.setText('10');  ctrl.jFreqW.setText('2');  panel_bst_dynamics('OnAxisChange','freq');  drawnow;
    panel_bst_dynamics('OnMeasurement','Solen');  drawnow;       % Function = stream
    st = getappdata(0,'DynamicsTarget');  nOcc0 = 0;
    gC0 = find(arrayfun(@(k) strcmp(st.T.Groups(k).Function,'stream') && strcmp(st.T.Groups(k).bandName,'alpha'), 1:st.T.nGroups), 1);
    if ~isempty(gC0), nOcc0 = numel(st.T.Groups(gC0).times); end
    panel_time('SetCurrentTime', 22.0);  panel_bst_dynamics('OnAxisChange','time');  drawnow;
    panel_bst_dynamics('OnSaveCursor');  drawnow;
    st = getappdata(0,'DynamicsTarget');
    gC1 = find(arrayfun(@(k) strcmp(st.T.Groups(k).Function,'stream') && strcmp(st.T.Groups(k).bandName,'alpha'), 1:st.T.nGroups), 1);
    nOcc1 = 0; if ~isempty(gC1), nOcc1 = numel(st.T.Groups(gC1).times); end
    ok4 = ~isempty(gC1) && (nOcc1 == nOcc0 + 1);
    fprintf('T4 save cursor: group=%d occ %d->%d => %s\n', ~isempty(gC1), nOcc0, nOcc1, PF{ok4+1});
    pass = pass && ok4;
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```matlab
rehash; test_detect_save
```
Expected: FAIL at T4 — `OnSaveCursor` is an unknown verb.

- [ ] **Step 3: Implement `OnSaveCursor`**

Add to `toolbox/gui/panel_bst_dynamics.m`:

```matlab
%% ===== SAVE CURSOR: commit the live 4-D cursor (st.nav) as ONE atom =====
function OnSaveCursor() %#ok<DEFNU>
    [~, st] = i_cs();  if isempty(st) || isempty(st.nav), return; end
    lt = bst_atom('Get', st.nav, 'time');     lf = bst_atom('Get', st.nav, 'freq');
    ls = bst_atom('Get', st.nav, 'source');   lk = bst_atom('Get', st.nav, 'scale'); %#ok<NASGU>
    if ~isfinite(lt.center)
        java_dialog('warning', 'Move the time cursor first (no cursor time).', 'Save cursor');  return;
    end
    band = st.curBand;  bandName = i_field(st, 'curBandName', '');
    op   = i_field(st, 'curOp', 'Total');
    switch op
        case 'Irrot', Func = 'potential';
        case 'Solen', Func = 'stream';
        otherwise,    Func = 'magnitude';
    end
    % measured descriptor at the cursor (operator scalar at the seed, if localized + a field is present)
    strength = NaN;  charge = NaN;
    if isfinite(ls.center) && ~isempty(st.hFig) && ishandle(st.hFig)
        St = getappdata(st.hFig, 'HelmholtzState');
        if ~isempty(St) && isfield(St,'Cache')
            [TimeVec, iT] = bst_memory('GetTimeVector', St.srcDS, St.srcResult, 'CurrentTimeIndex'); %#ok<ASGLU>
            if ~isempty(St.Cache) && isKey(St.Cache, iT)
                Ht = St.Cache(iT);
                switch op
                    case 'Irrot', sc = Ht.Phi;   case 'Solen', sc = Ht.Psi;   otherwise, sc = Ht.Fmag;
                end
                if ls.center>=1 && ls.center<=numel(sc), strength = sc(ls.center);  charge = sign(strength); end
            end
        end
    end
    % find-or-create the (band, Function) group; append ONE occurrence
    g = i_find_group(st.T, bandName, Func);
    if g < 1
        G = bst_dynamics('NewGroup', strtrim(sprintf('%s %s', i_disp_band(bandName, band), Func)));
        G.type='simple';  G.band=band;  G.bandName=bandName;  G.Function=Func;
        G.color=i_op_color(op);  G.SurfaceFile=st.T.SurfaceFile;  G.DataFile=st.T.DataFile;  G.ResultsFile=i_first_results(st.T);
        st.T = bst_dynamics('AddGroup', st.T, G);  g = numel(st.T.Groups);
    end
    G = st.T.Groups(g);
    G.times(1,end+1) = lt.center;   G.epochs(end+1) = 1;
    if isfinite(ls.center)
        v = double(ls.center);  G.vertices(end+1) = v;
        if ~isempty(ls.pos), G.pos(end+1,:) = ls.pos; else, G.pos(end+1,:) = [NaN NaN NaN]; end
        G.hemi(end+1) = 1 + (~isempty(ls.pos) && ls.pos(2)<0);
    else
        G.vertices(end+1) = NaN;  G.pos(end+1,:) = [NaN NaN NaN];  G.hemi(end+1) = NaN;
    end
    G.strength(end+1) = strength;   G.charge(end+1) = charge;   G.type='simple';
    if ~isempty(band), G = bst_atom('Set', G, 'freq', 1, lf); end   % numeric freq index
    st.T.Groups(g) = G;  st.T.nGroups = numel(st.T.Groups);
    i_apply(st);
    if ~isempty(st.file), try, bst_dynamics('Save', st.file, st.T); catch, end; end %#ok<CTCH>
    bst_progress('text', sprintf('Saved cursor atom (%s) at %.3f s', Func, lt.center));
end
```

Add the button in the Actions row (after Save detection):
```matlab
    gui_component('button', jAct, 'br hfill', 'Save cursor', [], 'Commit the current 4-D cursor as one atom', @(h,e)bst_call(@OnSaveCursor));
```

- [ ] **Step 4: Run test to verify it passes**

Run:
```matlab
rehash; test_detect_save
```
Expected: `T4 save cursor: group=1 occ N->N+1 => PASS`, suite PASS.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m dev/test_detect_save.m
git commit -m "feat(dynamics): Save cursor — commit the live 4-D cursor as one atom"
```

---

### Task 4: Load atom into the navigator

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (`OnLoadAtom`; an Atoms-menu "Load into navigator" item)
- Modify: `dev/test_detect_save.m` (T5)

**Interfaces:**
- Consumes: `st.occMap` (the selected occurrence → `(g, o)`); `bst_atom('Get', st.T.Groups(g), axis, o)`; the navigator block handles (`jTimeC/jTimeW/jFreqC/jFreqW/jSrcC/jSrcW/jScaleC/jScaleW`); `i_drive` / `OnAxisChange`.
- Produces: `OnLoadAtom` — fills the four blocks + `st.nav` from the selected saved atom and drives the viewers to it.

- [ ] **Step 1: Write the failing test**

In `dev/test_detect_save.m`, before cleanup, add T5:

```matlab
    % ---------- T5: Load atom into navigator round-trips the saved coords ----------
    st = getappdata(0,'DynamicsTarget');
    gL = find(arrayfun(@(k) strcmp(st.T.Groups(k).Function,'stream'), 1:st.T.nGroups), 1);   % the saved-cursor group
    % select occurrence 1 of that group: drive the tree to its stack then the list row
    st.occMap = [gL 1 0];  setappdata(0,'DynamicsTarget', st);   % emulate a selected occurrence
    panel_bst_dynamics('OnLoadAtom');  drawnow;
    st = getappdata(0,'DynamicsTarget');
    lt = bst_atom('Get', st.nav, 'time');  lf = bst_atom('Get', st.nav, 'freq');
    gt = bst_atom('Get', st.T.Groups(gL), 'time', 1);  gf = bst_atom('Get', st.T.Groups(gL), 'freq');
    ftxt = str2double(char(ctrl.jFreqC.getText()));
    ok5b = (abs(lt.center-gt.center)<1e-6) && (abs(lf.center-gf.center)<1e-6) && (abs(ftxt-gf.center)<1e-6);
    fprintf('T5b load atom: navTime=%g navFreq=%g freqField=%g => %s\n', lt.center, lf.center, ftxt, PF{ok5b+1});
    pass = pass && ok5b;
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```matlab
rehash; test_detect_save
```
Expected: FAIL at T5b — `OnLoadAtom` is an unknown verb.

- [ ] **Step 3: Implement `OnLoadAtom`**

Add to `toolbox/gui/panel_bst_dynamics.m`:

```matlab
%% ===== LOAD ATOM: fill the navigator blocks + st.nav from the selected saved atom =====
function OnLoadAtom() %#ok<DEFNU>
    [ctrl, st] = i_cs();
    if isempty(ctrl) || isempty(st) || isempty(st.occMap), java_dialog('warning','Select a saved atom first.','Load into navigator');  return; end
    g = st.occMap(1,1);  o = st.occMap(1,2);
    if (g < 1) || (g > numel(st.T.Groups)), return; end
    G = st.T.Groups(g);
    for axis = {'time','freq','source','scale'}
        ax = axis{1};
        loc = bst_atom('Get', G, ax, o);
        st.nav = bst_atom('Set', st.nav, ax, 1, loc);
        i_fill_block(ctrl, ax, loc);
        i_drive(ax, loc);                          % drive the viewers to the atom
    end
    setappdata(0, 'DynamicsTarget', st);
    bst_progress('text', 'Loaded atom into the navigator');
end

% Write a Localization's center/window into the axis block's fields (source window in mm).
function i_fill_block(ctrl, axis, loc)
    switch axis
        case 'time',   jC=ctrl.jTimeC;  jW=ctrl.jTimeW;   wv=loc.extent;
        case 'freq',   jC=ctrl.jFreqC;  jW=ctrl.jFreqW;   wv=loc.extent;
        case 'source', jC=ctrl.jSrcC;   jW=ctrl.jSrcW;    wv=loc.extent*1000;   % metres -> mm
        case 'scale',  jC=ctrl.jScaleC; jW=ctrl.jScaleW;  wv=loc.extent;
        otherwise, return;
    end
    if isfinite(loc.center), jC.setText(num2str(loc.center)); else, jC.setText(''); end
    if isfinite(wv),         jW.setText(num2str(wv));         else, jW.setText(''); end
    if strcmp(axis,'freq') && isfield(ctrl,'jFreqBand') && ~isempty(loc.label)
        try, ctrl.jFreqBand.setSelectedItem(loc.label); catch, end %#ok<CTCH>
    end
end
```

Add an Atoms-menu item (in `CreatePanel`, the Atoms menu, after "Capture region -> active atom"):
```matlab
    gui_component('MenuItem', jMenuAtoms, [], 'Load into navigator', IconLoader.ICON_EVT_TYPE, [], @(h,e)bst_call(@OnLoadAtom));
```

- [ ] **Step 4: Run the tests**

Run:
```matlab
rehash; test_detect_save
test_dynamics_atoms
test_nav_panel
```
Expected: `test_detect_save` `==== SUITE: PASS ====` (T1–T5b); `test_dynamics_atoms` 8/8; `test_nav_panel` T1–T6 — all `==== SUITE: PASS ====`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m dev/test_detect_save.m
git commit -m "feat(dynamics): Load atom into the navigator (round-trip saved atom -> blocks + st.nav)"
```

---

## Self-Review

**1. Spec coverage:**
- Three buffers + writer set → Tasks 1-4 (Detect=events, Save detection/cursor=writers; Record/Capture unchanged).
- Detect → Brainstorm events (no st.T write) + time-series render + tree mirror → Task 1.
- Save detection (events→atoms, numeric freq from the detected band) → Task 2.
- Clear detection → Task 2.
- Save cursor (one atom from st.nav + descriptor) → Task 3.
- Load atom into navigator (explicit, round-trip) → Task 4.
- T5 split in the atom suite → Task 2 Step 5.
- Out of scope (source-feature detector, full Scale, non-freq atlases) — not in any task. Correct.

**2. Placeholder scan:** none — every code step is complete; the two ⚠️ verify-while-implementing items (the time-series event-render refresh in Tasks 1-2) name the exact fallback (`bst_figures('ReloadFigures', [], 0)`) and are confirmations of an existing API, not placeholders.

**3. Type consistency:** the verbs `OnDetect`/`OnSaveDetection`/`OnClearDetection`/`OnSaveCursor`/`OnLoadAtom` match between definition, buttons/menu, and tests. The detection event labels (`<band> (lo-hi Hz)`, `<band>_peak/...`) are produced by `i_detect_events` and consumed by `OnSaveDetection`/`OnClearDetection`/`BuildTree` with the same regexes. `bst_atom('Set'/'Get', …, 'freq', …)` numeric Localization is written in Save detection + Save cursor and read in the tests. `st.nav`/`st.occMap`/`st.curBand`/`st.curOp` are the same fields the navigator (Phase 3) and OnRecord/OnDetect use. `i_fill_block` source window mm↔metres matches Phase 3's `i_read_block` (÷1000) convention.
