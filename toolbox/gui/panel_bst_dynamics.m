function varargout = panel_bst_dynamics( varargin )
% PANEL_BST_DYNAMICS: Record-style panel for the spatiotemporal atom system (bst_dynamics).
%
% The atom-table component (one bordered section of the future Dynamics panel). The LEFT is a
% tree whose top-level band atom (e.g. "alpha (8-13 Hz)") is a STACK (ICON_DATA_LIST) that
% expands to its time-window occurrences. Selecting a window lists, on the RIGHT, every
% single-time atom whose time falls in that window -- FLAT, sorted by time, with a phase column
% (peak / trough / rising / falling). Selecting an occurrence highlights its marker on the
% cortex and jumps the recording time. A File menu opens/saves the dynamics_* table and an
% Atoms menu adds/renames/deletes/colors/sorts band atoms. Docked as a tools tab; opened by
% view_dynamics. The temporal / spatial / frequency / eigenmode axes fold in as sibling
% sections in later increments.
%
% USAGE:  bstPanel = panel_bst_dynamics('CreatePanel')
%                    panel_bst_dynamics('SetTarget', hFig, T)
%
% Authors: Diellor Basha, 2026

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

eval(macro_method);
end


%% ===== CREATE PANEL =====
function bstPanelNew = CreatePanel() %#ok<DEFNU>
    panelName = 'Dynamics';
    import java.awt.*;
    import javax.swing.*;
    import org.brainstorm.list.*;
    import org.brainstorm.icon.*;
    fontSize = java_scaled('value', 11);

    jPanelNew = gui_component('Panel');

    % ===== ATOMS section (own bordered/titled component, like Events in Record) =====
    jPanelAtoms = gui_component('Panel');
    jPanelAtoms.setBorder(BorderFactory.createCompoundBorder( ...
        BorderFactory.createEmptyBorder(0,7,7,7), java_scaled('titledborder', 'Atoms')));

    % --- menu bar: File, Atoms (ICON_MENU + per-item icons, like Record) ---
    jMenuBar = gui_component('MenuBar', jPanelAtoms, BorderLayout.NORTH);
    jMenuBar.setPreferredSize(java_scaled('dimension', 20, 20));
    jMenuFile = gui_component('Menu', jMenuBar, [], 'File', IconLoader.ICON_MENU, [], [], 11);
    gui_component('MenuItem', jMenuFile, [], 'Open dynamics table...', IconLoader.ICON_FOLDER_OPEN, [], @(h,e)bst_call(@FileOpen));
    jMenuFile.addSeparator();
    gui_component('MenuItem', jMenuFile, [], 'Save',       IconLoader.ICON_SAVE, [], @(h,e)bst_call(@FileSave));
    gui_component('MenuItem', jMenuFile, [], 'Save as...', IconLoader.ICON_SAVE, [], @(h,e)bst_call(@FileSaveAs));
    jMenuAtoms = gui_component('Menu', jMenuBar, [], 'Atoms', IconLoader.ICON_MENU, [], [], 11);
    gui_component('MenuItem', jMenuAtoms, [], 'Add group',    IconLoader.ICON_EVT_TYPE_ADD,    [], @(h,e)bst_call(@AtomAddGroup));
    gui_component('MenuItem', jMenuAtoms, [], 'Rename group', IconLoader.ICON_EDIT,            [], @(h,e)bst_call(@AtomRenameGroup));
    gui_component('MenuItem', jMenuAtoms, [], 'Delete group', IconLoader.ICON_EVT_TYPE_DEL,    [], @(h,e)bst_call(@AtomDeleteGroup));
    gui_component('MenuItem', jMenuAtoms, [], 'Set color',    IconLoader.ICON_COLOR_SELECTION, [], @(h,e)bst_call(@AtomSetColor));
    jMenuAtoms.addSeparator();
    jMenuPhases = gui_component('Menu', jMenuAtoms, [], 'Show phases', IconLoader.ICON_EVT_TYPE, [], []);
    phaseNames  = {'peak','trough','rising','falling'};
    jPhaseItems = javaArray('javax.swing.JCheckBoxMenuItem', 4);
    for ip = 1:4
        jit = gui_component('checkboxmenuitem', jMenuPhases, [], phaseNames{ip}, [], [], @(h,e)bst_call(@()OnTogglePhase(ip)));
        jit.setSelected(true);
        jPhaseItems(ip) = jit;
    end
    jMenuAtoms.addSeparator();
    jMenuSort = gui_component('Menu', jMenuAtoms, [], 'Sort groups', IconLoader.ICON_EVT_TYPE, [], []);
    gui_component('MenuItem', jMenuSort, [], 'By name', IconLoader.ICON_EVT_TYPE, [], @(h,e)bst_call(@()AtomSort('name')));
    gui_component('MenuItem', jMenuSort, [], 'By time', IconLoader.ICON_EVT_TYPE, [], @(h,e)bst_call(@()AtomSort('time')));
    jMenuAtoms.addSeparator();
    gui_component('MenuItem', jMenuAtoms, [], 'Record at cursor', IconLoader.ICON_EVT_TYPE_ADD, [], @(h,e)bst_call(@OnRecord));
    gui_component('MenuItem', jMenuAtoms, [], 'Capture region -> active atom', IconLoader.ICON_SCOUT_NEW, [], @(h,e)bst_call(@OnCaptureRegion));

    % --- split: band stack TREE (left) | flat per-window atom list (right) ---
    jTree = java_create('javax.swing.JTree');
    jTree.setRootVisible(0);  jTree.setShowsRootHandles(1);
    jTree.getSelectionModel().setSelectionMode(javax.swing.tree.TreeSelectionModel.SINGLE_TREE_SELECTION);
    jTree.setFont(Font('Monospaced', Font.PLAIN, fontSize));
    rend = javax.swing.tree.DefaultTreeCellRenderer();   % band = data-list STACK icon; windows = leaf
    rend.setClosedIcon(IconLoader.ICON_DATA_LIST);  rend.setOpenIcon(IconLoader.ICON_DATA_LIST);  rend.setLeafIcon(IconLoader.ICON_DATA);
    jTree.setCellRenderer(rend);
    java_setcb(jTree, 'ValueChangedCallback', @(h,e)TreeSel_Callback());
    jScrollTree = JScrollPane(jTree);  jScrollTree.setBorder([]);

    jListOccur = JList();
    jListOccur.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    jListOccur.setCellRenderer(BstStringListRenderer(fontSize));
    java_setcb(jListOccur, 'ValueChangedCallback', @(h,e)OccurSel_Callback());
    jScrollOccur = JScrollPane(jListOccur);  jScrollOccur.setBorder([]);

    jSplit = JSplitPane(JSplitPane.HORIZONTAL_SPLIT, jScrollTree, jScrollOccur);
    jSplit.setResizeWeight(0.5);  jSplit.setDividerSize(java_scaled('value', 4));  jSplit.setBorder([]);
    jSplit.setPreferredSize(java_scaled('dimension', 360, 420));
    jPanelAtoms.add(jSplit, BorderLayout.CENTER);

    % ===== CONTROL area: the 4-axis (center,extent) navigator =====
    jCtrl = JPanel();  jCtrl.setLayout(BoxLayout(jCtrl, BoxLayout.Y_AXIS));
    BW = java_scaled('value', 30);  BH = java_scaled('value', 22);

    % TIME block (no preset yet)
    [jTimeC, jTimeW] = i_axis_block(jCtrl, 'time', 'Time', 'center', char(177), []);

    % FREQUENCY block + band-atlas preset combobox (right slot)
    bnames = i_bands();  bandItems = [bnames(:,1); {'custom'}];
    jFreqBand = gui_component('combobox', [], [], [], {bandItems}, [], [], []);
    jFreqBand.setSelectedItem('custom');
    java_setcb(jFreqBand, 'ActionPerformedCallback', @(h,e)bst_call(@OnFreqPreset));
    [jFreqC, jFreqW] = i_axis_block(jCtrl, 'freq', 'Frequency', 'center', char(177), jFreqBand);

    % SOURCE block + Region tool (right slot) -- the seed/radius picker
    jRegionTool = gui_component('toggle', [], '', 'Region', {Insets(0,0,0,0), Dimension(java_scaled('value',54),BH)}, 'Heat-disk tool: click a cortex vertex to seed (center), scroll to grow the radius (window)', @(h,e)bst_call(@()bst_geodesic_tool('Toggle', ctrl_region_state())));
    [jSrcC, jSrcW] = i_axis_block(jCtrl, 'source', 'Source', 'center', 'radius', jRegionTool);

    % SCALE block (basic: window -> heat smoothing; center reserved for Phase 5)
    [jScaleC, jScaleW] = i_axis_block(jCtrl, 'scale', 'Scale', 'center', char(177), []);

    % MEASUREMENT row (descriptor, not an axis) + actions
    jMeas = gui_river([2 2], [0 7 2 7], 'Measurement');
    jMeasPot = gui_component('toggle', jMeas, '', char(934), {Insets(0,0,0,0), Dimension(BW,BH)}, 'Potential \Phi (divergence: sources / sinks)', @(h,e)bst_call(@()OnMeasurement('Irrot')));
    jMeasStr = gui_component('toggle', jMeas, '', char(936), {Insets(0,0,0,0), Dimension(BW,BH)}, 'Stream \Psi (curl: vortices)', @(h,e)bst_call(@()OnMeasurement('Solen')));
    gui_component('label', jMeas, 'tab', '  Peaks:', [], [], [], []);
    jPeaks = gui_component('text', jMeas, '', '3', {Dimension(java_scaled('value',26), BH)}, 'Extrema kept per sign', []);
    jCtrl.add(jMeas);

    % ACTIONS row (kept: Detect / Record / Capture)
    jAct = gui_river([2 2], [0 7 2 7], 'Actions');
    gui_component('button', jAct, 'hfill', 'Detect windows', [], 'Run the band-power detector (refphase) on the selected band: writes the band-window stack + phase markers', @(h,e)bst_call(@OnDetect));
    gui_component('button', jAct, 'br hfill', 'Record at cursor', [], 'Store the shaped field''s extrema at the cursor as atoms', @(h,e)bst_call(@OnRecord));
    gui_component('button', jAct, 'br hfill', 'Capture region -> active atom', [], 'Snapshot the Region tool''s heat-disk into the selected atom', @(h,e)bst_call(@OnCaptureRegion));
    jCtrl.add(jAct);

    jPanelNew.add(jCtrl, BorderLayout.NORTH);

    jPanelNew.add(jPanelAtoms, BorderLayout.CENTER);
    bstPanelNew = BstPanel(panelName, jPanelNew, struct( ...
        'jTree',jTree, 'jListOccur',jListOccur, 'jMenuFile',jMenuFile, 'jMenuAtoms',jMenuAtoms, ...
        'jTimeC',jTimeC, 'jTimeW',jTimeW, 'jFreqC',jFreqC, 'jFreqW',jFreqW, 'jFreqBand',jFreqBand, ...
        'jSrcC',jSrcC, 'jSrcW',jSrcW, 'jRegionTool',jRegionTool, 'jScaleC',jScaleC, 'jScaleW',jScaleW, ...
        'jMeasPot',jMeasPot, 'jMeasStr',jMeasStr, 'jPeaks',jPeaks, 'jPhaseItems',jPhaseItems));
end


% Standard EEG/MEG bands (delta/theta/alpha/beta/gamma): {name, [lo hi], greek-label}
function b = i_bands()
    b = {'delta',[2 4],  char(948); ...
         'theta',[4 8],  char(952); ...
         'alpha',[8 13], char(945); ...
         'beta', [13 30],char(946); ...
         'gamma',[30 60],char(947)};
end


%% ===== UNIFORM AXIS BLOCK BUILDER =====
% Symmetric row: [center field] [window field] [right selector slot].
%   axis      'time'|'freq'|'source'|'scale'
%   cLabel    left label for center; wLabel for window
%   rightSel  the axis selector component already created (combobox/toggle) or [] (none)
% Returns the center/window text fields.
function [jC, jW] = i_axis_block(jCtrl, axis, title, cLabel, wLabel, rightSel)
    import java.awt.*;
    BH = java_scaled('value', 22);  FW = java_scaled('value', 52);
    jB = gui_river([2 2], [0 7 2 7], title);
    gui_component('label', jB, '', cLabel, [], [], [], []);
    jC = gui_component('text', jB, '', '', {Dimension(FW,BH)}, ['Center (' axis ')'], []);
    gui_component('label', jB, 'tab', wLabel, [], [], [], []);
    jW = gui_component('text', jB, '', '', {Dimension(FW,BH)}, ['Window/extent (' axis ')'], []);
    if ~isempty(rightSel), jB.add('tab', rightSel); end
    java_setcb(jC, 'ActionPerformedCallback', @(h,e)bst_call(@()OnAxisChange(axis)));
    java_setcb(jW, 'ActionPerformedCallback', @(h,e)bst_call(@()OnAxisChange(axis)));
    jCtrl.add(jB);
end


%% ===== READ a block's (center, extent) into a Localization =====
function loc = i_read_block(ctrl, axis)
    loc = bst_atom('NewLoc', axis);
    switch axis
        case 'time',   jC = ctrl.jTimeC;  jW = ctrl.jTimeW;
        case 'freq',   jC = ctrl.jFreqC;  jW = ctrl.jFreqW;
        case 'source', jC = ctrl.jSrcC;   jW = ctrl.jSrcW;
        case 'scale',  jC = ctrl.jScaleC; jW = ctrl.jScaleW;
        otherwise, return;
    end
    c = str2double(char(jC.getText()));  w = str2double(char(jW.getText()));
    if ~isnan(c), loc.center = c; end
    if ~isnan(w), loc.extent = abs(w); else, loc.extent = 0; end
    if strcmp(axis,'source') && isfinite(loc.extent), loc.extent = loc.extent / 1000; end   % jSrcW is mm -> metres
    if isfinite(loc.center), if loc.extent>0, loc.state='window'; else, loc.state='point'; end; end
end


%% ===== AXIS CHANGE: write the cursor atom + drive the engine =====
function OnAxisChange(axis) %#ok<DEFNU>
    [ctrl, st] = i_cs();
    if isempty(ctrl) || isempty(st), return; end
    loc = i_read_block(ctrl, axis);
    st.nav = bst_atom('Set', st.nav, axis, 1, loc);
    setappdata(0, 'DynamicsTarget', st);
    i_drive(axis, loc);
end


%% ===== FREQ PRESET: the band combobox fills the freq fields, THEN drives =====
% Only the combobox calls this (not the field edits), so typing a custom value is never
% overwritten by the selected band.
function OnFreqPreset() %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'Dynamics');  if isempty(ctrl), return; end
    i_freq_preset(ctrl);
    OnAxisChange('freq');
end


%% ===== per-axis engine driver (reuses existing engines; keeps legacy coords) =====
function i_drive(axis, loc)
    [ctrl, st] = i_cs();  if isempty(st), return; end
    switch axis
        case 'time'
            if isfinite(loc.center), try, panel_time('SetCurrentTime', loc.center); catch, end; end %#ok<CTCH>
        case 'freq'
            if isfinite(loc.center) && (loc.extent>0)
                lo = loc.center - loc.extent;  hi = loc.center + loc.extent;
                panel_filter('SetFilters', 1, hi, 1, lo, 0, [], 0, 1);
                st.curBand = [lo hi];  st.curBandName = i_freq_name(ctrl);
            else
                panel_filter('SetFilters', 0, [], 0, [], 0, [], 0, 0);
                st.curBand = [];  st.curBandName = '';
            end
        case 'source'
            % center/window are populated by the Region tool (Task 2 syncs them); nothing to drive here
        case 'scale'
            if ~isempty(st.hFig) && ishandle(st.hFig) && ~isempty(st.Lambda) && (loc.extent>0)
                params = struct('t', loc.extent);
                try, view_helmholtz('SetSmoothing', st.hFig, 1, 'heat', params); catch, end %#ok<CTCH>
                st.curScale = struct('on',1,'name','heat','params',params);
            elseif ~isempty(st.hFig) && ishandle(st.hFig)
                try, view_helmholtz('SetSmoothing', st.hFig, 0, 'heat', struct('t',1)); catch, end %#ok<CTCH>
                st.curScale = struct('on',0,'name','heat','params',[]);
            end
    end
    setappdata(0, 'DynamicsTarget', st);
end


%% ===== MEASUREMENT (operator descriptor; not an axis) =====
function OnMeasurement(which) %#ok<DEFNU>
    [ctrl, st] = i_cs();
    if isempty(ctrl) || isempty(st) || ~ishandle(st.hFig), return; end
    if strcmp(which, 'Irrot')
        if ctrl.jMeasPot.isSelected(), ctrl.jMeasStr.setSelected(false); name = 'Irrot'; else, name = 'Total'; end
    else
        if ctrl.jMeasStr.isSelected(), ctrl.jMeasPot.setSelected(false); name = 'Solen'; else, name = 'Total'; end
    end
    view_helmholtz('SetComponent', st.hFig, name);
    st.curOp = name;  setappdata(0, 'DynamicsTarget', st);
end


%% ===== frequency-atlas preset: a chosen band fills the freq center/window fields =====
function i_freq_preset(ctrl)
    if ~isfield(ctrl,'jFreqBand') || isempty(ctrl.jFreqBand), return; end
    sel = char(ctrl.jFreqBand.getSelectedItem());
    b = i_bands();  k = find(strcmpi(b(:,1), sel), 1);
    if isempty(k), return; end                                 % 'custom' -> leave fields as typed
    lo = b{k,2}(1);  hi = b{k,2}(2);
    ctrl.jFreqC.setText(num2str((lo+hi)/2));  ctrl.jFreqW.setText(num2str((hi-lo)/2));
end
function nm = i_freq_name(ctrl)
    nm = '';
    if isfield(ctrl,'jFreqBand') && ~isempty(ctrl.jFreqBand), nm = char(ctrl.jFreqBand.getSelectedItem()); end
    if strcmpi(nm,'custom'), nm = ''; end
end


%% ===== DETECT: refphase time-window skeleton for the selected band =====
% The FIRST atom creation: runs process_evt_refphase on the recording for the current
% band and writes the temporal skeleton -- the band-window extended stack + the 4
% phase-marker trains (time + phase, NO source). Space is added later by Record.
function OnDetect() %#ok<DEFNU>
    [~, st] = i_cs();
    if isempty(st), return; end
    band = i_field(st, 'curBand', []);  bandName = i_field(st, 'curBandName', '');
    if isempty(band)
        java_dialog('warning', ['Select a frequency band first (' char(948) '/' char(952) '/' char(945) '/' char(946) '/' char(947) ').'], 'Detect windows');
        return;
    end
    DataFile = st.T.DataFile;
    if isempty(DataFile)
        for g = 1:numel(st.T.Groups), if ~isempty(st.T.Groups(g).DataFile), DataFile = st.T.Groups(g).DataFile; break; end; end
    end
    if isempty(DataFile)
        java_dialog('warning', 'No recording (DataFile) is associated with this table.', 'Detect windows');  return;
    end
    bst_progress('start', 'Detect windows', sprintf('Running refphase on the %s band...', bandName));
    ChannelMat = in_bst_channel(bst_get('ChannelFileForStudy', DataFile));
    iMEG = channel_find(ChannelMat.Channel, 'MEG');
    [F, TimeVector] = i_load_meg(DataFile, ChannelMat, iMEG);
    if isempty(F)
        bst_progress('stop');  java_dialog('error', 'Could not read the recording.', 'Detect windows');  return;
    end
    OPTIONS = process_evt_refphase('Compute');  OPTIONS.freqRange = band;
    [evt, markers] = process_evt_refphase('Compute', F, TimeVector, OPTIONS);
    bst_progress('stop');
    if isempty(evt)
        java_dialog('msgbox', sprintf('No %s windows detected (try a different band or threshold).', bandName), 'Detect windows');  return;
    end
    % Stage the detection as a Brainstorm EVENT group on the recording (NOT atoms).
    % The events auto-render on the time series + mirror in the tree; Save converts them to atoms.
    i_detect_events(bandName, band, evt, markers);
    i_apply(st);                                                  % refresh tree (mirrors the detection events)
    bst_progress('text', sprintf('Detected %d %s windows (events; not yet saved)', size(evt,2), bandName));
end

% Load MEG sensor data (matrix [nMEG x nT] + time) for a data block or raw file.
function [F, TimeVector] = i_load_meg(DataFile, ChannelMat, iMEG)
    F = [];  TimeVector = [];
    DataMat = in_bst_data(DataFile, 'F', 'Time');
    if isstruct(DataMat.F)                                  % raw link
        IO = db_template('ImportOptions');
        IO.ImportMode='Time';  IO.UseCtfComp=1;  IO.UseSsp=1;  IO.EventsMode='ignore';  IO.DisplayMessages=0;  IO.RemoveBaseline='no';
        [F, TimeVector] = in_fread(DataMat.F, ChannelMat, 1, [], iMEG, IO);
    else                                                    % imported data block
        F = DataMat.F(iMEG, :);  TimeVector = DataMat.Time;
    end
end

% Drop the band's time-window group(s) (extended top-level, this bandName) + their phase
% children, leaving recorded spatial groups (simple, with a Function) untouched.
function T = i_remove_band(T, bandName)
    if isempty(T.Groups), return; end
    keep = true(1, numel(T.Groups));  winLabels = {};
    for g = 1:numel(T.Groups)
        G = T.Groups(g);  bn = G.bandName;  if isempty(bn), bn = ''; end
        if isempty(G.parent) && strcmp(bn, bandName) && (size(G.times,1) == 2)
            winLabels{end+1} = G.label;  keep(g) = false; %#ok<AGROW>
        end
    end
    for g = 1:numel(T.Groups)
        if any(strcmp(T.Groups(g).parent, winLabels)), keep(g) = false; end
    end
    T.Groups = T.Groups(keep);  T.nGroups = numel(T.Groups);
end


% Install the refphase detection as a Brainstorm event group on the loaded recording
% (extended band-window + 4 phase-marker simple events). Replaces prior same-label events.
function i_detect_events(bandName, band, evt, markers)
    global GlobalData;
    % Detection events are render-only preview, not a permanent edit to the recording file.
    % SetEvents sets Measures.isModified=1, and recording unload auto-saves under nogui --
    % so capture the dataset's clean/dirty state and restore CLEAN afterwards (preserving any
    % genuine prior user event edits).
    iDS = panel_record('GetCurrentDataset');                 % the loaded recording's dataset
    wasMod = ~isempty(iDS) && ~isempty(GlobalData.DataSet(iDS).Measures.sFile) && GlobalData.DataSet(iDS).Measures.isModified;
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
    panel_record('UpdateEventsList');                            % refresh the Record-panel event list
    panel_record('ReplotEvents');                               % re-plot the event markers on the open time series
    if ~isempty(iDS) && ~wasMod, GlobalData.DataSet(iDS).Measures.isModified = 0; end   % staged events are render-only, don't persist
end


%% ===== shared panel/target accessor =====
function [ctrl, st] = i_cs()
    ctrl = bst_get('PanelControls', 'Dynamics');
    st   = getappdata(0, 'DynamicsTarget');
end


%% ===== SYNC the Source block fields from the geodesic tool state =====
function SyncSource() %#ok<DEFNU>
    [ctrl, st] = i_cs();
    if isempty(ctrl) || isempty(st) || ~isfield(ctrl,'jSrcC'), return; end
    gs = bst_geodesic_tool('GetState');
    if isempty(gs), return; end
    ctrl.jSrcC.setText(num2str(double(gs.seed)));
    ctrl.jSrcW.setText(num2str(round(gs.radius*1000)));      % radius in mm
    loc = bst_atom('NewLoc', 'source');
    loc.center = double(gs.seed);  loc.extent = gs.radius;  loc.pos = gs.pos;  loc.state = 'window';
    st.nav = bst_atom('Set', st.nav, 'source', 1, loc);
    setappdata(0, 'DynamicsTarget', st);
end


%% ===== RECORD: shaped field's extrema at the cursor -> atoms =====
% Reads the linked Helmholtz decomposition at the current time, detects extrema of the
% scalar selected by the operator, and appends them to the (band, Function) group --
% tagged with the panel's (time, band, scale, operator) coordinates. Auto-saves.
function OnRecord() %#ok<DEFNU>
    [ctrl, st] = i_cs();
    if isempty(ctrl) || isempty(st) || ~ishandle(st.hFig), return; end
    St = getappdata(st.hFig, 'HelmholtzState');
    if isempty(St)
        java_dialog('warning', 'Record needs the linked Helmholtz source view (open via a Dirac result).', 'Record atoms');
        return;
    end
    view_helmholtz('UpdateFrame', st.hFig);                        % make sure the cursor frame is current
    St = getappdata(st.hFig, 'HelmholtzState');
    [TimeVec, iT] = bst_memory('GetTimeVector', St.srcDS, St.srcResult, 'CurrentTimeIndex');
    if isempty(St.Cache) || ~isKey(St.Cache, iT)
        java_dialog('warning', 'No decomposition at the current time.', 'Record atoms');  return;
    end
    Ht = St.Cache(iT);  tCur = TimeVec(iT);
    % scalar field + Function from the current operator
    op = i_field(st, 'curOp', 'Total');
    switch op
        case 'Irrot', Scal = Ht.Phi;  Func = 'potential';  signed = true;
        case 'Solen', Scal = Ht.Psi;  Func = 'stream';     signed = true;
        otherwise,    Scal = Ht.Fmag; Func = 'magnitude';  signed = false;
    end
    Surf = getappdata(st.hFig, 'DynamicsSurf');
    if isempty(Surf), Surf = in_tess_bst(st.T.SurfaceFile, 0);  setappdata(st.hFig, 'DynamicsSurf', Surf); end
    ex = bst_dynamics('Extrema', Scal, Surf.VertConn, i_peaks(ctrl), signed);
    if isempty(ex.iVertex)
        java_dialog('msgbox', 'No extrema in the current field.', 'Record atoms');  return;
    end
    v    = double(ex.iVertex(:)');
    pos  = Surf.Vertices(v, :);
    hemi = 1 + (Surf.Vertices(v,2) < 0);                          % SCS Y>0 = left
    band = i_field(st, 'curBand', []);  bandName = i_field(st, 'curBandName', '');
    % find-or-create the (band, Function) group, then append the cursor-time occurrences
    g = i_find_group(st.T, bandName, Func);
    if g < 1
        G = bst_dynamics('NewGroup', strtrim(sprintf('%s %s', i_disp_band(bandName, band), Func)));
        G.type='simple';  G.band=band;  G.bandName=bandName;  G.Function=Func;  G.scaleName=i_scale_name(st);
        G.color = i_op_color(op);  G.SurfaceFile = st.T.SurfaceFile;  G.DataFile = st.T.DataFile;  G.ResultsFile = i_first_results(st.T);
        st.T = bst_dynamics('AddGroup', st.T, G);  g = numel(st.T.Groups);
    end
    G = st.T.Groups(g);  nNew = numel(v);
    G.times    = [G.times,    repmat(tCur,1,nNew)];   G.epochs   = [G.epochs,   ones(1,nNew)];
    G.vertices = [G.vertices, v];                     G.pos      = [G.pos;      pos];
    G.hemi     = [G.hemi,     hemi];                  G.strength = [G.strength, ex.value(:)'];
    G.charge   = [G.charge,   ex.charge(:)'];         G.type     = 'simple';
    if ~isempty(G.region), G.region(end+1:numel(G.vertices)) = {[]}; end
    st.T.Groups(g) = G;  st.T.nGroups = numel(st.T.Groups);
    i_apply(st);                                                  % redraw markers + rebuild tree (stores st)
    if ~isempty(st.file), try, bst_dynamics('Save', st.file, st.T); catch, end; end %#ok<CTCH>  % auto-save
    bst_progress('text', sprintf('Recorded %d %s atom(s) at %.3f s', nNew, Func, tCur));
end

%% ===== CAPTURE: the selected Scout's geodesic region -> the active (selected) atom =====
% The active atom is the occurrence selected in the right-hand list. Snapshots the currently
% selected Scout's vertices into that occurrence (region + seed), localizing a time-only marker.
function OnCaptureRegion() %#ok<DEFNU>
    [ctrl, st] = i_cs();
    if isempty(ctrl) || isempty(st) || ~ishandle(st.hFig), return; end
    if isempty(st.occMap)
        java_dialog('warning', 'Select an atom in the list first.', 'Capture region');  return;
    end
    row = ctrl.jListOccur.getSelectedIndex() + 1;
    if (row < 1) || (row > size(st.occMap,1))
        java_dialog('warning', 'Select an atom in the list first.', 'Capture region');  return;
    end
    g = st.occMap(row,1);  o = st.occMap(row,2);
    if (size(st.T.Groups(g).times,1) ~= 1)
        java_dialog('warning', 'Select a single atom, not a time window.', 'Capture region');  return;
    end
    SurfaceFile = st.T.SurfaceFile;
    if isempty(SurfaceFile) && ~isempty(st.T.Groups), SurfaceFile = st.T.Groups(g).SurfaceFile; end
    % the geodesic region = the dynamics Region tool's current heat-disk (no scout)
    gs = bst_geodesic_tool('GetState');
    if isempty(gs) || isempty(gs.vertices)
        java_dialog('warning', 'Seed a region with the Region tool first.', 'Capture region');  return;
    end
    if ~isempty(gs.SurfaceFile) && ~isempty(SurfaceFile) && ~file_compare(gs.SurfaceFile, SurfaceFile)
        java_dialog('warning', 'The region is on a different surface than the atoms.', 'Capture region');  return;
    end
    seed = double(gs.seed);
    pos  = gs.pos;
    hemi = 1 + (pos(2) < 0);                                       % SCS Y>0 = left
    st.T.Groups(g) = bst_dynamics('AttachRegion', st.T.Groups(g), o, gs.vertices, seed, pos, hemi);
    i_apply(st);                                                   % redraw markers/regions + rebuild tree
    if ~isempty(st.file), try, bst_dynamics('Save', st.file, st.T); catch, end; end %#ok<CTCH>
    bst_progress('text', sprintf('Captured %d-vertex region into "%s"', numel(gs.vertices), st.T.Groups(g).label));
end

% Region-tool toggle state (1 when pressed) -> bst_geodesic_tool('Toggle', state)
function s = ctrl_region_state()
    s = 0;
    ctrl = bst_get('PanelControls', 'Dynamics');
    if ~isempty(ctrl) && isfield(ctrl,'jRegionTool') && ~isempty(ctrl.jRegionTool)
        s = double(ctrl.jRegionTool.isSelected());
    end
end

%% ===== SHOW-PHASES FILTER (display-only; never deletes atoms) =====
function OnTogglePhase(ip) %#ok<DEFNU>
    [ctrl, st] = i_cs();
    if isempty(ctrl) || isempty(st) || ~isfield(ctrl,'jPhaseItems') || isempty(ctrl.jPhaseItems), return; end
    sp = i_field(st, 'showPhase', [1 1 1 1]);
    sp(ip) = ctrl.jPhaseItems(ip).isSelected();
    st.showPhase = sp;  setappdata(0, 'DynamicsTarget', st);
    i_apply(st);                                                  % rebuild list + redraw cortex
end

% Phase name -> filter index (peak=1 trough=2 rising=3 falling=4; 0 = not a phase group -> always shown).
function k = i_phase_index(ph)
    if isempty(ph), k = 0; return; end
    switch lower(char(ph))
        case 'peak',    k = 1;
        case 'trough',  k = 2;
        case 'rising',  k = 3;
        case 'falling', k = 4;
        otherwise,      k = 0;
    end
end

% Readable phase type ('peak'/'trough'/'rising'/'falling') from a group's label
% suffix; '' for non-phase groups. phase itself is now a NUMERIC value (radians),
% so the human-readable type is parsed from the label ('<band>_peak' etc.).
function t = i_phase_type(G)
    t = '';
    if isempty(G.label), return; end
    tok = regexp(G.label, '_(peak|trough|rising|falling)$', 'tokens', 'once');
    if ~isempty(tok), t = tok{1}; end
end

function n = i_peaks(ctrl)
    n = 3;
    if isfield(ctrl,'jPeaks') && ~isempty(ctrl.jPeaks)
        x = str2double(char(ctrl.jPeaks.getText()));
        if ~isnan(x) && (x >= 1), n = round(x); end
    end
end
function val = i_field(st, name, default)
    if isfield(st, name) && ~isempty(st.(name)), val = st.(name); else, val = default; end
end
function g = i_find_group(T, bandName, Func)
    g = 0;
    for k = 1:numel(T.Groups)
        G = T.Groups(k);
        bn = G.bandName;  if isempty(bn), bn = ''; end
        fn = G.Function;  if isempty(fn), fn = ''; end
        qb = bandName;    if isempty(qb), qb = ''; end
        if isempty(G.parent) && strcmp(bn, qb) && strcmp(fn, Func), g = k;  return; end
    end
end
function s = i_disp_band(bandName, band)
    if ~isempty(bandName),  s = bandName;
    elseif ~isempty(band),  s = sprintf('%g-%g Hz', band(1), band(2));
    else,                   s = 'broadband';  end
end
function s = i_scale_name(st)
    sc = i_field(st, 'curScale', []);
    if isempty(sc) || ~isstruct(sc) || ~isfield(sc,'on') || ~sc.on, s = 'none';  else, s = sc.name;  end
end
function c = i_op_color(op)
    switch op
        case 'Irrot', c = [0.95 0.55 0.10];   % potential / divergence  -> orange
        case 'Solen', c = [0.55 0.20 0.85];   % stream / curl           -> purple
        otherwise,    c = [0.40 0.40 0.40];   % magnitude               -> gray
    end
end
function r = i_first_results(T)
    r = '';
    for k = 1:numel(T.Groups), if ~isempty(T.Groups(k).ResultsFile), r = T.Groups(k).ResultsFile;  return; end; end
end
function [rows, occMap] = i_group_atoms(T, g)
    rows = {};  occMap = zeros(0,3);
    G = T.Groups(g);
    [~, ord] = sort(G.times(1,:));
    for k = ord
        ch = '+';  if (k <= numel(G.charge))   && (G.charge(k) < 0), ch = '-'; end
        sv = 0;    if (k <= numel(G.strength)),  sv = G.strength(k); end
        vx = 0;    if (k <= numel(G.vertices)),  vx = G.vertices(k); end
        rows{end+1} = sprintf(' %8.3fs  %s  %8.3g  v%d', G.times(1,k), ch, sv, vx); %#ok<AGROW>
        occMap(end+1,:) = [g, k, 0]; %#ok<AGROW>
    end
end


%% ===== SET TARGET (called by view_dynamics) =====
function SetTarget(hFig, T) %#ok<DEFNU>
    file = '';
    if ~isempty(hFig) && ishandle(hFig), file = getappdata(hFig, 'DynamicsFile'); end
    setappdata(0, 'DynamicsTarget', struct('hFig',hFig, 'T',T, 'file',file, 'curGroup',0, ...
        'nodeList',{ {} }, 'nodeInfo',[], 'occMap',[], 'Lambda',[], 'showPhase',[1 1 1 1], ...
        'nav', bst_dynamics('NewGroup', 'cursor')));
    % the scale driver needs the source eigenspectrum (Lambda) to build heat-kernel params
    if ~isempty(hFig) && ishandle(hFig)
        St = getappdata(hFig, 'HelmholtzState');
        if ~isempty(St) && isfield(St,'Lambda') && ~isempty(St.Lambda)
            st = getappdata(0, 'DynamicsTarget');  st.Lambda = St.Lambda;  setappdata(0, 'DynamicsTarget', st);
        end
    end
    BuildTree();
end


%% ===== BUILD THE BAND-STACK TREE =====
% Top-level EXTENDED group = a STACK that expands to its time-window leaves (select a
% window -> its atoms on the right). Top-level SIMPLE group (e.g. a recorded band-Function
% group) = a stack with NO leaves; selecting it lists its atoms on the right.
function BuildTree()
    import javax.swing.tree.*;
    ctrl = bst_get('PanelControls', 'Dynamics');
    st   = getappdata(0, 'DynamicsTarget');
    if isempty(ctrl) || isempty(st), return; end
    T = st.T;
    parents = {T.Groups.parent};
    root = DefaultMutableTreeNode('Atoms');
    nodeList = {};  nodeInfo = struct('kind',{},'g',{},'w',{});
    for g = find(cellfun(@isempty, parents))         % top-level (band) groups
        G = T.Groups(g);
        nOcc = max(size(G.times, 2), numel(G.vertices));
        stackNode = DefaultMutableTreeNode(sprintf('%s  (%d)', G.label, nOcc));
        root.add(stackNode);
        nodeList{end+1} = stackNode;  nodeInfo(end+1) = struct('kind','stack','g',g,'w',0); %#ok<AGROW>
        if (size(G.times,1) == 2)                    % extended -> one leaf per time window
            for w = 1:size(G.times,2)
                leaf = DefaultMutableTreeNode(sprintf(' %.3f - %.3f s', G.times(1,w), G.times(2,w)));
                stackNode.add(leaf);
                nodeList{end+1} = leaf;  nodeInfo(end+1) = struct('kind','window','g',g,'w',w); %#ok<AGROW>
            end
        end
    end
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
    ctrl.jTree.setModel(DefaultTreeModel(root));
    for r = 0:(root.getChildCount()-1), ctrl.jTree.expandRow(r); end
    st.nodeList = nodeList;  st.nodeInfo = nodeInfo;  st.occMap = [];
    setappdata(0, 'DynamicsTarget', st);
    ctrl.jListOccur.setModel(javax.swing.DefaultListModel());
end


%% ===== TREE SELECTION =====
function TreeSel_Callback()
    ctrl = bst_get('PanelControls', 'Dynamics');
    st   = getappdata(0, 'DynamicsTarget');
    if isempty(ctrl) || isempty(st), return; end
    sel = ctrl.jTree.getLastSelectedPathComponent();
    info = [];
    if ~isempty(sel)
        for i = 1:numel(st.nodeList)
            if sel.equals(st.nodeList{i}); info = st.nodeInfo(i); break; end
        end
    end
    model = javax.swing.DefaultListModel();
    occMap = zeros(0,3);  % [groupIdx, occIdx, <unused>]  -> the atom to highlight per right-list row
    hSel = findobj(st.hFig, 'Tag', 'AtomSel');  if ~isempty(hSel), set(hSel, 'Visible', 'off'); end
    if ~isempty(info)
        st.curGroup = info.g;
        if strcmp(info.kind, 'stack') && (size(st.T.Groups(info.g).times,1) == 1) && ~isempty(st.T.Groups(info.g).vertices)
            % simple (recorded) group: list ALL its atoms flat on the right
            [rows, occMap] = i_group_atoms(st.T, info.g);
            for k = 1:numel(rows), model.addElement(rows{k}); end
        elseif strcmp(info.kind, 'window')
            [rows, occMap] = i_window_atoms(st.T, info.g, info.w, i_field(st,'showPhase',[1 1 1 1]));
            for k = 1:numel(rows), model.addElement(rows{k}); end
            i_jump(st.T.Groups(info.g).times(1, info.w));   % selecting a window jumps to its onset
        elseif strcmp(info.kind, 'atom')
            % single atom of a simple band group: list just it (and it is highlightable)
            G = st.T.Groups(info.g);
            if (info.w <= numel(G.vertices))
                model.addElement(sprintf(' %.3fs  %-8s  v%d', G.times(1,info.w), i_phase_type(G), G.vertices(info.w)));
                occMap(end+1,:) = [info.g, info.w, 0]; %#ok<AGROW>
            end
            i_jump(G.times(1, info.w));                     % and to the atom's time
        elseif strcmp(info.kind, 'detevt')
            evs = panel_record('GetEvents', [], 1);
            if (info.g <= numel(evs)) && ~isempty(evs(info.g).times)
                i_jump(evs(info.g).times(1,1));
            end
        end
    else
        st.curGroup = 0;
    end
    ctrl.jListOccur.setModel(model);
    st.occMap = occMap;  setappdata(0, 'DynamicsTarget', st);
end


%% ===== OCCURRENCE SELECTION -> highlight marker + jump time =====
function OccurSel_Callback()
    ctrl = bst_get('PanelControls', 'Dynamics');
    st   = getappdata(0, 'DynamicsTarget');
    if isempty(ctrl) || isempty(st) || isempty(st.occMap) || ~ishandle(st.hFig), return; end
    row = ctrl.jListOccur.getSelectedIndex() + 1;
    if (row < 1) || (row > size(st.occMap,1)), return; end
    g = st.occMap(row,1);  o = st.occMap(row,2);
    G = st.T.Groups(g);
    hSel = findobj(st.hFig, 'Tag', 'AtomSel');
    GroupsPosOff = getappdata(st.hFig, 'GroupsPosOff');
    if ~isempty(hSel) && ~isempty(GroupsPosOff) && (g <= numel(GroupsPosOff)) && (o <= size(GroupsPosOff{g},1))
        p = GroupsPosOff{g}(o,:);
        set(hSel, 'XData',p(1), 'YData',p(2), 'ZData',p(3), 'Visible','on');
    elseif ~isempty(hSel)
        set(hSel, 'Visible', 'off');
    end
    if (o <= size(G.times,2)), i_jump(G.times(1,o)); end
end


%% ===== FLAT, TIME-SORTED ATOMS WITHIN ONE WINDOW (across the band's phase children) =====
function [rows, occMap] = i_window_atoms(T, gBand, w, showPhase)
    if (nargin < 4) || isempty(showPhase), showPhase = [1 1 1 1]; end
    rows = {};  occMap = zeros(0,3);
    G = T.Groups(gBand);
    on = G.times(1,w);  off = G.times(2,w);
    children = find(strcmpi({T.Groups.parent}, G.label));
    times = [];  phases = {};  verts = [];  cc = [];  oo = [];
    for c = children(:)'
        Gc = T.Groups(c);
        pk = i_phase_index(i_phase_type(Gc));
        if (pk >= 1) && ~showPhase(pk), continue; end          % phase filtered out
        nO = size(Gc.times, 2);                 % iterate by occurrence (markers may have NO vertices)
        for o = 1:nO
            t = Gc.times(1,o);
            if (t >= on - 1e-9) && (t <= off + 1e-9)
                times(end+1)  = t;            %#ok<AGROW>
                phases{end+1} = i_phase_type(Gc); %#ok<AGROW>
                if (o <= numel(Gc.vertices)), verts(end+1) = Gc.vertices(o); else, verts(end+1) = NaN; end %#ok<AGROW>
                cc(end+1) = c;  oo(end+1) = o; %#ok<AGROW>
            end
        end
    end
    [~, ord] = sort(times);
    for k = ord
        if isnan(verts(k))                      % temporal marker (no source yet): time + phase
            rows{end+1} = sprintf(' %8.3fs  %-8s', times(k), phases{k}); %#ok<AGROW>
        else
            rows{end+1} = sprintf(' %8.3fs  %-8s  v%d', times(k), phases{k}, verts(k)); %#ok<AGROW>
        end
        occMap(end+1,:) = [cc(k), oo(k), 0]; %#ok<AGROW>
    end
end


%% ===== FILE menu =====
function FileOpen()
    [fn, pth] = uigetfile('dynamics_*.mat', 'Open dynamics table');
    if isequal(fn, 0), return; end
    view_dynamics(fullfile(pth, fn));
end
function FileSave()
    st = getappdata(0, 'DynamicsTarget');  if isempty(st), return; end
    if isempty(st.file), FileSaveAs();  return; end
    bst_dynamics('Save', st.file, st.T);
end
function FileSaveAs()
    st = getappdata(0, 'DynamicsTarget');  if isempty(st), return; end
    [fn, pth] = uiputfile('dynamics_*.mat', 'Save dynamics table');
    if isequal(fn, 0), return; end
    out = fullfile(pth, fn);
    bst_dynamics('Save', out, st.T);
    st.file = out;  setappdata(0, 'DynamicsTarget', st);
    if ~isempty(st.hFig) && ishandle(st.hFig), setappdata(st.hFig, 'DynamicsFile', out); end
end


%% ===== ATOMS menu (act on the selected band group, then refresh) =====
function AtomAddGroup()
    st = getappdata(0, 'DynamicsTarget');  if isempty(st), return; end
    name = java_dialog('input', 'Name for the new atom group:', 'Add group', [], '');
    if isempty(name), return; end
    G = bst_dynamics('NewGroup', name);  G.color = [0.6 0.6 0.6];
    if ~isempty(st.T.Groups), G.SurfaceFile = st.T.Groups(1).SurfaceFile; G.DataFile = st.T.Groups(1).DataFile; end
    st.T = bst_dynamics('AddGroup', st.T, G);
    i_apply(st);
end
function AtomRenameGroup()
    [st, g] = i_selected();  if g < 1, return; end
    old = st.T.Groups(g).label;
    name = java_dialog('input', 'New name:', 'Rename group', [], old);
    if isempty(name) || strcmp(name, old), return; end
    for c = 1:numel(st.T.Groups)
        if strcmp(st.T.Groups(c).parent, old), st.T.Groups(c).parent = name; end
    end
    st.T.Groups(g).label = name;
    i_apply(st);
end
function AtomDeleteGroup()
    [st, g] = i_selected();  if g < 1, return; end
    if ~java_dialog('confirm', sprintf('Delete band "%s" and its atoms?', st.T.Groups(g).label), 'Delete group'), return; end
    lbl = st.T.Groups(g).label;
    kill = strcmpi({st.T.Groups.parent}, lbl);  kill(g) = true;   % the band + its children
    st.T.Groups(kill) = [];  st.T.nGroups = numel(st.T.Groups);  st.curGroup = 0;
    i_apply(st);
end
function AtomSetColor()
    [st, g] = i_selected();  if g < 1, return; end
    c0 = st.T.Groups(g).color;  if isempty(c0), c0 = [0.6 0.6 0.6]; end
    c = uisetcolor(c0, 'Group color');
    if isscalar(c) && (c == 0), return; end   % cancelled
    st.T.Groups(g).color = c(:)';
    i_apply(st);
end
function AtomSort(mode)
    st = getappdata(0, 'DynamicsTarget');  if isempty(st) || isempty(st.T.Groups), return; end
    if strcmpi(mode, 'time'), key = cellfun(@i_firsttime, {st.T.Groups.times});
    else,                     key = lower({st.T.Groups.label}); end
    [~, ord] = sort(key);
    st.T.Groups = st.T.Groups(ord);  st.curGroup = 0;
    i_apply(st);
end


%% ===== helpers =====
function [st, g] = i_selected()
    g = 0;  st = getappdata(0, 'DynamicsTarget');
    if isempty(st), return; end
    g = st.curGroup;
    if g < 1, java_dialog('warning', 'Select a band atom in the tree first.', 'Atoms'); end
end
function i_apply(st)
    setappdata(0, 'DynamicsTarget', st);
    if ~isempty(st.hFig) && ishandle(st.hFig)
        try, view_dynamics('Redraw', st.hFig, st.T, i_field(st,'showPhase',[1 1 1 1])); catch, end %#ok<CTCH>
    end
    BuildTree();
end
function t0 = i_firsttime(times)
    if isempty(times), t0 = inf; else, t0 = times(1,1); end
end
function s = i_str(x)
    if isempty(x), s = '-'; else, s = char(x); end
end
function i_jump(t)   % drive the global time cursor (like Record's JumpToEvent); no-op if no time context
    if isempty(t), return; end
    try, panel_time('SetCurrentTime', t(1)); catch, end %#ok<CTCH>
end
