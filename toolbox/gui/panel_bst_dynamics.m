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

    % ===== MENU BAR (NORTH of the whole panel): File + the rich Atoms menu (scout-style) =====
    jMenuBar = gui_component('MenuBar', jPanelNew, BorderLayout.NORTH);
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
    % EAST close button: a glue pushes it to the right, then an 'x' that ends the whole session
    % (hide the panel + close the linked source figure), like the main GUI's close-all button.
    jMenuBar.add(javax.swing.Box.createHorizontalGlue());
    gui_component('button', jMenuBar, [], '', {IconLoader.ICON_DELETE, java_scaled('dimension', 24, 20)}, ...
        'Close the Dynamics session (hide this panel + close its source figure)', @(h,e)bst_call(@OnCloseSession));

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

    % ===== MAIN (CENTER): atoms list (CENTER, dominant) + action toolbar (EAST) + navigator (SOUTH) =====
    BH = java_scaled('value', 22);
    jPanelMain = gui_component('Panel');
    jPanelMain.setBorder(BorderFactory.createEmptyBorder(0,7,7,7));

    % vertical action toolbar (EAST), scout jToolbar2 style: Detect / Save / Show / Clear / Measure
    jToolbar2 = gui_component('Toolbar', jPanelMain, BorderLayout.EAST);
    jToolbar2.setOrientation(jToolbar2.VERTICAL);
    jToolbar2.setPreferredSize(java_scaled('dimension', 28, 20));
    jToolbar2.setBorder([]);
    gui_component('ToolbarButton', jToolbar2, [], '', IconLoader.ICON_EVT_TYPE_ADD, 'Detect windows: run the band-power detector on the selected band (preview events; not saved)', @(h,e)bst_call(@OnDetect));
    gui_component('ToolbarButton', jToolbar2, [], '', IconLoader.ICON_SAVE, 'Save: commit the staged detection to atoms, or pin the current source region as an atom', @(h,e)bst_call(@OnSave));
    jToolbar2.addSeparator();
    jShow = gui_component('ToolbarToggle', jToolbar2, [], '', IconLoader.ICON_SCOUT_ALL, 'Show all atom phases', @(h,e)bst_call(@OnShowAll));
    jShow.setSelected(true);
    gui_component('ToolbarButton', jToolbar2, [], '', IconLoader.ICON_EVT_TYPE_DEL, 'Clear: discard the staged detection preview (no save)', @(h,e)bst_call(@OnClearDetection));
    jToolbar2.addSeparator();
    gui_component('ToolbarButton', jToolbar2, [], '', IconLoader.ICON_PROPERTIES, 'Measure: choose the differential map (None / Divergence / Curl / Potential / Stream)', @(h,e)bst_call(@()OnMeasureMenu(h)));

    % atoms list = CENTER (primary, dominant)
    jPanelMain.add(jSplit, BorderLayout.CENTER);

    % ATOM section (SOUTH): the atom tool -- pick a filter, set its contextual params, localize on the
    % cortex, threshold + store. A dynamics atom = a thresholded localized eigfilter filter (Scout+Event).
    % Replaces the old (non-functional) Frequency/Source/Scale navigator.
    [atomKeys, atomDisp] = panel_eigenfilter_design('AtomKernels');
    jAtom = JPanel();  jAtom.setLayout(BoxLayout(jAtom, BoxLayout.Y_AXIS));
    jAtom.setBorder(java_scaled('titledborder', 'Atom'));
    jRowK = gui_river([0 0], [0 2 0 2]);
    gui_component('label', jRowK, [], 'Filter: ');
    jKernel = gui_component('combobox', jRowK, 'hfill', [], {atomDisp}, [], [], []);
    java_setcb(jKernel, 'ActionPerformedCallback', @(h,e)bst_call(@OnKernelChange));
    jAtom.add(jRowK);
    jAtomParams = gui_river([0 0], [0 2 0 2]);                  % contextual per-kernel sliders live here
    jAtom.add(jAtomParams);
    jRowA = gui_river([0 0], [0 2 0 2]);
    jLocalize = gui_component('toggle', jRowA, [], 'Localize', {Insets(0,0,0,0), Dimension(java_scaled('value',64),BH)}, 'Localize: ON = click a cortex vertex to seed the atom; OFF = clear the seed', @(h,e)bst_call(@OnLocalize));
    gui_component('label', jRowA, [], '  Thr ');
    jThresh = gui_component('text', jRowA, [], '0.5', {Dimension(java_scaled('value',40),BH)}, 'Level-set threshold (fraction of peak) for the stored Scout + Event', []);
    jStore = gui_component('button', jRowA, [], 'Store', {Dimension(java_scaled('value',56),BH)}, 'Store: threshold the atom into a Scout + Event group (carrying its kernel generator)', @(h,e)bst_call(@OnStore)); %#ok<NASGU>
    jAtom.add(jRowA);
    jPanelMain.add(jAtom, BorderLayout.SOUTH);
    % initial sliders for the default kernel (diffusion) with placeholder bounds; SetTarget rebuilds
    % them with the real spectrum once a surface is linked.
    panel_eigenfilter_design('BuildAtomSliders', jAtomParams, atomKeys{1}, i_atom_default_bounds(), []);

    jPanelNew.add(jPanelMain, BorderLayout.CENTER);
    bstPanelNew = BstPanel(panelName, jPanelNew, struct( ...
        'jTree',jTree, 'jListOccur',jListOccur, 'jMenuFile',jMenuFile, 'jMenuAtoms',jMenuAtoms, ...
        'jKernel',jKernel, 'jAtomParams',jAtomParams, 'jLocalize',jLocalize, 'jThresh',jThresh, ...
        'atomKeys',{atomKeys}, 'jShow',jShow, 'jPhaseItems',jPhaseItems));
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
    loc = bst_dynamics('NewLoc', axis);
    % Time has NO widget -- the global Brainstorm cursor drives it: read the live current time.
    if strcmp(axis, 'time')
        global GlobalData; %#ok<TLEV>
        if ~isempty(GlobalData) && ~isempty(GlobalData.UserTimeWindow) && ~isempty(GlobalData.UserTimeWindow.CurrentTime)
            loc.center = GlobalData.UserTimeWindow.CurrentTime;
        end
        loc.extent = 0;
        if isfinite(loc.center), loc.state = 'point'; end
        return;
    end
    switch axis
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
    st.nav = bst_dynamics('Set', st.nav, axis, 1, loc);
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
                st = i_freq_overlay(st, lo, hi);                 % ensure PSD + drive its freq selection
            else
                panel_filter('SetFilters', 0, [], 0, [], 0, [], 0, 0);
                st.curBand = [];  st.curBandName = '';
                st = i_freq_overlay_clear(st);                   % clear the spectrum band strip
            end
        case 'source'
            % center/window are populated by the Region tool (Task 2 syncs them); nothing to drive here
        case 'scale'
            % Smoothing (eigenmode low-pass) is deferred to the eigenvalue-axis work; the Scale
            % widget is inert for now. Record the request so it round-trips, drive nothing.
            if loc.extent > 0
                st.curScale = struct('on',1,'name','heat','params',struct('t',loc.extent));
            else
                st.curScale = struct('on',0,'name','heat','params',[]);
            end
    end
    setappdata(0, 'DynamicsTarget', st);
end


%% ===== MEASUREMENT (differential operator descriptor; not an axis) =====
% Measure toolbar button -> popup menu of the differential maps (None default + Div/Curl/Pot/Stream).
function OnMeasureMenu(jButton) %#ok<DEFNU>
    [~, st] = i_cs();
    curOp = 'none';  if ~isempty(st), curOp = i_field(st, 'curOp', 'none'); end
    items = {'none','None'; 'Divergence','Divergence'; 'Curl','Curl'; 'Potential','Potential'; 'Stream','Stream'};
    jPopup = java_create('javax.swing.JPopupMenu');
    grp    = java_create('javax.swing.ButtonGroup');
    for i = 1:size(items,1)
        op = items{i,1};
        ri = gui_component('RadioMenuItem', jPopup, [], items{i,2}, [], [], @(h,e)bst_call(@()OnMeasurement(op)));
        grp.add(ri);
        if strcmpi(curOp, op), ri.setSelected(true); end
    end
    jPopup.show(jButton, 0, jButton.getHeight());
end

% Apply a chosen differential map (or 'none' = native source display).
function OnMeasurement(name) %#ok<DEFNU>
    [~, st] = i_cs();
    if isempty(st) || isempty(st.hFig) || ~ishandle(st.hFig), return; end
    D = getappdata(st.hFig, 'DynamicsOverlay');
    if ~isempty(D), D.Op = name; setappdata(st.hFig, 'DynamicsOverlay', D); end
    st.curOp = name;  setappdata(0, 'DynamicsTarget', st);
    if strcmpi(name, 'none')
        % restore the native source map: reset the cortex colormap back to 'source' (the overlay
        % left it 'stat2'), re-assert 'source' as the figure's active colorbar (the lazy stat2
        % registration had taken it over), then repaint the native RMS-norm scalar + raw quivers.
        TI = getappdata(st.hFig, 'Surface');
        if ~isempty(D) && ~isempty(TI) && (D.iTess <= numel(TI))
            TI(D.iTess).ColormapType = 'source';  setappdata(st.hFig, 'Surface', TI);
        end
        bst_colormaps('AddColormapToFigure', st.hFig, 'source');
        panel_surface('UpdateSurfaceData', st.hFig);
    else
        view_dynamics('RefreshOverlay', st.hFig);         % free re-select from the per-frame cache
    end
end

% Save toolbar button: commit the staged detection to atoms if one is staged, else pin the cursor.
function OnSave() %#ok<DEFNU>
    if i_has_staged_detection()
        OnSaveDetection();
    else
        OnSaveCursor();
    end
end
function tf = i_has_staged_detection()
    tf = false;
    evs = panel_record('GetEvents', [], 1);
    if isempty(evs), return; end
    isDet = ~cellfun(@isempty, regexp({evs.label}, '_(peak|trough|rising|falling)$|\([0-9.]+-[0-9.]+ Hz\)$', 'once'));
    tf = any(isDet);
end

% Show toolbar toggle: show/hide ALL atom phases at once (syncs the per-phase menu checkboxes).
function OnShowAll() %#ok<DEFNU>
    [ctrl, st] = i_cs();
    if isempty(ctrl) || isempty(st), return; end
    on = true;
    if isfield(ctrl,'jShow') && ~isempty(ctrl.jShow), on = ctrl.jShow.isSelected(); end
    st.showPhase = double([on on on on]);  setappdata(0, 'DynamicsTarget', st);
    if isfield(ctrl,'jPhaseItems') && ~isempty(ctrl.jPhaseItems)
        for ip = 1:4, try, ctrl.jPhaseItems(ip).setSelected(on); catch, end; end %#ok<CTCH>
    end
    i_apply(st);
end

% ===== CLOSE the whole Dynamics session (the 'x' button) =====
% Closes the linked source figure, which (via its DeleteFcn) tears the panel down. If no live
% figure remains, tears down directly. Single teardown path lives in OnFigureDeleted.
function OnCloseSession() %#ok<DEFNU>
    st = getappdata(0, 'DynamicsTarget');
    if ~isempty(st) && isfield(st,'hFig') && ~isempty(st.hFig) && ishandle(st.hFig)
        close(st.hFig);                 % CloseRequestFcn -> bst_figures; DeleteFcn -> OnFigureDeleted
    else
        i_close_panel();                % no live figure: hide the panel + clear state directly
    end
end

% ===== figure DeleteFcn hook: closing the linked figure ends the session (no orphan panel) =====
function OnFigureDeleted(hFig) %#ok<DEFNU>
    st = getappdata(0, 'DynamicsTarget');
    if isempty(st) || ~isfield(st,'hFig') || ~isequal(st.hFig, hFig), return; end   % not our figure / already gone
    i_close_panel();
end

function i_close_panel()
    setappdata(0, 'DynamicsTarget', []);            % clear the session state
    try, gui_hide('Dynamics'); catch, end %#ok<CTCH> % remove the docked tab
end


%% ===== frequency-atlas preset: a chosen band fills the freq center/window fields =====
function i_freq_preset(ctrl)
    if ~isfield(ctrl,'jFreqBand') || isempty(ctrl.jFreqBand), return; end
    sel = char(ctrl.jFreqBand.getSelectedItem());
    if strcmpi(sel, 'none')                                     % broadband: clear fields -> no filter (i_drive off-path)
        ctrl.jFreqC.setText('');  ctrl.jFreqW.setText('');  return;
    end
    b = i_bands();  k = find(strcmpi(b(:,1), sel), 1);
    if isempty(k), return; end                                 % 'custom' -> leave fields as typed
    lo = b{k,2}(1);  hi = b{k,2}(2);
    ctrl.jFreqC.setText(num2str((lo+hi)/2));  ctrl.jFreqW.setText(num2str((hi-lo)/2));
end
function nm = i_freq_name(ctrl)
    nm = '';
    if isfield(ctrl,'jFreqBand') && ~isempty(ctrl.jFreqBand), nm = char(ctrl.jFreqBand.getSelectedItem()); end
    if any(strcmpi(nm,{'custom','none'})), nm = ''; end
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
    [~, st] = i_cs();                                            % i_apply rewrote the target
    st.detSel = []; setappdata(0, 'DynamicsTarget', st);        % clear any stale staged-window index
    i_focus_time(st, [evt(1,1), evt(2,1)]);                      % focus the FIRST detected window
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
    fLoc = bst_dynamics('NewLoc', 'freq');  fLoc.center = mean(band);  fLoc.extent = (band(2)-band(1))/2;  fLoc.label = bandName;
    % rebuild this band's saved groups
    st.T = i_remove_band(st.T, bandName);
    W = bst_dynamics('NewGroup', winLabel);
    W.type='extended';  W.times=evs(iWin).times;  W.epochs=ones(1,size(evs(iWin).times,2));  W.bandName=bandName;
    W.color=[0.5 0.5 0.5];  W.SurfaceFile=st.T.SurfaceFile;  W.DataFile=DataFile;
    W = bst_dynamics('Set', W, 'freq', 1, fLoc);                      % numeric freq tensor index
    st.T = bst_dynamics('AddGroup', st.T, W);
    phaseColors = struct('peak',[0.90 0.10 0.10],'trough',[0.10 0.25 0.90],'rising',[0.10 0.70 0.20],'falling',[0.95 0.55 0.10]);
    for ph = {'peak','trough','rising','falling'}
        lbl = [bandName '_' ph{1}];  ie = find(strcmpi({evs.label}, lbl), 1);
        if isempty(ie) || isempty(evs(ie).times), continue; end
        P = bst_dynamics('NewGroup', lbl);
        P.type='simple';  P.parent=winLabel;  P.phase=process_evt_refphase('PhaseValue', ph{1});
        P.times=evs(ie).times(1,:);  P.epochs=ones(1,size(evs(ie).times,2));  P.bandName=bandName;
        P.color=phaseColors.(ph{1});  P.SurfaceFile=st.T.SurfaceFile;  P.DataFile=DataFile;
        P = bst_dynamics('Set', P, 'freq', 1, fLoc);                  % numeric freq on each child
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
    global GlobalData;
    evs = panel_record('GetEvents', [], 1);
    if isempty(evs), return; end
    % Clearing mutates events (SetEvents sets isModified=1); preserve the recording's clean/dirty
    % state so discarding a render-only preview never dirties/persists the recording.
    iDS = panel_record('GetCurrentDataset');
    wasMod = ~isempty(iDS) && ~isempty(GlobalData.DataSet(iDS).Measures.sFile) && GlobalData.DataSet(iDS).Measures.isModified;
    isDet = ~cellfun(@isempty, regexp({evs.label}, '_(peak|trough|rising|falling)$|\([0-9.]+-[0-9.]+ Hz\)$', 'once'));
    keep = evs(~isDet);
    if isempty(keep), keep = repmat(db_template('event'), 1, 0); end   % empty event array, not []
    panel_record('SetEvents', keep);                             % replace the whole event set with the non-detection ones
    panel_record('UpdateEventsList');
    panel_record('ReplotEvents');                               % re-plot the (now cleared) markers; do NOT ReloadFigures (would reload events from disk)
    if ~isempty(iDS) && ~wasMod, GlobalData.DataSet(iDS).Measures.isModified = 0; end   % discarding preview must not dirty the recording
    [~, st] = i_cs();  if ~isempty(st), st.detSel = []; i_apply(st); end
end


%% ===== shared panel/target accessor =====
function [ctrl, st] = i_cs()
    ctrl = bst_get('PanelControls', 'Dynamics');
    st   = getappdata(0, 'DynamicsTarget');
end


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
function NotifySelection(hFig, axis, range) %#ok<DEFNU,INUSD>
    return;   % freq/time selection-sync retired with the Navigator strip (atom tool phase)
    if i_driving(), return; end %#ok<UNRCH>
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
    if isempty(hFig) || ~ishandle(hFig) || isempty(st.T) || ~isfield(st.T,'DataFile') || isempty(st.T.DataFile), return; end
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
        [sTfStudy, ~, iTf] = bst_get('TimefreqFile', TfInfo.FileName);
        if isempty(sTfStudy) || isempty(iTf) || (iTf == 0), continue; end
        sTf = sTfStudy.Timefreq(iTf);
        if ~isempty(sTf.DataFile) && file_compare(sTf.DataFile, DataFile)
            hSpec = h;  i_fix_spec_xlim(hSpec);  return;
        end
    end
    % 2) precomputed PSD timefreq file for this recording?
    TfFile = i_find_psd_file(DataFile);
    % 3) else compute one
    if isempty(TfFile), TfFile = i_compute_psd(DataFile); end
    if isempty(TfFile), return; end
    try
        hSpec = view_spectrum(TfFile, 'Spectrum');
    catch
        hSpec = [];  return;
    end
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
    else
        st.hSpec = [];                                           % null a stale (deleted) handle
    end
end


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
    % reflect the band name in the combobox; an exact preset match snaps fields+strip to the standard band (magnetic-to-band), 'custom' keeps the dragged values; overlay re-drive is i_driving-guarded
    try, ctrl.jFreqBand.setSelectedItem(nm); catch, end %#ok<CTCH>
end


%% ===== TIME FOCUS: drive the recording's Time Selection box (panel -> view) =====
function i_focus_time(st, win)
    if isempty(win) || (numel(win) < 2) || any(~isfinite(win(1:2))), return; end
    hRec = i_rec_figure(st);
    if isempty(hRec) || ~ishandle(hRec), return; end
    i_driving(true);
    try, figure_timeseries('SetTimeSelectionManual', hRec, [win(1) win(2)]); catch, end %#ok<CTCH>
    i_driving(false);
end


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


%% ===== SYNC the Source block fields from the geodesic tool state =====
function SyncSource() %#ok<DEFNU>
    % The geodesic tool (repurposed as the atom Localize seed-picker) calls this on every Draw.
    % The atom tool needs only the SEED vertex (its kernel + threshold define the spatial extent).
    [ctrl, st] = i_cs();
    if isempty(ctrl) || isempty(st) || ~isfield(ctrl,'jLocalize'), return; end
    gs = bst_geodesic_tool('GetState');
    if isempty(gs) || ~isfield(gs,'seed') || isempty(gs.seed), return; end
    st.atomSeed = double(gs.seed);
    setappdata(0, 'DynamicsTarget', st);
    i_atom_preview();
end


%% ===== SAVE CURSOR: commit the live 4-D cursor (st.nav) as ONE atom =====
function OnSaveCursor() %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(st) || isempty(st.nav), return; end
    st.nav = bst_dynamics('Set', st.nav, 'time', 1, i_read_block(ctrl, 'time'));   % time from the global cursor (no widget)
    lt = bst_dynamics('Get', st.nav, 'time');     lf = bst_dynamics('Get', st.nav, 'freq');
    ls = bst_dynamics('Get', st.nav, 'source');   lk = bst_dynamics('Get', st.nav, 'scale'); %#ok<NASGU>
    if ~isfinite(lt.center)
        java_dialog('warning', 'Move the time cursor first (no cursor time).', 'Save cursor');  return;
    end
    band = st.curBand;  bandName = i_field(st, 'curBandName', '');
    op   = i_field(st, 'curOp', 'none');
    switch op
        case 'Divergence', Func = 'divergence';
        case 'Curl',       Func = 'curl';
        case 'Potential',  Func = 'potential';
        case 'Stream',     Func = 'stream';
        otherwise,         Func = 'source';     % 'none' = the raw source map (no differential measured)
    end
    % measured descriptor at the cursor (operator scalar at the seed, if localized + cached + a
    % differential is selected). In 'none' the atom records only its location (strength NaN).
    strength = NaN;  charge = NaN;
    if ~strcmpi(op,'none') && isfinite(ls.center) && ~isempty(st.hFig) && ishandle(st.hFig)
        D = getappdata(st.hFig, 'DynamicsOverlay');
        if ~isempty(D) && isfield(D,'Cache')
            [~, iT] = bst_memory('GetTimeVector', D.srcDS, D.srcResult, 'CurrentTimeIndex');
            if ~isempty(D.Cache) && isKey(D.Cache, iT)
                sc = view_dynamics('PickScalar', D.Cache(iT), op);
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
    if ~isempty(band), G = bst_dynamics('Set', G, 'freq', 1, lf); end   % numeric freq index
    st.T.Groups(g) = G;  st.T.nGroups = numel(st.T.Groups);
    i_apply(st);
    if ~isempty(st.file), try, bst_dynamics('Save', st.file, st.T); catch, end; end %#ok<CTCH>
    bst_progress('text', sprintf('Saved cursor atom (%s) at %.3f s', Func, lt.center));
end


%% ===== LOAD ATOM: fill the navigator blocks + st.nav from the selected saved atom =====
function OnLoadAtom() %#ok<DEFNU>
    [ctrl, st] = i_cs();
    if isempty(ctrl) || isempty(st) || isempty(st.occMap), java_dialog('warning','Select a saved atom first.','Load into navigator');  return; end
    g = st.occMap(1,1);  o = st.occMap(1,2);
    if (g < 1) || (g > numel(st.T.Groups)), return; end
    G = st.T.Groups(g);
    for axis = {'time','freq','source','scale'}
        ax = axis{1};
        loc = bst_dynamics('Get', G, ax, o);
        st.nav = bst_dynamics('Set', st.nav, ax, 1, loc);
        i_fill_block(ctrl, ax, loc);
        i_drive(ax, loc);                          % drive the viewers to the atom
    end
    % final write: keep curBand/curScale that i_drive set during the loop; overlay the loaded cursor
    st2 = getappdata(0, 'DynamicsTarget');
    st2.nav = st.nav;
    setappdata(0, 'DynamicsTarget', st2);
    bst_progress('text', 'Loaded atom into the navigator');
end

% Write a Localization's center/window into the axis block's fields (source window in mm).
function i_fill_block(ctrl, axis, loc)
    switch axis
        case 'time',   return;   % no widget -- the global cursor drives time; nothing to fill
        case 'freq',   jC=ctrl.jFreqC;  jW=ctrl.jFreqW;   wv=loc.extent;
        case 'source', jC=ctrl.jSrcC;   jW=ctrl.jSrcW;    wv=loc.extent*1000;   % metres -> mm
        case 'scale',  jC=ctrl.jScaleC; jW=ctrl.jScaleW;  wv=loc.extent;
        otherwise, return;
    end
    % Set the freq band combobox FIRST: its ActionPerformed callback (OnFreqPreset) runs on the
    % Swing EDT and overwrites the center/window fields with the band midpoint, so drain the EDT
    % (drawnow) before writing the loaded coords -- otherwise the async preset clobbers them.
    if strcmp(axis,'freq') && isfield(ctrl,'jFreqBand') && ~isempty(loc.label)
        try, ctrl.jFreqBand.setSelectedItem(loc.label);  drawnow; catch, end %#ok<CTCH>
    end
    if isfinite(loc.center), jC.setText(num2str(loc.center)); else, jC.setText(''); end
    if isfinite(wv),         jW.setText(num2str(wv));         else, jW.setText(''); end
end


%% ===== RECORD: shaped field's extrema at the cursor -> atoms =====
% Reads the linked Helmholtz decomposition at the current time, detects extrema of the
% scalar selected by the operator, and appends them to the (band, Function) group --
% tagged with the panel's (time, band, scale, operator) coordinates. Auto-saves.
function OnRecord() %#ok<DEFNU>
    [ctrl, st] = i_cs();
    if isempty(ctrl) || isempty(st) || isempty(st.hFig) || ~ishandle(st.hFig), return; end
    if strcmpi(i_field(st, 'curOp', 'none'), 'none')
        java_dialog('warning', 'Choose a Measurement (Divergence / Curl / Potential / Stream) before recording extrema.', 'Record atoms');
        return;
    end
    D = getappdata(st.hFig, 'DynamicsOverlay');
    if isempty(D)
        java_dialog('warning', 'Record needs the linked dynamics source view (open via a Dirac result).', 'Record atoms');
        return;
    end
    view_dynamics('RefreshOverlay', st.hFig);                     % make sure the cursor frame is computed+cached
    D = getappdata(st.hFig, 'DynamicsOverlay');
    [TimeVec, iT] = bst_memory('GetTimeVector', D.srcDS, D.srcResult, 'CurrentTimeIndex');
    if isempty(D.Cache) || ~isKey(D.Cache, iT)
        java_dialog('warning', 'No decomposition at the current time.', 'Record atoms');  return;
    end
    Ht = D.Cache(iT);  tCur = TimeVec(iT);
    % scalar field + Function from the current operator (all signed)
    op   = i_field(st, 'curOp', 'Divergence');
    Scal = view_dynamics('PickScalar', Ht, op);
    switch op
        case 'Divergence', Func = 'divergence';
        case 'Curl',       Func = 'curl';
        case 'Potential',  Func = 'potential';
        case 'Stream',     Func = 'stream';
        otherwise,         Func = 'divergence';
    end
    signed = true;
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

%% ===== REGION TOOL toggle: ON = pick a region; OFF = clear the Source selection =====
function OnRegionTool() %#ok<DEFNU>
    [ctrl, st] = i_cs();
    state = ctrl_region_state();
    bst_geodesic_tool('Toggle', state);
    if ~state                                            % turning OFF -> clear the Source selection
        if ~isempty(ctrl) && isfield(ctrl,'jSrcC')
            ctrl.jSrcC.setText('');  ctrl.jSrcW.setText('');
        end
        if ~isempty(st)
            st.nav = bst_dynamics('Set', st.nav, 'source', 1, bst_dynamics('NewLoc','source'));   % unlocalized
            setappdata(0, 'DynamicsTarget', st);
        end
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
        case 'Divergence', c = [0.95 0.55 0.10];   % sources / sinks  -> orange
        case 'Curl',       c = [0.55 0.20 0.85];   % vorticity        -> purple
        case 'Potential',  c = [0.90 0.75 0.10];   % source potential -> amber
        case 'Stream',     c = [0.30 0.45 0.85];   % stream function  -> blue
        otherwise,         c = [0.40 0.40 0.40];   % gray
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
        'hSpec',[], 'focusTime',[], 'detSel',[], ...
        'nav', bst_dynamics('NewGroup', 'cursor')));
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
            st.detSel = [];
            % simple (recorded) group: list ALL its atoms flat on the right
            [rows, occMap] = i_group_atoms(st.T, info.g);
            for k = 1:numel(rows), model.addElement(rows{k}); end
        elseif strcmp(info.kind, 'window')
            [rows, occMap] = i_window_atoms(st.T, info.g, info.w, i_field(st,'showPhase',[1 1 1 1]));
            for k = 1:numel(rows), model.addElement(rows{k}); end
            i_jump(st.T.Groups(info.g).times(1, info.w));   % selecting a window jumps to its onset
            i_focus_time(st, st.T.Groups(info.g).times(:, info.w)');   % and focuses the [onset offset] box
            st.detSel = [];                                  % saved window -> no staged-edit target
        elseif strcmp(info.kind, 'atom')
            st.detSel = [];
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
            st.detSel = [];
        elseif strcmp(info.kind, 'detwinroot')
            st.detSel = [];
            evs = panel_record('GetEvents', [], 1);
            if (info.g <= numel(evs)) && ~isempty(evs(info.g).times)
                i_focus_time(st, evs(info.g).times(:,1)');
                st.detSel = [info.g, 1];
            end
        elseif strcmp(info.kind, 'detwin')
            st.detSel = [];
            evs = panel_record('GetEvents', [], 1);
            if (info.g <= numel(evs)) && (info.w <= size(evs(info.g).times,2))
                win = evs(info.g).times(:, info.w)';
                i_jump(win(1));
                i_focus_time(st, win);
                st.detSel = [info.g, info.w];   % remember for time-drag sync-back (Task 7)
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

%% ===== atom preview: normalize a realised field W[nGv x nT] for display =====
% Mirrors the atom designer's normalization. W rows are aligned with the eigenbasis support; Mass
% is the block lumped-mass [nGv x nGv]. One-signed mass-conserving kernels (heat/diffusion) -> per-
% frame unit-mass DENSITY (isSigned=false -> sequential colormap). Zero-mean/oscillatory kernels
% (mexhat/waves) -> the mass integral collapses, fall back to global PEAK (isSigned=true -> diverging).
function [W, isSigned] = i_atom_normalize(W, Mass) %#ok<DEFNU>
    mvec = full(sum(Mass, 2));                       % [nGv x 1] lumped vertex areas
    s    = mvec.' * W;                               % [1 x nT] signed mass per frame
    l1   = mvec.' * abs(W);                          % [1 x nT] absolute mass per frame
    r    = abs(s) ./ max(l1, eps);                   % one-signedness in [0,1]
    if ~isempty(r) && (min(r) > 0.1)
        W = W ./ s;  isSigned = false;               % unit-mass density (broadcast over columns)
    else
        pk = max(abs(W(:)));  if pk > 0, W = W / pk; end
        isSigned = true;                             % relative amplitude (density n/a)
    end
end

%% ===== atom tool: filter selector + contextual params + localize + store =====
function b = i_atom_default_bounds()
    b = struct('scaleMinMM',5, 'scaleMaxMM',120, 'rateMinMM2',25, 'rateMaxMM2',14400);
end

function k = i_atom_current_kernel(ctrl)
    keys = ctrl.atomKeys;
    idx  = max(1, min(numel(keys), double(ctrl.jKernel.getSelectedIndex()) + 1));
    k = keys{idx};
end

% (Re)build the contextual sliders for the selected filter; preview if a seed is already placed.
function OnKernelChange() %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl), return; end
    k = i_atom_current_kernel(ctrl);
    b = i_field(st, 'atomBounds', i_atom_default_bounds());
    panel_eigenfilter_design('BuildAtomSliders', ctrl.jAtomParams, k, b, @()bst_call(@OnParamSettle));
    if ~isempty(st) && ~isempty(i_field(st,'atomSeed',[])), i_atom_preview(); end
end

% A slider drag settled -> re-realise the preview.
function OnParamSettle() %#ok<DEFNU>
    i_atom_preview();
end

% Arm click-to-seed on the cortex (the repurposed geodesic tool); OFF clears the seed + preview.
function OnLocalize() %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl), return; end
    state = ctrl.jLocalize.isSelected();
    bst_geodesic_tool('Toggle', state);
    if ~state && ~isempty(st)
        st.atomSeed = [];  setappdata(0, 'DynamicsTarget', st);
        if ~isempty(st.hFig) && ishandle(st.hFig), view_dynamics('ClearAtomField', st.hFig); end
    end
end

% Build (once) the eigen-axes + physical-scale bounds for the linked surface; cache on st.
function st = i_atom_ensure_axes(st)
    if ~isempty(i_field(st,'atomAx',[])), return; end
    if isempty(st.hFig) || ~ishandle(st.hFig), return; end
    TI = getappdata(st.hFig, 'Surface');  if isempty(TI), return; end
    D  = getappdata(st.hFig, 'DynamicsOverlay');
    iTess = 1;  if ~isempty(D) && isfield(D,'iTess') && ~isempty(D.iTess), iTess = D.iTess; end
    SurfaceFile = TI(iTess).SurfaceFile;
    nFrames = 100;
    ax = bst_eigen('Axes', struct('SurfaceFile',SurfaceFile, 'Variant','Laplace-Beltrami', ...
                   'nModes',60, 'TimeWindow',[0 (nFrames-1)/100], 'SampleRate',100));
    lamAll = ax.Lambda{1}(:);  if numel(ax.Lambda)>1 && ~isempty(ax.Lambda{2}), lamAll = [lamAll; ax.Lambda{2}(:)]; end
    lmax = max(lamAll);  lminPos = min(lamAll(lamAll > 1e-9));
    mm = @(l) 2*pi ./ sqrt(l) * 1000;
    sMin = mm(lmax);  sMax = mm(lminPos);
    st.atomAx     = ax;
    st.atomBounds = struct('scaleMinMM',sMin, 'scaleMaxMM',sMax, 'rateMinMM2',sMin^2, 'rateMaxMM2',sMax^2);
    setappdata(0, 'DynamicsTarget', st);
end

% Realise the atom field for (kernel, slider values, seed) and normalize it for display.
function [W, gv, isSigned] = i_atom_realise(st, kernel, vals, seed)
    ax   = st.atomAx;
    lmax = max(ax.Lambda{1}(:));
    kp   = bst_eigfilter_controls('ToKernel', kernel, vals, lmax);
    [W, gv] = bst_eigenfilter('Atom', ax, kernel, kp, seed);
    blk = 1;
    for bI = 1:numel(ax.GlobalVertices)
        if any(ax.GlobalVertices{bI} == seed), blk = bI; break; end
    end
    [W, isSigned] = i_atom_normalize(W, ax.Mass{blk});
end

% Read the controls + seed, realise, and paint the live preview through the overlay.
function i_atom_preview() %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    seed = i_field(st, 'atomSeed', []);  if isempty(seed), return; end
    st = i_atom_ensure_axes(st);  if isempty(i_field(st,'atomAx',[])), return; end
    k    = i_atom_current_kernel(ctrl);
    vals = panel_eigenfilter_design('ReadAtomVals', ctrl.jAtomParams);
    [W, gv, isSigned] = i_atom_realise(st, k, vals, seed);
    if ~isempty(st.hFig) && ishandle(st.hFig)
        view_dynamics('SetAtomField', st.hFig, W, gv, isSigned);
    end
end

% Threshold the current atom into a stored Scout + Event group (carrying its kernel generator).
function OnStore() %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    seed = i_field(st, 'atomSeed', []);
    if isempty(seed)
        java_dialog('warning', 'Localize the atom first: toggle Localize and click a cortex vertex.', 'Atom tool');  return;
    end
    st = i_atom_ensure_axes(st);  if isempty(i_field(st,'atomAx',[])), return; end
    k    = i_atom_current_kernel(ctrl);
    vals = panel_eigenfilter_design('ReadAtomVals', ctrl.jAtomParams);
    lmax = max(st.atomAx.Lambda{1}(:));
    kp   = bst_eigfilter_controls('ToKernel', k, vals, lmax);
    thr  = str2double(char(ctrl.jThresh.getText()));
    if isnan(thr) || (thr <= 0) || (thr >= 1), thr = 0.5; end
    G    = bst_dynamics('AtomFromKernel', st.atomAx, k, kp, seed, thr);
    st.T = bst_dynamics('AddGroup', st.T, G);
    i_apply(st);
    if ~isempty(st.file), try, bst_dynamics('Save', st.file, st.T); catch, end; end %#ok<CTCH>
end
function s = i_str(x)
    if isempty(x), s = '-'; else, s = char(x); end
end
function i_jump(t)   % drive the global time cursor (like Record's JumpToEvent); no-op if no time context
    if isempty(t), return; end
    try, panel_time('SetCurrentTime', t(1)); catch, end %#ok<CTCH>
end
