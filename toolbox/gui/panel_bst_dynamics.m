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
    % Set operator: per-atom eigenbasis (mirrors the Scout panel's "Set function")
    jMenuOp = gui_component('Menu', jMenuAtoms, [], 'Set operator', IconLoader.ICON_PROPERTIES, [], []); %#ok<NASGU>
    bgOp = javax.swing.ButtonGroup();
    opDefs = {'Geometric','Laplace-Beltrami'; 'Connectomic','LB-Connectome'; 'Tangent (connection Laplacian)','Connection Laplacian'; 'Dirac','Dirac'};
    jOpItems = javaArray('javax.swing.JRadioButtonMenuItem', size(opDefs,1));
    for io = 1:size(opDefs,1)
        opv = opDefs{io,2};
        jit = gui_component('radiomenuitem', jMenuOp, [], opDefs{io,1}, [], [], @(h,e)bst_call(@()OnSetOperator(opv)));
        bgOp.add(jit);  jOpItems(io) = jit;
    end
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
    % EAST close button: a glue pushes it to the right, then an 'x' that ends the whole session
    % (hide the panel + close the linked source figure), like the main GUI's close-all button.
    jMenuBar.add(javax.swing.Box.createHorizontalGlue());
    gui_component('button', jMenuBar, [], '', {IconLoader.ICON_DELETE, java_scaled('dimension', 24, 20)}, ...
        'Close the Dynamics session (hide this panel + close its source figure)', @(h,e)bst_call(@OnCloseSession));

    % --- atom list (filterbank): a BstClusterList of atoms (coloured dot + label), panel_scout style ---
    jListAtoms = java_create('org.brainstorm.list.BstClusterList');
    jListAtoms.setCellRenderer(java_create('org.brainstorm.list.BstClusterListRenderer', 'I', fontSize));
    java_setcb(jListAtoms, 'ValueChangedCallback', @(h,e)bst_call(@AtomsListValueChanged_Callback,h,e));
    jScrollList = JScrollPane(jListAtoms);  jScrollList.setBorder(java_scaled('titledborder',''));
    jScrollList.setPreferredSize(java_scaled('dimension', 360, 300));

    % ===== MAIN (CENTER): atoms list (CENTER, dominant) + action toolbar (EAST) + Atom section (SOUTH) =====
    jPanelMain = gui_component('Panel');
    jPanelMain.setBorder(BorderFactory.createEmptyBorder(0,7,7,7));

    % vertical action toolbar (EAST), panel_scout jToolbar2 style
    TB_DIM = java_scaled('dimension', 25, 25);
    jToolbar2 = gui_component('Toolbar', jPanelMain, BorderLayout.EAST);
    jToolbar2.setOrientation(jToolbar2.VERTICAL);
    jToolbar2.setPreferredSize(java_scaled('dimension', 28, 20));
    jToolbar2.setBorder([]);
    % --- filterbank actions ---
    gui_component('ToolbarButton', jToolbar2, [], '', {IconLoader.ICON_SCOUT_NEW, TB_DIM}, ...
        '<HTML><B>Create atom</B>:<BR><BLOCKQUOTE> - Click to add a default (diffusion) atom to the filterbank<BR> - Edit its parameters in the Atom section below</BLOCKQUOTE>', @(h,e)bst_call(@OnCreateAtom));
    gui_component('ToolbarButton', jToolbar2, [], '', {IconLoader.ICON_SAVE, TB_DIM}, 'Save the filterbank (atom table) to disk', @(h,e)bst_call(@OnSaveFilterbank));
    jToolbar2.addSeparator();
    jLocalize = gui_component('ToolbarToggle', jToolbar2, [], '', {IconLoader.ICON_SCOUT_SEL, TB_DIM}, 'Localize: click a cortex vertex to re-seed the selected atom', @(h,e)bst_call(@OnLocalize));
    gui_component('ToolbarButton', jToolbar2, [], '', {IconLoader.ICON_PROPERTIES, TB_DIM}, 'Threshold: set the level-set threshold for the optional Scout+Event export', @(h,e)bst_call(@OnThresholdMenu));
    jApply = gui_component('ToolbarToggle', jToolbar2, [], '', {IconLoader.ICON_TS_DISPLAY, TB_DIM}, 'Apply: filter the REAL source through the selected atom over a 4 s window (Preview); OFF = impulse response (Design)', @(h,e)bst_call(@OnApply));
    gui_component('ToolbarButton', jToolbar2, [], '', {IconLoader.ICON_TIMEFREQ, TB_DIM}, 'Analyze: decompose the current window''s source through the frame -> spatial scalogram + residual', @(h,e)bst_call(@OnAnalyzeWindow));
    gui_component('ToolbarButton', jToolbar2, [], '', {IconLoader.ICON_SCOUT_SEL, TB_DIM}, 'Localize bands: localize each frame band into a marker atom -> a separate dynamics table', @(h,e)bst_call(@OnLocalizeBands));
    jToolbar2.addSeparator();
    % --- legacy detection / differential maps (untouched this step) ---
    jShow = gui_component('ToolbarToggle', jToolbar2, [], '', {IconLoader.ICON_SCOUT_ALL, TB_DIM}, 'Show all atom phases', @(h,e)bst_call(@OnShowAll));
    jShow.setSelected(true);
    jToolbar2.addSeparator();
    gui_component('ToolbarButton', jToolbar2, [], '', {IconLoader.ICON_PROPERTIES, TB_DIM}, 'Measure: choose the differential map (None / Divergence / Curl / Potential / Stream)', @(h,e)bst_call(@()OnMeasureMenu(h)));

    % atom list = CENTER (primary, dominant)
    jPanelMain.add(jScrollList, BorderLayout.CENTER);

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
    % Impulse-direction control (Task 7): a preset combobox for a quaternion (Dirac) operator, an angle
    % spinner for a tangent (connection-Laplacian) operator, hidden for a scalar operator. Built from
    % i_dir_control_spec(kind); (re)shown by i_build_dir_control, called from i_select_atom_load /
    % OnSetOperator so it always matches the SELECTED atom's operator.
    jRowD = gui_river([0 0], [0 2 0 2]);
    gui_component('label', jRowD, [], 'Direction: ');
    jDirCombo = gui_component('combobox', jRowD, 'hfill', [], {{'Normal','+X','+Y','+Z','Pick-on-surface'}}, [], [], []);
    java_setcb(jDirCombo, 'ActionPerformedCallback', @(h,e)bst_call(@OnPickSeedDir));
    jDirAngle = gui_component('spinner', jRowD, 'tab', []);
    jDirAngle.setModel(javax.swing.SpinnerNumberModel(0, 0, 359, 5));
    java_setcb(jDirAngle, 'StateChangedCallback', @(h,e)bst_call(@OnPickSeedDir));
    jRowD.setVisible(false);                                    % no atom selected yet -> hidden
    jAtom.add(jRowD);
    jRowI = gui_river([0 0], [0 2 0 2]);                        % selected-atom readout (kernel . seed . params)
    jAtomInfo = gui_component('label', jRowI, 'hfill', '');
    jAtom.add(jRowI);
    % initial sliders for the default kernel (diffusion) with placeholder bounds; SetTarget rebuilds
    % them with the real spectrum once a surface is linked.
    panel_eigenfilter_design('BuildAtomSliders', jAtomParams, atomKeys{1}, i_atom_default_bounds(), []);

    % Frame section (SOUTH, below Atom): bounds readout + tight-frame generate + coverage toggle
    jFrame = gui_river([0 0], [0 2 0 2], 'Frame');
    jFrameA = gui_component('label', jFrame, [], 'A —');
    jFrameB = gui_component('label', jFrame, 'tab', 'B —');
    jFrameT = gui_component('label', jFrame, 'tab', 'B/A —');
    jFrame.add('br', javax.swing.JLabel('N'));
    jFrameN = gui_component('spinner', jFrame, 'tab', []);
    jFrameN.setModel(javax.swing.SpinnerNumberModel(int32(6), int32(2), int32(24), int32(1)));
    gui_component('button', jFrame, 'tab', 'Design tight frame', [], 'Replace the bank with an itersine tight frame of N members', @(h,e)bst_call(@OnDesignFrame));
    jFrameShow = gui_component('checkbox', jFrame, 'br', 'Show coverage', [], 'Open/close the frame coverage response view', @(h,e)bst_call(@OnFrameShow));

    % wrap the Atom + Frame sections into a vertical SOUTH container
    jSouth = gui_component('Panel');  jSouth.setLayout(BoxLayout(jSouth, BoxLayout.Y_AXIS));
    jSouth.add(jAtom);  jSouth.add(jFrame);
    jPanelMain.add(jSouth, BorderLayout.SOUTH);

    jPanelNew.add(jPanelMain, BorderLayout.CENTER);
    bstPanelNew = BstPanel(panelName, jPanelNew, struct( ...
        'jListAtoms',jListAtoms, 'jMenuFile',jMenuFile, 'jMenuAtoms',jMenuAtoms, ...
        'jKernel',jKernel, 'jAtomParams',jAtomParams, 'jLocalize',jLocalize, 'jAtomInfo',jAtomInfo, ...
        'jApply',jApply, 'jOpItems',jOpItems, 'opVariants',{opDefs(:,2)'}, ...
        'atomKeys',{atomKeys}, 'jShow',jShow, 'jPhaseItems',jPhaseItems, ...
        'jDirRow',jRowD, 'jDirCombo',jDirCombo, 'jDirAngle',jDirAngle, ...
        'jFrameA',jFrameA, 'jFrameB',jFrameB, 'jFrameT',jFrameT, 'jFrameN',jFrameN, 'jFrameShow',jFrameShow));
end


%% ===== MEASUREMENT (differential operator descriptor; not an axis) =====
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


%% ===== shared panel/target accessor =====
function [ctrl, st] = i_cs()
    ctrl = bst_get('PanelControls', 'Dynamics');
    st   = getappdata(0, 'DynamicsTarget');
end


%% ===== NOTIFY SELECTION (view -> panel) =====
function NotifySelection(hFig, axis, range) %#ok<DEFNU,INUSD>
    % Retired no-op. The freq/time selection-sync was removed with the Navigator
    % strip; figure_spectrum / figure_timeseries still call this hook, so the entry
    % point is retained as a no-op to preserve their contract.
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
    ia = i_field(st, 'curAtom', 0);
    if (ia >= 1) && (ia <= numel(st.T.Groups))
        st.T.Groups(ia).vertices = st.atomSeed;
        ctrl.jAtomInfo.setText(i_atom_detail(st.T.Groups(ia)));
    end
    setappdata(0, 'DynamicsTarget', st);
    i_atom_preview();
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

function val = i_field(st, name, default)
    if isfield(st, name) && ~isempty(st.(name)), val = st.(name); else, val = default; end
end

%% ===== SET TARGET (called by view_dynamics) =====
function SetTarget(hFig, T) %#ok<DEFNU>
    file = '';
    if ~isempty(hFig) && ishandle(hFig), file = getappdata(hFig, 'DynamicsFile'); end
    setappdata(0, 'DynamicsTarget', struct('hFig',hFig, 'T',T, 'file',file, ...
        'curAtom',0, 'atomSeed',[], 'showPhase',[1 1 1 1], 'curOp','none'));
    setappdata(0, 'DynamicsApplyCache', []);                          % new session -> stale projection
    BuildTree();
    i_gate_operators(getappdata(0, 'DynamicsTarget'));
end


%% ===== BUILD THE BAND-STACK TREE =====
% Top-level EXTENDED group = a STACK that expands to its time-window leaves (select a
% window -> its atoms on the right). Top-level SIMPLE group (e.g. a recorded band-Function
% group) = a stack with NO leaves; selecting it lists its atoms on the right.
function BuildTree()
    UpdateAtomList();
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
    st.T.Groups(kill) = [];  st.T.nGroups = numel(st.T.Groups);  st.curAtom = 0;
    i_apply(st);
    i_frame_refresh();
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
    st.T.Groups = st.T.Groups(ord);  st.curAtom = 0;
    i_apply(st);
end


%% ===== helpers =====
function [st, g] = i_selected()
    g = 0;  st = getappdata(0, 'DynamicsTarget');
    if isempty(st), return; end
    g = i_field(st, 'curAtom', 0);
    if g < 1, java_dialog('warning', 'Select an atom in the list first.', 'Atoms'); end
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

%% ===== filterbank: an atom IS a filter (generator), not a thresholded marker =====
% Build a generator atomgroup: kernel + params + seed, with Threshold/region(Scout)/times(Event) UNSET.
function G = i_default_atom(kernelName, kp, seed, surfaceFile, label, operator) %#ok<DEFNU>
    if (nargin < 6) || isempty(operator), operator = 'Laplace-Beltrami'; end
    G = bst_dynamics('NewGroup', label);
    G.Operator     = operator;
    G.KernelName   = kernelName;
    G.KernelParams = kp;
    G.vertices     = seed;
    G.SurfaceFile  = surfaceFile;
end
% One-line summary of an atom's generator for the Atom-section readout.
function s = i_atom_detail(G) %#ok<DEFNU>
    seed = 0;  if ~isempty(G.vertices), seed = G.vertices(1); end
    op = '';  if isfield(G,'Operator') && ~isempty(G.Operator), op = [char(G.Operator) ' | ']; end
    s = sprintf('%s%s . vtx %d', op, G.KernelName, seed);
    if isstruct(G.KernelParams)
        f = fieldnames(G.KernelParams);
        for i = 1:numel(f)
            if strcmpi(f{i}, 'lmax'), continue; end
            v = G.KernelParams.(f{i});
            if isnumeric(v) && isscalar(v), s = sprintf('%s . %s=%.3g', s, f{i}, v); end
        end
    end
end

% Default atom operator for a freshly launched panel. Operators are canonical properties of the cortical
% SURFACE (find-or-create via bst_eigen): scalar->Laplace-Beltrami, tangent->Connection Laplacian,
% 3D vector->Dirac. An atom is a kernel on one of those operators + the time base; it does not depend on
% the dSPM. The panel is launched from an inverse kernel only to know the DATA kind, which sets a sensible
% DEFAULT: constrained (scalar) source -> Laplace-Beltrami; unconstrained (3-vector) source -> Dirac.
% (The user can switch operators afterward within what i_gate_mask allows.)
function op = i_launch_operator(st) %#ok<DEFNU>
    op = 'Laplace-Beltrami';                              % scalar default (constrained, or if unknown)
    src = i_src_resultfile_from_target(st);
    if ~isempty(src)
        try, R = in_bst_results(src, 0, 'nComponents');
            if isequal(R.nComponents, 3), op = 'Dirac'; end   % unconstrained 3-vector -> quaternion operator
        catch %#ok<CTCH>
        end
    end
end

% Per-variant eigen-axes over a 4 s window at the recording Fs; cached (variant|surface) across the session.
function ax = i_atom_axes(st, variant) %#ok<DEFNU>
    ax = [];
    % Dirac-dSPM: use the INVERSE's own eigenbasis (DiracEigenFile) so atoms filter the mode kernel's
    % coefficients directly, at the inverse's full mode count (not a fresh 60-mode canonical basis).
    if strcmp(variant, 'Dirac')
        D = getappdata(st.hFig, 'DynamicsOverlay');
        if ~isempty(D) && i_is_dirac_dspm(D)
            surf = i_atom_surface(st);
            key  = ['dspm|' variant '|' surf];
            Mc = getappdata(0,'DynamicsAtomAx');  if isempty(Mc)||~isa(Mc,'containers.Map'), Mc=containers.Map('KeyType','char','ValueType','any'); end
            if isKey(Mc,key), ax = Mc(key); return; end
            src = i_src_resultfile(D);
            R = in_bst_results(src, 0, 'DiracEigenFile','Eigenvalues','ModeHemisphere','SurfaceFile');
            E = [];  O = [];   % in_bst_eigen/in_bst_operator THROW on a missing/corrupt file -> catch and fall through
            try, E = in_bst_eigen(R.DiracEigenFile);  O = in_bst_operator(E.OperatorFile);  catch, E = [];  O = [];  end %#ok<CTCH>
            if ~isempty(E) && isstruct(E) && isfield(E,'Phi') && ~isempty(E.Phi) && isstruct(O) && isfield(O,'Mass') && numel(O.Mass)==2
                lam = double(R.Eigenvalues(:));  hemi = double(R.ModeHemisphere(:));
                ax = struct('Variant','Dirac','SurfaceFile',R.SurfaceFile, ...
                            'Phi',{E.Phi},'GlobalVertices',{E.GlobalVertices},'Mass',{O.Mass},'Operator',O);
                for h = 1:2
                    ord = find(hemi==h);  ls = sort(lam(ord),'ascend');  ax.Lambda{h} = ls;
                end
                Fs = 100; try, tv = bst_memory('GetTimeVector', D.srcDS, D.srcResult); if numel(tv)>1, Fs=1/median(diff(tv)); end, catch, end %#ok<CTCH>
                nF = max(2, round(4*Fs));  ax.nT = nF;  ax.tlag = (0:nF-1)/Fs;  ax.omega = (0:nF-1)*(Fs/nF);  ax.NFFT = nF;
                Mc(key) = ax;  setappdata(0,'DynamicsAtomAx', Mc);
                return;
            end
            % basis missing/malformed -> fall through to the canonical path
        end
    end
    surf = i_atom_surface(st);  if isempty(surf), return; end
    key = [variant '|' surf];
    M = getappdata(0, 'DynamicsAtomAx');
    if isempty(M) || ~isa(M, 'containers.Map'), M = containers.Map('KeyType','char','ValueType','any'); end
    if isKey(M, key), ax = M(key); return; end
    Fs = 100;
    D = getappdata(st.hFig, 'DynamicsOverlay');
    if ~isempty(D) && isfield(D,'srcDS') && isfield(D,'srcResult')
        try, tv = bst_memory('GetTimeVector', D.srcDS, D.srcResult);  if numel(tv) > 1, Fs = 1/median(diff(tv)); end, catch, end %#ok<CTCH>
    end
    nF = max(2, round(4*Fs));                                       % 4 s window at the recording sample rate
    bst_progress('start', 'Atom', sprintf('Building %s eigenbasis...', variant));
    ax = bst_eigen('Axes', struct('SurfaceFile',surf, 'Variant',variant, 'nModes',60, 'TimeWindow',[0 (nF-1)/Fs], 'SampleRate',Fs));
    bst_progress('stop');
    M(key) = ax;  setappdata(0, 'DynamicsAtomAx', M);
end
function surf = i_atom_surface(st)
    surf = '';
    if ~isfield(st,'hFig') || isempty(st.hFig) || ~ishandle(st.hFig), return; end
    TI = getappdata(st.hFig, 'Surface');  if isempty(TI), return; end
    D = getappdata(st.hFig, 'DynamicsOverlay');  iTess = 1;
    if ~isempty(D) && isfield(D,'iTess') && ~isempty(D.iTess), iTess = D.iTess; end
    surf = TI(iTess).SurfaceFile;
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

% Resolve the (kernel, kp) to realise for the SELECTED atom. Hand-pickable kernels come from the
% combobox + sliders (the live editor); a kernel NOT in the combobox (e.g. an itersine tight-frame
% member) uses the atom's STORED generator directly -- so the atom, not the combobox, is the source
% of truth for what gets previewed/applied.
function [kernel, kp] = i_selected_generator(st, ctrl, lmax)
    kernel = i_atom_current_kernel(ctrl);  kp = struct();
    ia = i_field(st, 'curAtom', 0);
    G = [];  if (ia >= 1) && (ia <= numel(st.T.Groups)), G = st.T.Groups(ia); end
    if ~isempty(G) && ~any(strcmp(ctrl.atomKeys, G.KernelName))     % non-combobox (itersine) -> stored generator
        kernel = G.KernelName;  kp = G.KernelParams;  if ~isstruct(kp), kp = struct(); end
    else                                                            % hand-pickable -> combobox + sliders (live)
        vals = panel_eigenfilter_design('ReadAtomVals', ctrl.jAtomParams);
        kp   = bst_eigfilter_controls('ToKernel', kernel, vals, lmax);
    end
    if ~isfield(kp,'lmax') || isempty(kp.lmax), kp.lmax = lmax; end
end

% (Re)build the contextual sliders for the selected filter; preview if a seed is already placed.
function OnKernelChange() %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl), return; end
    k = i_atom_current_kernel(ctrl);
    ax = []; if ~isempty(st), ax = i_atom_axes(st, i_atom_op(st)); end
    if ~isempty(ax), b = i_atom_bounds(ax); else, b = i_atom_default_bounds(); end
    panel_eigenfilter_design('BuildAtomSliders', ctrl.jAtomParams, k, b, @()bst_call(@OnParamSettle));
    i_atom_writeback();
    if ~isempty(st) && ~isempty(i_field(st,'atomSeed',[])), i_atom_preview(); end
end

% A slider drag settled -> write back to the selected atom + re-realise the preview.
function OnParamSettle() %#ok<DEFNU>
    i_atom_writeback();  i_atom_preview();
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

% Pure realise-core: run the atom on ax at (seed,dir); return the raw field W [C*n x nT], its global
% vertices gv, and the decoded full-surface ambient vectors V3 [nV x 3] ([] for scalar).
function [W, gv, V3] = i_atom_realise_core(ax, kernel, kp, seed, seedDir) %#ok<DEFNU>
    [W, gv] = bst_eigenfilter('Atom', ax, kernel, kp, seed, seedDir);
    [C, kind] = bst_eigenfilter('Fiber', ax);
    nV = 0; for h=1:numel(ax.GlobalVertices), nV = max(nV, max(ax.GlobalVertices{h}(:))); end
    V3 = [];
    switch kind
        case 'quaternion'
            n = numel(gv);  im = [W(2:4:end,1) W(3:4:end,1) W(4:4:end,1)];   % imag 3-vector at frame 1
            V3 = zeros(nV,3);  V3(gv,:) = im;
        case 'tangent'
            % complex (a+bi) in the operator frame -> ambient 3-vector a*e1 + b*e2 (Op.Frame per hemi)
            blk = 1; for h=1:numel(ax.GlobalVertices), if any(ax.GlobalVertices{h}==seed), blk=h; break; end, end
            Fr = ax.Operator.Frame{blk};  a = real(W(:,1));  b = imag(W(:,1));
            V3 = zeros(nV,3);  V3(gv,:) = a.*Fr.e1 + b.*Fr.e2;
    end
    if C == 1, V3 = []; end
end

% Realise the atom field on its operator's eigenbasis; reduce vector/complex bases to magnitude.
function [W, gv, isSigned, V3] = i_atom_realise(st, kernel, kp, seed, variant, seedDir)
    if (nargin < 5) || isempty(variant), variant = 'Laplace-Beltrami'; end
    W = [];  gv = [];  isSigned = false;  V3 = [];
    ax = i_atom_axes(st, variant);  if isempty(ax), return; end
    if ~isstruct(kp), kp = struct(); end
    if ~isfield(kp,'lmax') || isempty(kp.lmax), kp.lmax = max(ax.Lambda{1}(:)); end
    if (nargin < 6) || isempty(seedDir), seedDir = i_atom_default_dir(ax, seed); end
    try
        [Wraw, gv, V3] = i_atom_realise_core(ax, kernel, kp, seed, seedDir);
    catch %#ok<CTCH>
        W = [];  gv = [];  V3 = [];  return;                        % operator not realisable -> caller guards
    end
    nGv = numel(gv);
    W = i_paintable_scalar(Wraw, nGv);
    if size(W,1) ~= nGv, W = [];  V3 = [];  return; end              % unexpected shape -> not paintable (guarded)
    if any(strcmp(variant, {'Laplace-Beltrami','LB-Connectome'}))   % scalar basis: density/peak by kernel class
        [W, isSigned] = i_atom_normalize(W, ax.Mass{i_seed_block(ax, seed)});
    else                                                            % magnitude: peak-normalized, sequential
        pk = max(abs(W(:)));  if pk > 0, W = W / pk; end;  isSigned = false;
    end
end
function blk = i_seed_block(ax, seed)
    blk = 1;
    for bI = 1:numel(ax.GlobalVertices)
        if any(ax.GlobalVertices{bI} == seed), blk = bI; break; end
    end
end

% Default impulse direction per operator dimensionality: Dirac -> seed surface normal (unit 3-vector);
% tangent/scalar -> 1 (frame e1 / amplitude). See atom-operator-applicability.
function dir = i_atom_default_dir(ax, seedVert) %#ok<DEFNU>
    [~, kind] = bst_eigenfilter('Fiber', ax);
    if ~strcmp(kind, 'quaternion'), dir = 1; return; end
    dir = [0 0 1];                                       % fallback if normals are missing
    try
        S = in_tess_bst(ax.SurfaceFile, 0);
        if isfield(S,'VertNormals') && ~isempty(S.VertNormals) && seedVert <= size(S.VertNormals,1)
            n = S.VertNormals(seedVert, :);  if norm(n) > 0, dir = n / norm(n); end
        end
    catch %#ok<CTCH>
    end
    dir = dir(:)';
end

%% ===== impulse-direction control (Task 7): pure model + GUI wiring =====
% Direction-control model for an operator's Fiber kind: scalar -> hidden; tangent -> an angle field
% (theta -> complex(cos theta, sin theta) in the operator frame); quaternion -> a preset dropdown
% (ambient 3-vector dipole direction) with a Pick-on-surface option for an arbitrary target vertex.
function s = i_dir_control_spec(kind) %#ok<DEFNU>
    switch kind
        case 'tangent',    s = struct('show',true,'type','angle','presets',{{}});
        case 'quaternion', s = struct('show',true,'type','preset','presets',{{'Normal','+X','+Y','+Z','Pick-on-surface'}});
        otherwise,         s = struct('show',false,'type','none','presets',{{}});
    end
end
% Named preset -> ambient 3-vector; 'Normal' resolves to the seed normal by the caller (i_atom_default_dir).
function d = i_preset_dir(name) %#ok<DEFNU>
    switch name
        case 'Normal', d = [];                 % resolved to the seed normal by the caller
        case '+X', d = [1 0 0];  case '+Y', d = [0 1 0];  case '+Z', d = [0 0 1];
        otherwise, d = [];
    end
end

% (Re)build the direction control for the given operator Fiber kind: show/hide the row, and swap
% between the preset combobox (quaternion) and the angle spinner (tangent). Disarms any pending
% one-shot cortex pick (WaveletDesignerPick) left over from a previous atom/operator selection.
function i_build_dir_control(ctrl, kind, hFig) %#ok<DEFNU>
    if isempty(ctrl) || ~isfield(ctrl,'jDirRow') || isempty(ctrl.jDirRow), return; end
    if (nargin >= 3) && ~isempty(hFig) && ishandle(hFig) && isappdata(hFig, 'WaveletDesignerPick')
        rmappdata(hFig, 'WaveletDesignerPick');
    end
    spec = i_dir_control_spec(kind);
    ctrl.jDirRow.setVisible(spec.show);
    if ~spec.show, return; end
    isPreset = strcmp(spec.type, 'preset');
    if isfield(ctrl,'jDirCombo') && ~isempty(ctrl.jDirCombo)
        ctrl.jDirCombo.setVisible(isPreset);
        if isPreset                                            % refresh the model, suppress the change callback
            cb = java_getcb(ctrl.jDirCombo, 'ActionPerformedCallback');
            java_setcb(ctrl.jDirCombo, 'ActionPerformedCallback', []);
            ctrl.jDirCombo.setModel(javax.swing.DefaultComboBoxModel(spec.presets));
            ctrl.jDirCombo.setSelectedItem('Normal');
            java_setcb(ctrl.jDirCombo, 'ActionPerformedCallback', cb);
        end
    end
    if isfield(ctrl,'jDirAngle') && ~isempty(ctrl.jDirAngle)
        ctrl.jDirAngle.setVisible(strcmp(spec.type, 'angle'));
    end
    ctrl.jDirRow.revalidate();  ctrl.jDirRow.repaint();
end

% Direction-control callback. No-arg (combobox / angle-spinner change) resolves the direction from the
% control's current value; a numeric arg is a cortex vertex resolved by the one-shot Pick-on-surface
% hook (figure_3d's native WaveletDesignerPick: dir = unit(target - seed)). Either way, stores the
% result on st.atomSeedDir + the selected atom's SeedDir and re-previews (i_atom_preview).
function OnPickSeedDir(pickedVertex) %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    ia = i_field(st, 'curAtom', 0);  if (ia < 1) || (ia > numel(st.T.Groups)), return; end
    seed = i_field(st, 'atomSeed', []);  if isempty(seed), return; end
    ax = i_atom_axes(st, i_atom_op(st));  if isempty(ax), return; end
    dir = [];
    if (nargin >= 1) && ~isempty(pickedVertex)
        if ~isempty(st.hFig) && ishandle(st.hFig) && isappdata(st.hFig, 'WaveletDesignerPick')
            rmappdata(st.hFig, 'WaveletDesignerPick');                     % one-shot: disarm on first use
        end
        try, Surf = in_tess_bst(ax.SurfaceFile, 0); catch, return; end     %#ok<CTCH>
        v = Surf.Vertices(double(pickedVertex),:) - Surf.Vertices(seed,:);
        n = norm(v);  if n <= 0, return; end
        dir = v / n;
    else
        [~, kind] = bst_eigenfilter('Fiber', ax);
        spec = i_dir_control_spec(kind);
        switch spec.type
            case 'preset'
                if ~isfield(ctrl,'jDirCombo') || isempty(ctrl.jDirCombo), return; end
                name = char(ctrl.jDirCombo.getSelectedItem());
                if strcmp(name, 'Pick-on-surface')
                    if ~isempty(st.hFig) && ishandle(st.hFig)
                        setappdata(st.hFig, 'WaveletDesignerPick', @OnPickSeedDir);   % arm one-shot pick
                    end
                    if isfield(ctrl,'jAtomInfo') && ~isempty(ctrl.jAtomInfo)
                        ctrl.jAtomInfo.setText('Click a cortex vertex to set the impulse direction');
                    end
                    return;
                elseif strcmp(name, 'Normal')
                    dir = i_atom_default_dir(ax, seed);
                else
                    dir = i_preset_dir(name);
                end
            case 'angle'
                if ~isfield(ctrl,'jDirAngle') || isempty(ctrl.jDirAngle), return; end
                th  = double(ctrl.jDirAngle.getValue()) * pi / 180;
                dir = complex(cos(th), sin(th));
            otherwise
                return;
        end
    end
    if isempty(dir), return; end
    st.atomSeedDir = dir;
    st.T.Groups(ia).SeedDir = dir;
    setappdata(0, 'DynamicsTarget', st);
    i_atom_preview();
end

% Live preview dispatcher: Apply OFF -> impulse response (Design); Apply ON -> filtered real source (Preview).
% Every existing call site (param edit, operator change, atom select, seed) routes through here, so the
% displayed field always matches the current mode.
function i_atom_preview() %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'Dynamics');
    if ~isempty(ctrl) && isfield(ctrl,'jApply') && ~isempty(ctrl.jApply) && ctrl.jApply.isSelected()
        i_atom_apply();
    else
        i_atom_preview_impulse();
    end
end

% Read the controls + seed, realise on the atom's operator, and paint the impulse-response preview (Design).
function i_atom_preview_impulse()
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    seed = i_field(st, 'atomSeed', []);  if isempty(seed), return; end
    variant = i_atom_op(st);
    axk = i_atom_axes(st, variant);  lmax = 1;  if ~isempty(axk), lmax = max(axk.Lambda{1}(:)); end
    [k, kp] = i_selected_generator(st, ctrl, lmax);
    [W, gv, isSigned, V3] = i_atom_realise(st, k, kp, seed, variant, i_field(st,'atomSeedDir',[]));
    if isempty(W)                                                   % operator not realisable -> guard
        ctrl.jAtomInfo.setText(sprintf('%s: not realisable for this atom', variant));
        if ~isempty(st.hFig) && ishandle(st.hFig), view_dynamics('ClearAtomField', st.hFig); end
        return;
    end
    if ~isempty(st.hFig) && ishandle(st.hFig)
        view_dynamics('SetAtomField', st.hFig, W, gv, isSigned, V3);
    end
end
% The selected atom's operator (Variant), default Laplace-Beltrami.
function v = i_atom_op(st)
    v = 'Laplace-Beltrami';
    ia = i_field(st, 'curAtom', 0);
    if (ia >= 1) && (ia <= numel(st.T.Groups)) && isfield(st.T.Groups(ia),'Operator') && ~isempty(st.T.Groups(ia).Operator)
        v = st.T.Groups(ia).Operator;
    end
end
% Physical-scale bounds from an axes' spectrum (for the contextual sliders).
function b = i_atom_bounds(ax)
    lamAll = ax.Lambda{1}(:);  if numel(ax.Lambda) > 1 && ~isempty(ax.Lambda{2}), lamAll = [lamAll; ax.Lambda{2}(:)]; end
    lmax = max(lamAll);  lminPos = min(lamAll(lamAll > 1e-9));
    mm = @(l) 2*pi ./ sqrt(l) * 1000;  sMin = mm(lmax);  sMax = mm(lminPos);
    b = struct('scaleMinMM',sMin, 'scaleMaxMM',sMax, 'rateMinMM2',sMin^2, 'rateMaxMM2',sMax^2);
end

%% ===== Apply: filter the REAL source through the selected atom (Preview mode) =====
% Filter a scalar surface source field F[nV x nT] through the atom kernel on the operator's
% eigenbasis. Domain-aware: a static spatial kernel g(lambda) uses bst_eigenfilter('Analysis')
% (pure spatial filter); a dynamic ts/js kernel g(lambda,t|omega) uses the joint time-vertex
% transform bst_eigenwavelet('JTVAnalysis') (spatial AND temporal filtering -- the eigenwavelet).
% Scalar bases only; vector/Dirac/Tangent dynamic apply is a follow-up (caller guards).
function Ffilt = i_atom_filter_field(F, ax, variant, kernel, kp) %#ok<DEFNU>
    Ffilt = [];
    if ~any(strcmp(variant, {'Laplace-Beltrami','LB-Connectome'})), return; end
    meta = bst_eigfilter_kernel('info', kernel);
    dom  = '';  if isfield(meta,'domain') && ~isempty(meta.domain), dom = meta.domain; end
    if isempty(dom) || strcmpi(dom, 'static')                       % static g(lambda): pure spatial filter
        EigenMat    = struct();  EigenMat.Phi = ax.Phi;  EigenMat.Lambda = ax.Lambda;  EigenMat.Variant = variant;
        EigenMat.GlobalVertices = ax.GlobalVertices;
        if isfield(ax,'GlobalFaces'), EigenMat.GlobalFaces = ax.GlobalFaces; end
        OperatorMat = struct();  OperatorMat.Mass = ax.Mass;
        [Ffilt, ~, isError] = bst_eigenfilter('Analysis', F, EigenMat, OperatorMat, kernel, kp);
        if isError, Ffilt = []; end
    else                                                           % dynamic ts/js: joint time-vertex filter
        kernels = { struct('name', kernel, 'params', kp) };
        [W, ~, isError] = bst_eigenwavelet('JTVAnalysis', F, ax, kernels);
        if isError || isempty(W), Ffilt = []; else, Ffilt = W(:, :, 1); end
    end
end
% Windowed-source modal coefficients C{h} = Phi{h}' * B{h} * F(gv{h},:), cached by
% (srcResult, iWin, operator). Seed and kernel params do NOT invalidate it (Apply filters the
% whole field). Scalar operators only (F reduced to per-vertex magnitude).
function [C, gvAll] = i_apply_projection(st, ax, D, iWin, nV)
    C = {};  gvAll = [];
    key = sprintf('%s|%d-%d|%s', D.srcResult, iWin(1), iWin(end), i_atom_op(st));
    M = getappdata(0, 'DynamicsApplyCache');
    if ~isempty(M) && isstruct(M) && isfield(M,'key') && strcmp(M.key, key)
        C = M.C;  gvAll = M.gvAll;  return;
    end
    F = double(bst_memory('GetResultsValues', D.srcDS, D.srcResult, [], iWin, 0));
    Fr = i_paintable_scalar(F, nV);                       % scalar per-vertex magnitude field
    C = cell(1, numel(ax.Phi));  gvAll = [];
    for h = 1:numel(ax.Phi)
        if isempty(ax.Phi{h}), continue; end
        gv = ax.GlobalVertices{h}(:);
        C{h} = manifold_ft(ax.Phi{h}, ax.Mass{h}, Fr(gv,:));
        gvAll = [gvAll; gv]; %#ok<AGROW>
    end
    setappdata(0, 'DynamicsApplyCache', struct('key',key, 'C',{C}, 'gvAll',gvAll));
end
% Sample indices of a `secs`-long window from the cursor (round(secs*Fs) samples, clamped).
function iWin = i_cursor_window_core(tv, cursor, secs)
    iWin = [];
    if numel(tv) < 2, return; end
    Fs = 1 / median(diff(tv));  nF = max(1, round(secs * Fs));
    [~, i0] = min(abs(tv - cursor));
    iWin = i0 : min(i0 + nF - 1, numel(tv));
end
function iWin = i_cursor_window_test(tv, cursor, secs) %#ok<DEFNU>
    iWin = i_cursor_window_core(tv, cursor, secs);
end
function iWin = i_cursor_window(srcDS, srcResult, secs) %#ok<DEFNU>
    global GlobalData;
    tv  = bst_memory('GetTimeVector', srcDS, srcResult);
    cur = GlobalData.UserTimeWindow.CurrentTime;  if isempty(cur), cur = tv(1); end
    iWin = i_cursor_window_core(tv, cur, secs);
end

% Apply toggle: ON -> Preview (filtered real source); OFF -> Design (impulse response). i_atom_preview
% dispatches on the toggle, so both branches just re-run the preview.
function OnApply() %#ok<DEFNU>
    i_clear_dirac_sensor_overlay();
    i_atom_preview();
end

% Clear the Dirac filtered-sensor overlay; i_atom_preview redraws it if still in Dirac-Apply.
function i_clear_dirac_sensor_overlay()
    st = getappdata(0,'DynamicsTarget');
    if ~isempty(st)
        hRec = i_rec_figure(st);
        if ~isempty(hRec) && ishandle(hRec), try, view_dynamics('ClearFilteredSensors', hRec); catch, end, end %#ok<CTCH>
    end
end

% Preview: reconstruct the real source over a 4 s window at the cursor, filter it through the selected
% atom's operator with bst_eigenfilter('Analysis'), and paint the per-vertex magnitude on the cortex.
function i_atom_apply() %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    D = getappdata(st.hFig, 'DynamicsOverlay');
    if isempty(D) || ~isfield(D,'srcDS') || ~isfield(D,'srcResult') || isempty(D.srcResult)
        ctrl.jAtomInfo.setText('Apply: no real source linked (Design only)');  return;
    end
    seed = i_field(st, 'atomSeed', []);  if isempty(seed), ctrl.jAtomInfo.setText('Apply: place a seed first'); return; end
    variant = i_atom_op(st);
    ax = i_atom_axes(st, variant);  if isempty(ax), return; end
    lmax   = max(ax.Lambda{1}(:));
    [kernel, kp] = i_selected_generator(st, ctrl, lmax);
    nV     = i_overlay_nv(ax);
    % --- Dirac operator: filter in the Dirac eigenbasis -> cortex magnitude + filtered-sensor overlay ---
    if strcmp(variant, 'Dirac')
        % the Dirac cortex/sensor forward uses a STATIC spatial gain g(lambda); a dynamic (ts/js)
        % kernel is g(lambda,t|omega) -> not supported here (use a static kernel, e.g. heat/mexhat/itersine).
        meta = bst_eigfilter_kernel('info', kernel);
        dom = 'static'; if isfield(meta,'domain') && ~isempty(meta.domain), dom = meta.domain; end
        if any(strcmpi(dom, {'ts','js'}))
            ctrl.jAtomInfo.setText(sprintf('Dirac Apply: %s is dynamic; use a static kernel', kernel));
            return;
        end
        Leig = i_dirac_leadfield(st, ax);
        iWin = i_cursor_window(D.srcDS, D.srcResult, 4);  if isempty(iWin), ctrl.jAtomInfo.setText('Apply: no window'); return; end
        J = double(bst_memory('GetResultsValues', D.srcDS, D.srcResult, [], iWin, 0));   % [3nV x nWin] Dirac field
        g = bst_eigfilter_kernel(kernel, kp);
        bst_progress('start','Atom','Dirac filter -> sensors...');
        [Dfilt, Jfilt] = i_dirac_forward(ax, Leig, J, g);
        bst_progress('stop');
        % cortex: filtered magnitude
        Jmag = i_paintable_scalar(Jfilt, nV);
        pk = max(abs(Jmag(:))); if pk>0, Jmag = Jmag/pk; end
        if ~isempty(st.hFig) && ishandle(st.hFig), view_dynamics('SetFilteredField', st.hFig, Jmag, (1:nV)', iWin, false); end
        % sensor: overlay Dfilt vs raw (only if L_eig available)
        if ~isempty(Leig) && ~isempty(Dfilt)
            hRec = i_rec_figure(st);
            if isempty(hRec) || ~ishandle(hRec), try, hRec = view_timeseries(i_rec_datafile(st), 'MEG'); catch, hRec=[]; end, end %#ok<CTCH>
            tv = bst_memory('GetTimeVector', D.srcDS, D.srcResult);  tWin = tv(iWin);
            if ~isempty(hRec) && ishandle(hRec)
                view_dynamics('SetFilteredSensors', hRec, Dfilt, tWin);
                try, figure_timeseries('SetTimeSelectionManual', hRec, [tWin(1) tWin(end)]); catch, end %#ok<CTCH>  % mark the filtered window
                ctrl.jAtomInfo.setText(sprintf('Dirac | %s [Preview: cortex + %d-sensor overlay]', kernel, size(Dfilt,1)));
            else
                ctrl.jAtomInfo.setText(sprintf('Dirac | %s [cortex filtered; open the recording for the sensor overlay]', kernel));
            end
        else
            ctrl.jAtomInfo.setText('Dirac: cortex filtered (no Dirac-dSPM leadfield -> no sensor view)');
        end
        return;
    end
    % --- reduce the reconstructed field to the operator's expected row layout ---
    % The joint time-vertex filter (dynamic atoms) is scalar-only for now; vector/Dirac/Tangent
    % dynamic real-source Preview is a follow-up. Scalar operators act on the source magnitude.
    if ~any(strcmp(variant, {'Laplace-Beltrami','LB-Connectome'}))
        ctrl.jAtomInfo.setText(sprintf('%s: real-source Preview is scalar-only for now (use Geometric/Connectomic)', variant));
        return;
    end
    % --- domain-gated apply: static g(lambda) uses the cached projection (instant on param drag);
    %     dynamic ts/js kernels keep the joint time-vertex path (cache does not apply). ---
    iWin = i_cursor_window(D.srcDS, D.srcResult, 4);
    if isempty(iWin), ctrl.jAtomInfo.setText('Apply: no recording window'); return; end
    meta = bst_eigfilter_kernel('info', kernel);
    dom  = 'static';  if isfield(meta,'domain') && ~isempty(meta.domain), dom = meta.domain; end
    bst_progress('start', 'Atom', 'Filtering the real source...');
    if strcmpi(dom, 'static')
        [C, ~] = i_apply_projection(st, ax, D, iWin, nV);         % cached by (srcResult,iWin,operator)
        if isempty(C), bst_progress('stop'); ctrl.jAtomInfo.setText('Apply: projection failed'); return; end
        g = bst_eigfilter_kernel(kernel, kp);
        Ffilt = zeros(nV, numel(iWin));
        for h = 1:numel(ax.Phi)
            if isempty(ax.Phi{h}) || isempty(C{h}), continue; end
            gv = ax.GlobalVertices{h}(:);  hgain = g(ax.Lambda{h}(:));
            Ffilt(gv,:) = manifold_ift(ax.Phi{h}, hgain(:) .* C{h});
        end
    else
        F  = double(bst_memory('GetResultsValues', D.srcDS, D.srcResult, [], iWin, 0));
        Fr = i_paintable_scalar(F, nV);
        Ffilt = i_atom_filter_field(Fr, ax, variant, kernel, kp);   % dynamic: JTVAnalysis (unchanged)
    end
    bst_progress('stop');
    if isempty(Ffilt), ctrl.jAtomInfo.setText(sprintf('%s: filter not applicable to this source', variant)); return; end
    Ffilt = i_paintable_scalar(Ffilt, nV);
    if size(Ffilt,1) ~= nV, ctrl.jAtomInfo.setText('Apply: source/operator shape mismatch'); return; end
    pk = max(abs(Ffilt(:)));  if pk > 0, Ffilt = Ffilt / pk; end
    if ~isempty(st.hFig) && ishandle(st.hFig)
        view_dynamics('SetFilteredField', st.hFig, Ffilt, (1:nV)', iWin, false);
    end
    ctrl.jAtomInfo.setText(sprintf('%s | %s  [Preview: filtered source, %d-sample window]', variant, kernel, numel(iWin)));
end

% The recording DataTimeSeries figure for this session (matches st.T.DataFile), or [] if none open.
function hFig = i_rec_figure(st)
    hFig = [];  global GlobalData; %#ok<TLEV>
    recFile = i_rec_datafile(st);  if isempty(recFile), return; end
    for h = bst_figures('GetFiguresByType', {'DataTimeSeries'})'
        [~,~,iDS] = bst_figures('GetFigure', h);
        if ~isempty(iDS) && ~isempty(GlobalData.DataSet(iDS).DataFile) && file_compare(GlobalData.DataSet(iDS).DataFile, recFile)
            hFig = h;  return;
        end
    end
end
% The raw recording DataFile the linked Dirac source was computed from (st.T.DataFile is a link|... path).
function recFile = i_rec_datafile(st)
    recFile = '';
    if isempty(st) || ~isfield(st,'hFig') || isempty(st.hFig) || ~ishandle(st.hFig), return; end
    D = getappdata(st.hFig, 'DynamicsOverlay');
    if isempty(D) || ~isfield(D,'srcResult') || isempty(D.srcResult), return; end
    src = i_src_resultfile(D);  if isempty(src), return; end
    try, RR = in_bst_results(src, 0, 'DataFile'); recFile = RR.DataFile; catch, end %#ok<CTCH>
end

% Reduce a real/complex/vector field [k*nRows x nT] to a per-row magnitude [nRows x nT]
% (scalar passes through). nRows is the divisor for THIS call: the basis-support count nGv
% at the impulse site, the full-surface count nV at the apply site.
function s = i_paintable_scalar(F, nRows)
    if ~isreal(F), F = abs(F); end
    if size(F,1) == nRows, s = F; return; end
    if mod(size(F,1), nRows) == 0
        nc = size(F,1) / nRows;
        s = reshape(sqrt(sum(reshape(F, nc, nRows, []).^2, 1)), nRows, []);
    else
        s = F;                                                     % unexpected shape -> caller guards
    end
end

% Surface vertex count from the operator's eigenbasis support (global indices span the full surface).
function nV = i_overlay_nv(ax)
    nV = 0;
    for b = 1:numel(ax.GlobalVertices), nV = max(nV, max(double(ax.GlobalVertices{b}(:)))); end
end

%% ===== ANALYZE: decompose the current window's real source through the frame -> scalogram =====
% Build a scalogram TimefreqMat [3 x nT x M] {Global,LH,RH} from a bst_eigenwavelet('Scalogram') result.
function FileMat = i_scalogram_timefreq(scal, timeVec, surfaceFile, dataFile, comment)
    FileMat = db_template('timefreqmat');
    FileMat.TF        = scal.energy;                          % [3 x nT x M]
    FileMat.Time      = timeVec(:)';
    FileMat.Freqs     = scal.centers(:);                      % sqrt(lambda) scale centers
    FileMat.RowNames  = {'Global','LH','RH'};
    FileMat.Measure   = 'power';
    FileMat.Method    = 'framescalogram';
    FileMat.DataType  = 'matrix';
    FileMat.SurfaceFile = surfaceFile;
    FileMat.nAvg = 1;  FileMat.Leff = 1;
    if ~isempty(dataFile), FileMat.DataFile = file_short(dataFile); end
    FileMat.Comment = comment;
end

% Analyze the current 4 s window's source through the frame -> scalogram TimefreqMat + residual.
function OnAnalyzeWindow() %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    D = getappdata(st.hFig, 'DynamicsOverlay');
    if isempty(D) || ~isfield(D,'srcDS') || ~isfield(D,'srcResult') || isempty(D.srcResult)
        ctrl.jAtomInfo.setText('Analyze: no real source linked');  return;
    end
    variant = i_atom_op(st);
    if ~any(strcmp(variant, {'Laplace-Beltrami','LB-Connectome'}))
        ctrl.jAtomInfo.setText(sprintf('%s: Analyze is scalar-only for now (use Geometric/Connectomic)', variant));  return;
    end
    ax = i_atom_axes(st, variant);  if isempty(ax), return; end
    fr = i_frame_response(st, ax);
    if fr.nMembers < 1, ctrl.jAtomInfo.setText('Analyze: no static frame members'); return; end
    nV = i_overlay_nv(ax);
    iWin = i_cursor_window(D.srcDS, D.srcResult, 4);
    if isempty(iWin), ctrl.jAtomInfo.setText('Analyze: no recording window'); return; end
    bst_progress('start', 'Frame', 'Analyzing the source through the frame...');
    [C, ~] = i_apply_projection(st, ax, D, iWin, nV);          % reuse B's cache
    scal = bst_eigenwavelet('Scalogram', ax, fr.gCell, C);
    tv = bst_memory('GetTimeVector', D.srcDS, D.srcResult);  tv = tv(iWin);
    surf = ax.SurfaceFile;  srcFile = i_src_resultfile(D);   % D.srcResult is an INDEX; resolve the results filename
    FileMat = i_scalogram_timefreq(scal, tv, surf, srcFile, sprintf('Frame scalogram (window) | %d members', scal_nmembers(scal)));
    TfFile = i_save_scalogram(srcFile, FileMat, 'Frame scalogram (window)');   % find-or-replace in the source study
    bst_progress('stop');
    if ~isempty(TfFile), try, view_timefreq(TfFile, 'SingleSensor'); catch, end, end %#ok<CTCH>
    ctrl.jAtomInfo.setText(sprintf('Analyze: residual %.1f%% (frame completeness)', 100*scal.resScalar));
end
function n = scal_nmembers(scal), n = size(scal.energy, 3); end

% Save (find-or-replace by Comment) a timefreq FileMat into the source result's study; return its path.
% Resolve the results FILE name from the overlay's dataset/result INDICES (D.srcResult is an index).
function rf = i_src_resultfile(D)
    rf = '';  global GlobalData; %#ok<TLEV>
    try, rf = GlobalData.DataSet(D.srcDS).Results(D.srcResult).FileName; catch, end %#ok<CTCH>
end

% True when the linked source is a Dirac-dSPM inverse (carries the mode kernel + its eigenbasis file).
function tf = i_is_dirac_dspm(D) %#ok<DEFNU>
    tf = false;  src = i_src_resultfile(D);  if isempty(src), return; end
    try
        R = in_bst_results(src, 0, 'ImagingKernelMode', 'DiracEigenFile');
        tf = isfield(R,'ImagingKernelMode') && ~isempty(R.ImagingKernelMode) ...
          && isfield(R,'DiracEigenFile')   && ~isempty(R.DiracEigenFile);
    catch %#ok<CTCH>
    end
end

% Mode coefficients c = ImagingKernelMode * recordings(GoodChannel, iWin), split per hemisphere and
% ordered ascending-lambda (mirrors view_eigen_timeseries). Free vs projecting the reconstructed field.
function [cCell, meta] = i_mode_coeffs(st, D, iWin) %#ok<DEFNU>
    cCell = {[],[]};  meta = struct('Eigenvalues',{{[],[]}},'DiracEigenFile','');
    src = i_src_resultfile(D);  if isempty(src), return; end
    R = in_bst_results(src, 0, 'ImagingKernelMode','Eigenvalues','ModeHemisphere','GoodChannel','DiracEigenFile');
    if isempty(R.ImagingKernelMode), return; end
    key = sprintf('%s|%s|%d-%d', R.DiracEigenFile, src, iWin(1), iWin(end));
    M = getappdata(0,'DynamicsModeCoeffCache');
    if ~isempty(M) && isstruct(M) && strcmp(M.key,key), cCell = M.cCell; meta = M.meta; return; end
    gc = R.GoodChannel;  if isempty(gc), gc = 1:size(R.ImagingKernelMode,2); end
    d  = double(bst_memory('GetRecordingsValues', D.srcDS, gc, iWin, 0));   % [nGoodChan x nWin] RAW (no gradmag scaling; the Dirac kernel is built for raw d_raw)
    cAll = double(R.ImagingKernelMode) * d;                             % [nMode x nWin]
    lam  = double(R.Eigenvalues(:));  hemi = double(R.ModeHemisphere(:));
    for h = 1:2
        ord = find(hemi==h);  [ls, s] = sort(lam(ord),'ascend');  ord = ord(s);
        cCell{h} = cAll(ord,:);  meta.Eigenvalues{h} = ls;
    end
    meta.DiracEigenFile = R.DiracEigenFile;
    setappdata(0,'DynamicsModeCoeffCache', struct('key',key,'cCell',{cCell},'meta',meta));
end

% Which operators a source with nComponents supports (order = ctrl.opVariants):
%   {'Laplace-Beltrami','LB-Connectome','Connection Laplacian','Dirac'}.
% Scalar source (1) -> only the two scalar operators; vector (3) -> all; unknown -> permissive.
function m = i_gate_mask(nComponents) %#ok<DEFNU>
    % Order: {Laplace-Beltrami, LB-Connectome, Connection Laplacian, Dirac}. Applicability follows the
    % inverse type: constrained (scalar) supports the two scalar bases + the tangent operator (for the
    % gradient of the scalar map), but not Dirac. Unconstrained (vector) supports Dirac (native) + the
    % two scalar bases (on the norm |J|), but not the tangent-only Connection Laplacian.
    m = true(1,4);
    if     isequal(nComponents, 1), m = logical([1 1 1 0]);   % constrained: LBO, LB-Connectome, Connection Laplacian
    elseif isequal(nComponents, 3), m = logical([1 1 0 1]);   % unconstrained: LBO, LB-Connectome, Dirac
    end
end

% Grey out the Set-operator radio items the linked source can't support (reads nComponents).
function i_gate_operators(st)
    ctrl = bst_get('PanelControls', 'Dynamics');
    if isempty(ctrl) || ~isfield(ctrl,'jOpItems') || isempty(ctrl.jOpItems), return; end
    nComp = [];
    src = i_src_resultfile_from_target(st);
    if ~isempty(src)
        try, R = in_bst_results(src, 0, 'nComponents'); nComp = R.nComponents; catch, end %#ok<CTCH>
    end
    m = i_gate_mask(nComp);
    for k = 1:min(numel(ctrl.jOpItems), numel(m))
        try, ctrl.jOpItems(k).setEnabled(logical(m(k))); catch, end %#ok<CTCH>
    end
    % if the currently-selected op is now disabled, fall back to the first enabled one
    for k = 1:numel(ctrl.jOpItems)
        if ctrl.jOpItems(k).isSelected() && ~m(k)
            f = find(m, 1);  if ~isempty(f), ctrl.jOpItems(f).setSelected(1); OnSetOperator(ctrl.opVariants{f}); end
            break;
        end
    end
end

% Resolve the launched source's results filename from the target (link|... form ok), or '' if none.
function src = i_src_resultfile_from_target(st)
    src = '';
    if isempty(st) || ~isfield(st,'hFig') || isempty(st.hFig) || ~ishandle(st.hFig), return; end
    D = getappdata(st.hFig, 'DynamicsOverlay');
    if isempty(D) || ~isfield(D,'srcDS') || ~isfield(D,'srcResult') || isempty(D.srcResult), return; end
    src = i_src_resultfile(D);      % the C-era index->filename resolver
end

function TfFile = i_save_scalogram(srcFile, FileMat, tag)
    TfFile = '';
    [sStudy, iStudy] = bst_get('AnyFile', srcFile);
    if isempty(sStudy), return; end
    % reuse an existing same-tag file for this source, else make a new path
    old = '';  iOld = [];
    if isfield(sStudy,'Timefreq') && ~isempty(sStudy.Timefreq)
        for i=1:numel(sStudy.Timefreq)
            if strncmp(sStudy.Timefreq(i).Comment, tag, numel(tag)), old = sStudy.Timefreq(i).FileName; iOld = i; break; end
        end
    end
    if ~isempty(old)
        TfFile = file_fullpath(old);
    else
        TfFile = bst_process('GetNewFilename', bst_fileparts(file_fullpath(sStudy.FileName)), 'timefreq_framescalo');
    end
    bst_save(TfFile, FileMat, 'v6');
    if ~isempty(iOld)
        db_add_data(iStudy, file_short(TfFile), FileMat, iOld);     % replace the same tree slot in place
    else
        db_add_data(iStudy, file_short(TfFile), FileMat);
    end
    panel_protocols('UpdateNode', 'Study', iStudy);
end

% Localize each frame band into a marker atom (peak vertex / time window / level set) -> separate dynamicsmat.
function OnLocalizeBands() %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    D = getappdata(st.hFig, 'DynamicsOverlay');
    if isempty(D) || ~isfield(D,'srcResult') || isempty(D.srcResult)
        ctrl.jAtomInfo.setText('Localize bands: no real source linked');  return;
    end
    variant = i_atom_op(st);
    if ~any(strcmp(variant, {'Laplace-Beltrami','LB-Connectome'}))
        ctrl.jAtomInfo.setText(sprintf('%s: Localize bands is scalar-only for now', variant));  return;
    end
    ax = i_atom_axes(st, variant);  if isempty(ax), return; end
    fr = i_frame_response(st, ax);  if fr.nMembers < 1, ctrl.jAtomInfo.setText('Localize: no static frame members'); return; end
    nV = i_overlay_nv(ax);
    iWin = i_cursor_window(D.srcDS, D.srcResult, 4);  if isempty(iWin), return; end
    bst_progress('start', 'Frame', 'Localizing frame bands...');
    [C, ~] = i_apply_projection(st, ax, D, iWin, nV);
    scal = bst_eigenwavelet('Scalogram', ax, fr.gCell, C);
    axL = ax;  axL.Time = bst_memory('GetTimeVector', D.srcDS, D.srcResult);  axL.Time = axL.Time(iWin);  axL.tlag = axL.Time;
    axL.TimeFile = i_src_resultfile(D);   % JTVAtoms reads ax.TimeFile for the table's DataFile
    thr = i_field(st, 'atomThreshold', 0.5);
    T = bst_eigenwavelet('JTVAtoms', scal.W, axL, thr);
    T.SurfaceFile = ax.SurfaceFile;  T.DataFile = i_src_resultfile(D);   % index -> results filename
    T.Comment = sprintf('Frame bands (%s, %d members)', variant, fr.nMembers);
    bst_progress('stop');
    [sStudyL,~] = bst_get('AnyFile', T.DataFile);  studyDirL = pwd;
    if ~isempty(sStudyL), studyDirL = bst_fileparts(file_fullpath(sStudyL.FileName)); end
    out = bst_fullfile(studyDirL, sprintf('dynamics_framebands_%s.mat', datestr_safe()));
    bst_dynamics('Save', out, T);
    try, view_dynamics(out); catch, end %#ok<CTCH>
    ctrl.jAtomInfo.setText(sprintf('Localized %d frame-band atoms', numel(T.Groups)));
end
function s = datestr_safe(), s = sprintf('%09d', mod(round(now*1e5), 1e9)); end

%% ===== filterbank: create / list / select / save =====
% + Create atom: append a default DIFFUSION filter atom (no threshold) on a default seed, select it.
function OnCreateAtom() %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end %#ok<ASGLU>
    op = i_launch_operator(st);                                     % default = the launch source's operator
    ax = i_atom_axes(st, op);  if isempty(ax), return; end
    lmax = max(ax.Lambda{1}(:));  seed = ax.GlobalVertices{1}(1);
    S  = bst_eigfilter_controls('Sliders', 'diffusion', i_atom_bounds(ax));
    vals = [0 0 0];  for i = 1:3, if ~isempty(S(i).def), vals(i) = S(i).def; end, end
    kp = bst_eigfilter_controls('ToKernel', 'diffusion', vals, lmax);  kp.vals = vals;
    sdir = i_atom_default_dir(ax, seed);
    G  = i_default_atom('diffusion', kp, seed, ax.SurfaceFile, sprintf('atom%d', numel(st.T.Groups)+1), op);
    G.SeedDir = sdir;
    st.T = bst_dynamics('AddGroup', st.T, G);  setappdata(0,'DynamicsTarget', st);
    UpdateAtomList();  SetSelectedAtom(numel(st.T.Groups));
    i_frame_refresh();
end

% Rebuild the BstClusterList model from the atom table (one coloured row per atom).
function UpdateAtomList() %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'Dynamics');  if isempty(ctrl) || ~isfield(ctrl,'jListAtoms'), return; end
    st = getappdata(0, 'DynamicsTarget');
    listModel = javax.swing.DefaultListModel();
    pal = [31 119 180; 255 127 14; 44 160 44; 214 39 40; 148 103 189; 140 86 75; 227 119 194; 127 127 127] / 255;
    if ~isempty(st) && ~isempty(st.T) && ~isempty(st.T.Groups)
        for i = 1:numel(st.T.Groups)
            G = st.T.Groups(i);  c = pal(mod(i-1, size(pal,1))+1, :);
            listModel.addElement(org.brainstorm.list.BstListItem(char(G.KernelName), [], char(G.label), i, c(1), c(2), c(3)));
        end
    end
    ctrl.jListAtoms.setModel(listModel);
end

% Select an atom programmatically (callback-suppressed, panel_scout idiom) -> load it.
function SetSelectedAtom(iAtom) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'Dynamics');  if isempty(ctrl), return; end
    cb = java_getcb(ctrl.jListAtoms, 'ValueChangedCallback');
    java_setcb(ctrl.jListAtoms, 'ValueChangedCallback', []);
    ctrl.jListAtoms.setSelectedIndex(iAtom - 1);
    java_setcb(ctrl.jListAtoms, 'ValueChangedCallback', cb);
    st = getappdata(0,'DynamicsTarget');  if ~isempty(st), st.curAtom = iAtom; setappdata(0,'DynamicsTarget', st); end
    i_select_atom_load(iAtom);
end

% User clicked an atom in the list -> load it.
function AtomsListValueChanged_Callback(h, ev) %#ok<DEFNU,INUSL>
    if ev.getValueIsAdjusting(), return; end
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    iSel = double(ctrl.jListAtoms.getSelectedIndex()) + 1;
    if (iSel < 1) || (iSel > numel(st.T.Groups)), return; end
    st.curAtom = iSel;  setappdata(0,'DynamicsTarget', st);
    i_select_atom_load(iSel);
    i_frame_refresh();
end

% Load the selected atom's generator into the Atom section (combobox + sliders + seed) and preview.
function i_select_atom_load(iAtom)
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    if (iAtom < 1) || (iAtom > numel(st.T.Groups)), return; end
    G  = st.T.Groups(iAtom);
    op = 'Laplace-Beltrami';  if isfield(G,'Operator') && ~isempty(G.Operator), op = G.Operator; end
    ax = i_atom_axes(st, op);
    ik = find(strcmp(ctrl.atomKeys, G.KernelName), 1);
    if ~isempty(ik)                                                % hand-pickable kernel: combobox + sliders (live editor)
        ctrl.jKernel.setSelectedIndex(ik - 1);
        if ~isempty(ax), b = i_atom_bounds(ax); else, b = i_atom_default_bounds(); end
        panel_eigenfilter_design('BuildAtomSliders', ctrl.jAtomParams, G.KernelName, b, @()bst_call(@OnParamSettle));
        if isstruct(G.KernelParams) && isfield(G.KernelParams, 'vals') && ~isempty(G.KernelParams.vals)
            panel_eigenfilter_design('SetAtomVals', ctrl.jAtomParams, G.KernelParams.vals);
        end
    else                                                          % non-combobox (itersine tight-frame member): no editable sliders
        ctrl.jAtomParams.removeAll();
        for si = 1:3, ctrl.jAtomParams.putClientProperty(sprintf('atomslot_%d', si), []); end   % clear slot handles (BuildAtomSliders convention)
        ctrl.jAtomParams.revalidate();  ctrl.jAtomParams.repaint();
    end
    i_select_op_radio(op);                                         % check the matching operator radio
    st.atomSeed = G.vertices;
    st.atomSeedDir = [];  if isfield(G,'SeedDir') && ~isempty(G.SeedDir), st.atomSeedDir = G.SeedDir; end
    setappdata(0, 'DynamicsTarget', st);
    kind = 'scalar';  if ~isempty(ax), [~, kind] = bst_eigenfilter('Fiber', ax); end
    i_build_dir_control(ctrl, kind, st.hFig);                       % show/hide + refresh to match this atom's operator
    ctrl.jAtomInfo.setText(i_atom_detail(G));
    i_atom_preview();
    i_frame_refresh();
end

% Persist the edited params back to the selected atom (+ refresh the readout).
function i_atom_writeback()
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    ia = i_field(st, 'curAtom', 0);  if (ia < 1) || (ia > numel(st.T.Groups)), return; end
    % non-combobox kernels (itersine tight-frame members) are not slider-editable; the atom's stored
    % generator is authoritative -- never overwrite it from the (stale) combobox + sliders.
    if ~any(strcmp(ctrl.atomKeys, st.T.Groups(ia).KernelName)), return; end
    ax = i_atom_axes(st, i_atom_op(st));  if isempty(ax), return; end
    k = i_atom_current_kernel(ctrl);  vals = panel_eigenfilter_design('ReadAtomVals', ctrl.jAtomParams);
    lmax = max(ax.Lambda{1}(:));
    kp = bst_eigfilter_controls('ToKernel', k, vals, lmax);  kp.vals = vals;
    st.T.Groups(ia).KernelName = k;  st.T.Groups(ia).KernelParams = kp;
    setappdata(0, 'DynamicsTarget', st);
    ctrl.jAtomInfo.setText(i_atom_detail(st.T.Groups(ia)));
    i_frame_refresh();
end

% Save the filterbank (atom table) to disk.
function OnSaveFilterbank() %#ok<DEFNU>
    [~, st] = i_cs();  if isempty(st) || isempty(st.T), return; end
    f = i_field(st, 'file', '');
    if isempty(f)
        f = java_dialog('save', 'Save filterbank', '', 'dynamics_filterbank.mat');
        if isempty(f), return; end
        st.file = f;  setappdata(0, 'DynamicsTarget', st);
    end
    bst_dynamics('Save', f, st.T);
end

% Set the level-set threshold (for the optional Scout+Event export; not the atom's identity).
function OnThresholdMenu() %#ok<DEFNU>
    [~, st] = i_cs();  if isempty(st), return; end
    cur = i_field(st, 'atomThreshold', 0.5);
    v = java_dialog('input', 'Level-set threshold (0..1) for the optional Scout+Event export:', 'Atom threshold', [], num2str(cur));
    if isempty(v), return; end
    t = str2double(v);  if ~isnan(t) && (t > 0) && (t < 1), st.atomThreshold = t; setappdata(0, 'DynamicsTarget', st); end
end

% Set the selected atom's operator (eigenbasis), rebuild + re-preview on that basis.
function OnSetOperator(variant) %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end %#ok<ASGLU>
    ia = i_field(st, 'curAtom', 0);  if (ia < 1) || (ia > numel(st.T.Groups)), return; end
    st.T.Groups(ia).Operator = variant;  setappdata(0, 'DynamicsTarget', st);
    setappdata(0, 'DynamicsApplyCache', []);                        % operator changed -> stale projection
    i_clear_dirac_sensor_overlay();
    i_select_op_radio(variant);
    i_select_atom_load(ia);                                        % reload sliders/bounds on the new basis + preview
    i_frame_refresh();
end
% Check the operator radio matching the variant.
function i_select_op_radio(variant)
    ctrl = bst_get('PanelControls', 'Dynamics');  if isempty(ctrl) || ~isfield(ctrl,'jOpItems'), return; end
    k = find(strcmp(ctrl.opVariants, variant), 1);
    if ~isempty(k), ctrl.jOpItems(k).setSelected(1); end
end
function s = i_str(x)
    if isempty(x), s = '-'; else, s = char(x); end
end
function i_jump(t)   % drive the global time cursor (like Record's JumpToEvent); no-op if no time context
    if isempty(t), return; end
    try, panel_time('SetCurrentTime', t(1)); catch, end %#ok<CTCH>
end

%% ===== FRAME: coverage S(lambda) + bounds A/B/tightness over the CURRENT operator's atoms =====
% Frame response over the CURRENT operator's atoms: coverage S(lambda) + bounds A/B/tightness.
function fr = i_frame_response(st, ax)
    fr = struct('lam',[], 'S',[], 'A',NaN, 'B',NaN, 'tightness',NaN, 'nMembers',0, 'gCell',{{}});
    if isempty(ax) || ~isfield(ax,'Lambda') || isempty(ax.Lambda), return; end
    variant = i_atom_op(st);
    lamAll = ax.Lambda{1}(:);
    if numel(ax.Lambda) > 1 && ~isempty(ax.Lambda{2}), lamAll = [lamAll; ax.Lambda{2}(:)]; end
    lmax = max(lamAll);  lminPos = min(lamAll(lamAll > 1e-9));  if isempty(lminPos), lminPos = 0; end
    % gather g(lambda) for every atom on the current operator
    gCell = {};
    for i = 1:numel(st.T.Groups)
        G = st.T.Groups(i);
        op = 'Laplace-Beltrami'; if isfield(G,'Operator') && ~isempty(G.Operator), op = G.Operator; end
        if ~strcmp(op, variant), continue; end
        % frame coverage is a SPATIAL concept: a dynamic (ts/js) kernel is g(lambda,t|omega),
        % not a single g(lambda), so it contributes no spatial band -- skip it.
        meta = struct(); try, meta = bst_eigfilter_kernel('info', G.KernelName); catch, end %#ok<CTCH>
        dom = 'static'; if isfield(meta,'domain') && ~isempty(meta.domain), dom = meta.domain; end
        if any(strcmpi(dom, {'ts','js'})), continue; end
        kp = G.KernelParams;  if ~isstruct(kp), kp = struct(); end
        if ~isfield(kp,'lmax') || isempty(kp.lmax), kp.lmax = lmax; end
        try, g = bst_eigfilter_kernel(G.KernelName, kp); catch, continue; end
        if iscell(g), gCell = [gCell, g(:)']; else, gCell{end+1} = g; end %#ok<AGROW>
    end
    fr.nMembers = numel(gCell);  fr.gCell = gCell;
    if fr.nMembers == 0, return; end
    lam = linspace(max(lminPos,eps), lmax, 400)';
    S = zeros(size(lam));
    for m = 1:numel(gCell), y = gCell{m}(lam); S = S + real(y(:)).^2; end
    fr.lam = lam;  fr.S = S;  fr.A = min(S);  fr.B = max(S);
    if fr.nMembers >= 2 && fr.A > 0, fr.tightness = fr.B / fr.A; end   % undefined for a single band
end

% Recompute the frame response and update the Frame labels (+ the coverage view if shown).
function i_frame_refresh()
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    if ~isfield(ctrl,'jFrameA'), return; end
    ax = [];  try, ax = i_atom_axes(st, i_atom_op(st)); catch, end
    fr = i_frame_response(st, ax);
    if fr.nMembers == 0
        ctrl.jFrameA.setText('A —');  ctrl.jFrameB.setText('B —');  ctrl.jFrameT.setText('B/A —');
    else
        ctrl.jFrameA.setText(sprintf('A %.3g', fr.A));
        ctrl.jFrameB.setText(sprintf('B %.3g', fr.B));
        if isnan(fr.tightness), ctrl.jFrameT.setText('B/A —');
        else
            chk = ''; if abs(fr.tightness-1) < 0.05, chk = ' ✓'; end
            ctrl.jFrameT.setText(sprintf('B/A %.3g%s', fr.tightness, chk));
        end
    end
    if isfield(ctrl,'jFrameShow') && ~isempty(ctrl.jFrameShow) && ctrl.jFrameShow.isSelected() && ~isempty(ax)
        bnd = struct('A',fr.A,'B',fr.B,'tightness',fr.tightness);
        gstruct = struct('Kernels',{fr.gCell}, 'Active', max(1,i_field(st,'curAtom',1)), ...
            'OnSelect', @(j)bst_call(@()SetSelectedAtom(j)), 'Coverage', true, 'Bounds', bnd);
        lamMark = ax.Lambda{1}(:);
        try, view_eigfilter_response(gstruct, lamMark, sprintf('Frame coverage (%s)', i_atom_op(st))); catch, end
    end
end

% Show-coverage toggle: open/refresh or close the frame coverage response view.
function OnFrameShow() %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end %#ok<ASGLU>
    if isfield(ctrl,'jFrameShow') && ~isempty(ctrl.jFrameShow) && ctrl.jFrameShow.isSelected()
        i_frame_refresh();
    else
        try, view_eigfilter_response('close'); catch, end %#ok<CTCH>
    end
end

% Design tight frame: replace the bank with N itersine members spanning the current operator's spectrum.
function OnDesignFrame() %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    op = i_atom_op(st);
    ax = i_atom_axes(st, op);  if isempty(ax), return; end
    lamAll = ax.Lambda{1}(:);
    if numel(ax.Lambda) > 1 && ~isempty(ax.Lambda{2}), lamAll = [lamAll; ax.Lambda{2}(:)]; end
    lmax = max(lamAll);  lminPos = min(lamAll(lamAll > 1e-9));  if isempty(lminPos), lminPos = eps; end
    N = 6;  if isfield(ctrl,'jFrameN') && ~isempty(ctrl.jFrameN), N = double(ctrl.jFrameN.getValue()); end
    N = max(2, round(N));
    if ~isempty(st.T.Groups)
        if ~java_dialog('confirm', sprintf('Replace the current %d atom(s) with a %d-member itersine tight frame?', numel(st.T.Groups), N), 'Design tight frame')
            return;
        end
        st.T.Groups(:) = [];  st.T.nGroups = 0;
    end
    seed = ax.GlobalVertices{1}(1);
    for ii = 1:N
        kp = struct('member',ii, 'Nf',N, 'lmin',lminPos, 'lmax',lmax, 'vals',[]);
        G  = i_default_atom('itersine', kp, seed, ax.SurfaceFile, sprintf('itersine %d/%d', ii, N), op);
        st.T = bst_dynamics('AddGroup', st.T, G);
    end
    setappdata(0, 'DynamicsTarget', st);
    UpdateAtomList();
    SetSelectedAtom(1);              % selects member 1, loads it, and triggers i_frame_refresh
    bst_progress('text', sprintf('Designed %d-member itersine tight frame on %s', N, op));
end

%% ===== Dirac sensor forward (Task 3: filtered-SENSOR view for the Preview mode) =====
% Dirac sensor forward: filter J in the Dirac eigenbasis and forward to sensors.
%   ax   : i_atom_axes(st,'Dirac') (per-hemi Phi/Lambda/Mass/GlobalVertices; quaternion basis)
%   Leig : eigenbasis leadfield [nCh x 2K] (L-then-R), from i_dirac_leadfield
%   J    : source 3-vector field [3nV x nT]
%   g    : gain handle g(lambda)
% Returns Dfilt [nCh x nT], Jfilt [3nV x nT] (filtered 3-vector field), cfilt [2K x nT].
function [Dfilt, Jfilt, cfilt] = i_dirac_forward(ax, Leig, J, g) %#ok<DEFNU>
    nV = 0; for h=1:numel(ax.GlobalVertices), nV = max(nV, max(ax.GlobalVertices{h}(:))); end
    nT = size(J, 2);  Jfilt = zeros(3*nV, nT);  cfilt = [];
    for h = 1:numel(ax.Phi)
        Phi = ax.Phi{h};  if isempty(Phi), continue; end
        idx = ax.GlobalVertices{h}(:);  n = numel(idx);  Lam = ax.Lambda{h}(:);
        % embed the 3-vector source into the quaternion imag slots (w=0), per bst_eigenwavelet i_hemimap
        gIn = reshape([(idx-1)*3+1, (idx-1)*3+2, (idx-1)*3+3].', [], 1);            % global 3-vec rows
        lIn = reshape([(0:n-1)*4+2; (0:n-1)*4+3; (0:n-1)*4+4], [], 1);              % local quat imag slots
        U = zeros(4*n, nT);  U(lIn, :) = J(gIn, :);
        C  = manifold_ft(Phi, ax.Mass{h}, U);            % [K x nT] eigenmode coeffs
        Ch = g(Lam) .* C;                                % gain
        cfilt = [cfilt; Ch]; %#ok<AGROW>                 % stack L-then-R -> [2K x nT]
        Uf = manifold_ift(Phi, Ch);                      % [4n x nT] filtered quaternion field
        Jfilt(gIn, :) = Uf(lIn, :);                      % extract imag 3-vector
    end
    if ~isempty(Leig), Dfilt = Leig * cfilt;             % [nCh x nT]
    else,              Dfilt = [];                       % basis mismatch -> no sensor view (cortex still previews)
    end
end

% Dirac eigenbasis leadfield L_eig [nCh x 2K] via bst_dirac(HeadModel); cached per (headmodel,nModes,tau).
function Leig = i_dirac_leadfield(st, ax) %#ok<DEFNU>
    Leig = [];
    D = getappdata(st.hFig, 'DynamicsOverlay');  if isempty(D), return; end
    src = i_src_resultfile(D);  if isempty(src), return; end
    R = in_bst_results(src, 0, 'HeadModelFile');
    hmFile = R.HeadModelFile;
    if isempty(hmFile)
        sS = bst_get('AnyFile', src);
        if isfield(sS,'HeadModel') && ~isempty(sS.HeadModel), hmFile = sS.HeadModel(sS.iHeadModel).FileName; end
    end
    if isempty(hmFile), return; end
    % nModes/Tau MUST match the atom's Dirac eigenbasis (ax) so L_eig's eigenmode ordering aligns with
    % c_filt; pull Tau from the eigen node ax resolved to (ax is built any-Tau) rather than hardcoding it.
    K = size(ax.Lambda{1},1);  tau = 0.5;
    if isfield(ax,'EigenMat') && isstruct(ax.EigenMat) && isfield(ax.EigenMat,'Provenance') ...
            && isstruct(ax.EigenMat.Provenance) && isfield(ax.EigenMat.Provenance,'Tau') ...
            && ~isempty(ax.EigenMat.Provenance.Tau)
        tau = ax.EigenMat.Provenance.Tau;
    end
    key = sprintf('%s|%d|%g', hmFile, K, tau);
    M = getappdata(0, 'DynamicsDiracFwd');
    if ~isempty(M) && isstruct(M) && strcmp(M.key, key), Leig = M.Leig; return; end
    HeadModel = in_bst_headmodel(hmFile);
    CompHM = bst_dirac(HeadModel, 'nModes', K, 'Tau', tau);
    Leig = CompHM.Gain;                     % [nCh x 2K]
    % Guard (spec §4/§8): assert L_eig's Dirac spectrum matches the atom's ax eigenvalues (same
    % eigendecomposition, L-then-R). On mismatch, bail the sensor view (Leig=[]) -> cortex still previews.
    lamAx = ax.Lambda{1}(:);
    if numel(ax.Lambda) > 1 && ~isempty(ax.Lambda{2}), lamAx = [lamAx; ax.Lambda{2}(:)]; end
    okAlign = isfield(CompHM,'Eigenvalues') && ~isempty(CompHM.Eigenvalues) ...
        && (numel(CompHM.Eigenvalues) == numel(lamAx)) ...
        && (max(abs(sort(double(CompHM.Eigenvalues(:))) - sort(double(lamAx)))) <= 1e-6*max(abs(lamAx)) + 1e-12);
    if ~okAlign, Leig = [];  return; end    % basis mismatch -> no sensor view
    setappdata(0, 'DynamicsDiracFwd', struct('key',key, 'Leig',Leig));
end
