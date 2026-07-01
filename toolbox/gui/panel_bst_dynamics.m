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
    jRowI = gui_river([0 0], [0 2 0 2]);                        % selected-atom readout (kernel . seed . params)
    jAtomInfo = gui_component('label', jRowI, 'hfill', '');
    jAtom.add(jRowI);
    jPanelMain.add(jAtom, BorderLayout.SOUTH);
    % initial sliders for the default kernel (diffusion) with placeholder bounds; SetTarget rebuilds
    % them with the real spectrum once a surface is linked.
    panel_eigenfilter_design('BuildAtomSliders', jAtomParams, atomKeys{1}, i_atom_default_bounds(), []);

    jPanelNew.add(jPanelMain, BorderLayout.CENTER);
    bstPanelNew = BstPanel(panelName, jPanelNew, struct( ...
        'jListAtoms',jListAtoms, 'jMenuFile',jMenuFile, 'jMenuAtoms',jMenuAtoms, ...
        'jKernel',jKernel, 'jAtomParams',jAtomParams, 'jLocalize',jLocalize, 'jAtomInfo',jAtomInfo, ...
        'jApply',jApply, 'jOpItems',jOpItems, 'opVariants',{opDefs(:,2)'}, ...
        'atomKeys',{atomKeys}, 'jShow',jShow, 'jPhaseItems',jPhaseItems));
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
    BuildTree();
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

% Launch-derived default operator: Dirac if the source result is a Dirac inverse, else Laplace-Beltrami.
function op = i_launch_operator(st) %#ok<DEFNU>
    op = 'Laplace-Beltrami';
    cmt = '';
    if isfield(st,'srcComment')
        cmt = st.srcComment;
    elseif isfield(st,'hFig') && ~isempty(st.hFig) && ishandle(st.hFig)
        D = getappdata(st.hFig, 'DynamicsOverlay');
        if ~isempty(D) && isfield(D,'srcDS') && isfield(D,'srcResult')
            try
                gd = []; global GlobalData; gd = GlobalData; %#ok<TLEV>
                cmt = gd.DataSet(D.srcDS).Results(D.srcResult).Comment;
            catch %#ok<CTCH>
            end
        end
    end
    if ~isempty(cmt) && contains(lower(char(cmt)), 'dirac'), op = 'Dirac'; end
end

% Per-variant eigen-axes over a 4 s window at the recording Fs; cached (variant|surface) across the session.
function ax = i_atom_axes(st, variant) %#ok<DEFNU>
    ax = [];
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

% Realise the atom field on its operator's eigenbasis; reduce vector/complex bases to magnitude.
function [W, gv, isSigned] = i_atom_realise(st, kernel, vals, seed, variant)
    if (nargin < 5) || isempty(variant), variant = 'Laplace-Beltrami'; end
    W = [];  gv = [];  isSigned = false;
    ax = i_atom_axes(st, variant);  if isempty(ax), return; end
    lmax = max(ax.Lambda{1}(:));
    kp   = bst_eigfilter_controls('ToKernel', kernel, vals, lmax);
    try
        [W, gv] = bst_eigenfilter('Atom', ax, kernel, kp, seed);
    catch %#ok<CTCH>
        W = [];  gv = [];  return;                                  % operator not realisable -> caller guards
    end
    nGv = numel(gv);
    W = i_paintable_scalar(W, nGv);
    if size(W,1) ~= nGv, W = [];  return; end                       % unexpected shape -> not paintable (guarded)
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
    k    = i_atom_current_kernel(ctrl);
    vals = panel_eigenfilter_design('ReadAtomVals', ctrl.jAtomParams);
    [W, gv, isSigned] = i_atom_realise(st, k, vals, seed, variant);
    if isempty(W)                                                   % operator not realisable -> guard
        ctrl.jAtomInfo.setText(sprintf('%s: not realisable for this atom', variant));
        if ~isempty(st.hFig) && ishandle(st.hFig), view_dynamics('ClearAtomField', st.hFig); end
        return;
    end
    if ~isempty(st.hFig) && ishandle(st.hFig)
        view_dynamics('SetAtomField', st.hFig, W, gv, isSigned);
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
    i_atom_preview();
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
    kernel = i_atom_current_kernel(ctrl);
    vals   = panel_eigenfilter_design('ReadAtomVals', ctrl.jAtomParams);
    kp     = bst_eigfilter_controls('ToKernel', kernel, vals, lmax);
    nV     = i_overlay_nv(ax);
    % --- reconstruct the windowed real source (imaging kernel x recording, via bst_memory) ---
    iWin = i_cursor_window(D.srcDS, D.srcResult, 4);
    if isempty(iWin), ctrl.jAtomInfo.setText('Apply: no recording window'); return; end
    bst_progress('start', 'Atom', 'Filtering the real source...');
    F = double(bst_memory('GetResultsValues', D.srcDS, D.srcResult, [], iWin, 0));
    % --- reduce the reconstructed field to the operator's expected row layout ---
    % The joint time-vertex filter (dynamic atoms) is scalar-only for now; vector/Dirac/Tangent
    % dynamic real-source Preview is a follow-up. Scalar operators act on the source magnitude.
    if any(strcmp(variant, {'Laplace-Beltrami','LB-Connectome'}))
        Fr = i_paintable_scalar(F, nV);                             % scalar operator: per-vertex magnitude
    else
        bst_progress('stop');
        ctrl.jAtomInfo.setText(sprintf('%s: real-source Preview is scalar-only for now (use Geometric/Connectomic)', variant));
        return;
    end
    % --- filter through the atom, reduce to a paintable per-vertex magnitude ---
    Ffilt = i_atom_filter_field(Fr, ax, variant, kernel, kp);
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
    G  = i_default_atom('diffusion', kp, seed, ax.SurfaceFile, sprintf('atom%d', numel(st.T.Groups)+1), op);
    st.T = bst_dynamics('AddGroup', st.T, G);  setappdata(0,'DynamicsTarget', st);
    UpdateAtomList();  SetSelectedAtom(numel(st.T.Groups));
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
end

% Load the selected atom's generator into the Atom section (combobox + sliders + seed) and preview.
function i_select_atom_load(iAtom)
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    if (iAtom < 1) || (iAtom > numel(st.T.Groups)), return; end
    G  = st.T.Groups(iAtom);
    op = 'Laplace-Beltrami';  if isfield(G,'Operator') && ~isempty(G.Operator), op = G.Operator; end
    ax = i_atom_axes(st, op);
    ik = find(strcmp(ctrl.atomKeys, G.KernelName), 1);
    if ~isempty(ik), ctrl.jKernel.setSelectedIndex(ik - 1); end
    if ~isempty(ax), b = i_atom_bounds(ax); else, b = i_atom_default_bounds(); end
    panel_eigenfilter_design('BuildAtomSliders', ctrl.jAtomParams, G.KernelName, b, @()bst_call(@OnParamSettle));
    if isstruct(G.KernelParams) && isfield(G.KernelParams, 'vals')
        panel_eigenfilter_design('SetAtomVals', ctrl.jAtomParams, G.KernelParams.vals);
    end
    i_select_op_radio(op);                                         % check the matching operator radio
    st.atomSeed = G.vertices;  setappdata(0, 'DynamicsTarget', st);
    ctrl.jAtomInfo.setText(i_atom_detail(G));
    i_atom_preview();
end

% Persist the edited params back to the selected atom (+ refresh the readout).
function i_atom_writeback()
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    ia = i_field(st, 'curAtom', 0);  if (ia < 1) || (ia > numel(st.T.Groups)), return; end
    ax = i_atom_axes(st, i_atom_op(st));  if isempty(ax), return; end
    k = i_atom_current_kernel(ctrl);  vals = panel_eigenfilter_design('ReadAtomVals', ctrl.jAtomParams);
    lmax = max(ax.Lambda{1}(:));
    kp = bst_eigfilter_controls('ToKernel', k, vals, lmax);  kp.vals = vals;
    st.T.Groups(ia).KernelName = k;  st.T.Groups(ia).KernelParams = kp;
    setappdata(0, 'DynamicsTarget', st);
    ctrl.jAtomInfo.setText(i_atom_detail(st.T.Groups(ia)));
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
    i_select_op_radio(variant);
    i_select_atom_load(ia);                                        % reload sliders/bounds on the new basis + preview
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
