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
    jMenuSort = gui_component('Menu', jMenuAtoms, [], 'Sort groups', IconLoader.ICON_EVT_TYPE, [], []);
    gui_component('MenuItem', jMenuSort, [], 'By name', IconLoader.ICON_EVT_TYPE, [], @(h,e)bst_call(@()AtomSort('name')));
    gui_component('MenuItem', jMenuSort, [], 'By time', IconLoader.ICON_EVT_TYPE, [], @(h,e)bst_call(@()AtomSort('time')));

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

    jPanelNew.add(jPanelAtoms, BorderLayout.CENTER);
    bstPanelNew = BstPanel(panelName, jPanelNew, ...
        struct('jTree',jTree, 'jListOccur',jListOccur, 'jMenuFile',jMenuFile, 'jMenuAtoms',jMenuAtoms));
end


%% ===== SET TARGET (called by view_dynamics) =====
function SetTarget(hFig, T) %#ok<DEFNU>
    file = '';
    if ~isempty(hFig) && ishandle(hFig), file = getappdata(hFig, 'DynamicsFile'); end
    setappdata(0, 'DynamicsTarget', struct('hFig',hFig, 'T',T, 'file',file, 'curGroup',0, ...
        'nodeList',{ {} }, 'nodeInfo',[], 'occMap',[]));
    BuildTree();
end


%% ===== BUILD THE BAND-STACK TREE =====
% Top-level extended group = a STACK that expands to its time-window leaves.
% Top-level simple group = a stack that expands to its single-atom leaves.
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
        nWin = size(G.times, 2);
        stackNode = DefaultMutableTreeNode(sprintf('%s  (%d)', G.label, nWin));
        root.add(stackNode);
        nodeList{end+1} = stackNode;  nodeInfo(end+1) = struct('kind','stack','g',g,'w',0); %#ok<AGROW>
        isWin = (size(G.times,1) == 2);
        for w = 1:nWin
            if isWin
                lab  = sprintf(' %.3f - %.3f s', G.times(1,w), G.times(2,w));
                kind = 'window';
            elseif (w <= numel(G.vertices))
                lab  = sprintf(' %.3fs  v%d', G.times(1,w), G.vertices(w));
                kind = 'atom';
            else
                lab  = sprintf(' %.3fs', G.times(1,w));  kind = 'atom';
            end
            leaf = DefaultMutableTreeNode(lab);
            stackNode.add(leaf);
            nodeList{end+1} = leaf;  nodeInfo(end+1) = struct('kind',kind,'g',g,'w',w); %#ok<AGROW>
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
        if strcmp(info.kind, 'window')
            [rows, occMap] = i_window_atoms(st.T, info.g, info.w);
            for k = 1:numel(rows), model.addElement(rows{k}); end
            i_jump(st.T.Groups(info.g).times(1, info.w));   % selecting a window jumps to its onset
        elseif strcmp(info.kind, 'atom')
            % single atom of a simple band group: list just it (and it is highlightable)
            G = st.T.Groups(info.g);
            if (info.w <= numel(G.vertices))
                model.addElement(sprintf(' %.3fs  %-8s  v%d', G.times(1,info.w), i_str(G.phase), G.vertices(info.w)));
                occMap(end+1,:) = [info.g, info.w, 0]; %#ok<AGROW>
            end
            i_jump(G.times(1, info.w));                     % and to the atom's time
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
function [rows, occMap] = i_window_atoms(T, gBand, w)
    rows = {};  occMap = zeros(0,3);
    G = T.Groups(gBand);
    on = G.times(1,w);  off = G.times(2,w);
    children = find(strcmpi({T.Groups.parent}, G.label));
    times = [];  phases = {};  verts = [];  cc = [];  oo = [];
    for c = children(:)'
        Gc = T.Groups(c);
        for o = 1:numel(Gc.vertices)
            t = Gc.times(1,o);
            if (t >= on - 1e-9) && (t <= off + 1e-9)
                times(end+1)  = t;            %#ok<AGROW>
                phases{end+1} = i_str(Gc.phase); %#ok<AGROW>
                verts(end+1)  = Gc.vertices(o); %#ok<AGROW>
                cc(end+1) = c;  oo(end+1) = o; %#ok<AGROW>
            end
        end
    end
    [~, ord] = sort(times);
    for k = ord
        rows{end+1} = sprintf(' %8.3fs  %-8s  v%d', times(k), phases{k}, verts(k)); %#ok<AGROW>
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
        try, view_dynamics('Redraw', st.hFig, st.T); catch, end %#ok<CTCH>
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
