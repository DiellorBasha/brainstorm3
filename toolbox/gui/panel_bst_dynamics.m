function varargout = panel_bst_dynamics( varargin )
% PANEL_BST_DYNAMICS: Record-style panel for the spatiotemporal atom system (bst_dynamics).
%
% The atom-table component (increment 1): an Events-style panel for a dynamics table
% (db_template('dynamicsmat')) -- a File menu (open/save the dynamics_* table), an Atoms menu
% (add/rename/delete/color/sort groups), and a split pane with the group TREE on the left
% (nested: an extended "alpha (8-13 Hz)" window over simple "alpha_peak"/... children) and the
% selected group's OCCURRENCE list on the right. Selecting an occurrence highlights its marker
% on the cortex and jumps the recording time (like the Record panel's JumpToEvent). Opened by
% view_dynamics. The temporal / spatial / frequency / eigenmode axes fold in here in later
% increments.
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

    jPanelNew = gui_component('Panel');
    jPanelTop = gui_component('Panel');  jPanelNew.add(jPanelTop, BorderLayout.NORTH);

    % ===== MENU BAR: File, Atoms =====
    jMenuBar = gui_component('MenuBar', jPanelTop, BorderLayout.NORTH);
    jMenuFile = gui_component('Menu', jMenuBar, [], 'File', [], [], [], 11);
    gui_component('MenuItem', jMenuFile, [], 'Open dynamics table...', [], [], @(h,e)bst_call(@FileOpen));
    gui_component('MenuItem', jMenuFile, [], 'Save',                   [], [], @(h,e)bst_call(@FileSave));
    gui_component('MenuItem', jMenuFile, [], 'Save as...',             [], [], @(h,e)bst_call(@FileSaveAs));
    jMenuAtoms = gui_component('Menu', jMenuBar, [], 'Atoms', [], [], [], 11);
    gui_component('MenuItem', jMenuAtoms, [], 'Add group',    [], [], @(h,e)bst_call(@AtomAddGroup));
    gui_component('MenuItem', jMenuAtoms, [], 'Rename group', [], [], @(h,e)bst_call(@AtomRenameGroup));
    gui_component('MenuItem', jMenuAtoms, [], 'Delete group', [], [], @(h,e)bst_call(@AtomDeleteGroup));
    gui_component('MenuItem', jMenuAtoms, [], 'Set color...', [], [], @(h,e)bst_call(@AtomSetColor));
    jMenuAtoms.addSeparator();
    gui_component('MenuItem', jMenuAtoms, [], 'Sort by name', [], [], @(h,e)bst_call(@()AtomSort('name')));
    gui_component('MenuItem', jMenuAtoms, [], 'Sort by time', [], [], @(h,e)bst_call(@()AtomSort('time')));
    % Small toolbar (text buttons)
    jToolbar = gui_component('Toolbar', jMenuBar);
    gui_component('ToolbarButton', jToolbar, [], 'Open',    [], 'Open a dynamics table', @(h,e)bst_call(@FileOpen));
    gui_component('ToolbarButton', jToolbar, [], 'Save',    [], 'Save the dynamics table', @(h,e)bst_call(@FileSave));
    gui_component('ToolbarButton', jToolbar, [], '+ Group', [], 'Add an atom group', @(h,e)bst_call(@AtomAddGroup));

    % ===== SPLIT PANE: group tree | occurrence list =====
    jTree = java_create('javax.swing.JTree');
    jTree.setRootVisible(0);  jTree.setShowsRootHandles(1);
    jTree.getSelectionModel().setSelectionMode(javax.swing.tree.TreeSelectionModel.SINGLE_TREE_SELECTION);
    jTree.setFont(Font('Monospaced', Font.PLAIN, java_scaled('value', 11)));
    java_setcb(jTree, 'ValueChangedCallback', @(h,e)TreeSel_Callback());
    jScrollTree = JScrollPane(jTree);

    jListOccur = java_create('javax.swing.JList');
    jListOccur.setFont(Font('Monospaced', Font.PLAIN, java_scaled('value', 11)));
    java_setcb(jListOccur, 'ValueChangedCallback', @(h,e)OccurSel_Callback());
    jScrollOccur = JScrollPane(jListOccur);

    jSplit = JSplitPane(JSplitPane.HORIZONTAL_SPLIT, jScrollTree, jScrollOccur);
    jSplit.setResizeWeight(0.58);  jSplit.setDividerSize(java_scaled('value', 4));
    jSplit.setPreferredSize(java_scaled('dimension', 380, 440));
    jPanelNew.add(jSplit, BorderLayout.CENTER);

    bstPanelNew = BstPanel(panelName, jPanelNew, ...
        struct('jTree',jTree, 'jListOccur',jListOccur, 'jMenuFile',jMenuFile, 'jMenuAtoms',jMenuAtoms));
end


%% ===== SET TARGET (called by view_dynamics) =====
function SetTarget(hFig, T) %#ok<DEFNU>
    file = '';
    if ~isempty(hFig) && ishandle(hFig), file = getappdata(hFig, 'DynamicsFile'); end
    setappdata(0, 'DynamicsTarget', struct('hFig',hFig, 'T',T, 'file',file, 'curGroup',0, ...
        'nodeList',{ {} }, 'goList',[]));
    BuildTree();
end


%% ===== BUILD THE GROUP TREE (groups only, nested by parent) =====
function BuildTree()
    import javax.swing.tree.*;
    ctrl = bst_get('PanelControls', 'Dynamics');
    st   = getappdata(0, 'DynamicsTarget');
    if isempty(ctrl) || isempty(st), return; end
    T = st.T;
    root = DefaultMutableTreeNode('Atoms');
    nodeList = {};  goList = [];
    added = false(1, numel(T.Groups));
    parents = {T.Groups.parent};
    isTop = cellfun(@isempty, parents);
    for g = [find(isTop), find(~isTop)]
        i_addGroup(root, g);
    end
    ctrl.jTree.setModel(DefaultTreeModel(root));
    if root.getChildCount() > 0, ctrl.jTree.expandRow(0); end
    st.nodeList = nodeList;  st.goList = goList;
    setappdata(0, 'DynamicsTarget', st);
    % clear the occurrence list
    ctrl.jListOccur.setModel(javax.swing.DefaultListModel());

    function i_addGroup(parentNode, g)
        if added(g), return; end
        added(g) = true;
        G = T.Groups(g);
        nOcc = max(size(G.times,2), numel(G.vertices));
        node = DefaultMutableTreeNode(sprintf('%s   [%s, %d]', G.label, G.type, nOcc));
        parentNode.add(node);
        nodeList{end+1} = node;  goList(end+1,:) = g; %#ok<AGROW>
        for c = find(strcmpi(parents, G.label))
            i_addGroup(node, c);
        end
    end
end


%% ===== TREE SELECTION -> populate the occurrence list =====
function TreeSel_Callback()
    ctrl = bst_get('PanelControls', 'Dynamics');
    st   = getappdata(0, 'DynamicsTarget');
    if isempty(ctrl) || isempty(st), return; end
    sel = ctrl.jTree.getLastSelectedPathComponent();
    g = 0;
    if ~isempty(sel)
        for i = 1:numel(st.nodeList)
            if sel.equals(st.nodeList{i}); g = st.goList(i); break; end
        end
    end
    st.curGroup = g;  setappdata(0, 'DynamicsTarget', st);
    % Populate the occurrence list for this group
    model = javax.swing.DefaultListModel();
    if g >= 1
        G = st.T.Groups(g);
        nOcc = size(G.times, 2);
        for o = 1:nOcc
            if size(G.times,1) == 2
                model.addElement(sprintf('%8.3f - %8.3f s', G.times(1,o), G.times(2,o)));
            elseif (o <= numel(G.vertices))
                model.addElement(sprintf('%9.3fs  v%-6d %s', G.times(1,o), G.vertices(o), i_hemi(G.hemi(o))));
            else
                model.addElement(sprintf('%9.3fs', G.times(1,o)));
            end
        end
    end
    ctrl.jListOccur.setModel(model);
    % clear marker selection until an occurrence is picked
    hSel = findobj(st.hFig, 'Tag', 'AtomSel');  if ~isempty(hSel), set(hSel, 'Visible', 'off'); end
end


%% ===== OCCURRENCE SELECTION -> highlight marker + jump time =====
function OccurSel_Callback()
    ctrl = bst_get('PanelControls', 'Dynamics');
    st   = getappdata(0, 'DynamicsTarget');
    if isempty(ctrl) || isempty(st) || (st.curGroup < 1) || ~ishandle(st.hFig), return; end
    o = ctrl.jListOccur.getSelectedIndex() + 1;   % Java 0-based
    if (o < 1), return; end
    g = st.curGroup;  G = st.T.Groups(g);
    % marker (spatial occurrences only)
    hSel = findobj(st.hFig, 'Tag', 'AtomSel');
    GroupsPosOff = getappdata(st.hFig, 'GroupsPosOff');
    if ~isempty(hSel) && ~isempty(GroupsPosOff) && (g <= numel(GroupsPosOff)) && (o <= size(GroupsPosOff{g},1))
        p = GroupsPosOff{g}(o,:);
        set(hSel, 'XData',p(1), 'YData',p(2), 'ZData',p(3), 'Visible','on');
    elseif ~isempty(hSel)
        set(hSel, 'Visible', 'off');
    end
    % time cursor
    if (o <= size(G.times,2))
        try, panel_time('SetCurrentTime', G.times(1,o)); catch, end %#ok<CTCH>
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
    if isempty(st.file)
        FileSaveAs();  return;
    end
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


%% ===== ATOMS menu (edit the table, then refresh) =====
function AtomAddGroup()
    st = getappdata(0, 'DynamicsTarget');  if isempty(st), return; end
    name = java_dialog('input', 'Name for the new atom group:', 'Add group', [], '');
    if isempty(name), return; end
    G = bst_dynamics('NewGroup', name);
    G.color = [0.6 0.6 0.6];
    if ~isempty(st.T.Groups), G.SurfaceFile = st.T.Groups(1).SurfaceFile; G.DataFile = st.T.Groups(1).DataFile; end
    st.T = bst_dynamics('AddGroup', st.T, G);
    i_apply(st);
end
function AtomRenameGroup()
    [st, g] = i_selected();  if g < 1, return; end
    old = st.T.Groups(g).label;
    name = java_dialog('input', 'New name:', 'Rename group', [], old);
    if isempty(name) || strcmp(name, old), return; end
    % update children that reference the old parent label
    for c = 1:numel(st.T.Groups)
        if strcmp(st.T.Groups(c).parent, old), st.T.Groups(c).parent = name; end
    end
    st.T.Groups(g).label = name;
    i_apply(st);
end
function AtomDeleteGroup()
    [st, g] = i_selected();  if g < 1, return; end
    if ~java_dialog('confirm', sprintf('Delete group "%s" (and its occurrences)?', st.T.Groups(g).label), 'Delete group'), return; end
    lbl = st.T.Groups(g).label;
    % orphan any children (promote to top-level)
    for c = 1:numel(st.T.Groups)
        if strcmp(st.T.Groups(c).parent, lbl), st.T.Groups(c).parent = ''; end
    end
    st.T.Groups(g) = [];  st.T.nGroups = numel(st.T.Groups);
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
    if strcmpi(mode, 'time')
        key = cellfun(@(t) i_firsttime(t), {st.T.Groups.times});
    else
        key = lower({st.T.Groups.label});
    end
    [~, ord] = sort(key);
    st.T.Groups = st.T.Groups(ord);
    i_apply(st);
end


%% ===== helpers =====
function [st, g] = i_selected()
    g = 0;  st = getappdata(0, 'DynamicsTarget');
    if isempty(st), return; end
    g = st.curGroup;
    if g < 1, java_dialog('warning', 'Select a group in the tree first.', 'Atoms'); end
end
function i_apply(st)
    % store edited table, redraw cortex markers, rebuild tree
    setappdata(0, 'DynamicsTarget', st);
    if ~isempty(st.hFig) && ishandle(st.hFig)
        try, view_dynamics('Redraw', st.hFig, st.T); catch, end %#ok<CTCH>
    end
    BuildTree();
end
function t0 = i_firsttime(times)
    if isempty(times), t0 = inf; else, t0 = times(1,1); end
end
function s = i_hemi(h)
    if isempty(h) || (h < 1) || (h > 2), s = '?'; else, hc = 'LR'; s = hc(double(h)); end
end
