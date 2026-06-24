function varargout = panel_bst_dynamics( varargin )
% PANEL_BST_DYNAMICS: Record-style panel for the spatiotemporal atom system (bst_dynamics).
%
% The atom-table component (increment 1), built by REUSING the Record panel's Events-section
% UI components so it looks/behaves like Events: a File menu (open/save the dynamics_* table),
% an Atoms menu (add/rename/delete/color/sort groups), and a JSplitPane with the colored group
% list on the left (BstColorListRenderer; child phase groups indented under their window) and
% the selected group's occurrence list on the right (BstStringListRenderer). Selecting an
% occurrence highlights its marker on the cortex and jumps the recording time (like Record's
% JumpToEvent). Docked as a tools tab; opened by view_dynamics. The temporal / spatial /
% frequency / eigenmode axes fold in here in later increments.
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
    fontSize = java_scaled('value', 11);

    jPanelNew = gui_component('Panel');
    jPanelTop = gui_component('Panel');  jPanelNew.add(jPanelTop, BorderLayout.NORTH);

    % ===== MENU BAR: File, Atoms (+ small toolbar) =====
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
    jToolbar = gui_component('Toolbar', jMenuBar);
    gui_component('ToolbarButton', jToolbar, [], 'Open',    [], 'Open a dynamics table', @(h,e)bst_call(@FileOpen));
    gui_component('ToolbarButton', jToolbar, [], 'Save',    [], 'Save the dynamics table', @(h,e)bst_call(@FileSave));
    gui_component('ToolbarButton', jToolbar, [], '+ Group', [], 'Add an atom group', @(h,e)bst_call(@AtomAddGroup));

    % ===== SPLIT PANE: colored group list | occurrence list (Events-section components) =====
    jListGroups = JList();
    jListGroups.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    jListGroups.setCellRenderer(BstColorListRenderer(fontSize));
    java_setcb(jListGroups, 'ValueChangedCallback', @(h,e)GroupSel_Callback());
    jScrollGroups = JScrollPane(jListGroups);
    jScrollGroups.setBorder([]);

    jListOccur = JList();
    jListOccur.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    jListOccur.setCellRenderer(BstStringListRenderer(fontSize));
    java_setcb(jListOccur, 'ValueChangedCallback', @(h,e)OccurSel_Callback());
    jScrollOccur = JScrollPane(jListOccur);
    jScrollOccur.setBorder([]);

    jSplit = JSplitPane(JSplitPane.HORIZONTAL_SPLIT, jScrollGroups, jScrollOccur);
    jSplit.setResizeWeight(0.6);  jSplit.setDividerSize(java_scaled('value', 4));  jSplit.setBorder([]);
    jPanelNew.add(jSplit, BorderLayout.CENTER);

    bstPanelNew = BstPanel(panelName, jPanelNew, ...
        struct('jListGroups',jListGroups, 'jListOccur',jListOccur, 'jMenuFile',jMenuFile, 'jMenuAtoms',jMenuAtoms));
end


%% ===== SET TARGET (called by view_dynamics) =====
function SetTarget(hFig, T) %#ok<DEFNU>
    file = '';
    if ~isempty(hFig) && ishandle(hFig), file = getappdata(hFig, 'DynamicsFile'); end
    setappdata(0, 'DynamicsTarget', struct('hFig',hFig, 'T',T, 'file',file, 'curGroup',0, 'gIdx',[]));
    BuildGroupList();
end


%% ===== BUILD THE COLORED GROUP LIST (children indented under their window) =====
function BuildGroupList()
    import org.brainstorm.list.*;
    ctrl = bst_get('PanelControls', 'Dynamics');
    st   = getappdata(0, 'DynamicsTarget');
    if isempty(ctrl) || isempty(st), return; end
    T = st.T;
    % Depth-first display order (top-level groups, then their nested children)
    parents = {T.Groups.parent};
    order = [];  depth = [];  added = false(1, numel(T.Groups));
    function addG(g, d)
        if added(g), return; end
        added(g) = true;  order(end+1) = g;  depth(end+1) = d; %#ok<AGROW>
        for c = find(strcmpi(parents, T.Groups(g).label)), addG(c, d+1); end
    end
    for g = find(cellfun(@isempty, parents)), addG(g, 0); end
    for g = find(~added), addG(g, 0); end
    % Build the colored list model (disable the callback during the rebuild)
    bak = java_getcb(ctrl.jListGroups, 'ValueChangedCallback');
    java_setcb(ctrl.jListGroups, 'ValueChangedCallback', []);
    model = javax.swing.DefaultListModel();
    for k = 1:numel(order)
        g = order(k);  G = T.Groups(g);
        nOcc   = max(size(G.times,2), numel(G.vertices));
        indent = repmat('   ', 1, depth(k));
        item = BstListItem('', '', sprintf(' %s%s  (x%d)', indent, G.label, nOcc));
        col = G.color;  if isempty(col), col = [0.6 0.6 0.6]; end
        item.setColor(java.awt.Color(col(1), col(2), col(3)));
        model.addElement(item);
    end
    ctrl.jListGroups.setModel(model);
    ctrl.jListGroups.repaint();
    java_setcb(ctrl.jListGroups, 'ValueChangedCallback', bak);
    st.gIdx = order;  setappdata(0, 'DynamicsTarget', st);
    ctrl.jListOccur.setModel(javax.swing.DefaultListModel());
end


%% ===== GROUP SELECTION -> populate the occurrence list =====
function GroupSel_Callback()
    ctrl = bst_get('PanelControls', 'Dynamics');
    st   = getappdata(0, 'DynamicsTarget');
    if isempty(ctrl) || isempty(st), return; end
    row = ctrl.jListGroups.getSelectedIndex() + 1;   % Java 0-based
    g = 0;
    if (row >= 1) && (row <= numel(st.gIdx)), g = st.gIdx(row); end
    st.curGroup = g;  setappdata(0, 'DynamicsTarget', st);
    model = javax.swing.DefaultListModel();
    if g >= 1
        G = st.T.Groups(g);
        for o = 1:size(G.times, 2)
            if size(G.times,1) == 2
                model.addElement(sprintf(' %1.3f - %1.3f s', G.times(1,o), G.times(2,o)));
            elseif (o <= numel(G.vertices))
                model.addElement(sprintf(' %1.3fs   v%d  %s', G.times(1,o), G.vertices(o), i_hemi(G.hemi(o))));
            else
                model.addElement(sprintf(' %1.3fs', G.times(1,o)));
            end
        end
    end
    ctrl.jListOccur.setModel(model);
    hSel = findobj(st.hFig, 'Tag', 'AtomSel');  if ~isempty(hSel), set(hSel, 'Visible', 'off'); end
end


%% ===== OCCURRENCE SELECTION -> highlight marker + jump time =====
function OccurSel_Callback()
    ctrl = bst_get('PanelControls', 'Dynamics');
    st   = getappdata(0, 'DynamicsTarget');
    if isempty(ctrl) || isempty(st) || (st.curGroup < 1) || ~ishandle(st.hFig), return; end
    o = ctrl.jListOccur.getSelectedIndex() + 1;
    if (o < 1), return; end
    g = st.curGroup;  G = st.T.Groups(g);
    hSel = findobj(st.hFig, 'Tag', 'AtomSel');
    GroupsPosOff = getappdata(st.hFig, 'GroupsPosOff');
    if ~isempty(hSel) && ~isempty(GroupsPosOff) && (g <= numel(GroupsPosOff)) && (o <= size(GroupsPosOff{g},1))
        p = GroupsPosOff{g}(o,:);
        set(hSel, 'XData',p(1), 'YData',p(2), 'ZData',p(3), 'Visible','on');
    elseif ~isempty(hSel)
        set(hSel, 'Visible', 'off');
    end
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


%% ===== ATOMS menu (edit the table, then refresh) =====
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
    if ~java_dialog('confirm', sprintf('Delete group "%s" (and its occurrences)?', st.T.Groups(g).label), 'Delete group'), return; end
    lbl = st.T.Groups(g).label;
    for c = 1:numel(st.T.Groups)
        if strcmp(st.T.Groups(c).parent, lbl), st.T.Groups(c).parent = ''; end
    end
    st.T.Groups(g) = [];  st.T.nGroups = numel(st.T.Groups);  st.curGroup = 0;
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
    if g < 1, java_dialog('warning', 'Select a group first.', 'Atoms'); end
end
function i_apply(st)
    setappdata(0, 'DynamicsTarget', st);
    if ~isempty(st.hFig) && ishandle(st.hFig)
        try, view_dynamics('Redraw', st.hFig, st.T); catch, end %#ok<CTCH>
    end
    BuildGroupList();
end
function t0 = i_firsttime(times)
    if isempty(times), t0 = inf; else, t0 = times(1,1); end
end
function s = i_hemi(h)
    if isempty(h) || (h < 1) || (h > 2), s = '?'; else, hc = 'LR'; s = hc(double(h)); end
end
