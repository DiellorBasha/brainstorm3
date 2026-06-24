function varargout = panel_dynamics( varargin )
% PANEL_DYNAMICS: Tree of spatiotemporal atom groups (the joint Events+Scouts system).
%
% Milestone-1 skeleton panel. Shows a tree of atom groups (db_template('atomgroup')):
% top-level groups, their nested child groups (an extended "alpha (8-13 Hz)" window
% containing simple "alpha_peak"/"alpha_trough"/... children), and each group's
% occurrences as leaves. Selecting an occurrence highlights its marker on the
% cortex figure and jumps the recording time. Opened by view_dynamics().
%
% USAGE:  bstPanel = panel_dynamics('CreatePanel')
%                    panel_dynamics('SetTarget', hFig, T)   % populate from view_dynamics
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

    jPanelMain = java_create('javax.swing.JPanel');
    jPanelMain.setLayout(BorderLayout());
    jPanelMain.setPreferredSize(java_scaled('dimension', 300, 460));

    jTree = java_create('javax.swing.JTree');
    jTree.setRootVisible(0);
    jTree.setShowsRootHandles(1);
    jTree.getSelectionModel().setSelectionMode(javax.swing.tree.TreeSelectionModel.SINGLE_TREE_SELECTION);
    jTree.setFont(Font('Monospaced', Font.PLAIN, java_scaled('value', 11)));
    java_setcb(jTree, 'ValueChangedCallback', @(h,ev)TreeSel_Callback());
    jScroll = JScrollPane(jTree);
    jPanelMain.add(jScroll, BorderLayout.CENTER);

    bstPanelNew = BstPanel(panelName, jPanelMain, struct('jTree', jTree));
end


%% ===== SET TARGET (called by view_dynamics) =====
% Builds the group tree and stores the figure + node->(group,occurrence) map.
function SetTarget(hFig, T) %#ok<DEFNU>
    import javax.swing.tree.*;
    ctrl = bst_get('PanelControls', 'Dynamics');
    if isempty(ctrl), return; end

    root     = DefaultMutableTreeNode('Atoms');
    nodeList = {};
    goList   = zeros(0, 2);
    added    = false(1, numel(T.Groups));
    parents  = {T.Groups.parent};

    % Top-level groups first, then any orphans (parent set but not found)
    isTop = cellfun(@isempty, parents);
    for g = [find(isTop), find(~isTop)]
        i_addGroup(root, g);
    end
    ctrl.jTree.setModel(DefaultTreeModel(root));
    ctrl.jTree.expandRow(0);
    setappdata(0, 'DynamicsTarget', struct('hFig', hFig, 'T', T, 'nodeList', {nodeList}, 'goList', goList));

    % --- nested builder (captures nodeList/goList/added/T) ---
    function i_addGroup(parentNode, g)
        if added(g), return; end
        added(g) = true;
        G = T.Groups(g);
        nOcc  = max(size(G.times,2), numel(G.vertices));
        gNode = DefaultMutableTreeNode(sprintf('%s   [%s, %d]', G.label, G.type, nOcc));
        parentNode.add(gNode);
        nodeList{end+1} = gNode;  goList(end+1,:) = [g, 0]; %#ok<AGROW>
        % occurrence leaves (spatial occurrences only)
        for o = 1:numel(G.vertices)
            leaf = DefaultMutableTreeNode(sprintf('%9.3fs  v%-6d %s', G.times(1,o), G.vertices(o), i_hemi(G.hemi(o))));
            gNode.add(leaf);
            nodeList{end+1} = leaf;  goList(end+1,:) = [g, o]; %#ok<AGROW>
        end
        % nested child groups
        for c = find(strcmpi(parents, G.label))
            i_addGroup(gNode, c);
        end
    end
end


%% ===== TREE SELECTION CALLBACK =====
function TreeSel_Callback()
    ctrl = bst_get('PanelControls', 'Dynamics');
    st   = getappdata(0, 'DynamicsTarget');
    if isempty(ctrl) || isempty(st) || ~ishandle(st.hFig), return; end
    sel = ctrl.jTree.getLastSelectedPathComponent();
    if isempty(sel), return; end
    % Map the selected Java node back to (group, occurrence)
    g = 0;  o = 0;
    for i = 1:numel(st.nodeList)
        if sel.equals(st.nodeList{i})
            g = st.goList(i,1);  o = st.goList(i,2);  break;
        end
    end
    if (g < 1), return; end
    G = st.T.Groups(g);
    % Highlight the marker (only for spatial occurrences)
    hSel = findobj(st.hFig, 'Tag', 'AtomSel');
    GroupsPosOff = getappdata(st.hFig, 'GroupsPosOff');
    if (o >= 1) && ~isempty(hSel) && (g <= numel(GroupsPosOff)) && (o <= size(GroupsPosOff{g},1))
        p = GroupsPosOff{g}(o,:);
        set(hSel, 'XData', p(1), 'YData', p(2), 'ZData', p(3), 'Visible', 'on');
    elseif ~isempty(hSel)
        set(hSel, 'Visible', 'off');   % group node selected: clear the point marker
    end
    % Jump recording time (occurrence time, or the group's first onset)
    if (o >= 1) && (o <= size(G.times,2))
        t = G.times(1, o);
    elseif ~isempty(G.times)
        t = G.times(1, 1);
    else
        t = [];
    end
    if ~isempty(t)
        try, panel_time('SetCurrentTime', t); catch, end %#ok<CTCH>
    end
end


%% ===== HEMISPHERE CHAR =====
function s = i_hemi(h)
    if isempty(h) || (h < 1) || (h > 2)
        s = '?';
    else
        hc = 'LR';
        s = hc(double(h));
    end
end
