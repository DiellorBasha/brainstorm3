function varargout = panel_dynamics( varargin )
% PANEL_DYNAMICS: List/inspect spatiotemporal "atoms" (the joint Events+Scouts system).
%
% Milestone-1 skeleton panel. Lists the atoms of a loaded dynamics table
% (db_template('dynamicsmat')); selecting a row highlights the atom's marker on
% the cortex figure and jumps the recording time to the atom's time. Opened by
% view_dynamics(). The full interactive editor (add/remove/filter/stats) comes
% in later phases.
%
% USAGE:  bstPanel = panel_dynamics('CreatePanel')
%                    panel_dynamics('SetTarget', hFig, Atoms)   % populate from view_dynamics
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
    jPanelMain.setPreferredSize(java_scaled('dimension', 260, 420));

    % Header
    jLabel = JLabel('Atoms');
    jPanelMain.add(jLabel, BorderLayout.NORTH);

    % Atom list
    jList = java_create('javax.swing.JList');
    jList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    jList.setFont(Font('Monospaced', Font.PLAIN, java_scaled('value', 11)));
    java_setcb(jList, ...
        'ValueChangedCallback', @(h,ev)ListChanged_Callback(ev), ...
        'MouseClickedCallback', @(h,ev)ListChanged_Callback([]));
    jScroll = JScrollPane(jList);
    jPanelMain.add(jScroll, BorderLayout.CENTER);

    bstPanelNew = BstPanel(panelName, jPanelMain, struct('jList', jList, 'jLabel', jLabel));
end


%% ===== LIST SELECTION CALLBACK =====
function ListChanged_Callback(ev)
    if ~isempty(ev) && ev.getValueIsAdjusting()
        return;
    end
    ctrl = bst_get('PanelControls', 'Dynamics');
    if isempty(ctrl), return; end
    iAtom = ctrl.jList.getSelectedIndex() + 1;   % Java 0-based
    if (iAtom < 1), return; end
    HighlightAtom(iAtom);
end


%% ===== SET TARGET (called by view_dynamics) =====
% Stores the figure + atoms this panel drives, and populates the list.
function SetTarget(hFig, Atoms) %#ok<DEFNU>
    setappdata(0, 'DynamicsTarget', struct('hFig', hFig, 'Atoms', Atoms));
    ctrl = bst_get('PanelControls', 'Dynamics');
    if isempty(ctrl), return; end
    model = javax.swing.DefaultListModel();
    hemiChar = 'LR';
    for i = 1:numel(Atoms)
        a = Atoms(i);
        hc = '?';
        if ~isempty(a.hemi) && (a.hemi>=1) && (a.hemi<=2), hc = hemiChar(double(a.hemi)); end
        model.addElement(sprintf('%8.3fs  %-9s  v%-6d %c', a.time, a.category, a.vertex, hc));
    end
    ctrl.jList.setModel(model);
    ctrl.jLabel.setText(sprintf('Atoms (%d)', numel(Atoms)));
end


%% ===== HIGHLIGHT AN ATOM =====
function HighlightAtom(iAtom)
    st = getappdata(0, 'DynamicsTarget');
    if isempty(st) || ~ishandle(st.hFig) || (iAtom > numel(st.Atoms)), return; end
    a = st.Atoms(iAtom);
    % Move the selection marker on the cortex
    hSel   = findobj(st.hFig, 'Tag', 'AtomSel');
    posOff = getappdata(st.hFig, 'AtomOffsetPos');
    if ~isempty(hSel) && ~isempty(posOff) && (iAtom <= size(posOff,1))
        set(hSel, 'XData', posOff(iAtom,1), 'YData', posOff(iAtom,2), 'ZData', posOff(iAtom,3), 'Visible', 'on');
    end
    % Jump recording time (no-op if no time context is open)
    try, panel_time('SetCurrentTime', a.time); catch, end %#ok<CTCH>
end
