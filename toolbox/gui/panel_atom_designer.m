function varargout = panel_atom_designer(varargin)
% PANEL_ATOM_DESIGNER: docked controls panel for the atom designer (view_atom_designer).
%
% Hosts the Operator + Filter (kernel) selectors, the contextual per-kernel parameter sliders
% (panel_eigenfilter_design atom sliders, driven by bst_eigfilter_controls), a Connectome toggle,
% Save, and a status line. The designer (view_atom_designer) keeps the realise/overlay/save logic and
% passes its callback handles here via 'Configure'; this panel only HOSTS the controls + reads them
% back. The Swing callbacks relay to the designer handles (appdata 'AtomDesignerCB').
%
% API (macro_method):
%   bstPanel = panel_atom_designer('CreatePanel')
%   panel_atom_designer('Configure', cb, bounds, kernel0, variant0)   % cb = struct(Kernel/Operator/Param/Direction/Fibers/Save handles)
%   vals = panel_atom_designer('ReadVals')                            % [s1 s2 s3]
%   k    = panel_atom_designer('CurrentKernel')
%   v    = panel_atom_designer('CurrentOperator')                     % 'Laplace-Beltrami'|'LB-Connectome'|'Dirac'|'Dirac-Connectome'
%   d    = panel_atom_designer('CurrentDirection')                    % 'Normal'|'+X'|'+Y'|'+Z'|'Pick-on-surface'
%   panel_atom_designer('ShowDirection', tf)                          % show/hide the seed-direction control
%   panel_atom_designer('RebuildSliders', kernel, bounds)
%   panel_atom_designer('SetStatus', text)
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


function bstPanelNew = CreatePanel() %#ok<DEFNU>
    panelName = 'AtomDesigner';
    import java.awt.*;
    import javax.swing.*;
    [atomKeys, atomDisp] = panel_eigenfilter_design('AtomKernelsGrouped');
    BUTTON_WIDTH   = java_scaled('value', 56);
    DEFAULT_HEIGHT = java_scaled('value', 22);

    % options stack: one bordered sub-panel per group, pinned to the top (panel_surface convention)
    jPanelOptions = gui_component('Panel');
    jPanelOptions.setLayout(BoxLayout(jPanelOptions, BoxLayout.Y_AXIS));
    jPanelOptions.setBorder(BorderFactory.createEmptyBorder(7,7,7,7));

    % ===== FILTER DESIGN: operator + filter + contextual parameter sliders =====
    jDes = gui_river([1,1], [1,8,1,4], 'Filter design');
        % Operator: combobox over the SP1-routable operators (Connection Laplacian is out — unpersisted frame)
        gui_component('label', jDes, 'br', 'Operator:');
        opItems = {'Geometric','Connectomic','Dirac','Dirac (connectome)'};
        jOperator = gui_component('combobox', jDes, 'tab hfill', [], {opItems}, [], [], []);
        java_setcb(jOperator, 'ActionPerformedCallback', @(h,e) bst_call(@OnOperatorCb));
        % Filter: grouped combobox (Dynamic / Static headers; headers are not selectable)
        gui_component('label', jDes, 'br', 'Filter:');
        jKernel = gui_component('combobox', jDes, 'tab hfill', [], {atomDisp}, [], [], []);
        java_setcb(jKernel, 'ActionPerformedCallback', @(h,e) bst_call(@OnKernelCb));
        jParams = gui_river([1,1], [2,0,0,0]);                 % contextual sliders (own rows; removeAll-safe)
        jDes.add('br hfill', jParams);
        % Seed-direction control (quaternion/Dirac fibers): preset combobox, hidden for scalar operators.
        jDirRow = gui_river([0 0], [0 2 0 2]);
        gui_component('label', jDirRow, [], 'Direction: ');
        jDirCombo = gui_component('combobox', jDirRow, 'hfill', [], {{'Normal','+X','+Y','+Z','Pick-on-surface'}}, [], [], []);
        java_setcb(jDirCombo, 'ActionPerformedCallback', @(h,e) bst_call(@OnDirCb));
        jDirRow.setVisible(false);                               % no vector operator selected yet -> hidden
        jDes.add('br hfill', jDirRow);
    jPanelOptions.add(jDes);

    % ===== ACTIONS: connectome overlay + save (explicitly sized -> side by side) =====
    jAct = gui_river([1,1], [1,8,1,4]);
        jFibers = gui_component('toggle', jAct, 'br', 'Connectome', [], 'Overlay the connectome fibers, coloured by the atom', @(h,e) bst_call(@OnConnectomeCb));
        jFibers.setPreferredSize(Dimension(2*BUTTON_WIDTH, DEFAULT_HEIGHT));
        jSave = gui_component('button', jAct, [], 'Save', [], 'Save atom -> Scout + Event', @(h,e) bst_call(@OnSaveCb)); %#ok<NASGU>
        jSave.setPreferredSize(Dimension(BUTTON_WIDTH, DEFAULT_HEIGHT));
    jPanelOptions.add(jAct);

    % ===== STATUS (contained, line-wrapping text area) =====
    jStat = gui_river([1,1], [1,8,1,4], 'Status');
        jStatus = JTextArea(2, 1);
        jStatus.setEditable(0);  jStatus.setLineWrap(1);  jStatus.setWrapStyleWord(1);
        jStatus.setFont(Font('SansSerif', Font.PLAIN, java_scaled('value', 11)));
        jStatus.setBorder(BorderFactory.createEmptyBorder());
        jStat.add('hfill', jStatus);
    jPanelOptions.add(jStat);

    jPanelNew = gui_component('Panel');
    jPanelNew.add(jPanelOptions, BorderLayout.NORTH);
    ctrl = struct('jOperator',jOperator, 'jDirRow',jDirRow, 'jDirCombo',jDirCombo, ...
                  'jKernel',jKernel, 'jParams',jParams, ...
                  'jFibers',jFibers, 'jStatus',jStatus, 'atomKeys',{atomKeys});
    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);
end


%% ===== configure: wire the designer callbacks + populate the controls =====
function Configure(cb, bounds, kernel0, variant0) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'AtomDesigner');  if isempty(ctrl), return; end
    setappdata(0, 'AtomDesignerCB', cb);
    ctrl.jOperator.setSelectedItem(i_item_for_variant(variant0));
    iK = find(strcmp(ctrl.atomKeys, kernel0), 1);  if isempty(iK), iK = i_first_kernel_idx(ctrl.atomKeys); end
    ctrl.jKernel.setSelectedIndex(iK - 1);
    setappdata(0, 'AtomDesignerKIdx', iK);
    panel_eigenfilter_design('BuildAtomSliders', ctrl.jParams, kernel0, bounds, @()bst_call(@OnParamSettle));
end

%% ===== accessors =====
function vals = ReadVals() %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'AtomDesigner');
    if isempty(ctrl), vals = [0 0 0]; return; end
    vals = panel_eigenfilter_design('ReadAtomVals', ctrl.jParams);
end

function k = CurrentKernel() %#ok<DEFNU>
    k = '';
    ctrl = bst_get('PanelControls', 'AtomDesigner');  if isempty(ctrl), return; end
    idx = double(ctrl.jKernel.getSelectedIndex()) + 1;
    if (idx >= 1) && (idx <= numel(ctrl.atomKeys)) && ~isempty(ctrl.atomKeys{idx})
        k = ctrl.atomKeys{idx};
    else                                                                     % header row selected -> last valid
        last = getappdata(0, 'AtomDesignerKIdx');
        if ~isempty(last) && (last <= numel(ctrl.atomKeys)), k = ctrl.atomKeys{last}; end
    end
end

function v = CurrentOperator() %#ok<DEFNU>
    v = 'Laplace-Beltrami';
    ctrl = bst_get('PanelControls', 'AtomDesigner');  if isempty(ctrl), return; end
    v = i_variant_for_item(char(ctrl.jOperator.getSelectedItem()));
end

% Pure display->variant map (also used by Configure and the headless test).
function v = i_variant_for_item(name) %#ok<DEFNU>
    switch name
        case 'Connectomic',        v = 'LB-Connectome';
        case 'Dirac',              v = 'Dirac';
        case 'Dirac (connectome)', v = 'Dirac-Connectome';
        otherwise,                 v = 'Laplace-Beltrami';   % 'Geometric' / unknown
    end
end

% Pure variant->display map.
function nm = i_item_for_variant(variant) %#ok<DEFNU>
    switch variant
        case 'LB-Connectome',    nm = 'Connectomic';
        case 'Dirac',            nm = 'Dirac';
        case 'Dirac-Connectome', nm = 'Dirac (connectome)';
        otherwise,               nm = 'Geometric';
    end
end

% Selected seed-direction preset.
function v = CurrentDirection() %#ok<DEFNU>
    v = 'Normal';
    ctrl = bst_get('PanelControls', 'AtomDesigner');  if isempty(ctrl), return; end
    v = char(ctrl.jDirCombo.getSelectedItem());
end

% Show/hide the seed-direction control (vector fiber -> show; scalar -> hide).
function ShowDirection(tf) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'AtomDesigner');  if isempty(ctrl), return; end
    ctrl.jDirRow.setVisible(logical(tf));
end

% Relay the direction-combo change to the designer.
function OnDirCb()
    cb = getappdata(0, 'AtomDesignerCB');
    if ~isempty(cb) && isfield(cb,'Direction'), cb.Direction(); end
end

function ix = i_first_kernel_idx(atomKeys)
    ix = find(~cellfun('isempty', atomKeys), 1);  if isempty(ix), ix = 1; end
end

function RebuildSliders(kernel, bounds) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'AtomDesigner');  if isempty(ctrl), return; end
    panel_eigenfilter_design('BuildAtomSliders', ctrl.jParams, kernel, bounds, @()bst_call(@OnParamSettle));
end

function SetStatus(txt) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'AtomDesigner');  if isempty(ctrl), return; end
    ctrl.jStatus.setText(txt);
end

%% ===== internal Swing callbacks: relay to the designer handles =====
function OnOperatorCb()
    cb = getappdata(0, 'AtomDesignerCB');
    if ~isempty(cb) && isfield(cb,'Operator'), cb.Operator(); end
end
function OnKernelCb()
    ctrl = bst_get('PanelControls', 'AtomDesigner');  if isempty(ctrl), return; end
    idx = double(ctrl.jKernel.getSelectedIndex()) + 1;
    if (idx < 1) || (idx > numel(ctrl.atomKeys)) || isempty(ctrl.atomKeys{idx})   % header row -> revert, don't fire
        last = getappdata(0, 'AtomDesignerKIdx');
        if isempty(last), last = i_first_kernel_idx(ctrl.atomKeys); end
        ctrl.jKernel.setSelectedIndex(last - 1);
        return;
    end
    setappdata(0, 'AtomDesignerKIdx', idx);
    cb = getappdata(0, 'AtomDesignerCB');
    if ~isempty(cb) && isfield(cb,'Kernel'), cb.Kernel(); end
end
function OnParamSettle()
    cb = getappdata(0, 'AtomDesignerCB');
    if ~isempty(cb) && isfield(cb,'Param'), cb.Param(); end
end
function OnConnectomeCb()
    ctrl = bst_get('PanelControls', 'AtomDesigner');  if isempty(ctrl), return; end
    cb = getappdata(0, 'AtomDesignerCB');
    if ~isempty(cb) && isfield(cb,'Fibers'), cb.Fibers(double(ctrl.jFibers.isSelected())); end
end
function OnSaveCb()
    cb = getappdata(0, 'AtomDesignerCB');
    if ~isempty(cb) && isfield(cb,'Save'), cb.Save(); end
end
