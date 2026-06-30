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
%   panel_atom_designer('Configure', cb, bounds, kernel0, variant0)   % cb = struct(Kernel/Operator/Param/Fibers/Save handles)
%   vals = panel_atom_designer('ReadVals')                            % [s1 s2 s3]
%   k    = panel_atom_designer('CurrentKernel')
%   v    = panel_atom_designer('CurrentOperator')                     % 'LB-Connectome' | 'Laplace-Beltrami'
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
    [atomKeys, atomDisp] = panel_eigenfilter_design('AtomKernels');

    jPanelNew = gui_river([5,5], [10,15,12,10]);

    % --- Filter design: operator + kernel + contextual params ---
    jDes = gui_river([2,2], [0,10,12,10], 'Filter design');
        gui_component('label', jDes, '', 'Operator: ', [], [], [], []);
        jOperator = gui_component('combobox', jDes, 'tab hfill', [], {{'connectomic','geometric'}}, [], [], []);
        java_setcb(jOperator, 'ActionPerformedCallback', @(hh,ee) bst_call(@OnOperatorCb));
        gui_component('label', jDes, 'br', 'Filter: ', [], [], [], []);
        jKernel = gui_component('combobox', jDes, 'tab hfill', [], {atomDisp}, [], [], []);
        java_setcb(jKernel, 'ActionPerformedCallback', @(hh,ee) bst_call(@OnKernelCb));
        jParams = gui_river([2,2], [0,2,0,2]);
        jDes.add('br hfill', jParams);
    jPanelNew.add('br hfill', jDes);

    % --- actions: connectome overlay + save ---
    jAct = gui_river([2,2], [0,10,8,10]);
        jFibers = gui_component('toggle', jAct, '', 'Connectome', [], 'Overlay the connectome fibers, coloured by the atom', @(hh,ee) bst_call(@OnConnectomeCb));
        gui_component('button', jAct, 'tab right', 'Save', [], 'Save atom -> Scout + Event', @(hh,ee) bst_call(@OnSaveCb));
    jPanelNew.add('br hfill', jAct);

    % --- status line ---
    jStatus = gui_component('label', jPanelNew, 'br hfill', '', [], [], [], []);

    ctrl = struct('jOperator',jOperator, 'jKernel',jKernel, 'jParams',jParams, ...
                  'jFibers',jFibers, 'jStatus',jStatus, 'atomKeys',{atomKeys});
    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);
end


%% ===== configure: wire the designer callbacks + populate the controls =====
function Configure(cb, bounds, kernel0, variant0) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'AtomDesigner');  if isempty(ctrl), return; end
    setappdata(0, 'AtomDesignerCB', cb);
    ctrl.jOperator.setSelectedIndex(strcmpi(variant0, 'Laplace-Beltrami'));   % connectomic=0, geometric=1
    iK = find(strcmp(ctrl.atomKeys, kernel0), 1);  if isempty(iK), iK = 1; end
    ctrl.jKernel.setSelectedIndex(iK - 1);
    panel_eigenfilter_design('BuildAtomSliders', ctrl.jParams, kernel0, bounds, @()bst_call(@OnParamSettle));
end

%% ===== accessors =====
function vals = ReadVals() %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'AtomDesigner');
    if isempty(ctrl), vals = [0 0 0]; return; end
    vals = panel_eigenfilter_design('ReadAtomVals', ctrl.jParams);
end

function k = CurrentKernel() %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'AtomDesigner');
    if isempty(ctrl), k = ''; return; end
    idx = max(1, min(numel(ctrl.atomKeys), double(ctrl.jKernel.getSelectedIndex()) + 1));
    k = ctrl.atomKeys{idx};
end

function v = CurrentOperator() %#ok<DEFNU>
    v = 'Laplace-Beltrami';
    ctrl = bst_get('PanelControls', 'AtomDesigner');  if isempty(ctrl), return; end
    vmap = {'LB-Connectome','Laplace-Beltrami'};                              % connectomic=1, geometric=2
    v = vmap{max(1, min(2, double(ctrl.jOperator.getSelectedIndex()) + 1))};
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
