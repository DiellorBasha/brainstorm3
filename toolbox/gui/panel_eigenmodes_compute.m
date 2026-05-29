function varargout = panel_eigenmodes_compute(varargin)
% PANEL_EIGENMODES_COMPUTE: Options for computing Laplace-Beltrami eigenmodes (GUI).
%
% USAGE:  bstPanel = panel_eigenmodes_compute('CreatePanel')
%                s = panel_eigenmodes_compute('GetPanelContents')

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
%
% Authors: Diellor Basha, 2026

eval(macro_method);
end


%% ===== CREATE PANEL =====
function [bstPanelNew, panelName] = CreatePanel() %#ok<DEFNU>
    panelName = 'EigenmodesCompute';
    import java.awt.*;
    import javax.swing.*;
    MassList = {'barycentric', 'voronoi', 'galerkin'};

    jPanelNew = gui_river([5,5], [10,15,12,10], 'Compute eigenmodes');
        gui_component('label', jPanelNew, '', 'Number of eigenmodes per hemisphere: ', [], [], [], []);
        jTextN = gui_component('text', jPanelNew, 'tab', '300', [], [], [], []);
        gui_component('label', jPanelNew, 'br', 'Mass matrix type: ', [], [], [], []);
        jComboMass = gui_component('combobox', jPanelNew, 'tab', [], {MassList}, [], [], []);
        jComboMass.setSelectedIndex(0);   % barycentric
    % Validation buttons
    gui_component('button', jPanelNew, 'br right', 'Cancel', [], [], @ButtonCancel_Callback, []);
    gui_component('button', jPanelNew, [], 'OK', [], [], @ButtonOk_Callback, []);

    bst_mutex('create', panelName);
    ctrl = struct('jTextN', jTextN, 'jComboMass', jComboMass);
    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);

    function ButtonCancel_Callback(varargin)
        gui_hide(panelName);
    end
    function ButtonOk_Callback(varargin)
        bst_mutex('release', panelName);
    end
end


%% ===== GET PANEL CONTENTS =====
function s = GetPanelContents() %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'EigenmodesCompute');
    nModes = str2double(char(ctrl.jTextN.getText()));
    if isnan(nModes) || (nModes < 1)
        error('Number of eigenmodes must be a positive integer.');
    end
    s.nModes   = round(nModes);
    s.MassType = char(ctrl.jComboMass.getSelectedItem());
end
