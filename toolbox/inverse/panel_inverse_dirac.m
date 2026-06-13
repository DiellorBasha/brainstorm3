function varargout = panel_inverse_dirac(varargin)
% PANEL_INVERSE_DIRAC: Options dialog for the interactive Dirac eigenmode inverse.
%
% USAGE:  bstPanel = panel_inverse_dirac('CreatePanel')
%         s        = panel_inverse_dirac('GetPanelContents')
%
% Small companion dialog for process_inverse_dirac('ComputeInteractive'): collects
% the handful of options consumed by bst_inverse_dirac (measure, SNR, number of
% Dirac modes per hemisphere, tau, noise regularization, sensor types). The Dirac
% source space is always the unconstrained 3-vector cortical current, so there are
% no orientation / output-mode sections (unlike panel_inverse_2018).
%
% SEE ALSO: process_inverse_dirac, bst_inverse_dirac, panel_inverse_2018

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
    panelName = 'DiracInverseOptions';
    % Java initializations
    import java.awt.*;
    import javax.swing.*;

    % Create main panel
    jPanelNew = gui_river([5,5], [0,10,10,10]);

    % ==== PANEL: MEASURE ====
    jPanelMeasure = gui_river([1,1], [0,6,6,6], 'Measure');
        jGroupMeasure = ButtonGroup();
        jRadioCurrent = gui_component('radio', jPanelMeasure, [],   'Current density map', jGroupMeasure, '', [], []);
        jRadioDspm    = gui_component('radio', jPanelMeasure, 'br', 'dSPM',                jGroupMeasure, '', [], []);
        jRadioSloreta = gui_component('radio', jPanelMeasure, 'br', 'sLORETA',             jGroupMeasure, '', [], []);
        jRadioDspm.setSelected(1);
    jPanelNew.add('br hfill', jPanelMeasure);

    % ==== PANEL: DIRAC OPTIONS ====
    jPanelOpt = gui_river([3,3], [0,6,6,6], 'Dirac eigenmode options');
        % SNR
        gui_component('label', jPanelOpt, [], 'Signal-to-noise ratio (SNR): ', [], '', [], []);
        jTextSnr = gui_component('texttime', jPanelOpt, 'tab', '', [], '', [], []);
        gui_validate_text(jTextSnr, [], [], {0, 10000, 100}, '', 2, 3, []);
        % Modes per hemisphere
        gui_component('label', jPanelOpt, 'br', 'Dirac modes per hemisphere: ', [], '', [], []);
        jTextModes = gui_component('texttime', jPanelOpt, 'tab', '', [], '', [], []);
        gui_validate_text(jTextModes, [], [], {1, 100000, 1}, '', 0, 400, []);
        % Tau
        gui_component('label', jPanelOpt, 'br', 'Dirac tau (intrinsic-extrinsic mix): ', [], '', [], []);
        jTextTau = gui_component('texttime', jPanelOpt, 'tab', '', [], '', [], []);
        gui_validate_text(jTextTau, [], [], {0, 1, 100}, '', 2, 0.5, []);
        % Noise regularization
        gui_component('label', jPanelOpt, 'br', 'Noise covariance regularization (fraction): ', [], '', [], []);
        jTextNoiseReg = gui_component('texttime', jPanelOpt, 'tab', '', [], '', [], []);
        gui_validate_text(jTextNoiseReg, [], [], {0, 1, 1000}, '', 3, 0.1, []);
        % Sensor types
        gui_component('label', jPanelOpt, 'br', 'Sensor types or names: ', [], '', [], []);
        jTextSensors = gui_component('text', jPanelOpt, 'tab hfill', 'MEG', [], '', [], []);
    jPanelNew.add('br hfill', jPanelOpt);

    % ===== VALIDATION BUTTONS =====
    jPanelValid = gui_river([1,1], [0,6,6,6]);
        gui_component('label', jPanelValid, 'hfill', ' ');
        gui_component('Button', jPanelValid, 'right', 'Cancel', [], [], @ButtonCancel_Callback, []);
        gui_component('Button', jPanelValid, [],      'OK',     [], [], @ButtonOk_Callback,     []);
    jPanelNew.add('br right', jPanelValid);

    % ===== PANEL CREATION =====
    % Return a mutex to wait for panel close
    bst_mutex('create', panelName);
    % Controls struct
    ctrl = struct(...
        'jRadioCurrent', jRadioCurrent, ...
        'jRadioDspm',    jRadioDspm, ...
        'jRadioSloreta', jRadioSloreta, ...
        'jTextSnr',      jTextSnr, ...
        'jTextModes',    jTextModes, ...
        'jTextTau',      jTextTau, ...
        'jTextNoiseReg', jTextNoiseReg, ...
        'jTextSensors',  jTextSensors);
    % Create the BstPanel object that is returned by the function
    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);


%% =================================================================================
%  === LOCAL CALLBACKS  ============================================================
%  =================================================================================
    %% ===== BUTTON: CANCEL =====
    function ButtonCancel_Callback(varargin)
        gui_hide(panelName);
    end

    %% ===== BUTTON: OK =====
    function ButtonOk_Callback(varargin)
        bst_mutex('release', panelName);
    end
end


%% =================================================================================
%  === EXTERNAL CALLBACKS  =========================================================
%  =================================================================================
%% ===== GET PANEL CONTENTS =====
function s = GetPanelContents() %#ok<DEFNU>
    % Get panel controls handles
    ctrl = bst_get('PanelControls', 'DiracInverseOptions');
    if isempty(ctrl)
        s = [];
        return;
    end
    % Measure
    if ctrl.jRadioCurrent.isSelected()
        s.Measure = 'amplitude';
    elseif ctrl.jRadioSloreta.isSelected()
        s.Measure = 'sloreta';
    else
        s.Measure = 'dspm2018';
    end
    % Numeric options
    s.SNR         = str2double(char(ctrl.jTextSnr.getText()));
    s.nModes      = str2double(char(ctrl.jTextModes.getText()));
    s.Tau         = str2double(char(ctrl.jTextTau.getText()));
    s.NoiseReg    = str2double(char(ctrl.jTextNoiseReg.getText()));
    s.SensorTypes = strtrim(char(ctrl.jTextSensors.getText()));
end
