function varargout = panel_import_pet(varargin)
% PANEL_IMPORT_PET: Options for importing and processing a PET volume.
%
% USAGE: [bstPanelNew, panelName] = panel_import_pet('CreatePanel')
%        Options = gui_show_dialog('Import & process PET', @panel_import_pet, 1, [])
%
% Interactive front-end for process_import_pet(). A processing-mode selector:
%   - Recommended : import + full validated pipeline with default settings (no
%                   options to choose).
%   - Advanced    : reveals the pipeline knobs for customization.
%   - Import only : raw volume, no processing.
% The Advanced knobs are shown ONLY in Advanced mode; their initial values are
% the recommended defaults, so Recommended just reads those defaults.
%
% Returns an Options struct consumed by process_import_pet (fields: RunPipeline,
% Aggregate, Register, Reslice, PvcMethod, PvcFwhm, SuvrRef, DoProject).

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
% Authors: Diellor Basha, 2025-2026
%          Raymundo Cassani, 2025

eval(macro_method);
end

%% ===== CREATE PANEL =====
function [bstPanelNew, panelName] = CreatePanel()
    panelName = 'panel_import_pet';
    import java.awt.*
    import javax.swing.*

    % === MAIN LAYOUT ===
    jPanelMain = gui_river([5, 5], [0, 10, 10, 10]);

    % === PROCESSING MODE ===
    jPanelMode = gui_river([2, 2], [0, 10, 10, 10], 'Processing');
    jGroupMode = ButtonGroup();
    jRadioRec = gui_component('radio', jPanelMode, 'br', 'Recommended  (import + full pipeline)');
    jRadioAdv = gui_component('radio', jPanelMode, 'br', 'Advanced  (customize each step)');
    jRadioRaw = gui_component('radio', jPanelMode, 'br', 'Import only  (raw volume)');
    jGroupMode.add(jRadioRec);
    jGroupMode.add(jRadioAdv);
    jGroupMode.add(jRadioRaw);
    jRadioRec.setSelected(true);
    jPanelMain.add('br hfill', jPanelMode);

    % === ADVANCED OPTIONS (initial values = recommended defaults) ===
    jPanelAdv = gui_river([2, 2], [0, 10, 10, 10], 'Advanced options');
    % Aggregate frames
    gui_component('label', jPanelAdv, 'br', 'Aggregate frames: ');
    jComboAggregate = gui_component('combobox', jPanelAdv, 'tab', [], {{'Mean', 'Sum', 'Median', 'Max', 'Min', 'First', 'Last'}});
    % Register to MRI
    gui_component('label', jPanelAdv, 'br', 'Register to MRI: ');
    jComboRegister = gui_component('combobox', jPanelAdv, 'tab', [], {{'SPM', 'MNI', 'Ignore'}});
    % Reslice
    jCheckReslice = gui_component('checkbox', jPanelAdv, 'br', 'Reslice on the MRI grid');
    jCheckReslice.setSelected(true);
    % PVC method
    gui_component('label', jPanelAdv, 'br', 'Partial volume correction: ');
    jComboPvc = gui_component('combobox', jPanelAdv, 'tab', [], {{'GTM (Rousset)', 'Muller-Gartner', 'None'}});
    % PVC FWHM (0 = auto from scanner metadata)
    gui_component('label', jPanelAdv, 'br', 'PVC PSF FWHM (0 = auto): ');
    jTextFwhm = gui_component('text', jPanelAdv, 'tab', '0');
    jTextFwhm.setMaximumSize(java.awt.Dimension(50, 20));
    % SUVR reference (empty = tracer-aware default)
    gui_component('label', jPanelAdv, 'br', 'SUVR reference (empty = tracer): ');
    jTextRef = gui_component('text', jPanelAdv, 'tab', '');
    jTextRef.setMaximumSize(java.awt.Dimension(120, 20));
    % Project to surface
    jCheckProject = gui_component('checkbox', jPanelAdv, 'br', 'Project SUVR to cortical surface');
    jCheckProject.setSelected(false);
    jPanelMain.add('br hfill', jPanelAdv);
    % Advanced options are shown ONLY in Advanced mode
    jPanelAdv.setVisible(false);

    % Mode selection reveals/hides the Advanced options and resizes the dialog
    java_setcb(jRadioRec, 'ActionPerformedCallback', @(h, ev)UpdateMode());
    java_setcb(jRadioAdv, 'ActionPerformedCallback', @(h, ev)UpdateMode());
    java_setcb(jRadioRaw, 'ActionPerformedCallback', @(h, ev)UpdateMode());

    % === BUTTONS ===
    jPanelButtons = gui_river([2 0], [0 5 0 5]);
    gui_component('button', jPanelButtons, 'br right', 'Cancel', [], [], @ButtonCancel_Callback);
    gui_component('button', jPanelButtons, '', 'OK', [], [], @ButtonOk_Callback);
    jPanelMain.add('br right', jPanelButtons);
    % === Create mutex ===
    bst_mutex('create', panelName);
    % === Return panel object ===
    bstPanelNew = BstPanel(panelName, ...
        jPanelMain, ...
        struct('jRadioRec',       jRadioRec, ...
               'jRadioAdv',       jRadioAdv, ...
               'jRadioRaw',       jRadioRaw, ...
               'jComboAggregate', jComboAggregate, ...
               'jComboRegister',  jComboRegister, ...
               'jCheckReslice',   jCheckReslice, ...
               'jComboPvc',       jComboPvc, ...
               'jTextFwhm',       jTextFwhm, ...
               'jTextRef',        jTextRef, ...
               'jCheckProject',   jCheckProject));

%% =================================================================================
%  === INTERNAL CALLBACKS ==========================================================
%  =================================================================================
    function UpdateMode()
        jPanelAdv.setVisible(jRadioAdv.isSelected());
        % Resize the dialog to fit the new content
        jTop = jPanelMain.getTopLevelAncestor();
        if ~isempty(jTop)
            jPanelMain.revalidate();
            jTop.pack();
        end
    end
    function ButtonCancel_Callback(~, ~)
        gui_hide(panelName);
    end
    function ButtonOk_Callback(~, ~)
        bst_mutex('release', panelName);  % Triggers gui_show_dialog to call GetPanelContents
    end
end

%% =================================================================================
%  === EXTERNAL CALLBACKS ==========================================================
%  =================================================================================
%% ===== GET PANEL CONTENTS =====
function s = GetPanelContents()
    ctrl = bst_get('PanelControls', 'panel_import_pet');
    % Pipeline runs unless "Import only" is selected. In Recommended mode the
    % Advanced widgets are untouched, so they still hold the recommended defaults.
    s.RunPipeline = ~ctrl.jRadioRaw.isSelected();
    s.Aggregate   = lower(char(ctrl.jComboAggregate.getSelectedItem()));
    s.Register    = lower(char(ctrl.jComboRegister.getSelectedItem()));
    s.Reslice     = ctrl.jCheckReslice.isSelected();
    % PVC method: map display label -> internal value
    pvcDisp = lower(char(ctrl.jComboPvc.getSelectedItem()));
    if ~isempty(strfind(pvcDisp, 'gtm')) %#ok<STREMP>
        s.PvcMethod = 'gtm';
    elseif ~isempty(strfind(pvcDisp, 'none')) %#ok<STREMP>
        s.PvcMethod = 'none';
    else
        s.PvcMethod = 'mg';
    end
    % PVC FWHM: 0 (or blank) -> auto
    fwhm = str2double(char(ctrl.jTextFwhm.getText()));
    if isnan(fwhm)
        fwhm = 0;
    end
    s.PvcFwhm   = fwhm;
    s.SuvrRef   = strtrim(char(ctrl.jTextRef.getText()));
    s.DoProject = ctrl.jCheckProject.isSelected();
end
