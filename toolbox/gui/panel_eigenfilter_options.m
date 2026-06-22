function varargout = panel_eigenfilter_options(varargin)
% PANEL_EIGENFILTER_OPTIONS: Options for the eigen-domain spatial filter (bst_eigen 'filter').
%
% USAGE:  [bstPanel, panelName] = panel_eigenfilter_options('CreatePanel', EigenFile)
%                            s  = panel_eigenfilter_options('GetPanelContents')
%
% The spatial analogue of panel_timefreq_options. EigenFile is an eigen_ node; the panel
% reads its eigenvalues (Lambda) to drive the kernel sliders, embeds the shared
% panel_eigenfilter_design widget, and (Display toggle) previews the spectral response
% h(lambda). GetPanelContents returns an OPTIONS struct ready for bst_eigen.
% CreatePanel is arg-type aware: a char/file argument is the EigenFile path (now); a
% struct argument is reserved for a future process_eigenfilter (sProcess, sInputs) path
% wired via an 'editpref' option, so the panel attaches with no rewrite.

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
function [bstPanelNew, panelName] = CreatePanel(EigenFile) %#ok<DEFNU>
    panelName = 'EigenfilterOptions';
    import java.awt.*;
    import javax.swing.*;

    % Arg-type aware: char/file => EigenFile path (now); struct => future (sProcess,sInputs)
    if isstruct(EigenFile)
        bst_error('panel_eigenfilter_options: the (sProcess, sInputs) path is not implemented yet; pass an eigen_ file.', 'Eigenfilter options', 0);
        bstPanelNew = []; panelName = []; return;
    end
    if isempty(EigenFile) || ~ischar(EigenFile)
        bst_error('panel_eigenfilter_options: a valid eigen_ file is required.', 'Eigenfilter options', 0);
        bstPanelNew = []; panelName = []; return;
    end

    % Load the eigenbasis (the spatial axis)
    EigenMat = in_bst_eigen(EigenFile);
    Lambda = [];
    for h = 1:numel(EigenMat.Lambda)
        if ~isempty(EigenMat.Lambda{h})
            Lambda = EigenMat.Lambda{h}(:);
            break;
        end
    end
    if isempty(Lambda)
        bst_error('panel_eigenfilter_options: the eigen_ node has no eigenvalues.', 'Eigenfilter options', 0);
        bstPanelNew = []; panelName = []; return;
    end
    K     = numel(Lambda);
    nHemi = sum(~cellfun(@isempty, EigenMat.Lambda));

    % Display figure handle (shared across nested callbacks)
    hFigResp = [];

    % ===== MAIN PANEL =====
    jPanelNew = gui_river([5,5], [10,15,12,10]);

    % ===== EIGEN BASIS INFO =====
    jPanelInfo = gui_river([2,2], [0,10,10,10], 'Eigen basis');
        gui_component('label', jPanelInfo, '',   sprintf('Variant:  %s', EigenMat.Variant), [], [], [], []);
        gui_component('label', jPanelInfo, 'br', sprintf('Modes:  %d      Hemispheres:  %d', K, nHemi), [], [], [], []);
    jPanelNew.add('br hfill', jPanelInfo);

    % ===== FILTER DESIGN (reuse the shared design widget) =====
    jPanelDes = gui_river([2,2], [0,10,12,10], 'Filter design');
        [keys, displays] = panel_eigenfilter_design('Kernels');
        gui_component('label', jPanelDes, '', 'Kernel: ', [], [], [], []);
        jKernel = gui_component('combobox', jPanelDes, 'tab hfill', [], {displays}, [], [], []);
        iHeat = find(strcmp(keys, 'heat'), 1);
        if ~isempty(iHeat); jKernel.setSelectedIndex(iHeat-1); end
        jParams = gui_river([2,2], [0,2,0,2]);
        jPanelDes.add('br hfill', jParams);
        panel_eigenfilter_design('BuildSliders', jParams, ...
            panel_eigenfilter_design('CurrentKernel', jKernel, keys), Lambda, @() UpdateResponse());
        java_setcb(jKernel, 'ActionPerformedCallback', @(hh,ee) OnKernel());
    jPanelNew.add('br hfill', jPanelDes);

    % ===== DISPLAY TOGGLE =====
    jPanelDisp = gui_river([2,2], [0,10,8,10], 'Display');
        jToggleDisp = gui_component('toggle', jPanelDisp, '', 'Show spectral response', [], [], @(hh,ee) ToggleDisplay());
    jPanelNew.add('br hfill', jPanelDisp);

    % ===== OUTPUT COMMENT =====
    jPanelCom = gui_river([2,2], [0,10,8,10], 'Output');
        gui_component('label', jPanelCom, '', 'Comment: ', [], [], [], []);
        jTextComment = gui_component('text', jPanelCom, 'tab hfill', '', [], [], [], []);
    jPanelNew.add('br hfill', jPanelCom);

    % ===== OK / CANCEL =====
    gui_component('button', jPanelNew, 'br right', 'Cancel', [], [], @ButtonCancel_Callback, []);
    gui_component('button', jPanelNew, [],         'OK',     [], [], @ButtonOk_Callback, []);

    % ===== ASSEMBLE =====
    jPanelScroll = javax.swing.JScrollPane(jPanelNew);
    bst_mutex('create', panelName);
    ctrl = struct('jKernel',      jKernel, ...
                  'KernelKeys',   {keys}, ...
                  'jParams',      jParams, ...
                  'jToggleDisp',  jToggleDisp, ...
                  'jTextComment', jTextComment, ...
                  'Lambda',       Lambda, ...
                  'EigenFile',    EigenFile);
    bstPanelNew = BstPanel(panelName, jPanelScroll, ctrl);

    %% ===== NESTED CALLBACKS =====
    function OnKernel(varargin)
        key = panel_eigenfilter_design('CurrentKernel', jKernel, keys);
        panel_eigenfilter_design('BuildSliders', jParams, key, Lambda, @() UpdateResponse());
        UpdateResponse();
    end

    function ToggleDisplay(varargin)
        if jToggleDisp.isSelected()
            if isempty(hFigResp) || ~ishandle(hFigResp)
                hFigResp = figure('MenuBar','none', 'Toolbar','none', 'NumberTitle','off', ...
                                  'Name','Eigenfilter spectral response', 'Pointer','arrow');
            end
            UpdateResponse();
        else
            if ~isempty(hFigResp) && ishandle(hFigResp); close(hFigResp); end
            hFigResp = [];
        end
    end

    function UpdateResponse(varargin)
        if isempty(hFigResp) || ~ishandle(hFigResp); return; end
        key    = panel_eigenfilter_design('CurrentKernel', jKernel, keys);
        params = panel_eigenfilter_design('ReadParams', jParams, Lambda);
        hAxes  = findobj(hFigResp, 'Type', 'axes');
        if isempty(hAxes)
            hAxes = axes('Parent', hFigResp);
        else
            hAxes = hAxes(1);
        end
        panel_eigenfilter_design('DrawResponse', hAxes, key, params, Lambda);
    end

    function ButtonCancel_Callback(varargin)
        if ~isempty(hFigResp) && ishandle(hFigResp); close(hFigResp); end
        gui_hide(panelName);
    end

    function ButtonOk_Callback(varargin)
        if ~isempty(hFigResp) && ishandle(hFigResp); close(hFigResp); end
        bst_mutex('release', panelName);
    end
end


%% ===== GET PANEL CONTENTS =====
function s = GetPanelContents() %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'EigenfilterOptions');
    key    = panel_eigenfilter_design('CurrentKernel', ctrl.jKernel, ctrl.KernelKeys);
    params = panel_eigenfilter_design('ReadParams', ctrl.jParams, ctrl.Lambda);
    s = struct();
    s.Method       = 'filter';
    s.EigenFile    = ctrl.EigenFile;
    s.KernelName   = key;
    s.KernelParams = params;
    s.Comment      = char(ctrl.jTextComment.getText());
end
