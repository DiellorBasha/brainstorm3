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

    % Display figure handles + impulse-response state (shared across nested callbacks)
    hFigResp = [];     % spectral response h(lambda)
    hFigImp  = [];     % cortical impulse response (point-spread of a delta)
    ImpSeed  = [];     % seed vertex (global surface index) for the impulse
    ImpCache = [];     % lazily loaded {Op, Surf, supported, nVerts}
    ImpPick  = false;  % one-shot "pick seed on cortex" armed flag

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
            panel_eigenfilter_design('CurrentKernel', jKernel, keys), Lambda, @() OnSettle());
        java_setcb(jKernel, 'ActionPerformedCallback', @(hh,ee) OnKernel());
    jPanelNew.add('br hfill', jPanelDes);

    % ===== DISPLAY TOGGLES =====
    jPanelDisp = gui_river([2,2], [0,10,8,10], 'Display');
        jToggleDisp = gui_component('toggle', jPanelDisp, '', 'Show spectral response  h(\lambda)', [], [], @(hh,ee) ToggleDisplay());
        % Cortical impulse response: the filter atom = filtered delta (point-spread on the cortex)
        jToggleImp  = gui_component('toggle', jPanelDisp, 'br', 'Show cortical impulse response', [], [], @(hh,ee) ToggleImpulse());
        gui_component('label', jPanelDisp, 'br', 'Seed vertex: ', [], [], [], []);
        jTextSeed = gui_component('text', jPanelDisp, 'tab hfill', '', [], [], @(hh,ee) OnSeedText(), []);
        gui_component('button', jPanelDisp, '', 'Pick', [], [], @(hh,ee) ArmPick(), []);
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
                  'jToggleImp',   jToggleImp, ...
                  'jTextSeed',    jTextSeed, ...
                  'jTextComment', jTextComment, ...
                  'Lambda',       Lambda, ...
                  'EigenFile',    EigenFile);
    bstPanelNew = BstPanel(panelName, jPanelScroll, ctrl);

    %% ===== NESTED CALLBACKS =====
    function OnKernel(varargin)
        key = panel_eigenfilter_design('CurrentKernel', jKernel, keys);
        panel_eigenfilter_design('BuildSliders', jParams, key, Lambda, @() OnSettle());
        OnSettle();
    end

    function OnSettle(varargin)
        UpdateResponse();
        UpdateImpulse();
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

    %% ----- Cortical impulse response (filter atom = filtered delta) -----
    function ToggleImpulse(varargin)
        if jToggleImp.isSelected()
            if ~EnsureImpCache(); jToggleImp.setSelected(false); return; end
            if isempty(hFigImp) || ~ishandle(hFigImp)
                hFigImp = figure('MenuBar','none', 'Toolbar','figure', 'NumberTitle','off', ...
                                 'Name','Eigenfilter cortical impulse response', 'Color','w', ...
                                 'CloseRequestFcn', @(h,e) i_imp_closed());
                rotate3d(hFigImp, 'on');
            end
            UpdateImpulse();
        else
            if ~isempty(hFigImp) && ishandle(hFigImp); delete(hFigImp); end
            hFigImp = [];
        end
    end

    function ok = EnsureImpCache()
        ok = false;
        if ~isempty(ImpCache); ok = ImpCache.supported; end
        if ~isempty(ImpCache) && ~ImpCache.supported; i_imp_unsupported(); return; end
        if ~isempty(ImpCache); return; end
        try
            Op   = in_bst_operator(EigenMat.OperatorFile);
            Surf = in_tess_bst(EigenMat.ParentSurface, 0);
        catch err
            bst_error(['Cortical impulse response: cannot load operator/surface.' 10 err.message], 'Eigenfilter options', 0);
            return;
        end
        supported = strcmpi(EigenMat.Variant, 'Laplace-Beltrami');   % scalar, real, vertex-indexed
        ImpCache = struct('Op', Op, 'Surf', Surf, 'supported', supported, 'nVerts', size(Surf.Vertices,1));
        if ~supported; i_imp_unsupported(); return; end
        % Default seed = vertex nearest the first hemisphere's centroid
        if isempty(ImpSeed)
            gv = EigenMat.GlobalVertices{1};
            if isempty(gv); gv = (1:ImpCache.nVerts)'; end
            Vh = Surf.Vertices(gv,:);
            [~, ii] = min(sum((Vh - mean(Vh,1)).^2, 2));
            ImpSeed = gv(ii);
            jTextSeed.setText(num2str(ImpSeed));
        end
        ok = true;
    end

    function i_imp_unsupported()
        bst_error(sprintf(['Cortical impulse response currently supports the Laplace-Beltrami variant only.' 10 ...
            'This basis is ''%s'' (vector/complex variants such as Dirac and Connection Laplacian are coming next).'], ...
            EigenMat.Variant), 'Eigenfilter options', 0);
    end

    function i_imp_closed()
        if ~isempty(hFigImp) && ishandle(hFigImp); delete(hFigImp); end
        hFigImp = [];
        try
            jToggleImp.setSelected(false);
        catch
        end
    end

    function UpdateImpulse(varargin)
        if isempty(hFigImp) || ~ishandle(hFigImp); return; end
        if isempty(ImpCache) || ~ImpCache.supported || isempty(ImpSeed); return; end
        key    = panel_eigenfilter_design('CurrentKernel', jKernel, keys);
        params = panel_eigenfilter_design('ReadParams', jParams, Lambda);
        % Atom = filter applied to a unit delta at the seed vertex
        delta = zeros(ImpCache.nVerts, 1);
        delta(ImpSeed) = 1;
        [psf, msg, isErr] = bst_eigenfilter('Analysis', delta, EigenMat, ImpCache.Op, key, params);
        if isErr; bst_error(['Cortical impulse response: ' msg], 'Eigenfilter options', 0); return; end
        DrawImpulse(psf(:,1), key);
    end

    function DrawImpulse(psf, key)
        Surf = ImpCache.Surf;
        hAx  = findobj(hFigImp, 'Type', 'axes');
        if isempty(hAx); hAx = axes('Parent', hFigImp); else; hAx = hAx(1); end
        [az, el] = view(hAx);                       % preserve the user's camera across redraws
        cla(hAx);
        patch('Parent', hAx, 'Vertices', Surf.Vertices, 'Faces', Surf.Faces, ...
              'FaceVertexCData', double(psf(:)), 'FaceColor', 'interp', 'EdgeColor', 'none', ...
              'FaceLighting', 'gouraud', 'AmbientStrength', 0.8, 'DiffuseStrength', 0.4, ...
              'SpecularStrength', 0.05, 'ButtonDownFcn', @(h,e) i_axes_click());
        hold(hAx, 'on');
        sp = Surf.Vertices(ImpSeed, :);
        plot3(hAx, sp(1), sp(2), sp(3), 'k.', 'MarkerSize', 22, 'ButtonDownFcn', @(h,e) i_axes_click());
        hold(hAx, 'off');
        axis(hAx, 'equal', 'off', 'tight');
        if all([az el] == [0 90])                  % first draw: aim the camera at the seed
            view(hAx, sp - mean(Surf.Vertices, 1));
        else
            view(hAx, az, el);
        end
        camlight(hAx, 'headlight'); material(hAx, 'dull');
        m = max(abs(double(psf(:))));  if ~(m > 0); m = 1; end
        if min(psf(:)) < -1e-12*m                  % sign-changing (band-pass) -> symmetric diverging
            colormap(hAx, i_diverging()); caxis(hAx, [-m m]);
        else                                       % non-negative (low-pass) -> hot
            colormap(hAx, hot(256)); caxis(hAx, [0 m]);
        end
        colorbar(hAx);
        title(hAx, sprintf('Impulse response: %s   (seed vertex %d)', key, ImpSeed), 'Interpreter', 'none');
    end

    function OnSeedText(varargin)
        if ~EnsureImpCache(); return; end
        v = str2double(char(jTextSeed.getText()));
        if isnan(v) || v < 1 || v > ImpCache.nVerts
            jTextSeed.setText(num2str(ImpSeed)); return;
        end
        ImpSeed = round(v);
        UpdateImpulse();
    end

    function ArmPick(varargin)
        if ~jToggleImp.isSelected() || isempty(hFigImp) || ~ishandle(hFigImp)
            bst_error('Open the cortical impulse response first (Show cortical impulse response).', 'Eigenfilter options', 0);
            return;
        end
        ImpPick = true;
        rotate3d(hFigImp, 'off');
        set(hFigImp, 'Pointer', 'crosshair');
    end

    function i_axes_click(varargin)
        if ~ImpPick; return; end                   % only re-seed when picking is armed
        ImpPick = false;
        set(hFigImp, 'Pointer', 'arrow');
        hAx = findobj(hFigImp, 'Type', 'axes'); if isempty(hAx); return; end
        cp = get(hAx(1), 'CurrentPoint');          % 2x3 click ray; use the near point
        p  = cp(1, :);
        V  = ImpCache.Surf.Vertices;
        [~, iv] = min(sum((V - p).^2, 2));         % nearest vertex (ray near-point heuristic)
        ImpSeed = iv;
        jTextSeed.setText(num2str(iv));
        rotate3d(hFigImp, 'on');
        UpdateImpulse();
    end

    function cmap = i_diverging()
        % blue-white-red diverging colormap (no toolbox dependency)
        n = 128; t = linspace(0,1,n)';
        lo = [t, t, ones(n,1)];                     % blue -> white
        hi = [ones(n,1), flipud(t), flipud(t)];     % white -> red
        cmap = [lo; hi(2:end,:)];
    end

    function ButtonCancel_Callback(varargin)
        if ~isempty(hFigResp) && ishandle(hFigResp); close(hFigResp); end
        if ~isempty(hFigImp)  && ishandle(hFigImp);  delete(hFigImp);  end
        gui_hide(panelName);
    end

    function ButtonOk_Callback(varargin)
        if ~isempty(hFigResp) && ishandle(hFigResp); close(hFigResp); end
        if ~isempty(hFigImp)  && ishandle(hFigImp);  delete(hFigImp);  end
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
