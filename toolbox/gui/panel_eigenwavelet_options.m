function varargout = panel_eigenwavelet_options(varargin)
% PANEL_EIGENWAVELET_OPTIONS: Options for the eigen-domain spatial WAVELET (bst_eigen 'wavelet').
%
% USAGE:  [bstPanel, panelName] = panel_eigenwavelet_options('CreatePanel', EigenFile)
%                            s  = panel_eigenwavelet_options('GetPanelContents')
%
% The wavelet sibling of panel_eigenfilter_options. EigenFile is an eigen_ node; the panel
% designs a spectral-graph-wavelet FRAME (family + number of scales Nf) over the operator
% spectrum and previews the cortical wavelet ATOM (the filtered delta = "what the wavelet
% looks like") for a chosen scale and seed vertex.
%
% Vector ATOM + ORIENTATION (Dirac): for the quaternionic Dirac basis the atom is a 3-D
% dipole packet whose orientation is steered EXACTLY by right-quaternion multiplication
% (bst_eigenwavelet('Steer'); the full quaternion is carried end-to-end). On the atom
% figure: SCROLL rotates the orientation about the current axis, keys X/Y/Z pick the axis,
% R resets, arrow keys nudge. A 3-D arrow GLYPH at the seed shows the current dipole
% orientation. Scalar (Laplace-Beltrami) atoms show a magnitude envelope (no orientation).
% Complex (Connection Laplacian) and face-domain variants are deferred (pending the
% per-vertex tangent frame / face-centroid display).
%
% GetPanelContents returns an OPTIONS struct ready for bst_eigen (Method='wavelet').

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
    panelName = 'EigenwaveletOptions';
    import java.awt.*;
    import javax.swing.*;

    if isstruct(EigenFile)
        bst_error('panel_eigenwavelet_options: the (sProcess, sInputs) path is not implemented yet; pass an eigen_ file.', 'Eigenwavelet options', 0);
        bstPanelNew = []; panelName = []; return;
    end
    if isempty(EigenFile) || ~ischar(EigenFile)
        bst_error('panel_eigenwavelet_options: a valid eigen_ file is required.', 'Eigenwavelet options', 0);
        bstPanelNew = []; panelName = []; return;
    end

    % Load the eigenbasis (full: Phi/GlobalVertices needed for the atom)
    EigenMat = in_bst_eigen(EigenFile);
    lams = EigenMat.Lambda(~cellfun(@isempty, EigenMat.Lambda));
    if isempty(lams)
        bst_error('panel_eigenwavelet_options: the eigen_ node has no eigenvalues.', 'Eigenwavelet options', 0);
        bstPanelNew = []; panelName = []; return;
    end
    Lrange = [min(cellfun(@min, lams)), max(cellfun(@max, lams))];   % GLOBAL spectrum
    K      = sum(cellfun(@numel, lams));
    nHemi  = numel(lams);
    Variant = EigenMat.Variant;
    Kind    = i_variant_kind(Variant);   % 'scalar' | 'vector' | 'complex' | 'deferred'

    % ---- atom-display state (shared across nested callbacks) ----
    hFigAtom = [];
    AtomCache = [];                 % {Op, Surf, nVerts [, E1,E2,Nrm for complex]}
    AtomSeed  = [];                 % seed vertex (global index)
    AtomScale = 1;                  % which frame member (scale) is shown
    AtomBase  = [];                 % base atom (full field) [rows x 1 x M] for the current seed/frame
    AtomQuat  = [1;0;0;0];          % orientation quaternion (Dirac, 'vector'), identity = base dipole
    AtomPhase = 0;                  % orientation phase (Connection Laplacian, 'complex'), U(1)
    AtomAxis  = 'x';               % current steering axis (x|y|z) for 'vector'
    AtomPick  = false;             % one-shot "pick seed" armed
    Frame     = [];                 % current designed frame
    M         = 0;                  % number of frame members

    % ===== MAIN PANEL =====
    jPanelNew = gui_river([5,5], [10,15,12,10]);

    % ===== EIGEN BASIS INFO =====
    jPanelInfo = gui_river([2,2], [0,10,10,10], 'Eigen basis');
        gui_component('label', jPanelInfo, '',   sprintf('Variant:  %s', Variant), [], [], [], []);
        gui_component('label', jPanelInfo, 'br', sprintf('Modes:  %d      Hemispheres:  %d', K, nHemi), [], [], [], []);
    jPanelNew.add('br hfill', jPanelInfo);

    % ===== FRAME DESIGN =====
    jPanelDes = gui_river([2,2], [0,10,12,10], 'Wavelet frame');
        famKeys = {'itersine','mexhat','heat'};
        famDisp = {'itersine (tight frame)','mexhat (band-pass + scaling fn)','heat (diffusion scale-space)'};
        gui_component('label', jPanelDes, '', 'Family: ', [], [], [], []);
        jFamily = gui_component('combobox', jPanelDes, 'tab hfill', [], {famDisp}, [], [], []);
        gui_component('label', jPanelDes, 'br', 'Scales (Nf): ', [], [], [], []);
        [jNf, ~] = i_slider(jPanelDes, 'tab hfill', 2, 12, 6);
        jNfLabel = gui_component('label', jPanelDes, '', '6', [], [], [], []);
    jPanelNew.add('br hfill', jPanelDes);

    % ===== ATOM DISPLAY =====
    jPanelDisp = gui_river([2,2], [0,10,8,10], 'Wavelet atom');
        jToggleAtom = gui_component('toggle', jPanelDisp, '', 'Show wavelet atom', [], [], @(hh,ee) ToggleAtom(), []);
        gui_component('label', jPanelDisp, 'br', 'Scale: ', [], [], [], []);
        [jScale, ~] = i_slider(jPanelDisp, 'tab hfill', 1, 6, 1);
        jScaleLabel = gui_component('label', jPanelDisp, '', '1', [], [], [], []);
        gui_component('label', jPanelDisp, 'br', 'Seed vertex: ', [], [], [], []);
        jTextSeed = gui_component('text', jPanelDisp, 'tab hfill', '', [], [], @(hh,ee) OnSeedText(), []);
        gui_component('button', jPanelDisp, '', 'Pick', [], [], @(hh,ee) ArmPick(), []);
        if strcmp(Kind, 'vector')
            gui_component('label', jPanelDisp, 'br', '<HTML><i>Orient: scroll=rotate, X/Y/Z=axis, R=reset</i>', [], [], [], []);
        elseif strcmp(Kind, 'complex')
            gui_component('label', jPanelDisp, 'br', '<HTML><i>Orient: scroll=rotate phase, R=reset</i>', [], [], [], []);
        elseif strcmp(Kind, 'deferred')
            gui_component('label', jPanelDisp, 'br', sprintf('<HTML><i>Atom display for ''%s'' is deferred.</i>', Variant), [], [], [], []);
        end
    jPanelNew.add('br hfill', jPanelDisp);

    % ===== OUTPUT COMMENT =====
    jPanelCom = gui_river([2,2], [0,10,8,10], 'Output');
        gui_component('label', jPanelCom, '', 'Comment: ', [], [], [], []);
        jTextComment = gui_component('text', jPanelCom, 'tab hfill', '', [], [], [], []);
    jPanelNew.add('br hfill', jPanelCom);

    % ===== OK / CANCEL =====
    gui_component('button', jPanelNew, 'br right', 'Cancel', [], [], @ButtonCancel_Callback, []);
    gui_component('button', jPanelNew, [],         'OK',     [], [], @ButtonOk_Callback, []);

    % wire frame controls
    java_setcb(jFamily, 'ActionPerformedCallback', @(h,e) OnFrame());
    java_setcb(jNf,     'StateChangedCallback',    @(h,e) OnNf(h));
    java_setcb(jScale,  'StateChangedCallback',    @(h,e) OnScale(h));
    DesignFrame();   % initial frame

    % ===== ASSEMBLE =====
    jPanelScroll = javax.swing.JScrollPane(jPanelNew);
    bst_mutex('create', panelName);
    ctrl = struct('jFamily',      jFamily, 'FamilyKeys', {famKeys}, ...
                  'jNf',          jNf, ...
                  'jToggleAtom',  jToggleAtom, ...
                  'jScale',       jScale, ...
                  'jTextSeed',    jTextSeed, ...
                  'jTextComment', jTextComment, ...
                  'Lrange',       Lrange, ...
                  'EigenFile',    EigenFile);
    bstPanelNew = BstPanel(panelName, jPanelScroll, ctrl);

    %% ===== NESTED CALLBACKS =====
    function OnFrame(varargin)
        DesignFrame();
        RecomputeBase();
        UpdateAtom();
    end

    function OnNf(js)
        n = double(js.getValue());
        jNfLabel.setText(num2str(n));
        if ~js.getValueIsAdjusting()
            DesignFrame(); RecomputeBase(); UpdateAtom();
        end
    end

    function OnScale(js)
        AtomScale = max(1, min(M, double(js.getValue())));
        jScaleLabel.setText(num2str(AtomScale));
        if ~js.getValueIsAdjusting(); UpdateAtom(); end
    end

    function DesignFrame()
        fam = famKeys{max(1, jFamily.getSelectedIndex()+1)};
        Nf  = double(jNf.getValue());
        Frame = bst_eigenwavelet('Design', fam, Nf, Lrange);
        M = numel(Frame.g);
        % keep the Scale slider in range
        jScale.setMaximum(M);
        if AtomScale > M; AtomScale = M; jScale.setValue(M); jScaleLabel.setText(num2str(M)); end
    end

    function ToggleAtom(varargin)
        if jToggleAtom.isSelected()
            if strcmp(Kind, 'deferred')
                bst_error(sprintf(['Wavelet atom display for ''%s'' is deferred (needs the per-vertex tangent ' ...
                    'frame for Connection Laplacian, or face-centroid display for face variants).'], Variant), 'Eigenwavelet options', 0);
                jToggleAtom.setSelected(false); return;
            end
            if ~EnsureCache(); jToggleAtom.setSelected(false); return; end
            if isempty(hFigAtom) || ~ishandle(hFigAtom)
                hFigAtom = figure('MenuBar','none', 'Toolbar','figure', 'NumberTitle','off', ...
                                  'Name','Eigenwavelet cortical atom', 'Color','w', ...
                                  'CloseRequestFcn', @(h,e) i_atom_closed());
                % Default 3-D interactions (drag = rotate view) stay on; we only override
                % scroll (orientation steer) + keys (axis/reset). rotate3d as a MODE would
                % block these callbacks, so it is intentionally NOT used.
                set(hFigAtom, 'WindowScrollWheelFcn', @(h,e) i_atom_scroll(e), ...
                              'WindowKeyPressFcn',    @(h,e) i_atom_key(e));
            end
            RecomputeBase();
            UpdateAtom();
        else
            if ~isempty(hFigAtom) && ishandle(hFigAtom); delete(hFigAtom); end
            hFigAtom = [];
        end
    end

    function ok = EnsureCache()
        ok = false;
        if ~isempty(AtomCache); ok = true; end
        if ~isempty(AtomCache); return; end
        try
            Op   = in_bst_operator(EigenMat.OperatorFile);
            Surf = in_tess_bst(EigenMat.ParentSurface, 0);
        catch err
            bst_error(['Wavelet atom: cannot load operator/surface.' 10 err.message], 'Eigenwavelet options', 0);
            return;
        end
        AtomCache = struct('Op', Op, 'Surf', Surf, 'nVerts', size(Surf.Vertices,1));
        if strcmp(Kind, 'complex')
            % Global canonical tangent frame to decode complex eigenmodes -> 3-D tangent
            % vectors: field(v) = real(z)*e1(v) + imag(z)*e2(v). Read from the operator
            % (bst_operator_frame: stored frame, or nxr fallback for older files).
            nVg = AtomCache.nVerts;
            E1 = zeros(nVg,3); E2 = zeros(nVg,3); Nrm = zeros(nVg,3);
            for hh2 = 1:numel(EigenMat.GlobalVertices)
                gvh = EigenMat.GlobalVertices{hh2}; if isempty(gvh); continue; end
                Fr = bst_operator_frame(Op, hh2);
                E1(gvh,:) = Fr.e1; E2(gvh,:) = Fr.e2; Nrm(gvh,:) = Fr.normal;
            end
            AtomCache.E1 = E1; AtomCache.E2 = E2; AtomCache.Nrm = Nrm;
        end
        if isempty(AtomSeed)
            gv = EigenMat.GlobalVertices{1};
            if isempty(gv); gv = (1:AtomCache.nVerts)'; end
            Vh = Surf.Vertices(gv,:);
            [~, ii] = min(sum((Vh - mean(Vh,1)).^2, 2));
            AtomSeed = gv(ii);
            jTextSeed.setText(num2str(AtomSeed));
        end
        ok = true;
    end

    function RecomputeBase()
        if isempty(AtomCache) || isempty(AtomSeed) || isempty(Frame); return; end
        switch Kind
            case 'scalar',  seedDir = 1;
            case 'vector',  seedDir = [1;0;0];     % base dipole; orientation steers from here
            case 'complex', seedDir = 1;           % base phase (e1 axis); U(1) phase steers from here
            otherwise,      return;
        end
        [A, msg, isErr] = bst_eigenwavelet('Atom', EigenMat, AtomCache.Op, Frame, AtomSeed, seedDir);
        if isErr; bst_error(['Wavelet atom: ' msg], 'Eigenwavelet options', 0); AtomBase = []; return; end
        AtomBase = A;
    end

    function UpdateAtom(varargin)
        if isempty(hFigAtom) || ~ishandle(hFigAtom) || isempty(AtomBase); return; end
        m = max(1, min(M, AtomScale));
        switch Kind
            case 'scalar'
                psf = real(AtomBase(:, 1, m));
                DrawScalarAtom(psf);
            case 'vector'
                As = bst_eigenwavelet('Steer', AtomBase(:, 1, m), AtomQuat, EigenMat);  % exact right-quat steer
                V  = bst_eigenwavelet('ToVec', As, EigenMat);                            % [3nV x 1]
                DrawVectorAtom(V, sprintf('axis %s', upper(AtomAxis)));
            case 'complex'
                z  = exp(1i*AtomPhase) .* AtomBase(:, 1, m);    % U(1) phase steer (exact)
                v3 = real(z).*AtomCache.E1 + imag(z).*AtomCache.E2;   % decode -> 3-D tangent
                DrawVectorAtom(reshape(v3.', [], 1), sprintf('phase %.0f\\circ', AtomPhase*180/pi));
        end
    end

    function DrawScalarAtom(psf)
        Surf = AtomCache.Surf;
        hAx = i_axes();
        [az, el] = view(hAx); cla(hAx);
        patch('Parent', hAx, 'Vertices', Surf.Vertices, 'Faces', Surf.Faces, ...
              'FaceVertexCData', double(psf(:)), 'FaceColor', 'interp', 'EdgeColor', 'none', ...
              'FaceLighting', 'gouraud', 'AmbientStrength', 0.8, 'DiffuseStrength', 0.4, ...
              'SpecularStrength', 0.05, 'ButtonDownFcn', @(h,e) i_axes_click());
        hold(hAx, 'on');
        sp = Surf.Vertices(AtomSeed, :);
        plot3(hAx, sp(1), sp(2), sp(3), 'k.', 'MarkerSize', 22, 'ButtonDownFcn', @(h,e) i_axes_click());
        hold(hAx, 'off');
        i_finish_axes(hAx, az, el, sp);
        m_ = max(abs(double(psf(:)))); if ~(m_>0); m_=1; end
        if min(psf(:)) < -1e-12*m_; colormap(hAx, i_diverging()); caxis(hAx, [-m_ m_]);
        else; colormap(hAx, hot(256)); caxis(hAx, [0 m_]); end
        colorbar(hAx);
        title(hAx, sprintf('Atom: %s  scale %d/%d (\\surd\\lambda=%.1f)  seed %d', ...
            Frame.Family, AtomScale, M, Frame.Centers(AtomScale), AtomSeed), 'Interpreter','tex');
    end

    function DrawVectorAtom(V, orientLabel)
        Surf = AtomCache.Surf;
        V3  = reshape(V, 3, [])';                       % [nV x 3]
        mag = sqrt(sum(V3.^2, 2));
        hAx = i_axes();
        [az, el] = view(hAx); cla(hAx);
        patch('Parent', hAx, 'Vertices', Surf.Vertices, 'Faces', Surf.Faces, ...
              'FaceVertexCData', mag, 'FaceColor', 'interp', 'EdgeColor', 'none', ...
              'FaceLighting', 'gouraud', 'AmbientStrength', 0.85, 'DiffuseStrength', 0.35, ...
              'SpecularStrength', 0.05, 'ButtonDownFcn', @(h,e) i_axes_click());
        hold(hAx, 'on');
        sel = find(mag > 0.15*max(mag));               % the localized packet
        if ~isempty(sel)
            q = Surf.Vertices(sel, :); vv = V3(sel, :);
            sc = 0.010 / max(sqrt(sum(vv.^2,2)) + eps);
            quiver3(hAx, q(:,1),q(:,2),q(:,3), vv(:,1)*sc,vv(:,2)*sc,vv(:,3)*sc, 0, 'k', 'LineWidth', 0.6);
        end
        % orientation GLYPH at the seed (the steered dipole direction)
        sp = Surf.Vertices(AtomSeed, :);
        gdir = V3(AtomSeed, :);  ng = norm(gdir);
        if ng < eps; [~,ip]=max(mag); gdir = V3(ip,:); ng = norm(gdir)+eps; end
        gdir = gdir / ng;  gl = 0.018;
        quiver3(hAx, sp(1),sp(2),sp(3), gdir(1)*gl,gdir(2)*gl,gdir(3)*gl, 0, 'r', 'LineWidth', 3, 'MaxHeadSize', 2);
        plot3(hAx, sp(1), sp(2), sp(3), 'b.', 'MarkerSize', 24);
        hold(hAx, 'off');
        i_finish_axes(hAx, az, el, sp);
        colormap(hAx, parula(256)); m_ = max(mag); if ~(m_>0); m_=1; end; caxis(hAx, [0 m_]); colorbar(hAx);
        title(hAx, sprintf('%s atom: %s  scale %d/%d  seed %d   [%s]', ...
            Variant, Frame.Family, AtomScale, M, AtomSeed, orientLabel), 'Interpreter','tex');
    end

    function hAx = i_axes()
        hAx = findobj(hFigAtom, 'Type', 'axes');
        if isempty(hAx); hAx = axes('Parent', hFigAtom); else; hAx = hAx(1); end
    end

    function i_finish_axes(hAx, az, el, sp)
        axis(hAx, 'equal', 'off', 'tight');
        if all([az el] == [0 90]); view(hAx, sp - mean(AtomCache.Surf.Vertices, 1));
        else; view(hAx, az, el); end
        camlight(hAx, 'headlight'); material(hAx, 'dull');
    end

    function OnSeedText(varargin)
        if ~EnsureCache(); return; end
        v = str2double(char(jTextSeed.getText()));
        if isnan(v) || v < 1 || v > AtomCache.nVerts; jTextSeed.setText(num2str(AtomSeed)); return; end
        AtomSeed = round(v);
        RecomputeBase(); UpdateAtom();
    end

    function ArmPick(varargin)
        if ~jToggleAtom.isSelected() || isempty(hFigAtom) || ~ishandle(hFigAtom)
            bst_error('Open the wavelet atom first (Show wavelet atom).', 'Eigenwavelet options', 0); return;
        end
        AtomPick = true; set(hFigAtom, 'Pointer', 'crosshair');
    end

    function i_axes_click(varargin)
        if ~AtomPick; return; end
        AtomPick = false; set(hFigAtom, 'Pointer', 'arrow');
        hAx = findobj(hFigAtom, 'Type', 'axes'); if isempty(hAx); return; end
        cp = get(hAx(1), 'CurrentPoint'); p = cp(1, :);
        V = AtomCache.Surf.Vertices;
        [~, iv] = min(sum((V - p).^2, 2));
        AtomSeed = iv; jTextSeed.setText(num2str(iv));
        RecomputeBase(); UpdateAtom();
    end

    function i_atom_scroll(e)
        step = -double(e.VerticalScrollCount) * (pi/12);   % 15 deg per notch
        switch Kind
            case 'vector',  AtomQuat = i_qnorm(i_qmul(i_axisangle(AtomAxis, step), AtomQuat));
            case 'complex', AtomPhase = mod(AtomPhase + step, 2*pi);
            otherwise, return;
        end
        UpdateAtom();
    end

    function i_atom_key(e)
        if strcmp(Kind, 'complex')
            switch e.Key
                case 'r',          AtomPhase = 0;
                case 'leftarrow',  AtomPhase = mod(AtomPhase - pi/12, 2*pi);
                case 'rightarrow', AtomPhase = mod(AtomPhase + pi/12, 2*pi);
                otherwise, return;
            end
            UpdateAtom(); return;
        end
        if ~strcmp(Kind, 'vector'); return; end
        switch e.Key
            case 'x', AtomAxis = 'x';
            case 'y', AtomAxis = 'y';
            case 'z', AtomAxis = 'z';
            case 'r', AtomQuat = [1;0;0;0];
            case 'leftarrow',  AtomQuat = i_qnorm(i_qmul(i_axisangle(AtomAxis, -pi/12), AtomQuat));
            case 'rightarrow', AtomQuat = i_qnorm(i_qmul(i_axisangle(AtomAxis,  pi/12), AtomQuat));
            otherwise, return;
        end
        UpdateAtom();
    end

    function i_atom_closed()
        if ~isempty(hFigAtom) && ishandle(hFigAtom); delete(hFigAtom); end
        hFigAtom = [];
        try
            jToggleAtom.setSelected(false);
        catch
        end
    end

    function ButtonCancel_Callback(varargin)
        if ~isempty(hFigAtom) && ishandle(hFigAtom); delete(hFigAtom); end
        gui_hide(panelName);
    end

    function ButtonOk_Callback(varargin)
        if ~isempty(hFigAtom) && ishandle(hFigAtom); delete(hFigAtom); end
        bst_mutex('release', panelName);
    end
end


%% ===== GET PANEL CONTENTS =====
function s = GetPanelContents() %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'EigenwaveletOptions');
    fam  = ctrl.FamilyKeys{max(1, ctrl.jFamily.getSelectedIndex()+1)};
    s = struct();
    s.Method      = 'wavelet';
    s.EigenFile   = ctrl.EigenFile;
    s.KernelName  = fam;                 % bst_eigen 'wavelet' uses KernelName as the frame family
    s.Nf          = double(ctrl.jNf.getValue());
    s.Comment     = char(ctrl.jTextComment.getText());
end


%% ===== internal (file-scope) =====
function kind = i_variant_kind(Variant)
    switch Variant
        case 'Laplace-Beltrami',     kind = 'scalar';
        case 'Dirac',                kind = 'vector';    % quaternion, ambient, right-quat steer
        case 'Connection Laplacian', kind = 'complex';   % tangent, frame-decoded, U(1) phase steer
        otherwise,                   kind = 'deferred';  % face variants (face-centroid display TODO)
    end
end

function [js, jTitle] = i_slider(jParent, constraints, mn, mx, val)
    import javax.swing.*;
    jTitle = [];
    js = JSlider(mn, mx, val);
    js.setPreferredSize(java_scaled('dimension', 60, 22));
    jParent.add(constraints, js);
end

function cmap = i_diverging()
    n = 128; t = linspace(0,1,n)';
    lo = [t, t, ones(n,1)];  hi = [ones(n,1), flipud(t), flipud(t)];
    cmap = [lo; hi(2:end,:)];
end

% --- quaternion helpers (orientation steering) ---
function q = i_axisangle(axis, ang)
    switch axis
        case 'x', u = [1;0;0];
        case 'y', u = [0;1;0];
        otherwise, u = [0;0;1];
    end
    q = [cos(ang/2); sin(ang/2)*u];
end

function r = i_qmul(a, b)
    r = [ a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4); ...
          a(1)*b(2)+a(2)*b(1)+a(3)*b(4)-a(4)*b(3); ...
          a(1)*b(3)-a(2)*b(4)+a(3)*b(1)+a(4)*b(2); ...
          a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+a(4)*b(1) ];
end

function q = i_qnorm(q)
    q = q / max(norm(q), eps);
end
