function varargout = panel_eigenmodes(varargin)
% PANEL_EIGENMODES: Eigenmode scale lever — a stateful, surface-scoped,
% live-broadcast selection in eigenmode (spatial-frequency) space. Changing the
% selection live-coarsens/smooths a displayed source map via a non-destructive
% display-time band-limited reconstruction.
%
% USAGE:  bstPanel = panel_eigenmodes('CreatePanel')
%         W  = panel_eigenmodes('BuildWeights', shape, kLo, kHi, iCenter, K)
%         panel_eigenmodes('SetBand', kLo, kHi)
%         panel_eigenmodes('SetCurrentMode', k)
%         panel_eigenmodes('SetWindowShape', shape)
%         panel_eigenmodes('SetActive', isActive)
%         W  = panel_eigenmodes('GetWeights')
%         uF = panel_eigenmodes('ApplyToColumn', SurfaceFile, u)
%         panel_eigenmodes('UpdatePanel', hFig)
%
% SEE ALSO: bst_eigenmodes_filter, in_tess_eigenmodes, panel_freq, panel_surface

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


%% ===== PURE: clamped band from a center + half-width =====
function [lo, hi] = BandFromCenterWidth(center, width, K) %#ok<DEFNU>
    center = min(max(round(center), 1), K);
    width  = max(round(width), 0);
    lo = min(max(center - width, 1), K);
    hi = min(max(center + width, 1), K);
end


%% ===== PURE: build the [1 x K] weight vector from a window shape =====
function W = BuildWeights(shape, kLo, kHi, iCenter, K) %#ok<DEFNU>
    % Clamp all indices into [1,K]
    kLo     = min(max(round(kLo),     1), K);
    kHi     = min(max(round(kHi),     1), K);
    iCenter = min(max(round(iCenter), 1), K);
    if (kHi < kLo)
        [kLo, kHi] = deal(kHi, kLo);
    end
    W = zeros(1, K);
    switch lower(shape)
        case 'single'
            W(iCenter) = 1;
        case 'box'
            W(kLo:kHi) = 1;
        case 'tapered'
            % Tukey window over [kLo,kHi]: cosine shoulders taper to 0 at kLo and kHi;
            % flat interior (=1) over the central portion of the band.
            n = kHi - kLo + 1;
            if (n <= 2)
                W(kLo:kHi) = 1;
            else
                r = 0.5;    % total taper fraction (each shoulder spans r/2 = 25% of the window)
                t = linspace(0, 1, n);
                wt = ones(1, n);
                edge = (r/2);
                iL = t < edge;
                iR = t > (1 - edge);
                wt(iL) = 0.5 * (1 + cos(pi * (2*t(iL)/r - 1)));
                wt(iR) = 0.5 * (1 + cos(pi * (2*t(iR)/r - 2/r + 1)));
                W(kLo:kHi) = wt;
            end
        case 'gain'
            % Gaussian bell centered at iCenter; sigma = half the band width.
            sigma = max((kHi - kLo) / 2, 1);
            k = 1:K;
            W = exp(-0.5 * ((k - iCenter) / sigma).^2);
        otherwise
            error('panel_eigenmodes:BuildWeights: unknown shape "%s".', shape);
    end
end


%% ===== CREATE PANEL =====
function bstPanelNew = CreatePanel() %#ok<DEFNU>
    panelName = 'EigenModes';
    import java.awt.*;
    import javax.swing.*;

    jPanelNew = gui_river([2,2], [4,4,6,6], 'Spatial scale (eigenmodes)');

    % Active toggle + readout (same row)
    jCheckActive = gui_component('Checkbox', jPanelNew, '', 'Active', [], ...
        'Live-filter the displayed source map', @(h,ev)CheckActive_Callback());
    jLabelReadout = gui_component('Label', jPanelNew, 'hfill', '');
    jLabelReadout.setHorizontalAlignment(JLabel.RIGHT);

    % Center mode slider
    jLabelCenter = gui_component('Label', jPanelNew, 'br', 'Center mode');
    jSliderCenter = JSlider(1, 100, 1);
    jSliderCenter.setPaintLabels(1);
    java_setcb(jSliderCenter, 'StateChangedCallback',  @(h,ev)CenterPreview_Callback());
    java_setcb(jSliderCenter, 'MouseReleasedCallback', @(h,ev)Center_Callback());
    jPanelNew.add('br hfill', jSliderCenter);

    % Width text field
    jLabelWidth = gui_component('Label', jPanelNew, 'br', 'Width (+/- modes):');
    jTextWidth  = gui_component('Text', jPanelNew, '', '0', [], '', @(h,ev)Width_Callback());
    try
        jTextWidth.setColumns(4);
    catch
        % setColumns not available on this Swing wrapper — safe to skip
    end

    % Window shape radios
    jGroup = ButtonGroup();
    jLabelShape = gui_component('Label', jPanelNew, 'br', 'Shape:');
    jRadioBox   = gui_component('Radio', jPanelNew, '',  'Box',   jGroup, '', @(h,ev)Shape_Callback('box'));
    jRadioTaper = gui_component('Radio', jPanelNew, '',  'Taper', jGroup, '', @(h,ev)Shape_Callback('tapered'));
    jRadioGauss = gui_component('Radio', jPanelNew, '',  'Gauss', jGroup, '', @(h,ev)Shape_Callback('gain'));
    jRadioBox.setSelected(1);

    ctrl = struct('jPanelTop',      jPanelNew, ...
                  'jCheckActive',   jCheckActive, ...
                  'jLabelReadout',  jLabelReadout, ...
                  'jLabelCenter',   jLabelCenter, ...
                  'jSliderCenter',  jSliderCenter, ...
                  'jLabelWidth',    jLabelWidth, ...
                  'jTextWidth',     jTextWidth, ...
                  'jLabelShape',    jLabelShape, ...
                  'jRadioBox',      jRadioBox, ...
                  'jRadioTaper',    jRadioTaper, ...
                  'jRadioGauss',    jRadioGauss);
    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);
end


%% ===== SLIDER AXIS LABELS (with [ ] band markers) =====
function ApplySliderLabels(jSlider, K, lo, hi, width)
    import javax.swing.JLabel;
    tbl = java.util.Hashtable();
    % Window markers first
    if (width > 0)
        tbl.put(java.lang.Integer(lo), JLabel('['));
        tbl.put(java.lang.Integer(hi), JLabel(']'));
    end
    % Eigenmode axis labels (do not clobber a bracket position)
    pos = unique(round(linspace(1, K, 5)));
    for p = pos
        if (width > 0) && (p == lo || p == hi)
            continue;
        end
        tbl.put(java.lang.Integer(p), JLabel(num2str(p)));
    end
    jSlider.setLabelTable(tbl);
    jSlider.setPaintLabels(1);
end


%% ===== CALLBACKS (read controls -> state verbs) =====
function CheckActive_Callback()
    ctrl = bst_get('PanelControls', 'EigenModes');
    SetActive(ctrl.jCheckActive.isSelected());
end

function Center_Callback()
    ctrl = bst_get('PanelControls', 'EigenModes');
    st = GetState();
    c = ctrl.jSliderCenter.getValue();
    w = ReadWidth(ctrl, st.nModes);
    [lo, hi] = BandFromCenterWidth(c, w, st.nModes);
    SetBand(lo, hi);
end

function CenterPreview_Callback()
    ctrl = bst_get('PanelControls', 'EigenModes');
    if ~ctrl.jSliderCenter.getValueIsAdjusting(), return; end
    st = GetState();
    c = ctrl.jSliderCenter.getValue();
    w = ReadWidth(ctrl, st.nModes);
    [lo, hi] = BandFromCenterWidth(c, w, st.nModes);
    ApplySliderLabels(ctrl.jSliderCenter, st.nModes, lo, hi, w);
end

function Width_Callback()
    ctrl = bst_get('PanelControls', 'EigenModes');
    st = GetState();
    w = ReadWidth(ctrl, st.nModes);
    c = ctrl.jSliderCenter.getValue();
    [lo, hi] = BandFromCenterWidth(c, w, st.nModes);
    if (w == 0)
        SetWindowShape('single');
    else
        SetWindowShape(CurrentShape(ctrl));
    end
    SetBand(lo, hi);
end

function Shape_Callback(shape)
    ctrl = bst_get('PanelControls', 'EigenModes');
    w = ReadWidth(ctrl, GetState().nModes);
    if (w == 0)
        SetWindowShape('single');
    else
        SetWindowShape(shape);
    end
end

function w = ReadWidth(ctrl, K)
    s = char(ctrl.jTextWidth.getText());
    w = str2double(s);
    if isnan(w) || ~isreal(w)
        st = GetState(); w = round((st.Band(2) - st.Band(1)) / 2);
    end
    w = min(max(round(w), 0), max(K-1, 0));
end

function shape = CurrentShape(ctrl)
    if ctrl.jRadioTaper.isSelected()
        shape = 'tapered';
    elseif ctrl.jRadioGauss.isSelected()
        shape = 'gain';
    else
        shape = 'box';
    end
end


%% ===== PURE: panel context from front-figure facts =====
function c = ClassifyContext(facts) %#ok<DEFNU>
    c = struct('kind','none','selectEnabled',false,'activeEnabled',false);
    if facts.isEigenView
        c.kind='view';   c.selectEnabled=true;  c.activeEnabled=false;
    elseif facts.hasSourceModes
        c.kind='source'; c.selectEnabled=true;  c.activeEnabled=true;
    end
end


%% ===== UPDATE PANEL: populate/enable from the active figure's surface =====
function UpdatePanel(hFig) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'EigenModes');
    if isempty(ctrl), return; end
    if (nargin < 1) || isempty(hFig) || ~ishandle(hFig)
        hFig = bst_figures('GetCurrentFigure', '3D');
    end
    isEigenView = ~isempty(hFig) && ishandle(hFig) && ~isempty(getappdata(hFig, 'EigenView'));
    SurfaceFile = '';
    if isEigenView
        ev = getappdata(hFig, 'EigenView'); SurfaceFile = ev.SurfaceFile;
    elseif ~isempty(hFig) && ishandle(hFig)
        SurfaceFile = GetFigureSurfaceWithModes(hFig);
    end
    facts = struct('isEigenView', isEigenView, 'hasSourceModes', ~isEigenView && ~isempty(SurfaceFile));
    c = ClassifyContext(facts);
    % Selection controls + the Active toggle, gated by context
    SetSelectEnabled(ctrl, c.selectEnabled);
    ctrl.jCheckActive.setEnabled(c.activeEnabled);
    ctrl.jCheckActive.setVisible(c.activeEnabled);
    if ~strcmp(c.kind, 'source')
        SetActive(0);                          % no stale filtering off-source
    end
    if strcmp(c.kind, 'none')
        ctrl.jLabelReadout.setText('no eigenmode view');
        return;
    end
    [Eig, ~] = in_tess_eigenmodes(SurfaceFile);
    K = double(max(Eig.CompRank));
    st = GetState();
    if ~file_compare(st.SurfaceFile, SurfaceFile) || (st.nModes ~= K)
        ResetState(SurfaceFile, K);
        if isEigenView
            SetWindowShape('single'); SetCurrentMode(1);
        end
    end
    ctrl.jSliderCenter.setMaximum(K);
    RefreshControls();
end


%% ===== helpers =====
function SetSelectEnabled(ctrl, isOn)
    sel = {'jSliderCenter','jTextWidth','jRadioBox','jRadioTaper','jRadioGauss'};
    for i = 1:numel(sel)
        if isfield(ctrl, sel{i}) && isa(ctrl.(sel{i}), 'javax.swing.JComponent')
            ctrl.(sel{i}).setEnabled(logical(isOn));
        end
    end
end


%% ===== Reflect state back into the controls + readout =====
function RefreshControls()
    ctrl = bst_get('PanelControls', 'EigenModes');
    % Guard: panel may be unregistered, or a stale instance without the new controls.
    if isempty(ctrl) || ~isfield(ctrl, 'jSliderCenter'), return; end
    st = GetState();
    K  = max(st.nModes, 1);
    lo = st.Band(1); hi = st.Band(2);
    width = round((hi - lo) / 2);
    ctrl.jSliderCenter.setMaximum(K);
    ctrl.jSliderCenter.setValue(st.iCurrentMode);
    ApplySliderLabels(ctrl.jSliderCenter, K, lo, hi, width);
    ctrl.jTextWidth.setText(num2str(width));
    isBand = (width > 0) && ~strcmpi(st.WindowShape, 'single');
    ctrl.jRadioBox.setEnabled(isBand);
    ctrl.jRadioTaper.setEnabled(isBand);
    ctrl.jRadioGauss.setEnabled(isBand);
    % Sync the toggle + shape radios back from state (e.g. after ResetState)
    ctrl.jCheckActive.setSelected(logical(st.isActive));
    switch st.WindowShape
        case 'tapered', ctrl.jRadioTaper.setSelected(true);
        case 'gain',    ctrl.jRadioGauss.setSelected(true);
        otherwise,      ctrl.jRadioBox.setSelected(true);
    end
    nKeep = nnz(st.Weights > 1e-6);
    lamStr = '';
    if ~isempty(st.CacheEig) && isfield(st.CacheEig,'Values') && ~isempty(st.CacheEig.Values) ...
            && isfield(st.CacheEig,'CompRank')
        cr = st.CacheEig.CompRank(:); vals = st.CacheEig.Values(:);
        iLo = find(cr == lo, 1); iHi = find(cr == hi, 1);
        if ~isempty(iLo) && ~isempty(iHi)
            if (lo == hi)
                lamStr = sprintf('    lambda %.3g', vals(iLo));
            else
                lamStr = sprintf('    lambda %.3g - %.3g', vals(iLo), vals(iHi));
            end
        end
    end
    if (width == 0)
        ctrl.jLabelReadout.setText(sprintf('Mode %d / %d%s', st.iCurrentMode, K, lamStr));
    else
        ctrl.jLabelReadout.setText(sprintf('Keeping %d modes%s', nKeep, lamStr));
    end
end


%% ===== helpers =====
function SurfaceFile = GetFigureSurfaceWithModes(hFig)
    SurfaceFile = '';
    if isempty(hFig) || ~ishandle(hFig)
        return;
    end
    TessInfo = getappdata(hFig, 'Surface');
    if isempty(TessInfo)
        return;
    end
    for iTess = 1:numel(TessInfo)
        sf = TessInfo(iTess).SurfaceFile;
        if isempty(sf) || ~isfield(TessInfo(iTess), 'DataSource') ...
                || isempty(TessInfo(iTess).DataSource) ...
                || isempty(TessInfo(iTess).DataSource.FileName)
            continue;
        end
        [~, isComputed] = in_tess_eigenmodes(sf);
        if isComputed
            SurfaceFile = sf;
            return;
        end
    end
end


%% ===== STATE: lazy default =====
function st = GetState()
    global GlobalData;
    if ~isfield(GlobalData, 'UserModes') || isempty(GlobalData.UserModes) ...
            || ~isfield(GlobalData.UserModes, 'Weights')
        GlobalData.UserModes = struct(...
            'SurfaceFile',  '', ...
            'nModes',       0,  ...
            'iCurrentMode', 1,  ...
            'Weights',      [], ...
            'WindowShape',  'box', ...
            'Band',         [1 1], ...
            'BandSpan',     0, ...
            'isActive',     0,  ...
            'CacheSurfaceFile', '', ...
            'CacheEig',     [], ...
            'CacheMass',    []);
    end
    st = GlobalData.UserModes;
end

%% ===== STATE: reset for a surface with K modes =====
function ResetState(SurfaceFile, K) %#ok<DEFNU>
    global GlobalData;
    GetState();                          % ensure struct exists
    GlobalData.UserModes.SurfaceFile  = SurfaceFile;
    GlobalData.UserModes.nModes       = K;
    GlobalData.UserModes.Band         = [1, min(30, K)];
    GlobalData.UserModes.BandSpan     = GlobalData.UserModes.Band(2) - GlobalData.UserModes.Band(1);
    GlobalData.UserModes.iCurrentMode = round(mean(GlobalData.UserModes.Band));
    GlobalData.UserModes.WindowShape  = 'box';
    GlobalData.UserModes.isActive     = 0;
    RecomputeWeights();
end

%% ===== STATE: recompute Weights from shape/band/center =====
function RecomputeWeights()
    global GlobalData;
    st = GlobalData.UserModes;
    GlobalData.UserModes.Weights = BuildWeights(st.WindowShape, ...
        st.Band(1), st.Band(2), st.iCurrentMode, st.nModes);
end

%% ===== STATE: set band (clamped); recentre iCurrentMode to band midpoint =====
function SetBand(kLo, kHi) %#ok<DEFNU>
    global GlobalData;
    st = GetState();
    K  = st.nModes;
    kLo = min(max(round(kLo), 1), K);
    kHi = min(max(round(kHi), 1), K);
    if (kHi < kLo), [kLo, kHi] = deal(kHi, kLo); end
    GlobalData.UserModes.Band         = [kLo, kHi];
    GlobalData.UserModes.BandSpan     = kHi - kLo;
    GlobalData.UserModes.iCurrentMode = round((kLo + kHi) / 2);
    RecomputeWeights();
    NotifyChanged();
end

%% ===== STATE: move center (coupled — slide band, preserve width) =====
function SetCurrentMode(k) %#ok<DEFNU>
    global GlobalData;
    st = GetState();
    K  = st.nModes;
    k  = min(max(round(k), 1), K);
    halfW = round(st.BandSpan / 2);
    kLo = min(max(k - halfW, 1), K);
    kHi = min(max(k + halfW, 1), K);
    GlobalData.UserModes.iCurrentMode = k;
    GlobalData.UserModes.Band         = [kLo, kHi];
    RecomputeWeights();
    NotifyChanged();
end

%% ===== STATE: window shape ('single' collapses band to center) =====
function SetWindowShape(shape) %#ok<DEFNU>
    global GlobalData;
    GetState();
    GlobalData.UserModes.WindowShape = lower(shape);
    if strcmpi(shape, 'single')
        c = GlobalData.UserModes.iCurrentMode;
        GlobalData.UserModes.Band = [c, c];
    end
    RecomputeWeights();
    NotifyChanged();
end

%% ===== STATE: active toggle =====
function SetActive(isActive) %#ok<DEFNU>
    global GlobalData;
    GetState();
    isActive = logical(isActive);
    if isActive == GlobalData.UserModes.isActive
        return;                      % no change -> no broadcast
    end
    GlobalData.UserModes.isActive = isActive;
    NotifyChanged();
end

%% ===== STATE: read canonical weights =====
% Precondition: ResetState must have been called; returns [] on uninitialized state.
function W = GetWeights() %#ok<DEFNU>
    st = GetState();
    W = st.Weights;
end

%% ===== STATE: read the current center mode index =====
function k = GetCurrentMode() %#ok<DEFNU>
    st = GetState();
    k = st.iCurrentMode;
end

%% ===== QUERY: is the lever actively filtering this surface? =====
function tf = IsActive(SurfaceFile) %#ok<DEFNU>
    st = GetState();
    tf = logical(st.isActive) && ~isempty(SurfaceFile) ...
         && file_compare(st.SurfaceFile, SurfaceFile);
end

%% ===== Broadcast: sync the panel (if shown) + repaint affected figures =====
% Guarded so the state verbs stay callable headlessly (unit tests run without a
% running Brainstorm and without toolbox/core on the path).
function NotifyChanged()
    global GlobalData;
    % Only query PanelControls when the Brainstorm GUI is fully initialised
    % (GlobalData.Program.GUI must exist — it does not in headless unit tests).
    if exist('bst_get', 'file') ...
            && ~isempty(GlobalData) ...
            && isfield(GlobalData, 'Program') ...
            && isfield(GlobalData.Program, 'GUI')
        ctrl = bst_get('PanelControls', 'EigenModes');
        if ~isempty(ctrl)
            RefreshControls();
        end
    end
    if exist('bst_figures', 'file')
        bst_figures('FireModesChanged');
    end
end


%% ===== CACHE: store eigenmodes + mass matrix for a surface =====
function SetCache(SurfaceFile, Eig, MassMatrix) %#ok<DEFNU>
    global GlobalData;
    GetState();
    GlobalData.UserModes.CacheSurfaceFile = SurfaceFile;
    GlobalData.UserModes.CacheEig         = Eig;
    GlobalData.UserModes.CacheMass        = MassMatrix;
end

%% ===== CACHE: ensure eigenmodes + mass are loaded for a surface =====
function isOk = EnsureCache(SurfaceFile)
    global GlobalData;
    st = GetState();
    isOk = false;
    if ~isempty(st.CacheEig) && ~isempty(st.CacheMass) ...
            && file_compare(st.CacheSurfaceFile, SurfaceFile)
        isOk = true;
        return;
    end
    [Eig, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if ~isComputed || isempty(Eig) || ~isfield(Eig, 'Vectors') || isempty(Eig.Vectors)
        return;
    end
    % Prefer the mass matrix stored with the eigenmodes; recompute only if absent.
    if isfield(Eig, 'MassMatrix') && ~isempty(Eig.MassMatrix)
        M = Eig.MassMatrix;
    else
        sSurf = in_tess_bst(SurfaceFile, 0);
        [~, M] = tess_laplacian(sSurf.Vertices, sSurf.Faces, 'MassType', Eig.MassType);
    end
    SetCache(SurfaceFile, Eig, M);
    isOk = true;
end

%% ===== APPLY: filter a displayed source column (guarded, non-destructive) =====
function uF = ApplyToColumn(SurfaceFile, u) %#ok<DEFNU>
    global GlobalData;
    uF = u;                                   % default: unchanged
    st = GetState();
    % Guards: lever off, surface mismatch, empty column
    if ~st.isActive || isempty(u) || isempty(SurfaceFile) ...
            || ~file_compare(st.SurfaceFile, SurfaceFile)
        return;
    end
    if ~EnsureCache(SurfaceFile)
        return;
    end
    Eig = GlobalData.UserModes.CacheEig;
    M   = GlobalData.UserModes.CacheMass;
    % Scalar-field guard: only filter when the column matches the mesh vertex
    % count (skip unconstrained/volume/mismatched maps rather than mis-filter).
    if (size(u,1) ~= size(Eig.Vectors,1)) || (size(u,2) ~= 1)
        return;
    end
    CompRank = Eig.CompRank(:);
    Kpaired  = max(CompRank);
    W = GlobalData.UserModes.Weights;
    if isempty(W) || (numel(W) ~= Kpaired)
        return;
    end
    wRaw = W(CompRank);                         % expand paired -> raw columns
    uF = bst_eigenmodes_filter(Eig, u, M, 'custom', 'TransferFn', @(l) wRaw(:));
end
