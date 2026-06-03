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


%% ===== PURE: representative eigenvalue per paired rank =====
function lam = PairedLambda(Eig) %#ok<DEFNU>
    % Assumes CompRank is contiguous 1..K (guaranteed by tess_eigenmodes /
    % in_tess_eigenmodes); a gap would yield mean([]) = NaN for that rank.
    cr  = Eig.CompRank(:);
    val = Eig.Values(:);
    K   = max(cr);
    lam = zeros(K,1);
    for k = 1:K
        lam(k) = mean(val(cr == k));   % L/R near-symmetric -> shared representative lambda
    end
end

%% ===== PURE: paired weights from a spectral kernel =====
function [W, g] = KernelPairedWeights(Eig, name, params) %#ok<DEFNU>
    g   = bst_eigfilter_kernel(name, params);
    lam = PairedLambda(Eig);
    W   = bst_eigfilter_evaluate(g, lam)';    % 1 x Kpaired
end

%% ===== PURE: slider bounds for a kernel parameter (log heuristic) =====
function [smin, smax] = KernelSliderRange(pdef, prange) %#ok<DEFNU>
    if isempty(pdef) || ~isfinite(pdef) || pdef <= 0; pdef = 1; end
    smin = pdef / 100;
    smax = pdef * 100;
    if numel(prange) == 2
        if isfinite(prange(1)); smin = max(smin, prange(1) + eps); end
        if isfinite(prange(2)); smax = min(smax, prange(2)); end
    end
    if smax <= smin; smax = smin * 100; end
end

%% ===== PURE: integer slider position (0..1000) <-> log-interpolated value =====
function v = SliderToValue(islide, smin, smax) %#ok<DEFNU>
    t = min(max(islide,0),1000) / 1000;
    v = smin * (smax/smin)^t;       % geometric (log) interpolation
end
function islide = ValueToSlider(v, smin, smax) %#ok<DEFNU>
    v = min(max(v, smin), smax);
    islide = round(1000 * log(v/smin) / log(smax/smin));
end

%% ===== PURE: delta (impulse) point-spread = Phi*diag(g)*Phi'*M*e_i =====
function ps = DeltaPointSpread(Eig, MassMatrix, g, iVertex) %#ok<DEFNU>
    nV = size(Eig.Vectors, 1);
    iVertex = min(max(round(iVertex), 1), nV);
    e = zeros(nV, 1); e(iVertex) = 1;
    ps = bst_eigenmodes_filter(Eig, e, MassMatrix, 'custom', 'TransferFn', g);
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

    jRadioKernel = gui_component('Radio', jPanelNew, '',  'Kernel', jGroup, '', @(h,ev)Shape_Callback('kernel'));
    % Kernel dropdown (populated from the eigfilter registry)
    jLabelKernel = gui_component('Label', jPanelNew, 'br', 'Kernel:'); %#ok<NASGU>
    kernelNames  = bst_eigfilter_kernel('list');
    jComboKernel = gui_component('ComboBox', jPanelNew, 'hfill', [], {kernelNames}, '', @(h,ev)KernelChanged_Callback());
    % Two generic parameter sliders (relabeled per kernel; integer 0..1000 -> log value)
    jLabelP1 = gui_component('Label', jPanelNew, 'br', 'p1');
    jSliderP1 = JSlider(0, 1000, 500);
    java_setcb(jSliderP1, 'StateChangedCallback', @(h,ev)ParamPreview_Callback(1), 'MouseReleasedCallback', @(h,ev)Param_Callback(1));
    jPanelNew.add('hfill', jSliderP1);
    jReadoutP1 = gui_component('Label', jPanelNew, '', '');
    jLabelP2 = gui_component('Label', jPanelNew, 'br', 'p2');
    jSliderP2 = JSlider(0, 1000, 500);
    java_setcb(jSliderP2, 'StateChangedCallback', @(h,ev)ParamPreview_Callback(2), 'MouseReleasedCallback', @(h,ev)Param_Callback(2));
    jPanelNew.add('hfill', jSliderP2);
    jReadoutP2 = gui_component('Label', jPanelNew, '', '');
    % Display mode + delta vertex (eigenmode-view figures)
    jGroupDisp = ButtonGroup();
    jLabelDisp = gui_component('Label', jPanelNew, 'br', 'Show:'); %#ok<NASGU>
    jRadioSynth = gui_component('Radio', jPanelNew, '', 'Synthesis', jGroupDisp, '', @(h,ev)Display_Callback('synthesis'));
    jRadioDelta = gui_component('Radio', jPanelNew, '', 'Delta',     jGroupDisp, '', @(h,ev)Display_Callback('delta'));
    jRadioSynth.setSelected(1);
    jLabelVtx = gui_component('Label', jPanelNew, 'br', 'Delta vertex:'); %#ok<NASGU>
    jTextVtx  = gui_component('Text', jPanelNew, '', '1', [], '', @(h,ev)Vertex_Callback());
    jTogglePick = gui_component('Toggle', jPanelNew, '', 'Pick', [], 'Click the cortex to place the delta', @(h,ev)Pick_Callback());

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
                  'jRadioGauss',    jRadioGauss, ...
                  'jRadioKernel',   jRadioKernel, ...
                  'jComboKernel',   jComboKernel, ...
                  'jLabelP1',       jLabelP1, ...
                  'jSliderP1',      jSliderP1, ...
                  'jReadoutP1',     jReadoutP1, ...
                  'jLabelP2',       jLabelP2, ...
                  'jSliderP2',      jSliderP2, ...
                  'jReadoutP2',     jReadoutP2, ...
                  'jRadioSynth',    jRadioSynth, ...
                  'jRadioDelta',    jRadioDelta, ...
                  'jTextVtx',       jTextVtx, ...
                  'jTogglePick',    jTogglePick);
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
    if strcmpi(shape, 'kernel')
        SetWeightMode('kernel');
        return;
    end
    SetWeightMode('window');     % any window shape leaves kernel mode
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


%% ===== KERNEL CALLBACKS (read controls -> kernel state verbs) =====
function KernelChanged_Callback()
    ctrl = bst_get('PanelControls', 'EigenModes');
    name = char(ctrl.jComboKernel.getSelectedItem());
    SetKernelName(name);
    RefreshKernelControls();
    UpdateResponsePlot();
end

function Param_Callback(idx)
    ctrl = bst_get('PanelControls', 'EigenModes');
    v = ReadParamSlider(ctrl, idx);
    if ~isempty(v); SetKernelParam(idx, v); end
    UpdateResponsePlot();
end

function ParamPreview_Callback(idx)
    ctrl = bst_get('PanelControls', 'EigenModes');
    sl = ctrl.(sprintf('jSliderP%d', idx));
    if ~sl.getValueIsAdjusting(); return; end
    v = ReadParamSlider(ctrl, idx);
    if ~isempty(v); ctrl.(sprintf('jReadoutP%d', idx)).setText(sprintf('%.4g', v)); end
end

function Display_Callback(mode)
    SetDisplayMode(mode);
end

function Vertex_Callback()
    ctrl = bst_get('PanelControls', 'EigenModes');
    v = str2double(char(ctrl.jTextVtx.getText()));
    if ~isnan(v); SetDeltaVertex(v); end
end

function Pick_Callback()
    global GlobalData;
    GetState();
    ctrl = bst_get('PanelControls', 'EigenModes');
    if ctrl.jTogglePick.isSelected()
        hFig = bst_figures('GetCurrentFigure', '3D');
        if isempty(hFig) || ~ishandle(hFig); ctrl.jTogglePick.setSelected(0); return; end
        setappdata(hFig, 'EigfilterPickBak', get(hFig, 'WindowButtonDownFcn'));
        set(hFig, 'WindowButtonDownFcn', @(h,e) PickVertex(h));
        GlobalData.UserModes.PickFig = hFig;        % remember which figure we hooked
    else
        % Restore the SAME figure we hooked (not whatever is current now)
        hFig = [];
        if isfield(GlobalData.UserModes, 'PickFig'); hFig = GlobalData.UserModes.PickFig; end
        if ~isempty(hFig) && ishandle(hFig) && ~isempty(getappdata(hFig,'EigfilterPickBak'))
            set(hFig, 'WindowButtonDownFcn', getappdata(hFig,'EigfilterPickBak'));
        end
        GlobalData.UserModes.PickFig = [];
    end
end

function PickVertex(hFig)
    st = GetState();
    if isempty(st.CacheEig); return; end
    % Use Brainstorm's ray/surface vertex picker (rejects MRI clicks with AcceptMri=0)
    iv = panel_coordinates('SelectPoint', hFig, 0);
    if isempty(iv); return; end
    ctrl = bst_get('PanelControls', 'EigenModes');
    ctrl.jTextVtx.setText(num2str(iv));
    SetDeltaVertex(iv);
end


%% ===== KERNEL HELPERS: param-slider read + relabel + response plot =====
%% ===== PURE: slider bounds for a kernel param; ok=false if not scalar-adjustable =====
% Robust to params that lack a 'range' field or carry a non-scalar default (e.g.
% ideal's 2-element 'band'): those are flagged non-adjustable so the slider logic
% disables them instead of crashing.
function [smin, smax, ok] = ParamRange(pinfo) %#ok<DEFNU>
    smin = 0; smax = 1;
    ok = isfield(pinfo,'default') && isnumeric(pinfo.default) && isscalar(pinfo.default);
    if ~ok; return; end
    prange = [];
    if isfield(pinfo,'range'); prange = pinfo.range; end
    [smin, smax] = KernelSliderRange(pinfo.default, prange);
end

function v = ReadParamSlider(ctrl, idx)
    v = [];
    st = GetState();
    meta = bst_eigfilter_kernel('info', st.KernelName);
    fn = fieldnames(meta.params);
    if idx > numel(fn); return; end
    [smin, smax, ok] = ParamRange(meta.params.(fn{idx}));
    if ~ok; return; end                          % non-scalar param: not slider-adjustable
    sl = ctrl.(sprintf('jSliderP%d', idx));
    v = SliderToValue(sl.getValue(), smin, smax);
end

function RefreshKernelControls()
    ctrl = bst_get('PanelControls', 'EigenModes');
    if isempty(ctrl) || ~isfield(ctrl,'jComboKernel'); return; end
    st = GetState();
    meta = bst_eigfilter_kernel('info', st.KernelName);
    fn = fieldnames(meta.params);
    for idx = 1:2
        L  = ctrl.(sprintf('jLabelP%d', idx));
        S  = ctrl.(sprintf('jSliderP%d', idx));
        RO = ctrl.(sprintf('jReadoutP%d', idx));
        if idx <= numel(fn)
            [smin, smax, ok] = ParamRange(meta.params.(fn{idx}));
            L.setEnabled(true); L.setText(fn{idx});
            S.setEnabled(ok); RO.setEnabled(ok);
            if ok
                S.setValue(ValueToSlider(st.KernelParams.(fn{idx}), smin, smax));
                RO.setText(sprintf('%.4g', st.KernelParams.(fn{idx})));
            else
                RO.setText('(fixed)');           % non-scalar param uses its default
            end
        else
            L.setEnabled(false); L.setText(sprintf('p%d', idx));
            S.setEnabled(false); RO.setEnabled(false); RO.setText('');
        end
    end
end

function UpdateResponsePlot()
    st = GetState();
    if ~strcmpi(st.WeightMode,'kernel') || isempty(st.CacheEig); return; end
    g = st.KernelFn;
    if isempty(g); g = bst_eigfilter_kernel(st.KernelName, st.KernelParams); end
    meta = bst_eigfilter_kernel('info', st.KernelName);
    view_eigfilter_response(g, st.CacheEig.Values(:), meta.display);
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
    sel = {'jSliderCenter','jTextWidth','jRadioBox','jRadioTaper','jRadioGauss', ...
           'jRadioKernel','jComboKernel','jSliderP1','jSliderP2','jRadioSynth','jRadioDelta','jTextVtx','jTogglePick'};
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
    if isfield(ctrl,'jRadioKernel')
        ctrl.jRadioKernel.setSelected(strcmpi(st.WeightMode,'kernel'));
        ctrl.jRadioSynth.setSelected(strcmpi(st.DisplayMode,'synthesis'));
        ctrl.jRadioDelta.setSelected(strcmpi(st.DisplayMode,'delta'));
        ctrl.jTextVtx.setText(num2str(st.DeltaVertex));
        RefreshKernelControls();
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
            'CacheMass',    [], ...
            'WeightMode',   'window', ...    % 'window' | 'kernel'
            'KernelName',   'heat', ...
            'KernelParams', struct('t',0.01), ...
            'KernelFn',     [], ...
            'DisplayMode',  'synthesis', ... % 'synthesis' | 'delta'
            'DeltaVertex',  1);
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
    GlobalData.UserModes.WeightMode   = 'window';     % fresh surface starts in window mode
    GlobalData.UserModes.DisplayMode  = 'synthesis';  % (not stale kernel/delta from a prior surface)
    GlobalData.UserModes.KernelFn     = [];
    RecomputeWeights();
end

%% ===== STATE: recompute Weights from shape/band/center =====
function RecomputeWeights()
    global GlobalData;
    st = GlobalData.UserModes;
    if strcmpi(st.WeightMode, 'kernel')
        % Kernel weighting: W(k)=g(lambda_paired_k); store the raw-lambda handle too.
        if ~EnsureCache(st.SurfaceFile)
            GlobalData.UserModes.KernelFn = [];   % no stale handle if cache unavailable
            return;
        end
        Eig = GlobalData.UserModes.CacheEig;
        [W, g] = KernelPairedWeights(Eig, st.KernelName, st.KernelParams);
        GlobalData.UserModes.Weights  = W;
        GlobalData.UserModes.KernelFn = g;
    else
        GlobalData.UserModes.Weights = BuildWeights(st.WindowShape, ...
            st.Band(1), st.Band(2), st.iCurrentMode, st.nModes);
        GlobalData.UserModes.KernelFn = [];
    end
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
        GlobalData.UserModes.Band     = [c, c];
        GlobalData.UserModes.BandSpan = 0;   % zero width so the slider stays single
    end
    RecomputeWeights();
    NotifyChanged();
end

%% ===== STATE: weight mode ('window' | 'kernel') =====
function SetWeightMode(mode) %#ok<DEFNU>
    global GlobalData;
    GetState();
    GlobalData.UserModes.WeightMode = lower(mode);
    if ~strcmpi(mode, 'kernel')
        view_eigfilter_response('close');   % no kernel response to show in window mode
    end
    RecomputeWeights();
    NotifyChanged();
end

%% ===== STATE: select a spectral kernel (resets params to defaults) =====
function SetKernelName(name) %#ok<DEFNU>
    global GlobalData;
    GetState();
    meta = bst_eigfilter_kernel('info', name);
    p = struct();
    fn = fieldnames(meta.params);
    for i = 1:numel(fn); p.(fn{i}) = meta.params.(fn{i}).default; end
    GlobalData.UserModes.KernelName   = name;
    GlobalData.UserModes.KernelParams = p;
    GlobalData.UserModes.WeightMode   = 'kernel';
    RecomputeWeights();
    NotifyChanged();
end

%% ===== STATE: set one kernel parameter by index =====
function SetKernelParam(idx, value) %#ok<DEFNU>
    global GlobalData;
    GetState();
    fn = fieldnames(GlobalData.UserModes.KernelParams);
    if idx > numel(fn); return; end
    pNew = GlobalData.UserModes.KernelParams;
    pNew.(fn{idx}) = value;
    % Validate before committing: some kernels have cross-parameter constraints
    % (e.g. dog requires t1 < t2). On violation, keep the prior valid parameters and
    % snap the slider back instead of throwing an uncaught error.
    try
        bst_eigfilter_kernel(GlobalData.UserModes.KernelName, pNew);
    catch
        RefreshKernelControls();   % revert the slider to the last valid value (GUI only)
        return;
    end
    GlobalData.UserModes.KernelParams = pNew;
    RecomputeWeights();
    NotifyChanged();
end

%% ===== STATE: display mode for eigenmode-view ('synthesis' | 'delta') =====
function SetDisplayMode(mode) %#ok<DEFNU>
    global GlobalData; GetState();
    GlobalData.UserModes.DisplayMode = lower(mode);
    NotifyChanged();
end

%% ===== STATE: delta (impulse) source vertex =====
function SetDeltaVertex(v) %#ok<DEFNU>
    global GlobalData; GetState();
    GlobalData.UserModes.DeltaVertex = max(round(v),1);
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
    if strcmpi(st.WeightMode, 'kernel') && ~isempty(st.KernelFn)
        % Kernel: apply g at each RAW mode's lambda (exact), not per paired rank.
        uF = bst_eigenmodes_filter(Eig, u, M, 'custom', 'TransferFn', st.KernelFn);
    else
        wRaw = W(CompRank);                    % expand paired -> raw columns
        uF = bst_eigenmodes_filter(Eig, u, M, 'custom', 'TransferFn', @(l) wRaw(:));
    end
end


%% ===== DISPLAY: column for an eigenmode-view figure (synthesis or delta) =====
function col = GetDisplayColumn(SurfaceFile, PairedGrid) %#ok<DEFNU>
    global GlobalData;
    st = GetState();
    W  = st.Weights;
    if strcmpi(st.DisplayMode, 'delta')
        if ~EnsureCache(SurfaceFile)
            col = PairedGrid * W(:); return;   % fallback if cache missing
        end
        Eig = GlobalData.UserModes.CacheEig;
        M   = GlobalData.UserModes.CacheMass;
        if strcmpi(st.WeightMode, 'kernel') && ~isempty(st.KernelFn)
            tf = st.KernelFn;
        else
            wRaw = W(Eig.CompRank(:));
            tf = @(l) wRaw(:);
        end
        col = DeltaPointSpread(Eig, M, tf, st.DeltaVertex);
    else
        col = PairedGrid * W(:);
    end
end
