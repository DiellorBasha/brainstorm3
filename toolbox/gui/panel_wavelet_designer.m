function varargout = panel_wavelet_designer(varargin)
% PANEL_FILTER_DESIGNER: Transient control panel for the spatial filterbank designer.
%
% Owns the design state and the live-preview controller for a filterbank design
% session. Shown by view_wavelet_designer via gui_show(...,'BrainstormTab',...) and
% torn down by gui_hide('WaveletDesigner'). Java control handles are stored in the
% panel's sControls (bst_get('PanelControls','WaveletDesigner')); the mutable MATLAB
% session state lives in the preview figure's appdata ('WaveletDesignerState'), which
% the orchestrator (view_wavelet_designer) creates and links via ctxFn callbacks.
%
% Workflow: the user designs ONE wavelet (kernel + scale + direction + chirality) and
% sees it live on the cortex; optionally they then tile the spectrum into a small bank
% of that wavelet at several scales. The Dirac structure exposes three filter axes:
% SCALE g(lambda) (shared eigfilter library, driven here in mode-index space), DIRECTION
% (the ambient 3D seed vector, set by azimuth/elevation sliders), and CHIRALITY (helicity
% projector). Source vectors are full 3D ambient fields (never normal/tangent components).
%
% Dispatched subfunctions (panel_wavelet_designer('Name', args...)):
%   CreatePanel(EigenMat, EigenFile, hFig, ctxFn) -> bstPanel
%   SetSeedVertex(panelName, iVertex)   seed a delta at a clicked vertex
%   GetState(panelName) -> S            the current session-state struct
%   OnSelectTile(panelName, j)          make tile j the displayed tile
%   OnSave(panelName) / OnCancel(panelName)
%
% SEE ALSO: view_wavelet_designer, bst_filterbank_tiles, bst_dirac_eigenmodes_filter,
%           bst_eigfilter_kernel, view_eigfilter_response
%
% Authors: Diellor Basha, 2026

    eval(macro_method);
end


%% ===== CREATE PANEL =====
function bstPanelNew = CreatePanel(EigenMat, EigenFile, hFig, ctxFn, Frame) %#ok<DEFNU>
    import java.awt.*;
    import javax.swing.*;
    panelName = 'WaveletDesigner';
    isDirac = strcmpi(EigenMat.Variant, 'Dirac');

    % Responsive container (mirrors panel_surface / panel_scout): a BorderLayout root
    % whose NORTH holds a BoxLayout column of titled river sub-panels. Each river resizes
    % to the tab width, and every slider is capped by a trailing label so 'hfill' is
    % bounded by the panel edge instead of overflowing.
    jPanelNew = gui_component('Panel');
    jOpt = JPanel(); jOpt.setLayout(BoxLayout(jOpt, BoxLayout.Y_AXIS));

    % ===== SECTION 1: INPUT (the seed / source) =====
    % The seed is the source the wavelet is built from. For Laplace-Beltrami a vertex
    % delta suffices; the Dirac operator acts on full 3D vector fields, so the input is
    % richer: a delta gets a 3D ambient DIRECTION, and (optionally) a CHIRALITY.
    jSec1 = gui_river([2 2], [2 8 3 6], '1. Input');
    gui_component('label', jSec1, 'br', '<HTML><I>Click a vertex on the cortex to place the wavelet.</I>');

    jAz = []; jEl = []; jAzVal = []; jElVal = [];
    jChirNone = []; jChirPlus = []; jChirMinus = [];
    if isDirac
        % seed direction: azimuth + elevation -> full ambient 3D unit vector
        gui_component('label', jSec1, 'br', '<HTML><I>Seed direction (local frame)</I>');
        [jAz, jAzVal] = i_labeled_slider(jSec1, '<HTML>In-plane angle: 0&deg;',  '0',  '360', 0, 360, 0);
        [jEl, jElVal] = i_labeled_slider(jSec1, '<HTML>Tilt to normal: 90&deg;', '-90','90', -90, 90, 90);
        % chirality (helicity of the reconstructed vector field; None = real field)
        % radios on their own row so '- left' is not clipped at the panel edge
        gui_component('label', jSec1, 'br', 'Chirality:');
        jChirNone  = gui_component('radio', jSec1, 'br tab', 'None');
        jChirPlus  = gui_component('radio', jSec1, '', '+ right');
        jChirMinus = gui_component('radio', jSec1, '', '- left');
        jChirNone.setSelected(true);
        grpCh = ButtonGroup(); grpCh.add(jChirNone); grpCh.add(jChirPlus); grpCh.add(jChirMinus);
    end
    jOpt.add(jSec1);

    % ===== SECTION 2: FILTER KERNEL =====
    [keys, displays] = i_kernel_list();
    jSec2 = gui_river([2 2], [2 8 3 6], '2. Filter kernel');
    gui_component('label', jSec2, 'br', 'Kernel:');
    jKernel = gui_component('combobox', jSec2, 'br hfill', [], {displays}, [], [], []);
    jParams = gui_river([2 2], [0 2 0 2]);        % scale sliders, rebuilt per kernel
    jSec2.add('br hfill', jParams);
    jOpt.add(jSec2);

    % ===== SECTION 3: SPECTRUM TILING (optional) =====
    jSec3 = gui_river([2 2], [2 8 3 6], '3. Spectrum tiling (optional)');
    [jTiles, jTilesVal] = i_labeled_slider(jSec3, 'Tiles: 1 (single wavelet)', '1', '12', 1, 12, 1);
    jChiSplit = gui_component('checkbox', jSec3, 'br', 'Cross with both chiralities');
    if ~isDirac; jChiSplit.setEnabled(false); end
    jActiveTile = gui_component('label', jSec3, 'br', '');
    jOpt.add(jSec3);

    % ===== ACTIONS =====
    jSec4 = gui_river([2 2], [2 8 6 6]);
    jSave   = gui_component('button', jSec4, 'br right', 'Save bank');
    jCancel = gui_component('button', jSec4, '', 'Cancel');
    jOpt.add(jSec4);

    jPanelNew.add(jOpt, BorderLayout.NORTH);

    % --- collect Java handles ---
    ctrl = struct('jKernel',jKernel, 'KernelKeys',{keys}, 'jParams',jParams, ...
                  'jAz',jAz, 'jEl',jEl, 'jAzVal',jAzVal, 'jElVal',jElVal, ...
                  'jChirNone',jChirNone, 'jChirPlus',jChirPlus, 'jChirMinus',jChirMinus, ...
                  'jTiles',jTiles, 'jTilesVal',jTilesVal, 'jChiSplit',jChiSplit, ...
                  'jActiveTile',jActiveTile, ...
                  'jSave',jSave, 'jCancel',jCancel, 'hFig',hFig);

    % --- session state on the preview figure ---
    [~, iSubject] = bst_get('EigenFile', EigenFile);
    if isempty(iSubject); iSubject = []; end
    S = struct('EigenMat',EigenMat, 'EigenFile',file_short(EigenFile), 'iSubject',iSubject, ...
               'Variant',EigenMat.Variant, 'isDirac',isDirac, 'ctxFn',ctxFn, ...
               'Op',load(file_fullpath(EigenMat.OperatorFile)), ...
               'Lambda',double(EigenMat.Lambda{1}(:)), ...
               'FrameU',Frame.U, 'FrameV',Frame.V, 'FrameN',Frame.N, ...
               'SeedCoeffs',[], 'ActiveTile',1, 'iVertex',[], 'Tiles',[], 'ParamNames',{{}});
    setappdata(hFig, 'WaveletDesignerState', S);

    % --- build the initial scale sliders for the default kernel ---
    BuildParamWidgets(ctrl, hFig);

    % --- wire callbacks ---
    % JSlider: StateChangedCallback fires continuously during a drag; we update the cheap
    % text label on every tick and run the expensive Refresh/Reseed only when the drag
    % settles (~getValueIsAdjusting). MouseReleasedCallback on a JSlider is unreliable.
    java_setcb(jKernel,   'ActionPerformedCallback', @(h,e) OnKernelChanged(panelName));
    java_setcb(jChiSplit, 'ActionPerformedCallback', @(h,e) Refresh(panelName));
    java_setcb(jTiles,    'StateChangedCallback',    @(h,e) OnTilesSlider(panelName, h));
    if isDirac
        java_setcb(jAz, 'StateChangedCallback', @(h,e) OnDirSlider(panelName, h));
        java_setcb(jEl, 'StateChangedCallback', @(h,e) OnDirSlider(panelName, h));
        java_setcb(jChirNone,  'ActionPerformedCallback', @(h,e) Refresh(panelName));
        java_setcb(jChirPlus,  'ActionPerformedCallback', @(h,e) Refresh(panelName));
        java_setcb(jChirMinus, 'ActionPerformedCallback', @(h,e) Refresh(panelName));
    end
    java_setcb(jSave,   'ActionPerformedCallback', @(h,e) OnSave(panelName));
    java_setcb(jCancel, 'ActionPerformedCallback', @(h,e) OnCancel(panelName));

    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);
end


%% ===== STATE ACCESS =====
function [S, ctrl] = GetState(panelName) %#ok<DEFNU>
    if (nargin < 1) || isempty(panelName); panelName = 'WaveletDesigner'; end
    ctrl = bst_get('PanelControls', panelName);
    S = [];
    if isempty(ctrl) || ~isfield(ctrl,'hFig') || ~ishandle(ctrl.hFig); return; end
    S = getappdata(ctrl.hFig, 'WaveletDesignerState');
end

function SetState(ctrl, S)
    setappdata(ctrl.hFig, 'WaveletDesignerState', S);
end


%% ===== KERNEL LIST (curated, friendly names) =====
function [keys, displays] = i_kernel_list()
    keys = {'mexhat','dog','heat','inverse_heat','tikhonov'};
    displays = cell(1, numel(keys));
    for i = 1:numel(keys)
        try
            m = bst_eigfilter_kernel('info', keys{i});
            displays{i} = m.display;
        catch
            displays{i} = keys{i};
        end
    end
end

function key = i_current_kernel(ctrl)
    idx = ctrl.jKernel.getSelectedIndex() + 1;          % java 0-based
    idx = max(1, min(numel(ctrl.KernelKeys), idx));
    key = ctrl.KernelKeys{idx};
end


%% ===== SCALE SLIDERS (mode-index space) =====
function BuildParamWidgets(ctrl, hFig)
    import javax.swing.*;
    key  = i_current_kernel(ctrl);
    meta = bst_eigfilter_kernel('info', key);
    S = getappdata(hFig, 'WaveletDesignerState');
    K = numel(S.Lambda);
    ctrl.jParams.removeAll();
    pf = fieldnames(meta.params);
    names = {};
    nP = numel(pf);
    for i = 1:nP
        nm = pf{i};
        % Stagger default mode positions so multi-scale kernels (e.g. dog: t1,t2) start
        % at DISTINCT scales rather than all at K/2 (which would make t1==t2).
        defMode = max(1, min(K, round(K * i/(nP+1))));
        % Title row "<label>: mode N" + a slider row capped by coarse ... fine.
        [js, jTitle] = i_labeled_slider(ctrl.jParams, ...
            sprintf('%s: mode %d', i_param_label(nm), defMode), 'coarse', 'fine', 1, K, defMode);
        % store handles as client properties so ReadParams / labels can retrieve them
        ctrl.jParams.putClientProperty(['slider_' nm], js);
        ctrl.jParams.putClientProperty(['title_' nm], jTitle);
        java_setcb(js, 'StateChangedCallback', @(h,e) OnParamSlider('WaveletDesigner', nm, h));
        names{end+1} = nm; %#ok<AGROW>
    end
    ctrl.jParams.revalidate(); ctrl.jParams.repaint();
    S.ParamNames = names;
    setappdata(hFig, 'WaveletDesignerState', S);
end

function params = ReadParams(S, ctrl)
    params = struct();
    K = numel(S.Lambda);
    for i = 1:numel(S.ParamNames)
        nm = S.ParamNames{i};
        js = ctrl.jParams.getClientProperty(['slider_' nm]);
        if isempty(js); continue; end
        k = max(1, min(K, double(js.getValue())));
        params.(nm) = i_param_value(nm, S.Lambda(k));
    end
    % dog (difference of Gaussians) requires t1 < t2; the two sliders are independent,
    % so order them and guarantee a strict separation regardless of slider positions.
    if isfield(params,'t1') && isfield(params,'t2')
        lo = min(params.t1, params.t2);
        hi = max(params.t1, params.t2);
        if hi <= lo * (1 + 1e-3); hi = lo * 1.5; end
        params.t1 = lo; params.t2 = hi;
    end
end

function v = i_param_value(name, lamk)
    % Map a mode-index's eigenvalue to the kernel's scale parameter.
    if strcmpi(name, 'beta')
        v = max(lamk, eps);            % tikhonov acts at lambda ~ beta
    else
        v = 1 ./ max(lamk, eps);       % heat/mexhat/dog: peak/cutoff at lambda = 1/t
    end
end

function lab = i_param_label(name)
    switch lower(name)
        case 't',    lab = 'Scale';
        case 't1',   lab = 'Band edge 1';      % dog passband ends (mode space); math
        case 't2',   lab = 'Band edge 2';      % orders them so t1 < t2 automatically
        case 'beta', lab = 'Scale';
        otherwise,   lab = name;
    end
end

function i_param_label_update(panelName, name)
    ctrl = bst_get('PanelControls', panelName);
    if isempty(ctrl); return; end
    S = getappdata(ctrl.hFig, 'WaveletDesignerState');
    js = ctrl.jParams.getClientProperty(['slider_' name]);
    jt = ctrl.jParams.getClientProperty(['title_' name]);
    if isempty(js) || isempty(jt); return; end
    k = max(1, min(numel(S.Lambda), double(js.getValue())));
    jt.setText(sprintf('%s: mode %d', i_param_label(name), k));
end


%% ===== SINGLE WAVELET DESIGN (sections 1-2; never touches tiling) =====
% Reads the Input + Filter-kernel controls into ONE designed wavelet. This is the
% canonical design/exploration object. It carries no tiling fields.
function wavelet = BuildWavelet(S, ctrl)
    wavelet = struct('Kernel', i_current_kernel(ctrl), ...
                     'Params', ReadParams(S, ctrl), ...
                     'Direction', SeedDirection(S, ctrl), ...
                     'Chirality', 0, 'Axis', [0 0 1]);
    if S.isDirac
        if ctrl.jChirPlus.isSelected()
            wavelet.Chirality = +1;
        elseif ctrl.jChirMinus.isSelected()
            wavelet.Chirality = -1;
        else
            wavelet.Chirality = 0;
        end
    end
end

%% ===== SPECTRUM TILING OPTIONS (section 3; consumed by bst_filterbank_tiles) =====
% Reads only the tiling controls. Separate from the wavelet design so the two cannot
% corrupt one another: a single wavelet is rendered directly; tiling is an explicit,
% opt-in transform of a finished wavelet.
function opts = BuildTiling(S, ctrl)
    lam = S.Lambda;
    N = double(ctrl.jTiles.getValue());
    if N < 1; N = 1; end
    % Span the MEANINGFUL spectrum: exclude the near-zero degenerate kernel modes (the
    % constant-quaternion null space, lambda ~ 1e-20) so the coarsest tile isn't a
    % degenerate DC-only filter. Floor the lower bound at 1e-4 of the max eigenvalue.
    hi  = max(lam);
    pos = lam(lam > hi * 1e-4);
    if isempty(pos); lo = max(eps, hi*1e-3); else; lo = min(pos); end
    opts = struct('N', N, 'Spacing', 'geometric', 'LambdaRange', [lo hi], 'Chiralities', []);
    if S.isDirac && ctrl.jChiSplit.isSelected()
        opts.Chiralities = [1 -1];
    end
end

% A tiling is active only when the user opts in (more than one spectral tile, or a
% chirality cross). Otherwise the design stays a single wavelet.
function tf = i_is_tiling(opts)
    tf = (opts.N > 1) || ~isempty(opts.Chiralities);
end

%% ===== PURE: embed a local-frame direction into an ambient 3-vector =====
function d = EmbedDirection(phiDeg, thetaDeg, U, V, N) %#ok<DEFNU>
    d = cosd(thetaDeg) * (cosd(phiDeg)*U(:).' + sind(phiDeg)*V(:).') + sind(thetaDeg)*N(:).';
    nrm = norm(d); if nrm > 0; d = d / nrm; end
end

function d = SeedDirection(S, ctrl)
    d = [1 0 0];
    if ~S.isDirac || isempty(ctrl.jAz) || isempty(S.iVertex); return; end
    v = S.iVertex;
    if v > size(S.FrameU,1); return; end
    phi = double(ctrl.jAz.getValue());   % in-plane angle (deg)
    th  = double(ctrl.jEl.getValue());   % tilt toward normal (deg)
    d = EmbedDirection(phi, th, S.FrameU(v,:), S.FrameV(v,:), S.FrameN(v,:));
end


%% ===== FILTER THE SEED FOR ONE TILE (uses cached coeffs) =====
function J = ComputeField(S, tile)
    args = {'TransferFn', bst_eigfilter_kernel(tile.Kernel, tile.Params), 'Coeffs', S.SeedCoeffs};
    if S.isDirac && (tile.Chirality ~= 0)
        args = [args, {'Chirality', struct('Axis', tile.Axis, 'Sign', tile.Chirality)}];
    end
    Jc = bst_dirac_eigenmodes_filter(S.EigenMat, S.Op.Mass, [], 'custom', args{:});
    J = real(Jc);
end


%% ===== RECOMPUTE + PUSH PREVIEW =====
function Refresh(panelName)
    [S, ctrl] = GetState(panelName);
    if isempty(S) || isempty(S.SeedCoeffs); return; end
    % SINGLE-WAVELET path by default; the tiling module is invoked ONLY on opt-in.
    wavelet = BuildWavelet(S, ctrl);
    opts    = BuildTiling(S, ctrl);
    if i_is_tiling(opts)
        Tiles = bst_filterbank_tiles(wavelet, opts);   % spectrum tiling module
    else
        Tiles = wavelet;                                % one designed wavelet, as-is
    end
    S.ActiveTile = min(max(1, S.ActiveTile), numel(Tiles));
    S.Tiles = Tiles;
    SetState(ctrl, S);
    J = ComputeField(S, Tiles(S.ActiveTile));
    S.ctxFn.PushField(J);
    % active-tile indicator
    if numel(Tiles) > 1
        ctrl.jActiveTile.setText(sprintf('<HTML>Showing tile %d / %d', S.ActiveTile, numel(Tiles)));
    else
        ctrl.jActiveTile.setText('');
    end
    % spectrum strip: all tiles, active highlighted, clickable
    kernels = arrayfun(@(t) bst_eigfilter_kernel(t.Kernel, t.Params), Tiles, 'UniformOutput', false);
    bank = struct('Kernels', {kernels}, 'Active', S.ActiveTile, ...
                  'OnSelect', @(j) OnSelectTile(panelName, j));
    view_eigfilter_response(bank, S.Lambda, sprintf('%s - %d tile(s)', wavelet.Kernel, numel(Tiles)));
end


%% ===== SEEDING =====
function SetSeedVertex(panelName, iVertex) %#ok<DEFNU>
    [S, ctrl] = GetState(panelName);
    if isempty(S); return; end
    S.iVertex = iVertex;                 % SeedDirection reads the vertex's local frame
    nV = double(max(cellfun(@(x) max(x(:)), S.EigenMat.GlobalVertices)));
    d = SeedDirection(S, ctrl);
    if isfield(S.ctxFn,'DrawSeed') && ~isempty(S.ctxFn.DrawSeed)
        S.ctxFn.DrawSeed(iVertex, d);    % cyan seed-vector marker that tracks the sliders
    end
    Jdelta = zeros(3*nV, 1);
    Jdelta(3*(iVertex-1) + (1:3)) = d(:);
    [~,~,c] = bst_dirac_eigenmodes_filter(S.EigenMat, S.Op.Mass, Jdelta, 'custom', ...
                  'TransferFn', @(l) ones(size(l)), 'ReturnCoeffs', true);
    S.SeedCoeffs = c;
    SetState(ctrl, S);
    Refresh(panelName);
end

function ReseedAndRefresh(panelName)
    [S, ctrl] = GetState(panelName); %#ok<ASGLU>
    if ~isempty(S) && ~isempty(S.iVertex)
        SetSeedVertex(panelName, S.iVertex);    % re-project the delta with the new direction
    else
        Refresh(panelName);
    end
end


%% ===== TILE SELECTION =====
function OnSelectTile(panelName, j) %#ok<DEFNU>
    [S, ctrl] = GetState(panelName);
    if isempty(S); return; end
    S.ActiveTile = j;
    SetState(ctrl, S);
    Refresh(panelName);
end


%% ===== GATED SLIDER HANDLERS (label live; recompute on settle) =====
function OnParamSlider(panelName, nm, js)
    i_param_label_update('WaveletDesigner', nm);
    if ~js.getValueIsAdjusting(); Refresh(panelName); end
end

function OnTilesSlider(panelName, js)
    [S, ctrl] = GetState(panelName); %#ok<ASGLU>
    if ~isempty(ctrl)
        n = double(ctrl.jTiles.getValue());
        if n <= 1; ctrl.jTilesVal.setText('Tiles: 1 (single wavelet)');
        else;      ctrl.jTilesVal.setText(sprintf('Tiles: %d', n)); end
    end
    if ~js.getValueIsAdjusting(); Refresh(panelName); end
end

function OnDirSlider(panelName, js)
    [S, ctrl] = GetState(panelName); %#ok<ASGLU>
    if ~isempty(ctrl) && ~isempty(ctrl.jAz)
        ctrl.jAzVal.setText(sprintf('<HTML>In-plane angle: %d&deg;', double(ctrl.jAz.getValue())));
        ctrl.jElVal.setText(sprintf('<HTML>Tilt to normal: %d&deg;', double(ctrl.jEl.getValue())));
    end
    if ~js.getValueIsAdjusting(); ReseedAndRefresh(panelName); end
end


%% ===== ON KERNEL CHANGE =====
function OnKernelChanged(panelName)
    ctrl = bst_get('PanelControls', panelName);
    if isempty(ctrl); return; end
    BuildParamWidgets(ctrl, ctrl.hFig);
    Refresh(panelName);
end


%% ===== LOAD A SAVED BANK INTO THE WIDGETS =====
function LoadBank(panelName, loadBank) %#ok<DEFNU>
    [S, ctrl] = GetState(panelName);
    if isempty(S) || ~isfield(loadBank,'Tiling') || isempty(loadBank.Tiling); return; end
    % saved Tiling = struct('Wavelet', <single design>, 'Opts', <tiling options>)
    wavelet = loadBank.Tiling.Wavelet;
    opts    = loadBank.Tiling.Opts;
    % --- restore the single-wavelet design (sections 1-2) ---
    for k = 1:numel(ctrl.KernelKeys)
        if strcmpi(ctrl.KernelKeys{k}, wavelet.Kernel); ctrl.jKernel.setSelectedIndex(k-1); break; end
    end
    BuildParamWidgets(ctrl, ctrl.hFig);
    [S, ctrl] = GetState(panelName);
    % scale sliders: closest mode index to each saved param value
    if isfield(wavelet,'Params') && isstruct(wavelet.Params)
        for i = 1:numel(S.ParamNames)
            nm = S.ParamNames{i};
            if ~isfield(wavelet.Params, nm); continue; end
            js = ctrl.jParams.getClientProperty(['slider_' nm]);
            if isempty(js); continue; end
            target = wavelet.Params.(nm);
            if strcmpi(nm,'beta'); lamWanted = target; else; lamWanted = 1/max(target,eps); end
            [~, k] = min(abs(S.Lambda - lamWanted));
            js.setValue(k);
        end
    end
    if S.isDirac && isfield(wavelet,'Chirality')
        switch wavelet.Chirality
            case 1,  ctrl.jChirPlus.setSelected(true);
            case -1, ctrl.jChirMinus.setSelected(true);
            otherwise, ctrl.jChirNone.setSelected(true);
        end
    end
    % --- restore the tiling options (section 3) ---
    if isstruct(opts)
        if isfield(opts,'N'); ctrl.jTiles.setValue(max(1, opts.N)); end
        if S.isDirac && isfield(opts,'Chiralities') && ~isempty(opts.Chiralities)
            ctrl.jChiSplit.setSelected(true);
        end
    end
    iVertex = [];
    if isfield(loadBank,'Provenance') && isstruct(loadBank.Provenance) && isfield(loadBank.Provenance,'DesignVertex')
        iVertex = loadBank.Provenance.DesignVertex;
    end
    if ~isempty(iVertex); SetSeedVertex(panelName, iVertex); end
end


%% ===== SAVE / CANCEL =====
function OnSave(panelName) %#ok<DEFNU>
    [S, ctrl] = GetState(panelName);
    if isempty(S); return; end
    if isempty(S.SeedCoeffs)
        java_dialog('warning', 'Seed the filter first (click a cortical vertex or pick a source map).', 'Save filterbank');
        return;
    end
    wavelet = BuildWavelet(S, ctrl);
    opts    = BuildTiling(S, ctrl);
    if i_is_tiling(opts)
        Tiles = bst_filterbank_tiles(wavelet, opts);   % materialize the bank via the module
    else
        Tiles = wavelet;                                % single wavelet
    end
    fb = db_template('waveletmat');
    fb.ParentEigen = S.EigenFile;
    fb.Variant     = S.Variant;
    fb.Tiles       = Tiles;                              % materialized bank (1 or N recipes)
    fb.Tiling      = struct('Wavelet', wavelet, 'Opts', opts);   % regenerable: design + tiling
    fb.Provenance  = struct('DesignVertex', S.iVertex, 'ComputeDate', datestr(now,'yyyy-mm-dd HH:MM:SS'));
    Comment = sprintf('%s filterbank (%d tile(s))', wavelet.Kernel, numel(Tiles));
    db_add_wavelet(S.iSubject, S.EigenFile, fb, Comment);
    S.ctxFn.Close();
end

function OnCancel(panelName) %#ok<DEFNU>
    S = GetState(panelName);
    if ~isempty(S) && isfield(S,'ctxFn') && ~isempty(S.ctxFn); S.ctxFn.Close(); end
end


%% ===== JSLIDER HELPER =====
function js = i_slider(jParent, constraints, mn, mx, val)
    import javax.swing.*;
    js = JSlider(mn, mx, val);
    % Small preferred width: 'hfill' stretches the slider to the available row width, so a
    % large preferred size is what pushed the right edge past the narrow tools-tab panel.
    js.setPreferredSize(java_scaled('dimension', 40, 22));
    jParent.add(constraints, js);
end

% Title row ("Title: value") + a slider row capped by two end labels (loLabel ... hiLabel)
% so 'hfill' is bounded by the panel edge. Returns the slider and the (value-bearing) title.
function [js, jTitle] = i_labeled_slider(jParent, titleText, loLabel, hiLabel, mn, mx, val)
    jTitle = gui_component('label', jParent, 'br', titleText);
    gui_component('label', jParent, 'br', loLabel);     % left cap
    js = i_slider(jParent, 'hfill', mn, mx, val);        % fills between the caps
    gui_component('label', jParent, '', hiLabel);        % right cap (bounds the fill)
end
