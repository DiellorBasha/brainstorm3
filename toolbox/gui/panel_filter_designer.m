function varargout = panel_filter_designer(varargin)
% PANEL_FILTER_DESIGNER: Transient control panel for the spatial filterbank designer.
%
% Owns the design state and the live-preview controller for a filterbank design
% session. Shown by view_filter_designer via gui_show(...,'BrainstormTab',...) and
% torn down by gui_hide('FilterDesigner'). Java control handles are stored in the
% panel's sControls (bst_get('PanelControls','FilterDesigner')); the mutable MATLAB
% session state lives in the preview figure's appdata ('FilterDesignerState'), which
% the orchestrator (view_filter_designer) creates and links via ctxFn callbacks.
%
% Workflow: the user designs ONE wavelet (kernel + scale + direction + chirality) and
% sees it live on the cortex; optionally they then tile the spectrum into a small bank
% of that wavelet at several scales. The Dirac structure exposes three filter axes:
% SCALE g(lambda) (shared eigfilter library, driven here in mode-index space), DIRECTION
% (the ambient 3D seed vector, set by azimuth/elevation sliders), and CHIRALITY (helicity
% projector). Source vectors are full 3D ambient fields (never normal/tangent components).
%
% Dispatched subfunctions (panel_filter_designer('Name', args...)):
%   CreatePanel(EigenMat, EigenFile, hFig, ctxFn) -> bstPanel
%   SetSeedVertex(panelName, iVertex)   seed a delta at a clicked vertex
%   SetSeedSource(panelName, J)         seed the active loaded source frame [3nV x 1]
%   GetState(panelName) -> S            the current session-state struct
%   OnSelectTile(panelName, j)          make tile j the displayed tile
%   OnSave(panelName) / OnCancel(panelName)
%
% SEE ALSO: view_filter_designer, bst_filterbank_tiles, bst_dirac_eigenmodes_filter,
%           bst_eigfilter_kernel, view_eigfilter_response
%
% Authors: Diellor Basha, 2026

    eval(macro_method);
end


%% ===== CREATE PANEL =====
function bstPanelNew = CreatePanel(EigenMat, EigenFile, hFig, ctxFn) %#ok<DEFNU>
    import java.awt.*;
    import javax.swing.*;
    panelName = 'FilterDesigner';
    isDirac = strcmpi(EigenMat.Variant, 'Dirac');

    jPanel = gui_river([4 4], [3 8 8 8]);

    % ===== WAVELET DESIGN =====
    gui_component('label', jPanel, '', '<HTML><B>Design a single wavelet</B>');
    gui_component('label', jPanel, 'br', ['Operator: ' EigenMat.Variant]);

    % --- input mode ---
    jPanel.add('br', JLabel('Input:'));
    jInputDelta  = gui_component('radio', jPanel, 'tab',  'Delta (click a vertex)');
    jInputSource = gui_component('radio', jPanel, 'br tab', 'Active source map');
    jInputDelta.setSelected(true);
    grpIn = ButtonGroup(); grpIn.add(jInputDelta); grpIn.add(jInputSource);

    % --- kernel dropdown (curated, friendly display names) ---
    [keys, displays] = i_kernel_list();
    jPanel.add('br', JLabel('Kernel:'));
    jKernel = gui_component('combobox', jPanel, 'tab', [], {displays}, [], [], []);

    % --- scale parameter sliders (rebuilt per kernel, in mode-index space) ---
    jParams = gui_river([2 2], [0 4 0 4]);
    jPanel.add('br hfill', jParams);

    % --- seed direction (Dirac, delta mode): azimuth + elevation -> ambient 3D vector ---
    jAz = []; jEl = []; jAzVal = []; jElVal = [];
    if isDirac
        gui_component('label', jPanel, 'br', '<HTML><I>Seed direction (ambient 3D)</I>');
        jPanel.add('br', JLabel('Azimuth:'));
        jAz    = i_slider(jPanel, 'tab hfill', 0, 360, 0);
        jAzVal = gui_component('label', jPanel, '', '0&deg;');
        jPanel.add('br', JLabel('Elevation:'));
        jEl    = i_slider(jPanel, 'tab hfill', -90, 90, 0);
        jElVal = gui_component('label', jPanel, '', '0&deg;');
    end

    % --- chirality (Dirac) ---
    jChirNone = []; jChirPlus = []; jChirMinus = [];
    if isDirac
        jPanel.add('br', JLabel('Chirality:'));
        jChirNone  = gui_component('radio', jPanel, 'tab', 'None');
        jChirPlus  = gui_component('radio', jPanel, '', '+ right');
        jChirMinus = gui_component('radio', jPanel, '', '- left');
        jChirNone.setSelected(true);
        grpCh = ButtonGroup(); grpCh.add(jChirNone); grpCh.add(jChirPlus); grpCh.add(jChirMinus);
    end

    % ===== SPECTRUM TILING (optional) + ACTIONS =====
    gui_component('label', jPanel, 'br', ' ');
    gui_component('label', jPanel, 'br', '<HTML><B>Spectrum tiling (optional)</B>');
    jPanel.add('br', JLabel('Tiles:'));
    jTiles    = i_slider(jPanel, 'tab hfill', 1, 12, 1);
    jTilesVal = gui_component('label', jPanel, '', '1 (single wavelet)');
    jChiSplit = gui_component('checkbox', jPanel, 'br', 'Cross with both chiralities');
    if ~isDirac; jChiSplit.setEnabled(false); end
    jActiveTile = gui_component('label', jPanel, 'br', '');

    jSave   = gui_component('button', jPanel, 'br right', 'Save bank');
    jCancel = gui_component('button', jPanel, '', 'Cancel');

    % --- collect Java handles ---
    ctrl = struct('jKernel',jKernel, 'KernelKeys',{keys}, 'jParams',jParams, ...
                  'jAz',jAz, 'jEl',jEl, 'jAzVal',jAzVal, 'jElVal',jElVal, ...
                  'jChirNone',jChirNone, 'jChirPlus',jChirPlus, 'jChirMinus',jChirMinus, ...
                  'jTiles',jTiles, 'jTilesVal',jTilesVal, 'jChiSplit',jChiSplit, ...
                  'jActiveTile',jActiveTile, 'jInputDelta',jInputDelta, 'jInputSource',jInputSource, ...
                  'jSave',jSave, 'jCancel',jCancel, 'hFig',hFig);

    % --- session state on the preview figure ---
    [~, iSubject] = bst_get('EigenFile', EigenFile);
    if isempty(iSubject); iSubject = []; end
    S = struct('EigenMat',EigenMat, 'EigenFile',file_short(EigenFile), 'iSubject',iSubject, ...
               'Variant',EigenMat.Variant, 'isDirac',isDirac, 'ctxFn',ctxFn, ...
               'Op',load(file_fullpath(EigenMat.OperatorFile)), ...
               'Lambda',double(EigenMat.Lambda{1}(:)), ...
               'SeedCoeffs',[], 'ActiveTile',1, 'iVertex',[], 'Tiles',[], 'ParamNames',{{}});
    setappdata(hFig, 'FilterDesignerState', S);

    % --- build the initial scale sliders for the default kernel ---
    BuildParamWidgets(ctrl, hFig);

    % --- wire callbacks ---
    java_setcb(jKernel,   'ActionPerformedCallback', @(h,e) OnKernelChanged(panelName));
    java_setcb(jChiSplit, 'ActionPerformedCallback', @(h,e) Refresh(panelName));
    java_setcb(jTiles,    'StateChangedCallback',    @(h,e) i_tiles_label(panelName));
    java_setcb(jTiles,    'MouseReleasedCallback',   @(h,e) Refresh(panelName));
    if isDirac
        java_setcb(jAz, 'StateChangedCallback',  @(h,e) i_dir_label(panelName));
        java_setcb(jEl, 'StateChangedCallback',  @(h,e) i_dir_label(panelName));
        java_setcb(jAz, 'MouseReleasedCallback', @(h,e) ReseedAndRefresh(panelName));
        java_setcb(jEl, 'MouseReleasedCallback', @(h,e) ReseedAndRefresh(panelName));
        java_setcb(jChirNone,  'ActionPerformedCallback', @(h,e) Refresh(panelName));
        java_setcb(jChirPlus,  'ActionPerformedCallback', @(h,e) Refresh(panelName));
        java_setcb(jChirMinus, 'ActionPerformedCallback', @(h,e) Refresh(panelName));
    end
    java_setcb(jSave,   'ActionPerformedCallback', @(h,e) OnSave(panelName));
    java_setcb(jCancel, 'ActionPerformedCallback', @(h,e) OnCancel(panelName));

    bstPanelNew = BstPanel(panelName, jPanel, ctrl);
end


%% ===== STATE ACCESS =====
function [S, ctrl] = GetState(panelName) %#ok<DEFNU>
    if (nargin < 1) || isempty(panelName); panelName = 'FilterDesigner'; end
    ctrl = bst_get('PanelControls', panelName);
    S = [];
    if isempty(ctrl) || ~isfield(ctrl,'hFig') || ~ishandle(ctrl.hFig); return; end
    S = getappdata(ctrl.hFig, 'FilterDesignerState');
end

function SetState(ctrl, S)
    setappdata(ctrl.hFig, 'FilterDesignerState', S);
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
    S = getappdata(hFig, 'FilterDesignerState');
    K = numel(S.Lambda);
    ctrl.jParams.removeAll();
    pf = fieldnames(meta.params);
    names = {};
    for i = 1:numel(pf)
        nm = pf{i};
        ctrl.jParams.add('br', JLabel([i_param_label(nm) ':']));
        js = i_slider(ctrl.jParams, 'tab hfill', 1, K, round(K/2));
        jv = gui_component('label', ctrl.jParams, '', '');
        % store handles as client properties so ReadParams can retrieve them
        ctrl.jParams.putClientProperty(['slider_' nm], js);
        ctrl.jParams.putClientProperty(['vlabel_' nm], jv);
        java_setcb(js, 'StateChangedCallback',  @(h,e) i_param_label_update(ctrl.hFig, nm));
        java_setcb(js, 'MouseReleasedCallback', @(h,e) Refresh('FilterDesigner'));
        names{end+1} = nm; %#ok<AGROW>
    end
    ctrl.jParams.revalidate(); ctrl.jParams.repaint();
    S.ParamNames = names;
    setappdata(hFig, 'FilterDesignerState', S);
    for i = 1:numel(names); i_param_label_update(hFig, names{i}); end
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
        case 't',    lab = 'Scale (coarse-fine)';
        case 't1',   lab = 'Coarse edge';
        case 't2',   lab = 'Fine edge';
        case 'beta', lab = 'Scale (coarse-fine)';
        otherwise,   lab = name;
    end
end

function i_param_label_update(hFig, name)
    ctrl = bst_get('PanelControls', 'FilterDesigner');
    if isempty(ctrl); return; end
    S = getappdata(hFig, 'FilterDesignerState');
    js = ctrl.jParams.getClientProperty(['slider_' name]);
    jv = ctrl.jParams.getClientProperty(['vlabel_' name]);
    if isempty(js) || isempty(jv); return; end
    k = max(1, min(numel(S.Lambda), double(js.getValue())));
    jv.setText(sprintf('mode %d', k));
end


%% ===== READ WIDGETS -> base design =====
function base = BuildDesign(S, ctrl)
    key = i_current_kernel(ctrl);
    params = ReadParams(S, ctrl);
    lam = S.Lambda;
    N = double(ctrl.jTiles.getValue());
    if N < 1; N = 1; end
    base = struct('Kernel',key, 'Params',params, ...
        'Direction', SeedDirection(S, ctrl), 'Chirality',0, 'Axis',[0 0 1], ...
        'N',N, 'Spacing','geometric', ...
        'LambdaRange',[max(eps,min(lam(lam>0))) max(lam)], 'Chiralities',[]);
    if S.isDirac
        if ctrl.jChirPlus.isSelected()
            base.Chirality = +1;
        elseif ctrl.jChirMinus.isSelected()
            base.Chirality = -1;
        else
            base.Chirality = 0;
        end
        if ctrl.jChiSplit.isSelected(); base.Chiralities = [1 -1]; end
    end
end

function d = SeedDirection(S, ctrl)
    d = [1 0 0];
    if ~S.isDirac || isempty(ctrl.jAz); return; end
    az = double(ctrl.jAz.getValue());   % degrees
    el = double(ctrl.jEl.getValue());
    d = [cosd(el)*cosd(az), cosd(el)*sind(az), sind(el)];
    nrm = norm(d); if nrm > 0; d = d / nrm; end
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
    base  = BuildDesign(S, ctrl);
    Tiles = bst_filterbank_tiles(base);
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
    view_eigfilter_response(bank, S.Lambda, sprintf('%s - %d tile(s)', i_current_kernel(ctrl), numel(Tiles)));
end


%% ===== SEEDING =====
function SetSeedVertex(panelName, iVertex) %#ok<DEFNU>
    [S, ctrl] = GetState(panelName);
    if isempty(S); return; end
    nV = double(max(cellfun(@(x) max(x(:)), S.EigenMat.GlobalVertices)));
    d = SeedDirection(S, ctrl);
    Jdelta = zeros(3*nV, 1);
    Jdelta(3*(iVertex-1) + (1:3)) = d(:);
    [~,~,c] = bst_dirac_eigenmodes_filter(S.EigenMat, S.Op.Mass, Jdelta, 'custom', ...
                  'TransferFn', @(l) ones(size(l)), 'ReturnCoeffs', true);
    S.SeedCoeffs = c;  S.iVertex = iVertex;
    SetState(ctrl, S);
    Refresh(panelName);
end

function SetSeedSource(panelName, J) %#ok<DEFNU>
    [S, ctrl] = GetState(panelName);
    if isempty(S); return; end
    [~,~,c] = bst_dirac_eigenmodes_filter(S.EigenMat, S.Op.Mass, J(:), 'custom', ...
                  'TransferFn', @(l) ones(size(l)), 'ReturnCoeffs', true);
    S.SeedCoeffs = c;  S.iVertex = [];
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


%% ===== SMALL LIVE LABELS =====
function i_tiles_label(panelName)
    [S, ctrl] = GetState(panelName); %#ok<ASGLU>
    if isempty(ctrl); return; end
    n = double(ctrl.jTiles.getValue());
    if n <= 1; ctrl.jTilesVal.setText('1 (single wavelet)');
    else;      ctrl.jTilesVal.setText(sprintf('%d tiles', n)); end
end

function i_dir_label(panelName)
    [S, ctrl] = GetState(panelName); %#ok<ASGLU>
    if isempty(ctrl) || isempty(ctrl.jAz); return; end
    ctrl.jAzVal.setText(sprintf('%d&deg;', double(ctrl.jAz.getValue())));
    ctrl.jElVal.setText(sprintf('%d&deg;', double(ctrl.jEl.getValue())));
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
    base = loadBank.Tiling;
    % kernel
    for k = 1:numel(ctrl.KernelKeys)
        if strcmpi(ctrl.KernelKeys{k}, base.Kernel); ctrl.jKernel.setSelectedIndex(k-1); break; end
    end
    BuildParamWidgets(ctrl, ctrl.hFig);
    [S, ctrl] = GetState(panelName);
    % scale sliders: closest mode index to each saved param value
    if isfield(base,'Params') && isstruct(base.Params)
        for i = 1:numel(S.ParamNames)
            nm = S.ParamNames{i};
            if ~isfield(base.Params, nm); continue; end
            js = ctrl.jParams.getClientProperty(['slider_' nm]);
            if isempty(js); continue; end
            target = base.Params.(nm);
            if strcmpi(nm,'beta'); lamWanted = target; else; lamWanted = 1/max(target,eps); end
            [~, k] = min(abs(S.Lambda - lamWanted));
            js.setValue(k);
        end
    end
    if isfield(base,'N'); ctrl.jTiles.setValue(max(1, base.N)); end
    if S.isDirac && isfield(base,'Chirality')
        switch base.Chirality
            case 1,  ctrl.jChirPlus.setSelected(true);
            case -1, ctrl.jChirMinus.setSelected(true);
            otherwise, ctrl.jChirNone.setSelected(true);
        end
        if isfield(base,'Chiralities') && ~isempty(base.Chiralities); ctrl.jChiSplit.setSelected(true); end
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
    base  = BuildDesign(S, ctrl);
    Tiles = bst_filterbank_tiles(base);
    fb = db_template('filterbankmat');
    fb.ParentEigen = S.EigenFile;
    fb.Variant     = S.Variant;
    fb.Tiles       = Tiles;
    fb.Tiling      = base;
    fb.Provenance  = struct('DesignVertex', S.iVertex, 'ComputeDate', datestr(now,'yyyy-mm-dd HH:MM:SS'));
    Comment = sprintf('%s filterbank (%d tile(s))', base.Kernel, numel(Tiles));
    db_add_filterbank(S.iSubject, S.EigenFile, fb, Comment);
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
    js.setPreferredSize(java_scaled('dimension', 130, 22));
    jParent.add(constraints, js);
end
