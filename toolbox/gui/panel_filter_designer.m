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
% The Dirac structure exposes three filter axes (see bst_dirac, bst_dirac_eigenmodes_filter):
% SCALE g(lambda) (shared eigfilter library), DIRECTION (seed quaternion), CHIRALITY
% (helicity projector). The panel reads one base design, tiles the spectrum into a bank
% (bst_filterbank_tiles), previews the selected tile live, and saves a recipe bank.
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

    jPanel = gui_river([5 5], [3 8 8 8]);

    % --- operator (read-only) ---
    gui_component('label', jPanel, '', ['<HTML><B>Operator:</B> ' EigenMat.Variant]);

    % --- input mode ---
    jPanel.add('br', JLabel('Input:'));
    jInputDelta  = gui_component('radio', jPanel, 'tab',  'Delta (click vertex)');
    jInputSource = gui_component('radio', jPanel, 'br tab', 'Active source map');
    jInputDelta.setSelected(true);
    grpIn = ButtonGroup(); grpIn.add(jInputDelta); grpIn.add(jInputSource);

    % --- kernel dropdown (from the registry) ---
    names = bst_eigfilter_kernel('list');
    jPanel.add('br', JLabel('Kernel:'));
    jKernel = gui_component('combobox', jPanel, 'tab', [], {names}, [], [], []);

    % --- params sub-panel (rebuilt when the kernel changes) ---
    jParams = gui_river([2 2], [0 4 0 4]);
    jPanel.add('br hfill', jParams);

    % --- direction + chirality (Dirac only) ---
    jDir = []; jChir = [];
    if isDirac
        jPanel.add('br', JLabel('Direction:'));
        jDir = gui_component('combobox', jPanel, 'tab', [], {{'Surface normal','Tangent','Ambient X'}}, [], [], []);
        jPanel.add('br', JLabel('Chirality:'));
        jChir = gui_component('combobox', jPanel, 'tab', [], {{'None','+ (right)','- (left)'}}, [], [], []);
    end

    % --- tiling ---
    jPanel.add('br', JLabel('Tiles:'));
    jTiles = gui_component('text', jPanel, 'tab', '4');
    jChiSplit = gui_component('checkbox', jPanel, 'br', 'Cross with both chiralities');
    if ~isDirac; jChiSplit.setEnabled(false); end

    % --- save / cancel ---
    jSave   = gui_component('button', jPanel, 'br right', 'Save bank');
    jCancel = gui_component('button', jPanel, '', 'Cancel');

    % --- collect Java handles into sControls (hFig links to the appdata state) ---
    ctrl = struct('jKernel',jKernel, 'jParams',jParams, 'jDir',jDir, 'jChir',jChir, ...
                  'jTiles',jTiles, 'jChiSplit',jChiSplit, 'jInputDelta',jInputDelta, ...
                  'jInputSource',jInputSource, 'jSave',jSave, 'jCancel',jCancel, 'hFig',hFig);

    % --- initialise the session state on the preview figure ---
    [~, iSubject] = bst_get('EigenFile', EigenFile);
    if isempty(iSubject); iSubject = []; end
    VertNormals = [];
    try
        sSurf = in_tess_bst(EigenMat.ParentSurface);
        if isfield(sSurf,'VertNormals'); VertNormals = sSurf.VertNormals; end
    catch
    end
    S = struct('EigenMat',EigenMat, 'EigenFile',file_short(EigenFile), 'iSubject',iSubject, ...
               'Variant',EigenMat.Variant, 'isDirac',isDirac, 'ctxFn',ctxFn, ...
               'Op',load(file_fullpath(EigenMat.OperatorFile)), 'VertNormals',VertNormals, ...
               'SeedCoeffs',[], 'ActiveTile',1, 'iVertex',[], 'Tiles',[], 'ParamFields',struct());
    setappdata(hFig, 'FilterDesignerState', S);

    % --- build the initial parameter widgets for the default kernel ---
    BuildParamWidgets(ctrl, hFig);

    % --- wire callbacks ---
    java_setcb(jKernel,   'ActionPerformedCallback', @(h,e) OnKernelChanged(panelName));
    java_setcb(jTiles,    'ActionPerformedCallback', @(h,e) Refresh(panelName));
    java_setcb(jChiSplit, 'ActionPerformedCallback', @(h,e) Refresh(panelName));
    if isDirac
        java_setcb(jDir,  'ActionPerformedCallback', @(h,e) ReseedAndRefresh(panelName));
        java_setcb(jChir, 'ActionPerformedCallback', @(h,e) Refresh(panelName));
    end
    java_setcb(jSave,     'ActionPerformedCallback', @(h,e) OnSave(panelName));
    java_setcb(jCancel,   'ActionPerformedCallback', @(h,e) OnCancel(panelName));

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


%% ===== PARAM WIDGETS (auto-built from kernel meta) =====
function BuildParamWidgets(ctrl, hFig)
    kname = char(ctrl.jKernel.getSelectedItem());
    meta  = bst_eigfilter_kernel('info', kname);
    ctrl.jParams.removeAll();
    pf = struct();
    fn = fieldnames(meta.params);
    for i = 1:numel(fn)
        d = meta.params.(fn{i});
        ctrl.jParams.add('br', javax.swing.JLabel([fn{i} ':']));
        jv = gui_component('text', ctrl.jParams, 'tab', num2str(d.default));
        java_setcb(jv, 'ActionPerformedCallback', @(h,e) Refresh('FilterDesigner'));
        pf.(fn{i}) = jv;
    end
    ctrl.jParams.revalidate(); ctrl.jParams.repaint();
    S = getappdata(hFig, 'FilterDesignerState');
    S.ParamFields = pf;
    setappdata(hFig, 'FilterDesignerState', S);
end

function params = ReadParams(S)
    params = struct();
    fn = fieldnames(S.ParamFields);
    for i = 1:numel(fn)
        v = str2double(char(S.ParamFields.(fn{i}).getText()));
        if ~isnan(v); params.(fn{i}) = v; end
    end
end


%% ===== READ WIDGETS -> base design =====
function base = BuildDesign(S, ctrl)
    kname = char(ctrl.jKernel.getSelectedItem());
    params = ReadParams(S);
    lam = double(S.EigenMat.Lambda{1}(:));
    N = round(str2double(char(ctrl.jTiles.getText())));
    if isnan(N) || N < 1; N = 1; end
    base = struct('Kernel',kname, 'Params',params, ...
        'Direction',[1 0 0], 'Chirality',0, 'Axis',[0 0 1], ...
        'N',N, 'Spacing','geometric', ...
        'LambdaRange',[max(eps,min(lam)) max(lam)], 'Chiralities',[]);
    if S.isDirac
        switch char(ctrl.jChir.getSelectedItem())
            case '+ (right)', base.Chirality = +1;
            case '- (left)',  base.Chirality = -1;
            otherwise,        base.Chirality = 0;
        end
        if ctrl.jChiSplit.isSelected(); base.Chiralities = [1 -1]; end
    end
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
    if isempty(S) || isempty(S.SeedCoeffs); return; end       % no seed yet
    base  = BuildDesign(S, ctrl);
    Tiles = bst_filterbank_tiles(base);
    S.ActiveTile = min(max(1, S.ActiveTile), numel(Tiles));
    S.Tiles = Tiles;
    SetState(ctrl, S);
    % push the active tile's field to the preview figure
    J = ComputeField(S, Tiles(S.ActiveTile));
    S.ctxFn.PushField(J);
    % spectrum strip with all tiles, active highlighted, clickable
    kernels = arrayfun(@(t) bst_eigfilter_kernel(t.Kernel, t.Params), Tiles, 'UniformOutput', false);
    bank = struct('Kernels', {kernels}, 'Active', S.ActiveTile, ...
                  'OnSelect', @(j) OnSelectTile(panelName, j));
    view_eigfilter_response(bank, S.EigenMat.Lambda{1}, sprintf('%s tiles (n=%d)', base.Kernel, numel(Tiles)));
end


%% ===== SEEDING =====
function SetSeedVertex(panelName, iVertex) %#ok<DEFNU>
    [S, ctrl] = GetState(panelName);
    if isempty(S); return; end
    nV = double(max(cellfun(@(x) max(x(:)), S.EigenMat.GlobalVertices)));
    dirVec = ResolveDirection(S, ctrl, iVertex);
    Jdelta = zeros(3*nV, 1);
    Jdelta(3*(iVertex-1) + (1:3)) = dirVec(:);
    [~,~,c] = bst_dirac_eigenmodes_filter(S.EigenMat, S.Op.Mass, Jdelta, 'custom', ...
                  'TransferFn', @(l) ones(size(l)), 'ReturnCoeffs', true);
    S.SeedCoeffs = c;  S.iVertex = iVertex;
    SetState(ctrl, S);
    Refresh(panelName);
end

function SetSeedSource(panelName, J) %#ok<DEFNU>
    % Seed from a loaded unconstrained source frame J [3nV x 1] (the source-map input).
    [S, ctrl] = GetState(panelName);
    if isempty(S); return; end
    [~,~,c] = bst_dirac_eigenmodes_filter(S.EigenMat, S.Op.Mass, J(:), 'custom', ...
                  'TransferFn', @(l) ones(size(l)), 'ReturnCoeffs', true);
    S.SeedCoeffs = c;  S.iVertex = [];
    SetState(ctrl, S);
    Refresh(panelName);
end

function ReseedAndRefresh(panelName)
    % Direction changed: if a delta seed exists, re-project it with the new direction.
    [S, ctrl] = GetState(panelName); %#ok<ASGLU>
    if ~isempty(S) && ~isempty(S.iVertex)
        SetSeedVertex(panelName, S.iVertex);
    else
        Refresh(panelName);
    end
end

function dirVec = ResolveDirection(S, ctrl, iVertex)
    dirVec = [1 0 0];
    if ~S.isDirac || isempty(ctrl.jDir); return; end
    switch char(ctrl.jDir.getSelectedItem())
        case 'Surface normal'
            if ~isempty(S.VertNormals) && (iVertex <= size(S.VertNormals,1))
                dirVec = S.VertNormals(iVertex, :);
            end
        case 'Tangent'
            if ~isempty(S.VertNormals) && (iVertex <= size(S.VertNormals,1))
                n = S.VertNormals(iVertex, :);
                ref = [1 0 0]; if abs(n*ref') > 0.9; ref = [0 1 0]; end
                dirVec = cross(n, ref);
            end
        otherwise   % 'Ambient X'
            dirVec = [1 0 0];
    end
    nrm = norm(dirVec); if nrm > 0; dirVec = dirVec / nrm; end
end


%% ===== TILE SELECTION =====
function OnSelectTile(panelName, j) %#ok<DEFNU>
    [S, ctrl] = GetState(panelName);
    if isempty(S); return; end
    S.ActiveTile = j;
    SetState(ctrl, S);
    Refresh(panelName);
end


%% ===== LOAD A SAVED BANK INTO THE WIDGETS =====
function LoadBank(panelName, loadBank) %#ok<DEFNU>
    [S, ctrl] = GetState(panelName);
    if isempty(S) || ~isfield(loadBank,'Tiling') || isempty(loadBank.Tiling); return; end
    base = loadBank.Tiling;
    % kernel
    items = ctrl.jKernel.getModel();
    for k = 0:items.getSize()-1
        if strcmpi(char(items.getElementAt(k)), base.Kernel); ctrl.jKernel.setSelectedIndex(k); break; end
    end
    BuildParamWidgets(ctrl, ctrl.hFig);
    [S, ctrl] = GetState(panelName);
    % params
    if isfield(base,'Params') && isstruct(base.Params)
        pf = fieldnames(S.ParamFields);
        for i = 1:numel(pf)
            if isfield(base.Params, pf{i}); S.ParamFields.(pf{i}).setText(num2str(base.Params.(pf{i}))); end
        end
    end
    % tiles
    if isfield(base,'N'); ctrl.jTiles.setText(num2str(base.N)); end
    % direction / chirality
    if S.isDirac
        if isfield(base,'Chirality')
            switch base.Chirality
                case 1,  ctrl.jChir.setSelectedItem('+ (right)');
                case -1, ctrl.jChir.setSelectedItem('- (left)');
                otherwise, ctrl.jChir.setSelectedItem('None');
            end
        end
        if isfield(base,'Chiralities') && ~isempty(base.Chiralities); ctrl.jChiSplit.setSelected(true); end
    end
    % re-seed at the saved design vertex (if any) -> triggers a refresh
    iVertex = [];
    if isfield(loadBank,'Provenance') && isstruct(loadBank.Provenance) && isfield(loadBank.Provenance,'DesignVertex')
        iVertex = loadBank.Provenance.DesignVertex;
    end
    if ~isempty(iVertex)
        SetSeedVertex(panelName, iVertex);
    end
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
    Comment = sprintf('%s filterbank (%d tiles)', base.Kernel, numel(Tiles));
    db_add_filterbank(S.iSubject, S.EigenFile, fb, Comment);
    S.ctxFn.Close();
end

function OnCancel(panelName) %#ok<DEFNU>
    S = GetState(panelName);
    if ~isempty(S) && isfield(S,'ctxFn') && ~isempty(S.ctxFn); S.ctxFn.Close(); end
end
