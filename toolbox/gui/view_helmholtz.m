function hFig = view_helmholtz(SrcResultsFile, varargin)
% VIEW_HELMHOLTZ: Helmholtz/Hodge component view of a Dirac source map. Opens the NATIVE
% unconstrained-source display and lets a panel choose which COMPONENT of the decomposition
% to show -- Total |J| / Irrotational grad(phi) / Solenoidal curl(psi). (No harmonic
% component: the cortex is genus-0, so the harmonic space is trivial.)
% Switching component swaps the cortex SCALAR colormap (to that component's potential); the
% quiver always shows the (smoothed) TOTAL source field. An optional Dirac-eigenmode smoothing
% low-passes the active frame before the decomposition. Active frame only (cached Cholesky
% factor; recomputed as the cursor moves).
%
% NOTE: singular-point detection (sources/sinks, vortices), the persistence gate, trajectory
% tracking and the time-derivative (velocity/acceleration) views were removed -- that
% machinery is being rebuilt from scratch on top of the atom system (bst_dynamics).
%
% USAGE:
%   hFig = view_helmholtz(SrcResultsFile)
%   view_helmholtz('SetComponent', hFig, name)             % 'Total'|'Irrot'|'Solen'
%   view_helmholtz('SetSmoothing', hFig, isOn, name, params)  % Dirac-eigenmode low-pass
%   view_helmholtz('Close', hFig)
%   view_helmholtz('UpdateFrame', hFig)
% Authors: Diellor Basha, 2026
    global GlobalData;
    hFig = [];
    VERBS = {'SetComponent','SetSmoothing','Close','UpdateFrame'};
    if (nargin >= 1) && ischar(SrcResultsFile) && any(strcmp(SrcResultsFile, VERBS))
        if ~strcmp(SrcResultsFile, 'Close') && ...
                (isempty(varargin) || isempty(varargin{1}) || ~all(ishandle(varargin{1})))
            return;
        end
        feval(SrcResultsFile, varargin{:});
        return;
    end

    [iDS, iResult] = bst_memory('GetDataSetResult', SrcResultsFile);
    if isempty(iResult)
        try, [iDS, iResult] = bst_memory('LoadResultsFileFull', SrcResultsFile); catch, iResult = []; end %#ok<CTCH>
    end
    if isempty(iResult)
        bst_error(['Could not load this source over time.' 10 ...
            'For a Dirac source, open the Helmholtz view on a recordings link (under a data block), not the shared kernel.'], ...
            'Helmholtz view', 0);
        return;
    end
    R = GlobalData.DataSet(iDS).Results(iResult);
    if isempty(R.nComponents) || (R.nComponents ~= 3)
        bst_error('Helmholtz view requires an unconstrained (3-component) source.', 'Helmholtz view', 0); return;
    end
    SurfaceFile = R.SurfaceFile;

    bst_progress('start', 'Helmholtz view', 'Loading operators...');
    Cov = bst_get_operator_node(SurfaceFile, 'Covariant');
    LBO   = bst_get_operator_node(SurfaceFile, 'Laplace-Beltrami');
    Surf  = in_tess_bst(SurfaceFile, 0);
    Mani  = tess_manifold(SurfaceFile);
    nV = size(Surf.Vertices,1);
    bst_progress('text', 'Factorizing the cotan operator...');
    Op = bst_helmholtz('Prepare', {Cov, LBO}, Mani, Surf, 'Domain','vertex');
    bst_progress('text', 'Loading Dirac eigenbasis...');
    EigenMat = tess_eigen(SurfaceFile, 'Dirac');
    OpMat    = in_bst_operator(EigenMat.OperatorFile);
    Lambda   = double(EigenMat.Lambda{1}(:));
    bst_progress('stop');

    [hFig, iDSf] = view_surface_data(SurfaceFile, SrcResultsFile, [], 'NewFigure');
    if isempty(hFig); return; end
    iTess = i_find_tess(hFig);
    % No amplitude threshold by default: the Data threshold slider drives the quiver +
    % colormap. Show the full field on open.
    TI = getappdata(hFig,'Surface');
    if iTess <= numel(TI); TI(iTess).DataThreshold = 0; setappdata(hFig,'Surface',TI); end
    try, panel_surface('UpdateSurfaceProperties'); catch, end %#ok<CTCH>

    St = struct('Op',Op, 'srcDS',iDSf, 'srcResult',iResult, 'Component','Total', ...
                'iTess',iTess, 'nV',nV, ...
                'EigenMat',EigenMat, 'Mass',{OpMat.Mass}, 'Lambda',Lambda, ...
                'Smooth',struct('on',false,'name','heat','params',struct()), ...
                'Cache',containers.Map('KeyType','double','ValueType','any'));
    setappdata(hFig, 'HelmholtzState', St);
    setappdata(hFig, 'CustomOverlayFcn', @(h) UpdateFrame(h));
    set(hFig, 'CloseRequestFcn', @(h,e) Close(h));
    UpdateFrame(hFig);
    i_attach_listener(hFig);   % re-gate the quiver on a hemisphere/resect toggle (Surfaces panel)

    gui_hide('Helmholtz');
    bstPanel = panel_helmholtz('CreatePanel', hFig, St.Lambda);
    gui_show(bstPanel, 'BrainstormTab', 'tools');
    try, gui_brainstorm('SetSelectedTab', 'Helmholtz', 0); catch, end %#ok<CTCH>
end

%% ===== active-frame decompose + per-component override (the time hook) =====
function UpdateFrame(hFig)
    if isempty(hFig) || ~ishandle(hFig); return; end
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    hAx = findobj(hFig,'-depth',1,'Tag','Axes3D'); if isempty(hAx); return; end
    hAx = hAx(1); %#ok<NASGU>
    TessInfo = getappdata(hFig,'Surface');
    if isempty(TessInfo) || (St.iTess > numel(TessInfo)) || ~ishandle(TessInfo(St.iTess).hPatch); return; end
    % Per-hemisphere visibility from the cortex patch alpha (Surfaces-panel resect / hemi
    % toggle): a hidden hemisphere also hides its quiver (mirrors view_manifold's sync).
    visV = i_visibleVerts(TessInfo(St.iTess).hPatch);
    setappdata(hFig, 'HelmholtzLastVis', visV);
    [~, iT] = bst_memory('GetTimeVector', St.srcDS, St.srcResult, 'CurrentTimeIndex');
    if isempty(iT) || iT < 1; iT = 1; end
    % Current unconstrained source vector at this frame
    Jt = double(bst_memory('GetResultsValues', St.srcDS, St.srcResult, [], iT, 0));
    if size(Jt,1) ~= 3*St.nV; return; end
    % Spatial eigenmode smoothing (low-pass the active frame before the decomposition)
    if St.Smooth.on
        g  = bst_eigfilter_kernel(St.Smooth.name, St.Smooth.params);
        Jt = real(bst_eigenfilter('Analysis', Jt, St.EigenMat, struct('Mass', {St.Mass}), g, struct()));
    end
    % Helmholtz-Hodge decomposition of the active frame (cores off)
    if isKey(St.Cache, iT)
        Ht = St.Cache(iT);
    else
        Ht = bst_helmholtz('Frame', St.Op, Jt, false);  St.Cache(iT) = Ht;
    end
    comp = i_component(Ht, St.Component);
    % --- cortex scalar + colormap ---
    TessInfo(St.iTess).Data = comp.Scal;
    if comp.Signed
        TessInfo(St.iTess).DataMinMax = i_minmax(comp.Scal);  TessInfo(St.iTess).ColormapType = 'stat2';
    else
        m = max(comp.Scal); if m <= 0; m = eps; end
        TessInfo(St.iTess).DataMinMax = [0 m];                TessInfo(St.iTess).ColormapType = 'source';
    end
    setappdata(hFig,'Surface',TessInfo);
    panel_surface('UpdateSurfaceColormap', hFig);
    % --- quiver is ALWAYS the (smoothed) total source field; component changes only the scalar ---
    Vq = Ht.Vtot;
    if numel(visV) == size(Vq,1); Vq(~visV,:) = 0; end   % zero-length arrows = hidden hemisphere
    setappdata(hFig, 'QuiverVectorOverride', Vq);
    try, figure_3d('SetShowSourceVectors', hFig, St.iTess, 1); catch, end %#ok<CTCH>
    try, figure_3d('PlotSourceVectors', hFig, St.iTess); catch, end %#ok<CTCH>
    i_readout(comp.Kind);
end

%% ===== panel actions =====
function SetComponent(hFig, name) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.Component = name; setappdata(hFig, 'HelmholtzState', St);
    if any(strcmp(name, {'Irrot','Solen'})); bst_colormaps('AddColormapToFigure', hFig, 'stat2');
    else; bst_colormaps('AddColormapToFigure', hFig, 'source'); end
    UpdateFrame(hFig);
end
function SetSmoothing(hFig, isOn, name, params) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.Smooth = struct('on',logical(isOn), 'name',name, 'params',params);
    St.Cache  = containers.Map('KeyType','double','ValueType','any');   % decompositions now stale
    setappdata(hFig, 'HelmholtzState', St);  UpdateFrame(hFig);
end

%% ===== close =====
function Close(hFig) %#ok<DEFNU>
    try, gui_hide('Helmholtz'); catch, end %#ok<CTCH>
    if ~isempty(hFig) && all(ishandle(hFig))
        try, lh = getappdata(hFig, 'HelmholtzPatchListener'); if ~isempty(lh); delete(lh); end; catch, end %#ok<CTCH>
        try, rmappdata(hFig, 'HelmholtzPatchListener'); catch, end %#ok<CTCH>
        try, rmappdata(hFig, 'CustomOverlayFcn'); catch, end %#ok<CTCH>
        try, rmappdata(hFig, 'QuiverVectorOverride'); catch, end %#ok<CTCH>
        try, set(hFig, 'CloseRequestFcn', ''); catch, end %#ok<CTCH>
    end
    try, bst_figures('DeleteFigure', hFig, []); catch, if ~isempty(hFig)&&all(ishandle(hFig)); delete(hFig); end; end %#ok<CTCH>
end

%% ===== per-hemisphere visibility (Surfaces-panel resect / hemi toggle) =====
function visV = i_visibleVerts(hPatch)
% Per-vertex visibility from the cortex patch's alpha MODE -- the same signal
% view_manifold/SyncScalar reads. Returns all-visible when nothing is resected.
    V  = get(hPatch, 'Vertices');
    F  = get(hPatch, 'Faces');
    A  = get(hPatch, 'FaceVertexAlphaData');
    EA = get(hPatch, 'EdgeAlpha');
    nV = size(V, 1);  nF = size(F, 1);
    if ischar(EA) && ~isempty(A)
        if numel(A) == nF                 % per-face ('flat')
            faceVis = A > 0;
        elseif numel(A) == nV             % per-vertex ('interp')
            vv = A(:) > 0;  faceVis = all(vv(F), 2);
        else
            faceVis = true(nF, 1);
        end
    elseif isnumeric(EA) && isscalar(EA) && (EA <= 0)
        faceVis = false(nF, 1);           % uniformly transparent
    else
        faceVis = true(nF, 1);            % opaque -> everything visible
    end
    visV = false(nV, 1);
    visV(unique(F(faceVis, :))) = true;
end

%% ===== re-gate the quiver when the cortex patch is redrawn (hemisphere/resect toggle) =====
function i_attach_listener(hFig)
    St = getappdata(hFig, 'HelmholtzState');  if isempty(St); return; end
    TessInfo = getappdata(hFig, 'Surface');
    if isempty(TessInfo) || (St.iTess > numel(TessInfo)) || ~ishandle(TessInfo(St.iTess).hPatch); return; end
    lh = addlistener(TessInfo(St.iTess).hPatch, 'MarkedClean', @(s,e) i_on_patch_clean(hFig));
    setappdata(hFig, 'HelmholtzPatchListener', lh);
end

function i_on_patch_clean(hFig)
    if isempty(hFig) || ~ishandle(hFig); return; end
    busy = getappdata(hFig, 'HelmholtzInSync');            % guard the re-entrant MarkedClean
    if ~isempty(busy) && busy; return; end                %   that UpdateFrame itself triggers
    St = getappdata(hFig, 'HelmholtzState');  if isempty(St); return; end
    TessInfo = getappdata(hFig, 'Surface');
    if isempty(TessInfo) || (St.iTess > numel(TessInfo)) || ~ishandle(TessInfo(St.iTess).hPatch); return; end
    visV = i_visibleVerts(TessInfo(St.iTess).hPatch);
    if isequal(visV, getappdata(hFig, 'HelmholtzLastVis')); return; end   % visibility unchanged -> no-op
    setappdata(hFig, 'HelmholtzInSync', true);
    try, UpdateFrame(hFig); catch, end %#ok<CTCH>
    setappdata(hFig, 'HelmholtzInSync', false);
end

%% ===== helpers =====
function c = i_component(Ht, name)
    % No harmonic component: the cortex is genus-0 (topological sphere), so H^1 is
    % trivial and the harmonic part is identically zero (numerical residual only).
    switch name
        case 'Irrot', c = struct('Vec',Ht.Virr, 'Scal',Ht.Phi,  'Signed',true,  'Kind','source');
        case 'Solen', c = struct('Vec',Ht.Vsol, 'Scal',Ht.Psi,  'Signed',true,  'Kind','vortex');
        otherwise,    c = struct('Vec',Ht.Vtot, 'Scal',Ht.Fmag, 'Signed',false, 'Kind','total');
    end
end
function iTess = i_find_tess(hFig)
    TessInfo = getappdata(hFig, 'Surface');
    iTess = find(arrayfun(@(t) ~isempty(t.DataSource) && strcmpi(t.DataSource.Type,'Source'), TessInfo), 1);
    if isempty(iTess); iTess = 1; end
end
function mm = i_minmax(scal)
    m = max(abs(scal));  if m == 0; m = eps; end
    mm = [-m, m];
end
function i_readout(kind)
    switch kind
        case 'vortex', txt = 'solenoidal: stream function \Psi';
        case 'source', txt = 'irrotational: potential \Phi';
        otherwise,     txt = 'total field |J|';
    end
    try, panel_helmholtz('SetReadout', txt); catch, end %#ok<CTCH>
end
