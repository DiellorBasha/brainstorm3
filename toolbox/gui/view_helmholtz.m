function hFig = view_helmholtz(SrcResultsFile, varargin)
% VIEW_HELMHOLTZ: Helmholtz/vorticity view of a Dirac source map. Opens the NATIVE
% unconstrained-source display (the real vector field + norm, exactly like "Display on
% cortex") and lets a control panel swap the scalar that colors the cortex:
%   Norm (|J|, native)  /  Curl.n (vorticity)  /  Divergence  /  Stream psi  /  Potential phi
% The signed Helmholtz scalars are decomposed on the ACTIVE FRAME ONLY (cached Cholesky
% factor; recomputed as the time cursor moves) and shown with a diverging colormap; the
% native source vectors and the optional vortex-core markers stay in sync.
%
% USAGE:
%   hFig = view_helmholtz(SrcResultsFile)            % launch on a 3-component Dirac source
%   view_helmholtz('SetScalar', hFig, scalarName)    % 'Norm'|'Curl'|'Div'|'Psi'|'Phi'
%   view_helmholtz('SetCores',  hFig, showCores)
%   view_helmholtz('Close', hFig)
%   view_helmholtz('UpdateFrame', hFig)              % the time-cursor hook
% Authors: Diellor Basha, 2026
    global GlobalData;
    hFig = [];
    % dispatch string-first calls; the frame/overlay updaters need a live figure
    if (nargin >= 1) && ischar(SrcResultsFile) && any(strcmp(SrcResultsFile, {'SetScalar','SetCores','Close','UpdateFrame'}))
        if any(strcmp(SrcResultsFile, {'SetScalar','SetCores','UpdateFrame'})) && ...
                (isempty(varargin) || isempty(varargin{1}) || ~all(ishandle(varargin{1})))
            return;
        end
        feval(SrcResultsFile, varargin{:});
        return;
    end

    % Resolve + load the source (the tree node may not be loaded). Use the FULL loader so
    % a Dirac kernel link resolves its imaging kernel + recordings; a bare shared kernel
    % (no recordings) errors -> tell the user to open a recordings link.
    [iDS, iResult] = bst_memory('GetDataSetResult', SrcResultsFile);
    if isempty(iResult)
        try
            [iDS, iResult] = bst_memory('LoadResultsFileFull', SrcResultsFile);
        catch
            iResult = [];
        end
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
    Dirac = i_op(SurfaceFile, 'Dirac');
    LBO   = i_op(SurfaceFile, 'Laplace-Beltrami');
    Surf  = in_tess_bst(SurfaceFile, 0);
    nV = size(Surf.Vertices,1);
    bst_progress('text', 'Factorizing the cotan operator...');
    Op = bst_dirac_helmholtz('Prepare', Dirac, LBO, Surf);   % ONE Cholesky factor, reused per frame
    bst_progress('stop');

    % open the NATIVE source display (vector field + norm) and show the vectors
    [hFig, iDSf] = view_surface_data(SurfaceFile, SrcResultsFile, [], 'NewFigure');
    if isempty(hFig); return; end
    iTess = i_find_tess(hFig);
    try, figure_3d('SetShowSourceVectors', hFig, iTess, 1); catch, end %#ok<CTCH>

    St = struct('Op',Op, 'srcDS',iDSf, 'srcResult',iResult, 'Scalar','Norm', ...
                'ShowCores',true, 'iTess',iTess, 'nV',nV, ...
                'Cache',containers.Map('KeyType','double','ValueType','any'));
    setappdata(hFig, 'HelmholtzState', St);
    setappdata(hFig, 'CustomOverlayFcn', @(h) UpdateFrame(h));   % rides the time cursor
    set(hFig, 'CloseRequestFcn', @(h,e) Close(h));
    UpdateFrame(hFig);

    % dock the control panel (clear any stale prior instance first)
    gui_hide('Helmholtz');
    bstPanel = panel_helmholtz('CreatePanel', hFig);
    gui_show(bstPanel, 'BrainstormTab', 'tools');
    try, gui_brainstorm('SetSelectedTab', 'Helmholtz', 0); catch, end %#ok<CTCH>
end

%% ===== active-frame override of the cortex scalar + vectors + cores (the time hook) =====
function UpdateFrame(hFig)
    if isempty(hFig) || ~ishandle(hFig); return; end
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    hAx = findobj(hFig,'-depth',1,'Tag','Axes3D'); if isempty(hAx); return; end
    hAx = hAx(1);
    TessInfo = getappdata(hFig,'Surface');
    if isempty(TessInfo) || (St.iTess > numel(TessInfo)) || ~ishandle(TessInfo(St.iTess).hPatch); return; end
    % current source frame
    [~, iT] = bst_memory('GetTimeVector', St.srcDS, St.srcResult, 'CurrentTimeIndex');
    if isempty(iT) || iT < 1; iT = 1; end
    Jt = double(bst_memory('GetResultsValues', St.srcDS, St.srcResult, [], iT, 0));   % [3nV x 1]
    if size(Jt,1) ~= 3*St.nV; return; end
    isSigned = ~strcmp(St.Scalar, 'Norm');
    % decompose this frame only when a signed scalar or the cores are needed (cache it)
    Ht = [];
    if isSigned || St.ShowCores
        if isKey(St.Cache, iT); Ht = St.Cache(iT);
        else; Ht = bst_dirac_helmholtz('Frame', St.Op, Jt); St.Cache(iT) = Ht; end
    end
    % --- cortex scalar: override with the signed Helmholtz field (Norm stays native) ---
    if isSigned
        scal = Ht.(St.Scalar);
        TessInfo(St.iTess).Data        = scal;
        TessInfo(St.iTess).DataMinMax  = i_minmax(scal);
        TessInfo(St.iTess).ColormapType = 'stat2';     % signed diverging
        setappdata(hFig, 'Surface', TessInfo);
        panel_surface('UpdateSurfaceColormap', hFig);  % recolor only (no overlay re-fire)
        TessInfo = getappdata(hFig,'Surface');
    end
    % --- native vectors: gate by their own magnitude regardless of the colored scalar ---
    setappdata(hFig, 'QuiverVectorOverride', reshape(Jt, 3, [])');
    try, figure_3d('PlotSourceVectors', hFig, St.iTess); catch, end %#ok<CTCH>
    % --- vortex-core markers ---
    delete(findobj(hAx,'Tag','HelmholtzCore'));
    if St.ShowCores && ~isempty(Ht)
        V = get(TessInfo(St.iTess).hPatch, 'Vertices');
        for k = 1:numel(Ht.Cores)
            v = Ht.Cores(k).iVertex; col = [1 0 0]; if Ht.Cores(k).charge < 0; col = [0 0 1]; end
            line('Parent',hAx,'XData',V(v,1),'YData',V(v,2),'ZData',V(v,3), 'Marker','o', ...
                'MarkerSize',9,'MarkerFaceColor',col,'MarkerEdgeColor','k','LineStyle','none', ...
                'Tag','HelmholtzCore','Clipping','off');
        end
        i_count_readout(Ht.Cores);
    else
        i_count_readout([]);
    end
end

%% ===== panel actions =====
function SetScalar(hFig, scalarName) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.Scalar = scalarName; setappdata(hFig, 'HelmholtzState', St);
    if strcmp(scalarName, 'Norm')
        % restore the native norm coloring (source colormap), then refresh overlays
        bst_colormaps('AddColormapToFigure', hFig, 'source');
        TessInfo = getappdata(hFig,'Surface');
        TessInfo(St.iTess).ColormapType = 'source';
        TessInfo(St.iTess).DataMinMax   = [];
        setappdata(hFig,'Surface',TessInfo);
        panel_surface('UpdateSurfaceData', hFig);      % recompute norm + recolor (fires UpdateFrame)
    else
        bst_colormaps('AddColormapToFigure', hFig, 'stat2');
        UpdateFrame(hFig);
    end
end
function SetCores(hFig, showCores) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.ShowCores = showCores; setappdata(hFig, 'HelmholtzState', St);
    UpdateFrame(hFig);
end

%% ===== close (remove hooks, delete figure) =====
function Close(hFig) %#ok<DEFNU>
    try, gui_hide('Helmholtz'); catch, end %#ok<CTCH>
    if ~isempty(hFig) && all(ishandle(hFig))
        try, rmappdata(hFig, 'CustomOverlayFcn'); catch, end %#ok<CTCH>
        try, rmappdata(hFig, 'QuiverVectorOverride'); catch, end %#ok<CTCH>
        try, set(hFig, 'CloseRequestFcn', ''); catch, end %#ok<CTCH>
    end
    try, bst_figures('DeleteFigure', hFig, []); catch, if ~isempty(hFig)&&all(ishandle(hFig)); delete(hFig); end; end %#ok<CTCH>
end

%% ===== helpers =====
function Op = i_op(SurfaceFile, variant)
    [sSubject,~,iSurf] = bst_get('SurfaceFile', SurfaceFile);
    Op = [];
    if ~isempty(iSurf) && isfield(sSubject.Surface(iSurf),'Operator')
        for k = 1:numel(sSubject.Surface(iSurf).Operator)
            S = load(file_fullpath(sSubject.Surface(iSurf).Operator(k).FileName));
            if strcmpi(S.Variant, variant); Op = S; break; end
        end
    end
    if isempty(Op); tess_operators(SurfaceFile, variant); Op = i_op(SurfaceFile, variant); end
end
function iTess = i_find_tess(hFig)
    TessInfo = getappdata(hFig, 'Surface');
    iTess = find(arrayfun(@(t) ~isempty(t.DataSource) && strcmpi(t.DataSource.Type,'Source'), TessInfo), 1);
    if isempty(iTess); iTess = 1; end
end
function mm = i_minmax(scal)
    m = max(abs(scal));  if m == 0; m = eps; end
    mm = [-m, m];                                      % signed -> symmetric diverging
end
function i_count_readout(cores)
    if isempty(cores); txt = 'no cores'; else
        nPos = sum([cores.charge] > 0); nNeg = sum([cores.charge] < 0);
        txt = sprintf('%d vortices (+), %d antivortices (-), net %+d', nPos, nNeg, nPos-nNeg);
    end
    try, panel_helmholtz('SetReadout', txt); catch, end %#ok<CTCH>
end
