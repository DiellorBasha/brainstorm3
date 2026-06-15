function hFig = view_helmholtz(SrcResultsFile, varargin)
% VIEW_HELMHOLTZ: Dedicated Helmholtz/vorticity view of a Dirac source map. Shows a
% switchable signed scalar (curl/div/psi/phi/|field|) on the cortex with field-quiver +
% vortex-core overlays. The decomposition runs ON THE ACTIVE FRAME ONLY: the cotan
% Poisson factor is built once (bst_dirac_helmholtz Prepare) and each time the cursor
% moves the displayed frame is decomposed on the fly (Frame) and cached. The whole-series
% version belongs in a batch process_ (not here).
%
% USAGE:
%   hFig = view_helmholtz(SrcResultsFile)            % launch on a 3-component source results
%   view_helmholtz('SetScalar', hFig, scalarName)    % dispatched from panel_helmholtz
%   view_helmholtz('SetLayers', hFig, showCores, showQuiver)
%   view_helmholtz('Close', hFig)
%   view_helmholtz('UpdateFrame', hFig)              % the time-cursor overlay hook
% Authors: Diellor Basha, 2026
    global GlobalData;
    hFig = [];
    % dispatch string-first calls (from the panel / the overlay hook)
    if (nargin >= 1) && ischar(SrcResultsFile) && any(strcmp(SrcResultsFile, {'SetScalar','SetLayers','Close','UpdateFrame'}))
        feval(SrcResultsFile, varargin{:});
        return;
    end

    [iDS, iResult] = bst_memory('GetDataSetResult', SrcResultsFile);
    if isempty(iResult)
        bst_error('Could not resolve the source results.', 'Helmholtz view', 0); return;
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
    % size sanity: the in-memory field (kernel or full) must match the surface
    if ~isempty(R.ImageGridAmp) && (size(R.ImageGridAmp,1) ~= 3*nV)
        bst_progress('stop'); bst_error('Field size mismatch with surface.', 'Helmholtz view', 0); return;
    end
    bst_progress('text', 'Factorizing the cotan operator...');
    Op = bst_dirac_helmholtz('Prepare', Dirac, LBO, Surf);   % ONE Cholesky factor, reused per frame

    % decompose the CURRENTLY-DISPLAYED frame only
    fullTime = bst_memory('GetTimeVector', iDS, iResult);
    nT = numel(fullTime);
    [~, iT0] = bst_memory('GetTimeVector', iDS, iResult, 'CurrentTimeIndex');
    if isempty(iT0) || iT0 < 1; iT0 = 1; end
    Jt0 = double(bst_memory('GetResultsValues', iDS, iResult, [], iT0, 0));   % [3nV x 1]
    Ht0 = bst_dirac_helmholtz('Frame', Op, Jt0);
    bst_progress('stop');

    % dedicated figure: a 1-component results acts as the time-aware, colormapped host.
    % Its ImageGridAmp is a placeholder (seeded non-zero so the native pipeline never warns
    % "all values null"); the displayed scalar is set per frame via TessInfo.Data below.
    Cache = containers.Map('KeyType','double', 'ValueType','any');
    Cache(iT0) = Ht0;
    tmpFile = i_make_scalar_results(SurfaceFile, repmat(Ht0.Curl, 1, nT), fullTime);
    [hFig, iDSp, iResp] = i_open(SurfaceFile, tmpFile);
    if isempty(hFig); return; end
    iTess = i_find_tess(hFig);
    St = struct('Op',Op, 'srcDS',iDS, 'srcResult',iResult, 'Scalar','Curl', ...
                'ShowCores',true, 'ShowQuiver',false, 'TmpFile',tmpFile, ...
                'iDS',iDSp, 'iResult',iResp, 'nV',nV, 'iTess',iTess, 'Cache',Cache);
    setappdata(hFig, 'HelmholtzState', St);

    % the active frame is decomposed + recolored here; this also rides the time cursor
    % (panel_surface UpdateSurfaceData fires CustomOverlayFcn on every cursor move)
    setappdata(hFig, 'CustomOverlayFcn', @(h) UpdateFrame(h));
    set(hFig, 'CloseRequestFcn', @(h,e) Close(h));
    UpdateFrame(hFig);

    % dock the control panel
    bstPanel = panel_helmholtz('CreatePanel', hFig);
    gui_show(bstPanel, 'BrainstormTab', 'tools');
    try, gui_brainstorm('SetSelectedTab', 'Helmholtz', 0); catch, end %#ok<CTCH>
end

%% ===== active-frame decompose + recolor + overlays (the time hook) =====
function UpdateFrame(hFig)
    if isempty(hFig) || ~ishandle(hFig); return; end
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    hAx = findobj(hFig,'-depth',1,'Tag','Axes3D'); if isempty(hAx); return; end
    hAx = hAx(1);
    TessInfo = getappdata(hFig,'Surface');
    if isempty(TessInfo) || (St.iTess > numel(TessInfo)) || ~ishandle(TessInfo(St.iTess).hPatch); return; end
    % current frame -> decompose on demand (cache visited frames)
    [~, iT] = bst_memory('GetTimeVector', St.srcDS, St.srcResult, 'CurrentTimeIndex');
    if isempty(iT) || iT < 1; iT = 1; end
    if isKey(St.Cache, iT)
        Ht = St.Cache(iT);
    else
        Jt = double(bst_memory('GetResultsValues', St.srcDS, St.srcResult, [], iT, 0));
        Ht = bst_dirac_helmholtz('Frame', St.Op, Jt);
        St.Cache(iT) = Ht;   % containers.Map is a handle: persists without setappdata
    end
    % recolor the cortex with the selected scalar for this frame
    scal = Ht.(St.Scalar);
    TessInfo(St.iTess).Data = scal;
    TessInfo(St.iTess).DataMinMax = i_minmax(scal, St.Scalar);
    setappdata(hFig, 'Surface', TessInfo);
    panel_surface('UpdateSurfaceColormap', hFig);     % recolor only -- does NOT re-fire this hook
    % overlays
    V = get(TessInfo(St.iTess).hPatch, 'Vertices');
    delete(findobj(hAx,'Tag','HelmholtzCore'));
    delete(findobj(hAx,'Tag','HelmholtzQuiver'));
    if St.ShowCores
        for k = 1:numel(Ht.Cores)
            v = Ht.Cores(k).iVertex; col = [1 0 0]; if Ht.Cores(k).charge < 0; col = [0 0 1]; end
            line('Parent',hAx,'XData',V(v,1),'YData',V(v,2),'ZData',V(v,3), 'Marker','o', ...
                'MarkerSize',10,'MarkerFaceColor',col,'MarkerEdgeColor','k','LineStyle','none', ...
                'Tag','HelmholtzCore','Clipping','off');
        end
        i_count_readout(Ht.Cores);
    end
    if St.ShowQuiver
        nVv = size(V,1); step = max(1, round(nVv/2000)); vi = (1:step:nVv)';
        bb = max(V,[],1)-min(V,[],1);  L = 0.05*norm(bb);
        Jt = double(bst_memory('GetResultsValues', St.srcDS, St.srcResult, [], iT, 0));
        d = [Jt(3*(vi-1)+1), Jt(3*(vi-1)+2), Jt(3*(vi-1)+3)];
        nd = sqrt(sum(d.^2,2)); d = d ./ max(nd,eps);
        P = V(vi,:); Q = P + L*d;
        line('Parent',hAx, 'XData',[P(:,1) Q(:,1) nan(numel(vi),1)]', ...
             'YData',[P(:,2) Q(:,2) nan(numel(vi),1)]', 'ZData',[P(:,3) Q(:,3) nan(numel(vi),1)]', ...
             'Color',[.2 .2 .2], 'Tag','HelmholtzQuiver', 'Clipping','off');
    end
end

%% ===== set the displayed scalar / layers (called by the panel) =====
function SetScalar(hFig, scalarName) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.Scalar = scalarName; setappdata(hFig, 'HelmholtzState', St);
    UpdateFrame(hFig);
end
function SetLayers(hFig, showCores, showQuiver) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.ShowCores = showCores; St.ShowQuiver = showQuiver;
    setappdata(hFig, 'HelmholtzState', St);
    UpdateFrame(hFig);
end

%% ===== close (delete figure first while registered, then the temp node) =====
function Close(hFig) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState');
    try, gui_hide('Helmholtz'); catch, end %#ok<CTCH>
    try, rmappdata(hFig, 'CustomOverlayFcn'); catch, end %#ok<CTCH>
    try, set(hFig, 'CloseRequestFcn', ''); catch, end %#ok<CTCH>
    try, bst_figures('DeleteFigure', hFig, []); catch, delete(hFig); end %#ok<CTCH>
    if ~isempty(St) && isfield(St,'TmpFile') && ~isempty(St.TmpFile)
        try, file_delete(file_fullpath(St.TmpFile), 1); db_reload_studies(-3); catch, end %#ok<CTCH>
    end
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
function tmpFile = i_make_scalar_results(SurfaceFile, scalar, Time)
    R = db_template('resultsmat');
    R.ImageGridAmp = scalar; R.ImagingKernel = []; R.nComponents = 1;
    R.Time = Time; R.HeadModelType = 'surface'; R.SurfaceFile = file_short(SurfaceFile);
    R.Comment = 'Helmholtz scalar';
    tmpFile = db_add(-3, R);
end
function [hFig, iDS, iResult] = i_open(SurfaceFile, tmpFile)
    [hFig, iDS] = view_surface_data(SurfaceFile, tmpFile, [], 'NewFigure');
    [iDS2, iResult] = bst_memory('GetDataSetResult', tmpFile);
    if ~isempty(iDS2); iDS = iDS2; end
end
function iTess = i_find_tess(hFig)
    TessInfo = getappdata(hFig, 'Surface');
    iTess = find(arrayfun(@(t) ~isempty(t.DataSource) && strcmpi(t.DataSource.Type,'Source'), TessInfo), 1);
    if isempty(iTess); iTess = 1; end
end
function mm = i_minmax(scal, scalarName)
    if strcmp(scalarName, 'Fmag')
        mm = [0, max(max(scal), eps)];                 % magnitude: one-sided
    else
        m = max(abs(scal));  if m == 0; m = eps; end
        mm = [-m, m];                                  % signed: symmetric diverging
    end
end
function i_count_readout(cores)
    nPos = sum([cores.charge] > 0); nNeg = sum([cores.charge] < 0);
    try, panel_helmholtz('SetReadout', sprintf('%d vortices (+), %d antivortices (-), net %+d', nPos, nNeg, nPos-nNeg)); catch, end %#ok<CTCH>
end
