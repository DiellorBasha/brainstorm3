function hFig = view_helmholtz(SrcResultsFile, varargin)
% VIEW_HELMHOLTZ: Dedicated Helmholtz/vorticity view of a Dirac source map. Shows a
% switchable signed scalar (curl/div/psi/phi/|field|) as a native 1-component results,
% with field-quiver + vortex-core overlays that follow the time cursor.
%
% USAGE:
%   hFig = view_helmholtz(SrcResultsFile)            % launch on a 3-component source results
%   view_helmholtz('SetScalar', hFig, scalarName)    % dispatched from panel_helmholtz
%   view_helmholtz('SetLayers', hFig, showCores, showQuiver)
%   view_helmholtz('Close', hFig)
%   view_helmholtz('RedrawOverlays', hFig)
% Authors: Diellor Basha, 2026
    global GlobalData;
    hFig = [];
    % dispatch string-first calls (from the panel / the overlay hook)
    if (nargin >= 1) && ischar(SrcResultsFile) && any(strcmp(SrcResultsFile, {'SetScalar','SetLayers','Close','RedrawOverlays'}))
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
    J  = double(bst_memory('GetResultsValues', iDS, iResult, [], [], 0));
    if size(J,1) ~= 3*nV
        bst_progress('stop'); bst_error('Field size mismatch with surface.', 'Helmholtz view', 0); return;
    end
    bst_progress('text', 'Helmholtz decomposition...');
    H = bst_dirac_helmholtz(Dirac, LBO, Surf, J);
    bst_progress('stop');

    % dedicated figure: show the default scalar (Curl) as a 1-comp results node
    tmpFile = i_make_scalar_results(SurfaceFile, H.Curl, R.Time);
    [hFig, iDSp, iResp] = i_open(SurfaceFile, tmpFile);
    if isempty(hFig); return; end
    St = struct('srcField',J, 'H',H, 'Scalar','Curl', 'ShowCores',true, 'ShowQuiver',false, ...
                'TmpFile',tmpFile, 'iDS',iDSp, 'iResult',iResp);
    setappdata(hFig, 'HelmholtzState', St);

    % overlays follow time via the CustomOverlayFcn hook (panel_surface UpdateSurfaceData)
    setappdata(hFig, 'CustomOverlayFcn', @(h) RedrawOverlays(h));
    set(hFig, 'CloseRequestFcn', @(h,e) Close(h));
    RedrawOverlays(hFig);

    % dock the control panel
    bstPanel = panel_helmholtz('CreatePanel', hFig);
    gui_show(bstPanel, 'BrainstormTab', 'tools');
    try, gui_brainstorm('SetSelectedTab', 'Helmholtz', 0); catch, end %#ok<CTCH>
end

%% ===== set the displayed scalar (called by the panel) =====
function SetScalar(hFig, scalarName) %#ok<DEFNU>
    global GlobalData;
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.Scalar = scalarName; setappdata(hFig, 'HelmholtzState', St);
    GlobalData.DataSet(St.iDS).Results(St.iResult).ImageGridAmp = St.H.(scalarName);
    TessInfo = getappdata(hFig,'Surface'); for k=1:numel(TessInfo); TessInfo(k).DataMinMax=[]; end
    setappdata(hFig,'Surface',TessInfo);
    panel_surface('UpdateSurfaceData', hFig);     % recolors + fires RedrawOverlays
    panel_surface('UpdateSurfaceColormap', hFig);
end

function SetLayers(hFig, showCores, showQuiver) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.ShowCores = showCores; St.ShowQuiver = showQuiver;
    setappdata(hFig, 'HelmholtzState', St);
    RedrawOverlays(hFig);
end

%% ===== overlay redraw (cores + quiver) at the current time =====
function RedrawOverlays(hFig)
    if isempty(hFig) || ~ishandle(hFig); return; end
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    hAx = findobj(hFig,'-depth',1,'Tag','Axes3D'); if isempty(hAx); return; end
    hAx = hAx(1);
    delete(findobj(hAx,'Tag','HelmholtzCore'));
    delete(findobj(hAx,'Tag','HelmholtzQuiver'));
    iT = i_current_time(St.iDS, St.iResult);
    TessInfo = getappdata(hFig,'Surface'); if isempty(TessInfo)||~ishandle(TessInfo(1).hPatch); return; end
    V = get(TessInfo(1).hPatch, 'Vertices');
    bb = max(V,[],1)-min(V,[],1);
    % cores
    if St.ShowCores
        cores = St.H.Cores{min(iT, numel(St.H.Cores))};
        for k = 1:numel(cores)
            v = cores(k).iVertex; col = [1 0 0]; if cores(k).charge < 0; col = [0 0 1]; end
            line('Parent',hAx,'XData',V(v,1),'YData',V(v,2),'ZData',V(v,3), 'Marker','o', ...
                'MarkerSize',10,'MarkerFaceColor',col,'MarkerEdgeColor','k','LineStyle','none', ...
                'Tag','HelmholtzCore','Clipping','off');
        end
        i_count_readout(St.H.Cores{min(iT,numel(St.H.Cores))});
    end
    % field quiver (subsampled), optional
    if St.ShowQuiver
        nV = size(V,1); step = max(1, round(nV/2000)); vi = (1:step:nV)';
        J = St.srcField; L = 0.05*norm(bb);
        d = [J(3*(vi-1)+1, iT), J(3*(vi-1)+2, iT), J(3*(vi-1)+3, iT)];
        nd = sqrt(sum(d.^2,2)); d = d ./ max(nd,eps);
        P = V(vi,:); Q = P + L*d;
        line('Parent',hAx, 'XData',[P(:,1) Q(:,1) nan(numel(vi),1)]', ...
             'YData',[P(:,2) Q(:,2) nan(numel(vi),1)]', 'ZData',[P(:,3) Q(:,3) nan(numel(vi),1)]', ...
             'Color',[.2 .2 .2], 'Tag','HelmholtzQuiver', 'Clipping','off');
    end
end

%% ===== close (remove hook, temp node, figure) =====
function Close(hFig) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState');
    try, gui_hide('Helmholtz'); catch, end %#ok<CTCH>
    try, rmappdata(hFig, 'CustomOverlayFcn'); catch, end %#ok<CTCH>
    if ~isempty(St) && isfield(St,'TmpFile') && ~isempty(St.TmpFile)
        try, file_delete(file_fullpath(St.TmpFile), 1); db_reload_studies(-3); catch, end %#ok<CTCH>
    end
    try, set(hFig, 'CloseRequestFcn', ''); catch, end %#ok<CTCH>
    try, bst_figures('DeleteFigure', hFig, []); catch, delete(hFig); end %#ok<CTCH>
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
function iT = i_current_time(iDS, iResult)
    [~, iT] = bst_memory('GetTimeVector', iDS, iResult, 'CurrentTimeIndex');
    if isempty(iT) || iT < 1; iT = 1; end
end
function i_count_readout(cores)
    nPos = sum([cores.charge] > 0); nNeg = sum([cores.charge] < 0);
    try, panel_helmholtz('SetReadout', sprintf('%d vortices (+), %d antivortices (-), net %+d', nPos, nNeg, nPos-nNeg)); catch, end %#ok<CTCH>
end
