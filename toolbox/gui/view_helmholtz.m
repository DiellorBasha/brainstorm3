function hFig = view_helmholtz(SrcResultsFile, varargin)
% VIEW_HELMHOLTZ: Helmholtz/Hodge component view of a Dirac source map. Opens the NATIVE
% unconstrained-source display and lets a panel choose which COMPONENT of the decomposition
% to show -- Total |J| / Irrotational grad(phi) / Solenoidal curl(psi) / Harmonic h. Each
% state swaps the quiver to that component's vector field and the cortex colormap to its
% scalar potential, with component-aware singular-point markers. Active frame only (cached
% Cholesky factor; recomputed as the cursor moves).
%
% USAGE:
%   hFig = view_helmholtz(SrcResultsFile)
%   view_helmholtz('SetComponent', hFig, name)   % 'Total'|'Irrot'|'Solen'|'Harm'
%   view_helmholtz('SetVectors', hFig, show)
%   view_helmholtz('SetMarkers', hFig, show)
%   view_helmholtz('Close', hFig)
%   view_helmholtz('UpdateFrame', hFig)
% Authors: Diellor Basha, 2026
    global GlobalData;
    hFig = [];
    if (nargin >= 1) && ischar(SrcResultsFile) && any(strcmp(SrcResultsFile, {'SetComponent','SetVectors','SetMarkers','Close','UpdateFrame'}))
        if any(strcmp(SrcResultsFile, {'SetComponent','SetVectors','SetMarkers','UpdateFrame'})) && ...
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
    Dirac = i_op(SurfaceFile, 'Dirac');
    LBO   = i_op(SurfaceFile, 'Laplace-Beltrami');
    Surf  = in_tess_bst(SurfaceFile, 0);
    nV = size(Surf.Vertices,1);
    bst_progress('text', 'Factorizing the cotan operator...');
    Op = bst_dirac_helmholtz('Prepare', Dirac, LBO, Surf);
    bst_progress('stop');

    [hFig, iDSf] = view_surface_data(SurfaceFile, SrcResultsFile, [], 'NewFigure');
    if isempty(hFig); return; end
    iTess = i_find_tess(hFig);

    St = struct('Op',Op, 'srcDS',iDSf, 'srcResult',iResult, 'Component','Total', ...
                'ShowVectors',true, 'ShowMarkers',true, 'iTess',iTess, 'nV',nV, ...
                'Cache',containers.Map('KeyType','double','ValueType','any'));
    setappdata(hFig, 'HelmholtzState', St);
    setappdata(hFig, 'CustomOverlayFcn', @(h) UpdateFrame(h));
    set(hFig, 'CloseRequestFcn', @(h,e) Close(h));
    UpdateFrame(hFig);

    gui_hide('Helmholtz');
    bstPanel = panel_helmholtz('CreatePanel', hFig);
    gui_show(bstPanel, 'BrainstormTab', 'tools');
    try, gui_brainstorm('SetSelectedTab', 'Helmholtz', 0); catch, end %#ok<CTCH>
end

%% ===== active-frame decompose + per-component override (the time hook) =====
function UpdateFrame(hFig)
    if isempty(hFig) || ~ishandle(hFig); return; end
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    hAx = findobj(hFig,'-depth',1,'Tag','Axes3D'); if isempty(hAx); return; end
    hAx = hAx(1);
    TessInfo = getappdata(hFig,'Surface');
    if isempty(TessInfo) || (St.iTess > numel(TessInfo)) || ~ishandle(TessInfo(St.iTess).hPatch); return; end
    [~, iT] = bst_memory('GetTimeVector', St.srcDS, St.srcResult, 'CurrentTimeIndex');
    if isempty(iT) || iT < 1; iT = 1; end
    Jt = double(bst_memory('GetResultsValues', St.srcDS, St.srcResult, [], iT, 0));
    if size(Jt,1) ~= 3*St.nV; return; end
    if isKey(St.Cache, iT); Ht = St.Cache(iT);
    else; Ht = bst_dirac_helmholtz('Frame', St.Op, Jt); St.Cache(iT) = Ht; end
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
    TessInfo = getappdata(hFig,'Surface');
    % --- component vector field as the native quiver ---
    if St.ShowVectors
        setappdata(hFig, 'QuiverVectorOverride', comp.Vec);
        try, figure_3d('SetShowSourceVectors', hFig, St.iTess, 1); catch, end %#ok<CTCH>
        try, figure_3d('PlotSourceVectors', hFig, St.iTess); catch, end %#ok<CTCH>
    else
        try, figure_3d('SetShowSourceVectors', hFig, St.iTess, 0); catch, end %#ok<CTCH>
    end
    % --- component-aware markers ---
    delete(findobj(hAx,'Tag','HelmholtzCore'));
    if St.ShowMarkers && ~isempty(comp.Markers)
        V = get(TessInfo(St.iTess).hPatch, 'Vertices');
        for k = 1:numel(comp.Markers)
            v = comp.Markers(k).iVertex; col = [1 0 0]; if comp.Markers(k).charge < 0; col = [0 0 1]; end
            line('Parent',hAx,'XData',V(v,1),'YData',V(v,2),'ZData',V(v,3), 'Marker','o', ...
                'MarkerSize',9,'MarkerFaceColor',col,'MarkerEdgeColor','k','LineStyle','none', ...
                'Tag','HelmholtzCore','Clipping','off');
        end
    end
    i_readout(comp, Ht);
end

%% ===== panel actions =====
function SetComponent(hFig, name) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.Component = name; setappdata(hFig, 'HelmholtzState', St);
    if any(strcmp(name, {'Irrot','Solen'})); bst_colormaps('AddColormapToFigure', hFig, 'stat2');
    else; bst_colormaps('AddColormapToFigure', hFig, 'source'); end
    UpdateFrame(hFig);
end
function SetVectors(hFig, show) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.ShowVectors = show; setappdata(hFig, 'HelmholtzState', St);  UpdateFrame(hFig);
end
function SetMarkers(hFig, show) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.ShowMarkers = show; setappdata(hFig, 'HelmholtzState', St);  UpdateFrame(hFig);
end

%% ===== close =====
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
function c = i_component(Ht, name)
    switch name
        case 'Irrot', c = struct('Vec',Ht.Virr, 'Scal',Ht.Phi,  'Signed',true,  'Markers',Ht.Sources, 'Kind','source', 'HarmFrac',Ht.HarmFrac);
        case 'Solen', c = struct('Vec',Ht.Vsol, 'Scal',Ht.Psi,  'Signed',true,  'Markers',Ht.Cores,   'Kind','vortex', 'HarmFrac',Ht.HarmFrac);
        case 'Harm',  c = struct('Vec',Ht.Vharm,'Scal',Ht.Hmag, 'Signed',false, 'Markers',[],          'Kind','harm',   'HarmFrac',Ht.HarmFrac);
        otherwise,    c = struct('Vec',Ht.Vtot, 'Scal',Ht.Fmag, 'Signed',false, 'Markers',[],          'Kind','total',  'HarmFrac',Ht.HarmFrac);
    end
end
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
    mm = [-m, m];
end
function i_readout(comp, Ht)
    switch comp.Kind
        case 'vortex'
            if isempty(comp.Markers); txt = '0 vortices, 0 antivortices';
            else
                np = sum([comp.Markers.charge] > 0); nn = sum([comp.Markers.charge] < 0);
                txt = sprintf('%d vortices (+), %d antivortices (-), net %+d', np, nn, np-nn);
            end
        case 'source'
            if isempty(comp.Markers); txt = '0 sources, 0 sinks';
            else
                np = sum([comp.Markers.charge] > 0); nn = sum([comp.Markers.charge] < 0);
                txt = sprintf('%d sources (+), %d sinks (-), net %+d', np, nn, np-nn);
            end
        case 'harm'
            txt = sprintf('harmonic energy: %.1f%% of |J|^2', 100*Ht.HarmFrac);
        otherwise
            txt = 'total field |J|';
    end
    try, panel_helmholtz('SetReadout', txt); catch, end %#ok<CTCH>
end
