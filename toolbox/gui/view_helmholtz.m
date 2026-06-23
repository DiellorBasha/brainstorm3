function hFig = view_helmholtz(SrcResultsFile, varargin)
% VIEW_HELMHOLTZ: Helmholtz/Hodge component view of a Dirac source map. Opens the NATIVE
% unconstrained-source display and lets a panel choose which COMPONENT of the decomposition
% to show -- Total |J| / Irrotational grad(phi) / Solenoidal curl(psi) / Harmonic h.
% Switching component swaps the cortex SCALAR colormap (to that component's potential) and
% the component-aware singular-point markers; the quiver always shows the (smoothed) TOTAL
% source field. An optional Dirac-eigenmode smoothing low-passes the active frame before the
% decomposition, and a magnitude gate prunes weak markers. Active frame only (cached Cholesky
% factor; recomputed as the cursor moves).
%
% USAGE:
%   hFig = view_helmholtz(SrcResultsFile)
%   view_helmholtz('SetComponent', hFig, name)             % 'Total'|'Irrot'|'Solen'|'Harm'
%   view_helmholtz('SetVectors', hFig, show)
%   view_helmholtz('SetMarkers', hFig, show)
%   view_helmholtz('SetSmoothing', hFig, isOn, name, params)  % Dirac-eigenmode low-pass
%   view_helmholtz('SetGate', hFig, frac)                  % marker magnitude gate (0..1)
%   view_helmholtz('Close', hFig)
%   view_helmholtz('UpdateFrame', hFig)
% Authors: Diellor Basha, 2026
    global GlobalData;
    hFig = [];
    if (nargin >= 1) && ischar(SrcResultsFile) && any(strcmp(SrcResultsFile, {'SetComponent','SetVectors','SetMarkers','SetSmoothing','SetGate','SetTrack','SetDeriv','Close','UpdateFrame'}))
        if any(strcmp(SrcResultsFile, {'SetComponent','SetVectors','SetMarkers','SetSmoothing','SetGate','SetTrack','SetDeriv','UpdateFrame'})) && ...
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
    bst_progress('text', 'Loading Dirac eigenbasis...');
    EigenMat = tess_eigen(SurfaceFile, 'Dirac');
    OpMat    = in_bst_operator(EigenMat.OperatorFile);
    Lambda   = double(EigenMat.Lambda{1}(:));
    bst_progress('stop');

    [hFig, iDSf] = view_surface_data(SurfaceFile, SrcResultsFile, [], 'NewFigure');
    if isempty(hFig); return; end
    iTess = i_find_tess(hFig);
    % No amplitude threshold by default: the Data threshold slider drives the quiver +
    % colormap (the override quiver is no longer hardcoded). Show the full field on open.
    TI = getappdata(hFig,'Surface');
    if iTess <= numel(TI); TI(iTess).DataThreshold = 0; setappdata(hFig,'Surface',TI); end
    try, panel_surface('UpdateSurfaceProperties'); catch, end %#ok<CTCH>

    St = struct('Op',Op, 'srcDS',iDSf, 'srcResult',iResult, 'Component','Total', ...
                'ShowVectors',true, 'ShowMarkers',true, 'iTess',iTess, 'nV',nV, ...
                'EigenMat',EigenMat, 'Mass',{OpMat.Mass}, 'Lambda',Lambda, ...
                'Smooth',struct('on',false,'name','heat','params',struct()), 'GateFrac',0, 'Deriv',0, ...
                'Cache',containers.Map('KeyType','double','ValueType','any'), ...
                'Track',false, 'LastIT',[], 'Tracks',[], 'Graph',[], 'Ctx',{cell(1,numel(Op.vH))}, ...
                'V2H',[], 'Vertices',Surf.Vertices, 'Faces',double(Surf.Faces), ...
                'TrackComp','', 'TrackSmooth',false);
    setappdata(hFig, 'HelmholtzState', St);
    setappdata(hFig, 'CustomOverlayFcn', @(h) UpdateFrame(h));
    set(hFig, 'CloseRequestFcn', @(h,e) Close(h));
    UpdateFrame(hFig);
    i_attach_listener(hFig);   % re-gate overlays on a hemisphere/resect toggle (Surfaces panel)

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
    hAx = hAx(1);
    TessInfo = getappdata(hFig,'Surface');
    if isempty(TessInfo) || (St.iTess > numel(TessInfo)) || ~ishandle(TessInfo(St.iTess).hPatch); return; end
    % Per-hemisphere visibility from the cortex patch alpha (Surfaces-panel resect / hemi
    % toggle). Every overlay below is gated by this, so a hidden hemisphere also hides its
    % quiver, markers, tracks and centroids (mirrors view_manifold's MarkedClean sync).
    visV = i_visibleVerts(TessInfo(St.iTess).hPatch);
    setappdata(hFig, 'HelmholtzLastVis', visV);
    [TimeVec, iT] = bst_memory('GetTimeVector', St.srcDS, St.srcResult, 'CurrentTimeIndex');
    if isempty(iT) || iT < 1; iT = 1; end
    k = St.Deriv;
    if iT <= k                                            % not enough history for this order
        i_blank_display(hFig, St, hAx, k);
        St.Tracks = []; St.LastIT = iT; setappdata(hFig, 'HelmholtzState', St);
        return;
    end
    % fetch the (k+1) consecutive frames, oldest -> newest (last col = current)
    F = zeros(3*St.nV, k+1);
    for j = 0:k
        fj = double(bst_memory('GetResultsValues', St.srcDS, St.srcResult, [], iT-(k-j), 0));
        if size(fj,1) ~= 3*St.nV; return; end
        F(:,j+1) = fj;
    end
    if k > 0, dt = TimeVec(iT) - TimeVec(iT-1); else, dt = 1; end
    Jt = bst_time_derivative(F, dt, k);                   % D^k J(t): Field / Velocity / Acceleration
    % spatial eigenmode smoothing (linear -> applying after the difference is equivalent)
    if St.Smooth.on
        g  = bst_eigfilter_kernel(St.Smooth.name, St.Smooth.params);
        % bst_eigenfilter reads the per-hemisphere mass Oper.Mass{h}; St.Mass is {B1,B2}.
        Jt = real(bst_eigenfilter('Analysis', Jt, St.EigenMat, struct('Mass', {St.Mass}), g, struct()));
    end
    % Detect singular points only when they are actually used (markers shown or tracking on),
    % so plain time-stepping with both off skips the per-frame persistence computation.
    needCores = St.ShowMarkers || St.Track;
    if isKey(St.Cache, iT)
        Ht = St.Cache(iT);
        if needCores && isempty(Ht.Cores)                 % cached without detection -> fill in now
            Ht = bst_dirac_helmholtz('Frame', St.Op, Jt, true);  St.Cache(iT) = Ht;
        end
    else
        Ht = bst_dirac_helmholtz('Frame', St.Op, Jt, needCores);  St.Cache(iT) = Ht;
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
    TessInfo = getappdata(hFig,'Surface');
    % --- quiver is ALWAYS the (smoothed) total source field, regardless of component;
    %     switching component changes only the scalar colormap + markers ---
    if St.ShowVectors
        Vq = Ht.Vtot;                                    % gate the quiver per hemisphere:
        if numel(visV) == size(Vq,1); Vq(~visV,:) = 0; end   % zero-length arrows = hidden
        setappdata(hFig, 'QuiverVectorOverride', Vq);
        try, figure_3d('SetShowSourceVectors', hFig, St.iTess, 1); catch, end %#ok<CTCH>
        try, figure_3d('PlotSourceVectors', hFig, St.iTess); catch, end %#ok<CTCH>
    else
        try, figure_3d('SetShowSourceVectors', hFig, St.iTess, 0); catch, end %#ok<CTCH>
    end
    % --- component markers, pruned PER HEMISPHERE by the persistence gate;
    %     each hemisphere's global core is always kept ---
    mk = comp.Markers;
    if ~isempty(mk) && St.GateFrac > 0
        mk = bst_persistence_gate(mk, St.GateFrac);
    end
    delete(findobj(hAx,'Tag','HelmholtzCore'));
    if St.ShowMarkers && ~isempty(mk)
        V = get(TessInfo(St.iTess).hPatch, 'Vertices');
        prA = [mk.persistence];
        mxf = max([prA(isfinite(prA)), eps]);              % frame max finite persistence
        for k = 1:numel(mk)
            v = mk(k).iVertex;
            if (v <= numel(visV)) && ~visV(v); continue; end   % skip markers on a hidden hemisphere
            col = [1 0 0]; if mk(k).charge < 0; col = [0 0 1]; end
            pk = mk(k).persistence; if ~isfinite(pk); pk = mxf; end   % globals -> max size
            sz = 5 + 11 * max(0, min(1, pk/mxf));          % marker size 5..16 by persistence
            line('Parent',hAx,'XData',V(v,1),'YData',V(v,2),'ZData',V(v,3), 'Marker','o', ...
                'MarkerSize',sz,'MarkerFaceColor',col,'MarkerEdgeColor','k','LineStyle','none', ...
                'Tag','HelmholtzCore','Clipping','off');
        end
    end
    ordName = {'', 'velocity ', 'acceleration '};  pfx = ordName{St.Deriv+1};
    if needCores
        i_readout(comp.Kind, mk, Ht, pfx);
    else
        try, panel_helmholtz('SetReadout', [pfx 'singular points hidden']); catch, end %#ok<CTCH>
    end
    if St.Track
        St = i_track_update(hFig, St, hAx, mk, iT, visV);
        setappdata(hFig, 'HelmholtzState', St);
    end
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
function SetDeriv(hFig, order) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.Deriv  = max(0, min(2, round(order)));
    St.Cache  = containers.Map('KeyType','double','ValueType','any');  % decompositions now stale
    St.Tracks = [];  St.LastIT = [];                                    % derivative field changed
    setappdata(hFig, 'HelmholtzState', St);
    UpdateFrame(hFig);
end
function SetTrack(hFig, isOn) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.Track = logical(isOn);
    if ~St.Track
        hAx = findobj(hFig,'-depth',1,'Tag','Axes3D');
        if ~isempty(hAx)
            delete(findobj(hAx(1),'Tag','HelmholtzTrack'));
            delete(findobj(hAx(1),'Tag','HelmholtzCentroid'));
        end
        St.Tracks = []; St.LastIT = [];
    end
    setappdata(hFig, 'HelmholtzState', St);
    UpdateFrame(hFig);
end
function SetSmoothing(hFig, isOn, name, params) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.Smooth = struct('on',logical(isOn), 'name',name, 'params',params);
    St.Cache  = containers.Map('KeyType','double','ValueType','any');   % decompositions now stale
    setappdata(hFig, 'HelmholtzState', St);  UpdateFrame(hFig);
end
function SetGate(hFig, frac) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.GateFrac = max(0, min(1, frac));  setappdata(hFig, 'HelmholtzState', St);  UpdateFrame(hFig);
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
% view_manifold/SyncScalar reads. Brainstorm hides a hemisphere by switching EdgeAlpha to
% per-element ('flat'/'interp') with zeroed FaceVertexAlphaData; a vertex is visible if any
% incident face is visible. Returns all-visible when nothing is resected.
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

%% ===== re-gate overlays when the cortex patch is redrawn (hemisphere/resect toggle) =====
function i_attach_listener(hFig)
% Mirror view_manifold: a MarkedClean listener re-gates the overlays whenever the cortex
% patch is redrawn, so a Surfaces-panel hemisphere toggle hides that hemisphere's quiver /
% markers / tracks / centroids immediately, without waiting for the time cursor to move.
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
            S = in_bst_operator(sSubject.Surface(iSurf).Operator(k).FileName);
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
function i_readout(kind, mk, Ht, pfx)
    if nargin < 4, pfx = ''; end
    switch kind
        case 'vortex'
            txt = i_count_str(mk, 'vortices', 'antivortices', Ht);
        case 'source'
            txt = i_count_str(mk, 'sources', 'sinks', []);
        case 'harm'
            txt = sprintf('harmonic energy: %.1f%% of |J|^2', 100*Ht.HarmFrac);
        otherwise
            txt = 'total field |J|';
    end
    try, panel_helmholtz('SetReadout', [pfx txt]); catch, end %#ok<CTCH>
end

function i_blank_display(hFig, St, hAx, k)
    TessInfo = getappdata(hFig,'Surface');
    TessInfo(St.iTess).Data         = zeros(St.nV,1);
    TessInfo(St.iTess).DataMinMax   = [0 1];
    TessInfo(St.iTess).ColormapType = 'source';
    setappdata(hFig,'Surface',TessInfo);
    panel_surface('UpdateSurfaceColormap', hFig);
    try, figure_3d('SetShowSourceVectors', hFig, St.iTess, 0); catch, end %#ok<CTCH>
    delete(findobj(hAx,'Tag','HelmholtzCore'));
    delete(findobj(hAx,'Tag','HelmholtzTrack'));
    delete(findobj(hAx,'Tag','HelmholtzCentroid'));
    nm = 'velocity'; if k >= 2, nm = 'acceleration'; end
    try, panel_helmholtz('SetReadout', sprintf('%s: needs %d earlier frame(s)', nm, k)); catch, end %#ok<CTCH>
end

function txt = i_count_str(mk, posName, negName, Ht)
    if isempty(mk)
        txt = sprintf('0 %s, 0 %s', posName, negName);
    else
        hh = [mk.hemi];  parts = {};
        for h = unique(hh)
            m = mk(hh==h);
            np = sum([m.charge] > 0);  nn = sum([m.charge] < 0);
            parts{end+1} = sprintf('H%d: %d(+), %d(-)', h, np, nn); %#ok<AGROW>
        end
        txt = strjoin(parts, '  |  ');
        if ~isempty(Ht)
            pr = [mk.persistence];  tp = max(pr(isfinite(pr)));  if isempty(tp), tp = 0; end
            txt = sprintf('%s  (top persistence %.2g)', txt, tp);
        end
    end
end

%% ===== trajectory tracking (accumulated over contiguous forward play) =====
function St = i_track_update(hFig, St, hAx, mk, iT, visV) %#ok<INUSL>
    MAXJUMP = 0.012;                                  % m, max core jump per frame
    V = St.Vertices;
    if isempty(St.Graph),  St.Graph = i_build_graph(V, St.Faces);  end
    if isempty(St.V2H)
        St.V2H = zeros(size(V,1),1);
        for h = 1:numel(St.Op.vH), St.V2H(St.Op.vH{h}) = h; end
    end
    % reset on toggle-on, non-contiguous time, backward step, or component/smoothing change
    reset = isempty(St.LastIT) || (iT ~= St.LastIT + 1) ...
            || ~strcmp(St.TrackComp, St.Component) || (St.TrackSmooth ~= St.Smooth.on);
    if reset
        St.Tracks = i_seed_tracks(mk, St.V2H);
    else
        heads = i_heads(St.Tracks, V);
        cur   = i_cores_struct(mk, St.V2H, V);
        match = bst_vortex_link_step(heads.s, cur, MAXJUMP);
        usedC = false(1, numel(cur));
        for hi = 1:numel(match)
            ti = heads.idx(hi);
            if match(hi) > 0
                c = match(hi);  usedC(c) = true;
                p = shortestpath(St.Graph, St.Tracks(ti).coreVerts(end), cur(c).iVertex);
                if numel(p) >= 2, St.Tracks(ti).path = [St.Tracks(ti).path, p(2:end)]; end
                St.Tracks(ti).coreVerts(end+1) = cur(c).iVertex;
                St.Tracks(ti).persist = cur(c).persistence;
            else
                St.Tracks(ti).open = false;
            end
        end
        born = i_seed_tracks(mk(~usedC), St.V2H);
        St.Tracks = [St.Tracks, born];
    end
    % draw accumulated polylines
    delete(findobj(hAx,'Tag','HelmholtzTrack'));
    for t = St.Tracks
        if numel(t.path) < 2, continue; end
        col = [1 0 0]; if t.chirality < 0, col = [0 0 1]; end
        xs = V(t.path,1); ys = V(t.path,2); zs = V(t.path,3);
        if numel(visV) >= max(t.path)                          % break the polyline at any
            hideP = ~visV(t.path); xs(hideP) = NaN; ys(hideP) = NaN; zs(hideP) = NaN;  % hidden-hemi vertex
        end
        line('Parent',hAx,'XData',xs,'YData',ys,'ZData',zs, ...
             'Color',col,'LineWidth',2,'Tag','HelmholtzTrack','Clipping','off');
    end
    % --- Karcher-mean centroids for the top-N longest-lived, strongest tracks ---
    delete(findobj(hAx,'Tag','HelmholtzCentroid'));
    openT = find([St.Tracks.open] & arrayfun(@(t)numel(t.coreVerts), St.Tracks) >= 3);
    if ~isempty(openT)
        [~, ord] = sort(arrayfun(@(i)St.Tracks(i).persist, openT), 'descend');
        openT = openT(ord(1:min(5, numel(ord))));            % cap to top-5
        for i = openT
            hh = St.Tracks(i).hemi;
            if (numel(visV) >= max(St.Op.vH{hh})) && ~any(visV(St.Op.vH{hh})); continue; end  % hidden hemi
            [ctx, loc] = i_hemi_context(St, hh);
            lv = loc(St.Tracks(i).coreVerts);  lv = lv(lv > 0);
            if numel(lv) < 3, continue; end
            try
                c = nxr.manifold.query.center(ctx, lv(:));
            catch
                continue;
            end
            if isscalar(c), p = St.Vertices(St.Op.vH{hh}(c), :); else, p = c(:)'; end
            col = [1 0 0]; if St.Tracks(i).chirality < 0, col = [0 0 1]; end
            line('Parent',hAx,'XData',p(1),'YData',p(2),'ZData',p(3),'Marker','d', ...
                 'MarkerSize',12,'MarkerFaceColor',col,'MarkerEdgeColor','k','LineStyle','none', ...
                 'Tag','HelmholtzCentroid','Clipping','off');
        end
    end
    St.LastIT = iT;  St.TrackComp = St.Component;  St.TrackSmooth = St.Smooth.on;
end

function [ctx, loc] = i_hemi_context(St, hh)
% Cache one nxr context per hemisphere (in root appdata, keyed by hemisphere -- survives
% the per-call St; rebuilt once per session) + a global->local vertex map.
    key = sprintf('HelmCtx%d', hh);
    ctx = getappdata(0, key);
    if isempty(ctx)
        vH = St.Op.vH{hh};  nVh = numel(vH);
        isV = false(size(St.Vertices,1),1); isV(vH) = true;
        fMask = all(isV(St.Faces), 2);
        mapV = zeros(size(St.Vertices,1),1); mapV(vH) = 1:nVh;
        Floc = mapV(St.Faces(fMask,:));  Vloc = St.Vertices(vH,:);
        ctx = nxr.manifold.context(Vloc, int32(Floc));
        setappdata(0, key, ctx);
    end
    loc = zeros(size(St.Vertices,1),1);
    loc(St.Op.vH{hh}) = 1:numel(St.Op.vH{hh});
end

function G = i_build_graph(V, F)
    E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
    E = unique(sort(E,2), 'rows');
    w = sqrt(sum((V(E(:,1),:) - V(E(:,2),:)).^2, 2));
    G = graph(E(:,1), E(:,2), w, size(V,1));
end

function T = i_seed_tracks(mk, V2H)
    T = i_empty_track();
    for k = 1:numel(mk)
        v = mk(k).iVertex;
        T(end+1) = struct('coreVerts',v, 'path',v, 'chirality',sign(mk(k).charge), ...
                          'hemi',V2H(v), 'open',true, 'persist',mk(k).persistence, ...
                          'centroid',[]); %#ok<AGROW>
    end
end
function t = i_empty_track()
    t = struct('coreVerts',{},'path',{},'chirality',{},'hemi',{},'open',{},'persist',{},'centroid',{});
end
function h = i_heads(Tracks, V)
    h.idx = find([Tracks.open]);
    s = struct('pos',{},'charge',{},'hemi',{});
    for i = h.idx
        v = Tracks(i).coreVerts(end);
        s(end+1) = struct('pos',V(v,:), 'charge',Tracks(i).chirality, 'hemi',Tracks(i).hemi); %#ok<AGROW>
    end
    h.s = s;
end
function cur = i_cores_struct(mk, V2H, V)
    cur = struct('pos',{},'charge',{},'hemi',{},'iVertex',{},'persistence',{});
    for k = 1:numel(mk)
        v = mk(k).iVertex;
        cur(end+1) = struct('pos',V(v,:), 'charge',mk(k).charge, 'hemi',V2H(v), ...
                            'iVertex',v, 'persistence',mk(k).persistence); %#ok<AGROW>
    end
end
