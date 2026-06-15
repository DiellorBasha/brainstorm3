# Helmholtz Components View Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (or subagent-driven-development) to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Turn the Helmholtz view into a component-selected view — Total `|J|`, Irrotational `∇φ`, Solenoidal `∇⊥ψ`, Harmonic `h` — where each state shows that component's own vector field (quiver) and its scalar potential (colormap), with component-aware singular-point markers.

**Architecture:** `bst_dirac_helmholtz` gains a per-face FEM gradient (in `Prepare`) so `Frame` can return the three component vector fields plus `|h|` and the harmonic-energy fraction (the harmonic part is the exact residual `J − ∇φ − ∇⊥ψ`). `view_helmholtz` overrides, per frame, the native source figure's quiver (`QuiverVectorOverride`) with the selected component's field and the cortex colormap with its scalar; markers are component-aware (vortex cores on Solenoidal, sources/sinks on Irrotational). `panel_helmholtz` is a 4-state radio + two overlay toggles + a per-component readout.

**Tech Stack:** MATLAB (Brainstorm), the stored first-order Dirac + cotan LBO (cached Cholesky), `figure_3d>PlotSourceVectors` (`QuiverVectorOverride`), `bst_colormaps` (`source`/`stat2`), the `CustomOverlayFcn` time hook. Tests via the MATLAB MCP `run_matlab_file` against Subject01 `cortex_pial_low` (Dirac + LBO operators exist).

**Repo / branch:** `brainstorm3` on `feat/helmholtz-view`.

**Spec:** `docs/superpowers/specs/2026-06-15-helmholtz-components-design.md`.

---

## File Structure

**Modify:**
- `toolbox/math/bst_dirac_helmholtz.m` — `Prepare` builds `Gx/Gy/Gz`; `Frame` returns component vector fields + `Hmag` + `HarmFrac` + `Sources`; `FindCores` reused for sources/sinks.
- `toolbox/gui/view_helmholtz.m` — component-based uniform override (quiver = component field, colormap = component scalar), component-aware markers/readout.
- `toolbox/gui/panel_helmholtz.m` — 4-state component radio + Show vectors + Show singular points + readout.
- `dev/tests/test_dirac_helmholtz.m` — decomposition-correctness checks.
- `dev/tests/test_helmholtz_view.m` — component-state checks.

---

## Task 1: Decomposition math — gradient operator + component fields

**Files:**
- Modify: `toolbox/math/bst_dirac_helmholtz.m`
- Test: `dev/tests/test_dirac_helmholtz.m`

- [ ] **Step 1: Add the failing decomposition test**

Append these checks inside `test_dirac_helmholtz()` just before the final `fprintf`/`error` lines (after the existing check (4) block):

```matlab
    % --- (5) Hodge decomposition: component fields, exact reconstruction, dominance ---
    Op2 = bst_dirac_helmholtz('Prepare', Dirac, LBO, Surf);
    Ht  = bst_dirac_helmholtz('Frame', Op2, J(:,1));
    nFail = nFail + chk('component fields are [nV x 3]', isequal(size(Ht.Virr),[size(V,1) 3]) && isequal(size(Ht.Vsol),[size(V,1) 3]));
    recon = Ht.Virr + Ht.Vsol + Ht.Vharm;
    nFail = nFail + chk('exact reconstruction Virr+Vsol+Vharm == J', max(abs(recon(:)-Ht.Vtot(:))) < 1e-9*max(abs(Ht.Vtot(:))));
    % re-decompose each component: irrotational is divergence-dominated, solenoidal curl-dominated
    Hirr = bst_dirac_helmholtz('Frame', Op2, reshape(Ht.Virr',[],1));
    Hsol = bst_dirac_helmholtz('Frame', Op2, reshape(Ht.Vsol',[],1));
    nFail = nFail + chk('irrotational is divergence-dominated', sum(Hirr.Div.^2) > sum(Hirr.Curl.^2));
    nFail = nFail + chk('solenoidal is curl-dominated',         sum(Hsol.Curl.^2) > sum(Hsol.Div.^2));
    % markers: vortex cores from psi (sign omega), sources/sinks from phi (sign div)
    nFail = nFail + chk('Cores + Sources are struct arrays', isstruct(Ht.Cores) && isstruct(Ht.Sources));
    nFail = nFail + chk('HarmFrac in [0,1]', isscalar(Ht.HarmFrac) && Ht.HarmFrac >= 0 && Ht.HarmFrac <= 1.0001);
```

- [ ] **Step 2: Run it to verify it fails**

Run (MCP `run_matlab_file`): `dev/tests/test_dirac_helmholtz.m`
Expected: FAIL — `Ht.Virr` (and `Vsol/Vharm/Vtot/Hmag/Sources/HarmFrac`) do not exist yet.

- [ ] **Step 3: Replace `Prepare` and `Frame` in `bst_dirac_helmholtz.m`**

Replace the whole `Prepare` function with this version (adds the per-face FEM gradient `Gx/Gy/Gz`):

```matlab
%% ===== PREPARE: build + cache the static operator pieces and the Cholesky factor =====
function Op = Prepare(DiracOp, LBO, Surf) %#ok<DEFNU>
    Vtx = Surf.Vertices;  Fcs = double(Surf.Faces);
    nVtot = size(Vtx,1);
    nH = numel(DiracOp.FirstOrder.Intrinsic);
    Op = struct();
    Op.nVtot    = nVtot;
    Op.VertConn = Surf.VertConn;
    [Op.D, Op.vH, Op.Nf, Op.Wfv, Op.M, Op.cholK, Op.free, Op.totMass, ...
     Op.Gx, Op.Gy, Op.Gz] = deal(cell(1,nH));
    for hh = 1:nH
        D  = DiracOp.FirstOrder.Intrinsic{hh};           % [4F x 4V]
        vH = double(DiracOp.GlobalVertices{hh}(:));
        nVh = numel(vH);
        % reconstruct the hemisphere's local faces in tess_operators' ordering
        isV = false(nVtot,1); isV(vH) = true;
        fMask = all(isV(Fcs), 2);
        mapV  = zeros(nVtot,1); mapV(vH) = 1:nVh;
        Floc  = mapV(Fcs(fMask, :));   Vloc = Vtx(vH, :);
        % face normals (oriented outward via the surface vertex normals)
        e1 = Vloc(Floc(:,2),:) - Vloc(Floc(:,1),:);
        e2 = Vloc(Floc(:,3),:) - Vloc(Floc(:,1),:);
        Nf = cross(e1, e2, 2);  Af = sqrt(sum(Nf.^2,2));   % Af = 2*area (un-normalized)
        Nf = Nf ./ max(Af, eps);
        vn = Surf.VertNormals(vH(Floc(:,1)), :);
        flip = sum(Nf .* vn, 2) < 0;  Nf(flip,:) = -Nf(flip,:);
        nFh = size(Floc,1);
        % face->vertex area-weighted incidence [nVh x nFh]
        I = [Floc(:,1);Floc(:,2);Floc(:,3)];  Jc = [1:nFh,1:nFh,1:nFh]';
        Wfv = sparse(I, Jc, repmat(Af,3,1), nVh, nFh);
        Wfv = spdiags(1./max(sum(Wfv,2),eps),0,nVh,nVh) * Wfv;
        % per-face FEM gradient of a per-vertex scalar: grad f|_face = sum_i f_i (n x e_i)/(2A),
        % e_i = edge OPPOSITE local vertex i. Three sparse [nFh x nVh] blocks.
        eO1 = Vloc(Floc(:,3),:) - Vloc(Floc(:,2),:);
        eO2 = Vloc(Floc(:,1),:) - Vloc(Floc(:,3),:);
        eO3 = Vloc(Floc(:,2),:) - Vloc(Floc(:,1),:);
        c1 = cross(Nf, eO1, 2) ./ Af;   c2 = cross(Nf, eO2, 2) ./ Af;   c3 = cross(Nf, eO3, 2) ./ Af;
        rows = [(1:nFh)';(1:nFh)';(1:nFh)'];  cols = [Floc(:,1);Floc(:,2);Floc(:,3)];
        Gx = sparse(rows, cols, [c1(:,1);c2(:,1);c3(:,1)], nFh, nVh);
        Gy = sparse(rows, cols, [c1(:,2);c2(:,2);c3(:,2)], nFh, nVh);
        Gz = sparse(rows, cols, [c1(:,3);c2(:,3);c3(:,3)], nFh, nVh);
        % LBO pieces + ONE Cholesky factor of the pinned (vertex 1 fixed) cotan stiffness
        K = LBO.Operator{hh};  M = LBO.Mass{hh};
        free = (2:size(K,1))';
        Op.D{hh}=D; Op.vH{hh}=vH; Op.Nf{hh}=Nf; Op.Wfv{hh}=Wfv; Op.M{hh}=M;
        Op.Gx{hh}=Gx; Op.Gy{hh}=Gy; Op.Gz{hh}=Gz;
        Op.cholK{hh}  = decomposition(K(free,free), 'chol');
        Op.free{hh}   = free;
        Op.totMass{hh}= sum(M(:));
    end
end
```

Replace the whole `Frame` function with this version (computes the component vector fields, `|h|`, `HarmFrac`, and the source/sink markers):

```matlab
%% ===== FRAME: decompose a single time frame using the cached factor =====
function Ht = Frame(Op, Jt) %#ok<DEFNU>
    nVtot = Op.nVtot;
    z1 = zeros(nVtot,1);  z3 = zeros(nVtot,3);
    Ht = struct('Curl',z1,'Div',z1,'Psi',z1,'Phi',z1,'Fmag',z1,'Hmag',z1, ...
                'Vtot',z3,'Virr',z3,'Vsol',z3,'Vharm',z3);
    harmNum = 0;  harmDen = 0;
    for hh = 1:numel(Op.D)
        vH = Op.vH{hh};  nVh = numel(vH);
        Jx = Jt(3*(vH-1)+1);  Jy = Jt(3*(vH-1)+2);  Jz = Jt(3*(vH-1)+3);
        psiQ = zeros(4*nVh,1);
        psiQ(2:4:end) = Jx; psiQ(3:4:end) = Jy; psiQ(4:4:end) = Jz;
        q = Op.D{hh} * psiQ;                               % [4F x 1]
        divF  = q(1:4:end);
        curlF = [q(2:4:end), q(3:4:end), q(4:4:end)];
        omF = sum(curlF .* Op.Nf{hh}, 2);
        omV = Op.Wfv{hh} * omF;   dvV = Op.Wfv{hh} * divF;
        psi = i_poisson(Op.cholK{hh}, Op.M{hh}, omV, Op.free{hh}, Op.totMass{hh});
        phi = i_poisson(Op.cholK{hh}, Op.M{hh}, dvV, Op.free{hh}, Op.totMass{hh});
        % component vector fields: grad(phi) (irrotational), skew-grad(psi) (solenoidal)
        gphi = [Op.Gx{hh}*phi, Op.Gy{hh}*phi, Op.Gz{hh}*phi];   % [nF x 3]
        gpsi = [Op.Gx{hh}*psi, Op.Gy{hh}*psi, Op.Gz{hh}*psi];
        skew = cross(Op.Nf{hh}, gpsi, 2);                       % n x grad(psi)
        Virr = Op.Wfv{hh} * gphi;                               % face->vertex [nVh x 3]
        Vsol = Op.Wfv{hh} * skew;
        Jv   = [Jx Jy Jz];
        Vharm = Jv - Virr - Vsol;                              % exact residual
        % assemble global
        Ht.Curl(vH)=omV; Ht.Div(vH)=dvV; Ht.Psi(vH)=psi; Ht.Phi(vH)=phi;
        Ht.Fmag(vH)=sqrt(Jx.^2+Jy.^2+Jz.^2);  Ht.Hmag(vH)=sqrt(sum(Vharm.^2,2));
        Ht.Vtot(vH,:)=Jv;  Ht.Virr(vH,:)=Virr;  Ht.Vsol(vH,:)=Vsol;  Ht.Vharm(vH,:)=Vharm;
        av = full(sum(Op.M{hh},2));                            % lumped vertex mass
        harmNum = harmNum + sum(av .* sum(Vharm.^2,2));
        harmDen = harmDen + sum(av .* sum(Jv.^2,2));
    end
    Ht.HarmFrac = harmNum / max(harmDen, eps);
    Ht.Cores    = FindCores(Ht.Psi, Op.VertConn, Ht.Curl);    % vortex cores (sign = vorticity)
    Ht.Sources  = FindCores(Ht.Phi, Op.VertConn, Ht.Div);     % sources/sinks (sign = divergence)
end
```

> `Decompose`, `PoissonSolve`, `i_poisson`, and `FindCores` are unchanged. `Decompose` still loops `Frame` and reads `Ht.Curl/Div/Psi/Phi/Fmag` + `Ht.Cores` (all still present), so the whole-series path and its test keep working.

- [ ] **Step 4: Run the test to verify it passes**

Run: `dev/tests/test_dirac_helmholtz.m`
Expected: PASS — all checks incl. (5). If `irrotational is divergence-dominated` / `solenoidal is curl-dominated` fail, the gradient sign is wrong; the reconstruction check must hold to ~1e-9 regardless (it is residual-defined).

- [ ] **Step 5: Commit**

```bash
git add toolbox/math/bst_dirac_helmholtz.m dev/tests/test_dirac_helmholtz.m
git commit -m "feat(helmholtz): gradient operator + component vector fields (irrot/solen/harm)"
```

---

## Task 2: view_helmholtz — component-based override + component-aware markers

**Files:**
- Modify: `toolbox/gui/view_helmholtz.m`

- [ ] **Step 1: Replace `view_helmholtz.m` with the component version**

Replace the whole file with:

```matlab
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
        try, [iDS, iResult] = bst_memory('LoadResultsFileFull', SrcResultsFile); catch, iResult = []; end
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
            np = sum([comp.Markers.charge] > 0); nn = sum([comp.Markers.charge] < 0);
            txt = sprintf('%d vortices (+), %d antivortices (-), net %+d', np, nn, np-nn);
        case 'source'
            np = sum([comp.Markers.charge] > 0); nn = sum([comp.Markers.charge] < 0);
            txt = sprintf('%d sources (+), %d sinks (-), net %+d', np, nn, np-nn);
        case 'harm'
            txt = sprintf('harmonic energy: %.1f%% of |J|^2', 100*Ht.HarmFrac);
        otherwise
            txt = 'total field |J|';
    end
    if isempty(comp.Markers) && strcmp(comp.Kind,'source'); txt = '0 sources, 0 sinks'; end
    if isempty(comp.Markers) && strcmp(comp.Kind,'vortex'); txt = '0 vortices, 0 antivortices'; end
    try, panel_helmholtz('SetReadout', txt); catch, end %#ok<CTCH>
end
```

- [ ] **Step 2: Parse check**

Run (MCP `evaluate_matlab_code`): `rehash; exist('view_helmholtz','file')` → expect `2`.

- [ ] **Step 3: Commit**

```bash
git add toolbox/gui/view_helmholtz.m
git commit -m "feat(helmholtz): component-based view (quiver=component field, colormap=potential)"
```

---

## Task 3: panel_helmholtz — 4-state radio + toggles + live component test

**Files:**
- Modify: `toolbox/gui/panel_helmholtz.m`
- Test: `dev/tests/test_helmholtz_view.m`

- [ ] **Step 1: Replace the panel body**

Replace `CreatePanel` and the `OnScalar`/`OnCores` callbacks in `panel_helmholtz.m` with:

```matlab
function bstPanelNew = CreatePanel(hFig) %#ok<DEFNU>
    import javax.swing.*;
    panelName = 'Helmholtz';
    jPanelNew = gui_component('Panel');
    jOpt = JPanel(); jOpt.setLayout(BoxLayout(jOpt, BoxLayout.Y_AXIS));
    jSec = gui_river([2 2], [2 8 3 6], 'Helmholtz / Hodge components');

    gui_component('label', jSec, 'br', 'Component:');
    names  = {'Total','Irrot','Solen','Harm'};
    labels = {'Total field |J|','Irrotational (grad phi)','Solenoidal (curl psi)','Harmonic (h)'};
    grp = ButtonGroup(); jRadio = javaArray('javax.swing.JRadioButton', numel(names));
    for i = 1:numel(names)
        jRadio(i) = gui_component('radio', jSec, 'br', labels{i});
        grp.add(jRadio(i));
        java_setcb(jRadio(i), 'ActionPerformedCallback', @(h,e) OnComponent(panelName, names{i}));
    end
    jRadio(1).setSelected(true);   % Total: native start
    jVec  = gui_component('checkbox', jSec, 'br', 'Show vectors');           jVec.setSelected(true);
    jMark = gui_component('checkbox', jSec, 'br', 'Show singular points');   jMark.setSelected(true);
    java_setcb(jVec,  'ActionPerformedCallback', @(h,e) OnVectors(panelName));
    java_setcb(jMark, 'ActionPerformedCallback', @(h,e) OnMarkers(panelName));
    jReadout = gui_component('label', jSec, 'br', '');
    jClose = gui_component('button', jSec, 'br', 'Close');
    java_setcb(jClose, 'ActionPerformedCallback', @(h,e) OnClose(panelName));

    jOpt.add(jSec); jPanelNew.add(jOpt, java.awt.BorderLayout.NORTH);
    ctrl = struct('hFig',hFig, 'jVec',jVec, 'jMark',jMark, 'jReadout',jReadout);
    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);
end

function OnComponent(panelName, name) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    view_helmholtz('SetComponent', ctrl.hFig, name);
end
function OnVectors(panelName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    view_helmholtz('SetVectors', ctrl.hFig, ctrl.jVec.isSelected());
end
function OnMarkers(panelName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    view_helmholtz('SetMarkers', ctrl.hFig, ctrl.jMark.isSelected());
end
```

Keep `SetReadout`, `OnClose`, and `i_valid` as they are (they already operate on `ctrl.hFig` / `ctrl.jReadout` / the `Helmholtz` panel). Update the file header comment to describe components.

- [ ] **Step 2: Replace the live test with the component version**

Replace the body of `dev/tests/test_helmholtz_view.m` with:

```matlab
function test_helmholtz_view()
% Live: open the Helmholtz components view on a synthetic 3-comp source. Default Total shows
% the field + |J|; each component sets the quiver override to ITS vector field and colors by
% ITS scalar (signed for Irrot/Solen); markers are component-aware; cursor move recomputes;
% close cleans up.
% Authors: Diellor Basha, 2026
    global GlobalData; %#ok<NUSED>
    nFail = 0;
    SurfaceFile = bst_get('Subject',1).Surface(5).FileName;
    nV = size(in_tess_bst(SurfaceFile,0).Vertices,1);
    R = db_template('resultsmat'); rng(0);
    R.ImageGridAmp = randn(3*nV,3)*1e-9; R.nComponents=3; R.Time=0:2;
    R.HeadModelType='surface'; R.SurfaceFile=file_short(SurfaceFile); R.Comment='SYN helmholtz src';
    srcFile = db_add(-3, R);

    hFig = view_helmholtz(srcFile); drawnow;
    nFail = nFail + chk('view opens', ishandle(hFig));
    St = getappdata(hFig,'HelmholtzState'); hAx = findobj(hFig,'-depth',1,'Tag','Axes3D'); hAx=hAx(1);
    nFail = nFail + chk('default component is Total', strcmp(St.Component,'Total'));
    nFail = nFail + chk('native source vectors shown', ~isempty(findobj(hFig,'Tag','SourceVectors')));

    [~, iT] = bst_memory('GetTimeVector', St.srcDS, St.srcResult, 'CurrentTimeIndex');
    Jt = double(bst_memory('GetResultsValues', St.srcDS, St.srcResult, [], iT, 0));
    Ht = bst_dirac_helmholtz('Frame', St.Op, Jt);

    % Total: quiver override == J, source colormap
    nFail = nFail + chk('Total quiver = J', isequal(getappdata(hFig,'QuiverVectorOverride'), Ht.Vtot));
    nFail = nFail + chk('Total colormap source', strcmpi(getappdata(hFig,'Colormap').Type,'source'));

    % Solenoidal: quiver = Vsol, colormap stat2, markers = vortex cores
    view_helmholtz('SetComponent', hFig, 'Solen'); drawnow;
    nFail = nFail + chk('Solen quiver = Vsol', isequal(getappdata(hFig,'QuiverVectorOverride'), Ht.Vsol));
    nFail = nFail + chk('Solen colormap stat2', strcmpi(getappdata(hFig,'Colormap').Type,'stat2'));
    nFail = nFail + chk('Solen markers = vortex cores', numel(findobj(hAx,'Tag','HelmholtzCore'))==numel(Ht.Cores));

    % Irrotational: quiver = Virr, markers = sources/sinks
    view_helmholtz('SetComponent', hFig, 'Irrot'); drawnow;
    nFail = nFail + chk('Irrot quiver = Virr', isequal(getappdata(hFig,'QuiverVectorOverride'), Ht.Virr));
    nFail = nFail + chk('Irrot markers = sources/sinks', numel(findobj(hAx,'Tag','HelmholtzCore'))==numel(Ht.Sources));

    % Harmonic: quiver = Vharm, no markers, source colormap
    view_helmholtz('SetComponent', hFig, 'Harm'); drawnow;
    nFail = nFail + chk('Harm quiver = Vharm', isequal(getappdata(hFig,'QuiverVectorOverride'), Ht.Vharm));
    nFail = nFail + chk('Harm has no markers', isempty(findobj(hAx,'Tag','HelmholtzCore')));

    % toggles
    view_helmholtz('SetVectors', hFig, false); drawnow;
    nFail = nFail + chk('vectors hide when off', isempty(findobj(hFig,'Tag','SourceVectors')));
    view_helmholtz('SetComponent', hFig, 'Solen'); view_helmholtz('SetMarkers', hFig, false); drawnow;
    nFail = nFail + chk('markers hide when off', isempty(findobj(hAx,'Tag','HelmholtzCore')));

    % time-following recompute
    nFail = nFail + chk('figure not static', ~isequal(getappdata(hFig,'isStatic'),1));
    panel_time('SetCurrentTime', 2.0); drawnow;
    [~, iT2] = bst_memory('GetTimeVector', St.srcDS, St.srcResult, 'CurrentTimeIndex');
    nFail = nFail + chk('cursor move decomposed a new frame', isKey(St.Cache, iT2) && (St.Cache.Count >= 2));

    % close + stale guard
    view_helmholtz('Close', hFig); drawnow;
    nFail = nFail + chk('panel closed', isempty(bst_get('PanelControls','Helmholtz')));
    nFail = nFail + chk('figure closed', ~ishandle(hFig));
    okGuard = true;
    try, view_helmholtz('SetComponent', hFig, 'Irrot'); catch; okGuard = false; end
    nFail = nFail + chk('stale dispatch ignored', okGuard);

    [~,iSt] = bst_get('AnyFile', srcFile); file_delete(file_fullpath(srcFile),1); db_reload_studies(iSt);
    fprintf('\n==== test_helmholtz_view: %d failed ====\n', nFail);
    if nFail > 0, error('test_helmholtz_view FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
```

- [ ] **Step 3: Run the live test**

Run: `dev/tests/test_helmholtz_view.m`
Expected: PASS — all checks. (If `Solen quiver = Vsol` fails because `QuiverVectorOverride` is single vs double, wrap both sides in `double(...)` — `Ht.Vsol` is double, the appdata should match.)

- [ ] **Step 4: Run the math test too (regression)**

Run: `dev/tests/test_dirac_helmholtz.m`  → Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_helmholtz.m dev/tests/test_helmholtz_view.m
git commit -m "feat(helmholtz): 4-component panel (Total/Irrot/Solen/Harm) + component-aware markers"
```

---

## Self-Review

**Spec coverage:**
- Three component vector fields `∇φ/∇⊥ψ/h` + Total → Task 1 (`Frame` returns `Virr/Vsol/Vharm/Vtot`). ✓
- Each state = component quiver + its potential colormap → Task 2 (`i_component`, `QuiverVectorOverride`, signed→stat2 / unsigned→source). ✓
- Component-aware markers (vortex cores on Solen, sources/sinks on Irrot, none on Total/Harm) → Task 1 (`Cores`,`Sources`) + Task 2 (`comp.Markers`, `i_readout`). ✓
- Exact-residual harmonic + energy readout → Task 1 (`Vharm = J − Virr − Vsol`, `HarmFrac`) + Task 2 (`i_readout` harm). ✓
- 4-state panel + two toggles + readout → Task 3. ✓
- Active-frame-only + cached factor preserved → Task 2 (cache + `Frame`). ✓

**Placeholder scan:** none — every step has complete code.

**Type consistency:** `Frame` returns `Vtot/Virr/Vsol/Vharm [nV×3]`, `Fmag/Phi/Psi/Hmag [nV×1]`, `Cores/Sources` struct arrays `(iVertex,charge)`, `HarmFrac` scalar — consumed identically by `i_component` (Task 2) and the test (Task 3). Component names `{'Total','Irrot','Solen','Harm'}` match between the panel radio (Task 3), `SetComponent`/`i_component` (Task 2). Dispatch names `{'SetComponent','SetVectors','SetMarkers','Close','UpdateFrame'}` match the shim and the panel callbacks. ✓

---

## Build order

1. Task 1 — gradient operator + component fields (math + correctness test).
2. Task 2 — component-based view override + markers/readout.
3. Task 3 — 4-state panel + live component test.
