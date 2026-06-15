# Helmholtz / Vorticity View Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A Helmholtz/Vorticity view launched from a Dirac source figure: a dedicated cortex figure with a switchable signed-scalar colormap (curl·n̂, divergence, stream function ψ, potential φ, |field|), a field quiver, and vortex-core markers, all following the time cursor.

**Architecture:** A pure `bst_dirac_helmholtz` helper turns the source field into div/curl (via the stored first-order Dirac), Poisson-solves the scalar LBO for ψ/φ, and detects cores as ψ extrema. `view_helmholtz` shows the chosen scalar as a native 1-component results node (native colormap + time) on a dedicated figure, and rides Brainstorm's per-figure time refresh — via a small `CustomOverlayFcn` hook in `panel_surface` — to redraw the quiver + cores. `panel_helmholtz` is the control panel.

**Tech Stack:** MATLAB (Brainstorm toolbox), the stored first-order Dirac (`OperatorMat.FirstOrder.Intrinsic`) + scalar LBO (`tess_operators 'Laplace-Beltrami'`), Java-Swing GUI, `bst_memory`/`panel_surface`/`view_surface_data`, low-level `line` glyph overlays. Tests via the MATLAB MCP against the live dev protocol (Subject01 `cortex_20484V` has Dirac + can make an LBO operator).

**Repo / branch:** `brainstorm3` on `feat/helmholtz-view`.

**Conventions:** Tests in `dev/tests/`, plain functions that `error()` on failure, run with the MATLAB MCP `run_matlab_file`. Panel/dispatch functions use `eval(macro_method)`. Commit after each task with the message in its final step.

---

## File Structure

**Create:**
- `toolbox/math/bst_dirac_helmholtz.m` — pure: field → div/curl/ψ/φ + cores (dispatched subfns `PoissonSolve`, `FindCores`, and the main entry).
- `toolbox/gui/view_helmholtz.m` — launcher + dedicated figure + overlays + register the time hook.
- `toolbox/gui/panel_helmholtz.m` — control panel (scalar radio, toggles, readout, Close).
- `dev/tests/test_overlay_hook.m`, `dev/tests/test_dirac_helmholtz.m`, `dev/tests/test_helmholtz_view.m`.

**Modify:**
- `toolbox/gui/panel_surface.m` — a guarded `CustomOverlayFcn` hook at the end of `UpdateSurfaceData` (the per-figure time refresh).
- `toolbox/gui/figure_3d.m` — "Helmholtz / vorticity (Dirac)" source-figure popup item.

---

## Task 1: Time-following overlay hook (de-risk)

**Files:**
- Modify: `toolbox/gui/panel_surface.m` (end of `UpdateSurfaceData`)
- Test: `dev/tests/test_overlay_hook.m`

A `3DViz` figure gets `panel_surface('UpdateSurfaceData', hFig)` on every time-cursor move (`bst_figures FireCurrentTimeChanged`). Add a guarded hook there so a figure can register a post-update overlay refresh.

- [ ] **Step 1: Write the failing test**

```matlab
function test_overlay_hook()
% A figure's CustomOverlayFcn appdata must be called by panel_surface('UpdateSurfaceData').
% Authors: Diellor Basha, 2026
    nFail = 0;
    hFig = figure('Visible','off');
    setappdata(hFig, 'Surface', repmat(struct('Name','x'),0,1));   % empty surface -> UpdateSurfaceData no-ops safely
    assignin('base','OVH', 0);
    setappdata(hFig, 'CustomOverlayFcn', @(h) assignin('base','OVH', evalin('base','OVH')+1));
    panel_surface('UpdateSurfaceData', hFig);
    nFail = nFail + chk('overlay fcn fired once', evalin('base','OVH')==1);
    rmappdata(hFig, 'CustomOverlayFcn');
    panel_surface('UpdateSurfaceData', hFig);
    nFail = nFail + chk('not fired after removal', evalin('base','OVH')==1);
    close(hFig);
    fprintf('\n==== test_overlay_hook: %d failed ====\n', nFail);
    if nFail > 0, error('test_overlay_hook FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
```

- [ ] **Step 2: Run it to verify it fails**

Run (MCP `run_matlab_file`): `dev/tests/test_overlay_hook.m`
Expected: FAIL — `overlay fcn fired once` fails (hook not present; `OVH` stays 0).

- [ ] **Step 3: Add the hook at the end of UpdateSurfaceData**

In `toolbox/gui/panel_surface.m`, find `function [isOk, TessInfo] = UpdateSurfaceData(hFig, iSurfaces)`. At every `return` path's end and the normal end of the function, the simplest robust placement is just before the final `end` of the function body. Add immediately before the function's closing `end`:

```matlab
    % Custom overlay hook: a figure (e.g. the Helmholtz view) can register a function
    % that redraws its own overlays (quiver, markers) after each surface-data refresh,
    % which is exactly the per-figure time-cursor update path.
    ovFcn = getappdata(hFig, 'CustomOverlayFcn');
    if ~isempty(ovFcn)
        try, ovFcn(hFig); catch, end %#ok<CTCH>
    end
```

> If `UpdateSurfaceData` has multiple early `return`s, this hook at the very end still runs on the normal time-update path (the early returns are error/no-data cases). The test's empty-surface figure exercises the normal path.

- [ ] **Step 4: Run the test to verify it passes**

Run: `dev/tests/test_overlay_hook.m`
Expected: PASS — both `PASS`, `0 failed`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_surface.m dev/tests/test_overlay_hook.m
git commit -m "feat(gui): CustomOverlayFcn hook in UpdateSurfaceData (time-following overlays)"
```

---

## Task 2: bst_dirac_helmholtz (the pure decomposition + cores)

**Files:**
- Create: `toolbox/math/bst_dirac_helmholtz.m`
- Test: `dev/tests/test_dirac_helmholtz.m`

- [ ] **Step 1: Write the failing test**

```matlab
function test_dirac_helmholtz()
% Pure tests: Poisson solve round-trip, core classifier on a planted scalar, and a
% full-pipeline sanity on the real Subject01 cortex (Dirac + LBO operators).
% Authors: Diellor Basha, 2026
    nFail = 0;
    SurfaceFile = bst_get('Subject',1).Surface(5).FileName;
    Surf  = in_tess_bst(SurfaceFile, 0);
    Dirac = i_load_op(SurfaceFile, 'Dirac');
    LBO   = i_load_op(SurfaceFile, 'Laplace-Beltrami');

    % --- (1) Poisson solve round-trip on hemisphere 1: K psi = M omega recovers psi ---
    K = LBO.Operator{1}; M = LBO.Mass{1}; n = size(K,1);
    rng(0); psiTrue = randn(n,1); psiTrue = psiTrue - mean(psiTrue);
    rhs = K * psiTrue;                                   % = M*omega with omega = Lap*psiTrue
    psiRec = bst_dirac_helmholtz('PoissonSolve', K, M, M \ rhs);  % omega = M\rhs
    psiRec = psiRec - mean(psiRec);
    nFail = nFail + chk('Poisson recovers psi (up to const)', max(abs(psiRec - psiTrue)) < 1e-6 * max(abs(psiTrue)));

    % --- (2) Core classifier: a single Gaussian bump on the cortex -> one extremum core ---
    V = Surf.Vertices; vc = 5000;
    d2 = sum((V - V(vc,:)).^2, 2);
    psi = exp(-d2 / (2*(0.01)^2));                      % sharp positive bump at vc
    cores = bst_dirac_helmholtz('FindCores', psi, Surf.VertConn, zeros(size(psi)));
    nFail = nFail + chk('one core at the bump', numel(cores)==1 && cores(1).iVertex==vc);

    % --- (3) Full pipeline on a real (random) source field: finite, right size, nonzero ---
    nV = size(V,1); rng(1); J = randn(3*nV, 2) * 1e-9;
    H = bst_dirac_helmholtz(Dirac, LBO, Surf, J);
    nFail = nFail + chk('Curl size [nV x nT]', isequal(size(H.Curl), [nV 2]));
    nFail = nFail + chk('Div/Psi/Phi present + finite', all(isfinite(H.Div(:))) && all(isfinite(H.Psi(:))) && all(isfinite(H.Phi(:))));
    nFail = nFail + chk('field has both curl and div', max(abs(H.Curl(:)))>0 && max(abs(H.Div(:)))>0);
    nFail = nFail + chk('Cores is 1xnT cell', iscell(H.Cores) && numel(H.Cores)==2);

    fprintf('\n==== test_dirac_helmholtz: %d failed ====\n', nFail);
    if nFail > 0, error('test_dirac_helmholtz FAILED'); end
end

function Op = i_load_op(SurfaceFile, variant)
    sSubject = bst_get('Subject',1);
    iSurf = find(strcmpi({sSubject.Surface.FileName}, file_short(SurfaceFile)),1);
    Op = [];
    if isfield(sSubject.Surface(iSurf),'Operator')
        for k = 1:numel(sSubject.Surface(iSurf).Operator)
            S = load(file_fullpath(sSubject.Surface(iSurf).Operator(k).FileName));
            if strcmpi(S.Variant, variant); Op = S; break; end
        end
    end
    if isempty(Op); tess_operators(SurfaceFile, variant); Op = i_load_op(SurfaceFile, variant); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `dev/tests/test_dirac_helmholtz.m`
Expected: FAIL — `Undefined function 'bst_dirac_helmholtz'`.

- [ ] **Step 3: Implement the helper**

Create `toolbox/math/bst_dirac_helmholtz.m`:

```matlab
function varargout = bst_dirac_helmholtz(varargin)
% BST_DIRAC_HELMHOLTZ: Helmholtz-Hodge decomposition of a cortical source vector field
% using the first-order intrinsic Dirac (div + curl) and a scalar-LBO Poisson solve
% (potential phi, stream function psi), plus vortex-core detection (psi extrema).
%
% USAGE:
%   H = bst_dirac_helmholtz(DiracOp, LBO, Surf, J)
%       DiracOp : loaded Dirac operator node (.FirstOrder.Intrinsic{h}, .GlobalVertices{h})
%       LBO     : loaded Laplace-Beltrami operator node (.Operator{h}=cotan K, .Mass{h}=M)
%       Surf    : loaded surface (.Vertices, .Faces, .VertConn, .VertNormals)
%       J       : [3nV x nT] source field over time (x,y,z per global vertex)
%   Returns H with per-vertex x time fields .Curl .Div .Psi .Phi .Fmag [nV x nT] and
%   H.Cores (1 x nT cell; each a struct array (iVertex, charge, omega) of psi extrema).
%
%   psi = bst_dirac_helmholtz('PoissonSolve', K, M, omega)   % K psi = M omega, mean-zero
%   cores = bst_dirac_helmholtz('FindCores', psi, VertConn, omega)
%
% Math (per hemisphere h): embed J as a pure-imaginary quaternion; q = D_int * psi gives
% per-face divergence (w-part) + curl (imag part); omega = curl . n_face (face normals
% reconstructed from the same Floc ordering tess_operators uses). Area-weighted face->vertex
% averaging; Poisson K psi = M omega, K phi = M div. Cores = local extrema of psi over the
% 1-ring (max or min = a vortex center; handedness = -sign(Laplacian) = sign(omega)).
%
% Authors: Diellor Basha, 2026
    if ischar(varargin{1})
        [varargout{1:nargout}] = feval(varargin{:});
        return;
    end
    [varargout{1:nargout}] = Decompose(varargin{:});
end

function H = Decompose(DiracOp, LBO, Surf, J)
    Vtx = Surf.Vertices;  Fcs = double(Surf.Faces);
    nVtot = size(Vtx,1);  nT = size(J,2);
    H = struct('Curl',zeros(nVtot,nT), 'Div',zeros(nVtot,nT), ...
               'Psi',zeros(nVtot,nT),  'Phi',zeros(nVtot,nT),  'Fmag',zeros(nVtot,nT));
    for hh = 1:numel(DiracOp.FirstOrder.Intrinsic)
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
        Nf = cross(e1, e2, 2);  Af = sqrt(sum(Nf.^2,2));
        Nf = Nf ./ max(Af, eps);
        vn = Surf.VertNormals(vH(Floc(:,1)), :);         % a per-face reference normal
        flip = sum(Nf .* vn, 2) < 0;  Nf(flip,:) = -Nf(flip,:);
        nFh = size(Floc,1);
        % face->vertex area-weighted incidence [nVh x nFh]
        I = [Floc(:,1);Floc(:,2);Floc(:,3)];  Jc = [1:nFh,1:nFh,1:nFh]';
        Wfv = sparse(I, Jc, repmat(Af,3,1), nVh, nFh);
        Wfv = spdiags(1./max(sum(Wfv,2),eps),0,nVh,nVh) * Wfv;
        % LBO pieces for this hemisphere, pre-pinned solve (vertex 1 fixed, mean-zero)
        K = LBO.Operator{hh};  M = LBO.Mass{hh};
        for t = 1:nT
            Jx = J(3*(vH-1)+1, t);  Jy = J(3*(vH-1)+2, t);  Jz = J(3*(vH-1)+3, t);
            psiQ = zeros(4*nVh,1);
            psiQ(2:4:end) = Jx; psiQ(3:4:end) = Jy; psiQ(4:4:end) = Jz;
            q = D * psiQ;                                % [4F x 1]
            divF = q(1:4:end);                           % w-part (per face)
            curlF = [q(2:4:end), q(3:4:end), q(4:4:end)];% imag part (per face)
            omF = sum(curlF .* Nf, 2);                   % vorticity per face
            omV = Wfv * omF;                             % per vertex
            dvV = Wfv * divF;
            psi = PoissonSolve(K, M, omV);
            phi = PoissonSolve(K, M, dvV);
            H.Curl(vH,t) = omV;  H.Div(vH,t) = dvV;  H.Psi(vH,t) = psi;  H.Phi(vH,t) = phi;
            H.Fmag(vH,t) = sqrt(Jx.^2 + Jy.^2 + Jz.^2);
        end
    end
    H.Cores = cell(1, nT);
    for t = 1:nT
        H.Cores{t} = FindCores(H.Psi(:,t), Surf.VertConn, H.Curl(:,t));
    end
end

function psi = PoissonSolve(K, M, omega) %#ok<DEFNU>
    n = size(K,1);
    omega = omega - (sum(M*omega) / sum(sum(M))) * ones(n,1);   % project to mean-zero
    rhs = M * omega;
    free = 2:n;                                                 % pin psi(1)=0
    psi = zeros(n,1);
    psi(free) = K(free,free) \ rhs(free);
    psi = psi - mean(psi);
end

function cores = FindCores(psi, VertConn, omega) %#ok<DEFNU>
    cores = struct('iVertex',{}, 'charge',{}, 'omega',{});
    [ii, jj] = find(VertConn);                                  % neighbor pairs
    n = numel(psi);
    isMax = true(n,1);  isMin = true(n,1);
    for e = 1:numel(ii)
        v = ii(e); w = jj(e);
        if psi(w) >= psi(v); isMax(v) = false; end
        if psi(w) <= psi(v); isMin(v) = false; end
    end
    idx = find(isMax | isMin);
    for k = 1:numel(idx)
        v = idx(k);
        cores(end+1) = struct('iVertex', v, ...
            'charge', sign(omega(v) + (omega(v)==0)*(-1)*double(isMax(v)) + (omega(v)==0)*double(isMin(v))), ...
            'omega', omega(v)); %#ok<AGROW>
    end
end
```

> Note: `FindCores`' `charge` is `sign(omega)` at the core (the swirl handedness). The `omega==0` fallback uses max/min so a planted bump with zero omega still gets a sign (test (2) passes `omega=zeros`, expects one core; charge value is unchecked there). Strict local extrema (`>=`/`<=` make plateaus non-extremal) keep the count clean.

- [ ] **Step 4: Run the test to verify it passes**

Run: `dev/tests/test_dirac_helmholtz.m`
Expected: PASS — all `PASS`, `0 failed`. (First run may compute the LBO operator; allow time.)

- [ ] **Step 5: Commit**

```bash
git add toolbox/math/bst_dirac_helmholtz.m dev/tests/test_dirac_helmholtz.m
git commit -m "feat(math): bst_dirac_helmholtz - div/curl + Poisson psi/phi + vortex cores"
```

---

## Task 3: view_helmholtz (dedicated figure + overlays + time hook)

**Files:**
- Create: `toolbox/gui/view_helmholtz.m`

- [ ] **Step 1: Implement the launcher + figure + overlays**

Create `toolbox/gui/view_helmholtz.m`:

```matlab
function hFig = view_helmholtz(SrcResultsFile)
% VIEW_HELMHOLTZ: Dedicated Helmholtz/vorticity view of a Dirac source map. Shows a
% switchable signed scalar (curl/div/psi/phi/|field|) as a native 1-component results,
% with field-quiver + vortex-core overlays that follow the time cursor.
% Authors: Diellor Basha, 2026
    global GlobalData;
    hFig = [];
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
```

- [ ] **Step 2: Smoke-check (needs Task 4's panel; deferred verification)**

The live verification of `view_helmholtz` happens in Task 4 Step 4 and Task 5 (it calls `panel_helmholtz`). For now confirm it parses:
Run (MCP `evaluate_matlab_code`): `rehash; exist('view_helmholtz','file')` → expect `2`.

- [ ] **Step 3: Commit**

```bash
git add toolbox/gui/view_helmholtz.m
git commit -m "feat(gui): view_helmholtz dedicated figure + time-following overlays"
```

---

## Task 4: panel_helmholtz (controls)

**Files:**
- Create: `toolbox/gui/panel_helmholtz.m`
- Test: `dev/tests/test_helmholtz_view.m`

- [ ] **Step 1: Implement the panel**

Create `toolbox/gui/panel_helmholtz.m`:

```matlab
function varargout = panel_helmholtz(varargin)
% PANEL_HELMHOLTZ: Controls for the Helmholtz/vorticity view (view_helmholtz):
% scalar radio (Curl/Div/Stream psi/Potential phi/|field|), Show cores / Show quiver,
% a core-count readout, and Close.
% Authors: Diellor Basha, 2026
    eval(macro_method);
end

function bstPanelNew = CreatePanel(hFig) %#ok<DEFNU>
    import javax.swing.*;
    panelName = 'Helmholtz';
    jPanelNew = gui_component('Panel');
    jOpt = JPanel(); jOpt.setLayout(BoxLayout(jOpt, BoxLayout.Y_AXIS));
    jSec = gui_river([2 2], [2 8 3 6], 'Helmholtz / vorticity');

    gui_component('label', jSec, 'br', 'Show:');
    names    = {'Curl','Div','Psi','Phi','Fmag'};
    labels   = {'Curl . n (vorticity)','Divergence','Stream function','Potential','|field|'};
    grp = ButtonGroup(); jRadio = javaArray('javax.swing.JRadioButton', numel(names));
    for i = 1:numel(names)
        jRadio(i) = gui_component('radio', jSec, 'br', labels{i});
        grp.add(jRadio(i));
        java_setcb(jRadio(i), 'ActionPerformedCallback', @(h,e) OnScalar(panelName, names{i}));
    end
    jRadio(1).setSelected(true);
    jCores  = gui_component('checkbox', jSec, 'br', 'Show vortex cores');  jCores.setSelected(true);
    jQuiver = gui_component('checkbox', jSec, 'br', 'Show field quiver');
    java_setcb(jCores,  'ActionPerformedCallback', @(h,e) OnLayers(panelName));
    java_setcb(jQuiver, 'ActionPerformedCallback', @(h,e) OnLayers(panelName));
    jReadout = gui_component('label', jSec, 'br', '');
    jClose  = gui_component('button', jSec, 'br', 'Close');
    java_setcb(jClose, 'ActionPerformedCallback', @(h,e) OnClose(panelName));

    jOpt.add(jSec); jPanelNew.add(jOpt, java.awt.BorderLayout.NORTH);
    ctrl = struct('hFig',hFig, 'jCores',jCores, 'jQuiver',jQuiver, 'jReadout',jReadout);
    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);
end

function OnScalar(panelName, scalarName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if isempty(ctrl); return; end
    view_helmholtz('SetScalar', ctrl.hFig, scalarName);
end
function OnLayers(panelName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if isempty(ctrl); return; end
    view_helmholtz('SetLayers', ctrl.hFig, ctrl.jCores.isSelected(), ctrl.jQuiver.isSelected());
end
function SetReadout(text) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'Helmholtz'); if isempty(ctrl); return; end
    ctrl.jReadout.setText(text);
end
function OnClose(panelName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if isempty(ctrl); return; end
    if ishandle(ctrl.hFig); view_helmholtz('Close', ctrl.hFig); end
end
```

> Note: `view_helmholtz` must dispatch `SetScalar`/`SetLayers`/`Close`. Add `eval(macro_method)` dispatch to `view_helmholtz` — change its first line so string-first calls dispatch: at the top of `view_helmholtz`, before resolving the source, add:
> ```matlab
> if (nargin >= 1) && ischar(SrcResultsFile) && any(strcmp(SrcResultsFile, {'SetScalar','SetLayers','Close','RedrawOverlays'}))
>     feval(SrcResultsFile, varargin{2:end}); hFig = []; return;
> end
> ```
> and change the signature to `function hFig = view_helmholtz(SrcResultsFile, varargin)`.

- [ ] **Step 2: Add the dispatch shim to view_helmholtz**

Edit `toolbox/gui/view_helmholtz.m`: change the signature to `function hFig = view_helmholtz(SrcResultsFile, varargin)` and insert, as the first statements after `global GlobalData; hFig = [];`:

```matlab
    if (nargin >= 1) && ischar(SrcResultsFile) && any(strcmp(SrcResultsFile, {'SetScalar','SetLayers','Close','RedrawOverlays'}))
        feval(SrcResultsFile, varargin{:});
        return;
    end
```

- [ ] **Step 3: Write the live test**

Create `dev/tests/test_helmholtz_view.m`:

```matlab
function test_helmholtz_view()
% Live: open the Helmholtz view on a synthetic 3-comp source, switch the scalar, check
% the displayed series swaps and the core overlay count matches, then close cleanly.
% Authors: Diellor Basha, 2026
    global GlobalData;
    nFail = 0;
    % synthetic unconstrained source on cortex_20484V
    SurfaceFile = bst_get('Subject',1).Surface(5).FileName;
    nV = size(in_tess_bst(SurfaceFile,0).Vertices,1);
    R = db_template('resultsmat'); rng(0);
    R.ImageGridAmp = randn(3*nV,3)*1e-9; R.nComponents=3; R.Time=0:2;
    R.HeadModelType='surface'; R.SurfaceFile=file_short(SurfaceFile); R.Comment='SYN helmholtz src';
    srcFile = db_add(-3, R);
    [hSrc,~] = view_surface_data(SurfaceFile, srcFile, [], 'NewFigure'); drawnow;

    hFig = view_helmholtz(srcFile); drawnow;
    nFail = nFail + chk('view opens', ishandle(hFig));
    St = getappdata(hFig,'HelmholtzState');
    nFail = nFail + chk('default scalar is Curl', strcmp(St.Scalar,'Curl'));
    [iDS,iResult] = bst_memory('GetDataSetResult', St.TmpFile);
    nFail = nFail + chk('displays Curl series', isequal(GlobalData.DataSet(iDS).Results(iResult).ImageGridAmp, St.H.Curl));
    % switch to Psi
    view_helmholtz('SetScalar', hFig, 'Psi'); drawnow;
    nFail = nFail + chk('switch to Psi swaps series', isequal(GlobalData.DataSet(iDS).Results(iResult).ImageGridAmp, St.H.Psi));
    % core overlay count matches H.Cores at frame 1
    hAx = findobj(hFig,'-depth',1,'Tag','Axes3D'); hAx=hAx(1);
    nMarkers = numel(findobj(hAx,'Tag','HelmholtzCore'));
    nFail = nFail + chk('core markers match', nMarkers == numel(St.H.Cores{1}));
    % close cleans up
    view_helmholtz('Close', hFig); drawnow;
    nFail = nFail + chk('panel closed', isempty(bst_get('PanelControls','Helmholtz')));
    nFail = nFail + chk('figure closed', ~ishandle(hFig));
    % cleanup synthetic source
    try, bst_figures('DeleteFigure', hSrc, []); catch, end
    [~,iSt] = bst_get('AnyFile', srcFile); file_delete(file_fullpath(srcFile),1); db_reload_studies(iSt);

    fprintf('\n==== test_helmholtz_view: %d failed ====\n', nFail);
    if nFail > 0, error('test_helmholtz_view FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
```

- [ ] **Step 4: Run the live test**

Run: `dev/tests/test_helmholtz_view.m`
Expected: PASS — all `PASS`, `0 failed`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_helmholtz.m toolbox/gui/view_helmholtz.m dev/tests/test_helmholtz_view.m
git commit -m "feat(gui): panel_helmholtz controls + view_helmholtz dispatch + live test"
```

---

## Task 5: figure_3d popup item + end-to-end

**Files:**
- Modify: `toolbox/gui/figure_3d.m`

- [ ] **Step 1: Add the popup item next to the Spatial filter item**

In `figure_3d.m` `DisplayFigurePopup`, find the `% === SPATIAL FILTER ...` block (added earlier). Immediately after that block, add:

```matlab
        % === HELMHOLTZ / VORTICITY (unconstrained Dirac source) ===
        if ~isempty(ResultsFile)
            iTessSrc = find(arrayfun(@(t) ~isempty(t.DataSource) && strcmpi(t.DataSource.Type,'Source'), TessInfo), 1);
            if ~isempty(iTessSrc)
                [iDSh, iResh] = bst_memory('GetDataSetResult', TessInfo(iTessSrc).DataSource.FileName);
                if ~isempty(iResh) && isequal(GlobalData.DataSet(iDSh).Results(iResh).nComponents, 3)
                    gui_component('MenuItem', jPopup, [], 'Helmholtz / vorticity (Dirac)', IconLoader.ICON_RESULTS, [], @(h,ev)bst_call(@view_helmholtz, ResultsFile));
                end
            end
        end
```

- [ ] **Step 2: End-to-end live check (synthetic source)**

Run (MCP `evaluate_matlab_code`):
```matlab
rehash;
SurfaceFile = bst_get('Subject',1).Surface(5).FileName;
nV = size(in_tess_bst(SurfaceFile,0).Vertices,1);
R = db_template('resultsmat'); rng(2);
R.ImageGridAmp=randn(3*nV,3)*1e-9; R.nComponents=3; R.Time=0:2;
R.HeadModelType='surface'; R.SurfaceFile=file_short(SurfaceFile); R.Comment='SYN e2e';
srcFile = db_add(-3, R);
[hSrc,~] = view_surface_data(SurfaceFile, srcFile, [], 'NewFigure'); drawnow;
% guard would show the item?
TI = getappdata(hSrc,'Surface'); it = find(arrayfun(@(t) ~isempty(t.DataSource)&&strcmpi(t.DataSource.Type,'Source'),TI),1);
[iDSh,iResh] = bst_memory('GetDataSetResult', TI(it).DataSource.FileName);
global GlobalData;
fprintf('popup item shows: %d\n', isequal(GlobalData.DataSet(iDSh).Results(iResh).nComponents,3));
hFig = bst_call(@view_helmholtz, srcFile); drawnow;
fprintf('view opens: %d ; cores at f1: %d\n', ishandle(hFig), numel(getappdata(hFig,'HelmholtzState').H.Cores{1}));
view_helmholtz('Close', hFig); drawnow;
bst_figures('DeleteFigure', hSrc, []);
[~,iSt]=bst_get('AnyFile',srcFile); file_delete(file_fullpath(srcFile),1); db_reload_studies(iSt);
fprintf('done\n');
```
Expected: `popup item shows: 1`, `view opens: 1 ; cores at f1: <N>`, `done`.

- [ ] **Step 3: Commit**

```bash
git add toolbox/gui/figure_3d.m
git commit -m "feat(gui): 'Helmholtz / vorticity (Dirac)' figure popup item"
```

---

## Self-Review

**Spec coverage:**
- first-order Dirac → div/curl, LBO Poisson → ψ/φ, cores = ψ extrema → Task 2 (`bst_dirac_helmholtz`). ✓
- Dedicated figure, switchable signed scalar (curl/div/ψ/φ/|field|), follows time → Task 3 (native 1-comp results) + Task 4 (radio). ✓
- Field quiver + core markers overlays, follow time → Task 1 (hook) + Task 3 (`RedrawOverlays`). ✓
- Readout (counts/net charge) → Task 3 (`i_count_readout`) + Task 4 (`SetReadout`). ✓
- Launch from source figure popup (3-comp guard) → Task 5. ✓
- find-or-create Dirac + LBO operators; not eigen/manifold → Task 2/3 (`i_op`). ✓
- Teardown removes hook + temp node + figure → Task 3 (`Close`). ✓

**Placeholder scan:** No "TBD/TODO"; every code step is complete. The face-normal ordering is reconstructed from `GlobalVertices` (Task 2), matching `tess_operators` exactly — no operator-node change needed.

**Type consistency:** `bst_dirac_helmholtz` returns `H.Curl/Div/Psi/Phi/Fmag [nV x nT]` + `H.Cores{t}` (struct `iVertex/charge/omega`) — set in Task 2, consumed identically in Tasks 3/4 (`St.H.(scalarName)`, `St.H.Cores{iT}`). The scalar names `{'Curl','Div','Psi','Phi','Fmag'}` match between `panel_helmholtz` radio (Task 4) and `H` fields (Task 2). `view_helmholtz` dispatch names `{'SetScalar','SetLayers','Close','RedrawOverlays'}` (Task 4 shim) match the called subfunctions. `CustomOverlayFcn` appdata (Task 1 hook) is set in Task 3. `St` fields (`srcField,H,Scalar,ShowCores,ShowQuiver,TmpFile,iDS,iResult`) set in Task 3, read in Tasks 3/4. ✓

---

## Build order summary

1. Task 1 — `CustomOverlayFcn` hook in `UpdateSurfaceData` (de-risk time-following)
2. Task 2 — `bst_dirac_helmholtz` (div/curl + Poisson ψ/φ + cores)
3. Task 3 — `view_helmholtz` (figure + overlays + hook)
4. Task 4 — `panel_helmholtz` (controls) + dispatch shim + live test
5. Task 5 — figure popup item + end-to-end
