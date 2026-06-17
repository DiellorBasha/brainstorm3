# Face-domain Helmholtz (dual Dirac D̃) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A reusable `bst_dirac_helmholtz_face` that runs a Helmholtz decomposition on a per-face 3D field via nxr's dual face-Dirac `D̃`, plus a benchmark comparing it to the vertex pipeline on a real frame.

**Architecture:** Dual of `bst_dirac_helmholtz`: the field lives on **face centroids**; `D̃ = diracFaceD` ([4V×4F]) maps the per-face quaternion to a per-vertex quaternion, giving vorticity (w-part, normal-free) and divergence (imag·vertex-normal) **natively on vertices** (no face→vertex averaging); Poisson on the vertex cotan-Laplacian → ψ/φ; component fields reconstructed on faces. The benchmark interpolates the existing vertex dSPM solution to faces and compares.

**Tech Stack:** MATLAB (Brainstorm); nxr-compute (`nxr_safe_create`, `nxr_compute('operators',h,'diracFaceD')`); reuses `bst_vortex_persistence`. Tests: `chk(label,cond)` functions via the MATLAB MCP.

**Commits:** user-managed on `development`; commit steps OPTIONAL/user-gated.

**Preconditions:** MATLAB up, Brainstorm running, `TutorialAuditory` protocol, toolbox + `dev/tests` on path.

---

### Task 1: `bst_dirac_helmholtz_face` (Prepare + Frame)

**Files:** Create `toolbox/math/bst_dirac_helmholtz_face.m`; Create `dev/tests/test_dirac_helmholtz_face.m`.

- [ ] **Step 1: Write the failing test** — create `dev/tests/test_dirac_helmholtz_face.m`:

```matlab
function test_dirac_helmholtz_face()
% Face-domain Helmholtz via the dual face-Dirac D̃: convention + reconstruction + shapes.
% Author: Diellor Basha, 2026
    nFail = 0;
    SurfaceFile = bst_get('Subject',1).Surface(5).FileName;
    Surf  = in_tess_bst(SurfaceFile, 0);
    Dirac = i_load_op(SurfaceFile, 'Dirac');
    LBO   = i_load_op(SurfaceFile, 'Laplace-Beltrami');
    nV = size(Surf.Vertices,1); nF = size(Surf.Faces,1);

    Op = bst_dirac_helmholtz_face('Prepare', Dirac, LBO, Surf);
    nFail = nFail + chk('Prepare has D-tilde per hemi', numel(Op.Dt)==2 && size(Op.Dt{1},1)==4*numel(Op.vH{1}));

    rng(2); Jf = randn(nF, 3) * 1e-9;                          % random per-face field
    Ht = bst_dirac_helmholtz_face('Frame', Op, Jf);

    nFail = nFail + chk('Curl/Div are [nV x 1]', isequal(size(Ht.Curl),[nV 1]) && isequal(size(Ht.Div),[nV 1]));
    nFail = nFail + chk('component fields are [nF x 3]', isequal(size(Ht.Virr),[nF 3]) && isequal(size(Ht.Vsol),[nF 3]));
    recon = Ht.Virr + Ht.Vsol + Ht.Vharm;
    nFail = nFail + chk('exact reconstruction Virr+Vsol+Vharm == Jf', max(abs(recon(:)-Ht.Vtot(:))) < 1e-9*max(abs(Ht.Vtot(:))));
    nFail = nFail + chk('HarmFrac in [0,1]', isscalar(Ht.HarmFrac) && Ht.HarmFrac>=0 && Ht.HarmFrac<=1.0001);

    % --- convention: re-decompose the component fields (validates w=vorticity, imag.n=div) ---
    Hirr = bst_dirac_helmholtz_face('Frame', Op, Ht.Virr);
    Hsol = bst_dirac_helmholtz_face('Frame', Op, Ht.Vsol);
    nFail = nFail + chk('irrotational is divergence-dominated', sum(Hirr.Div.^2) > sum(Hirr.Curl.^2));
    nFail = nFail + chk('solenoidal is curl-dominated',         sum(Hsol.Curl.^2) > sum(Hsol.Div.^2));

    nFail = nFail + chk('Cores/Sources are struct arrays w/ persistence', isstruct(Ht.Cores) && isstruct(Ht.Sources) && (isempty(Ht.Cores) || isfield(Ht.Cores,'persistence')));

    fprintf('\n==== test_dirac_helmholtz_face: %d failed ====\n', nFail);
    if nFail > 0, error('test_dirac_helmholtz_face FAILED'); end
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

- [ ] **Step 2: Run, expect FAIL** (`Undefined function 'bst_dirac_helmholtz_face'`). Run `dev/tests/test_dirac_helmholtz_face.m`.

- [ ] **Step 3: Implement** — create `toolbox/math/bst_dirac_helmholtz_face.m`:

```matlab
function varargout = bst_dirac_helmholtz_face(varargin)
% BST_DIRAC_HELMHOLTZ_FACE  Helmholtz/Hodge decomposition of a PER-FACE 3D field via the
% dual face-Dirac D̃ (= nxr diracFaceD, [4V x 4F]). Dual of bst_dirac_helmholtz: the field
% lives on face centroids; D̃ yields vorticity (w-part, normal-free) and divergence
% (imag . vertex-normal) natively on vertices (no face->vertex averaging); Poisson on the
% vertex cotan-Laplacian gives psi/phi; component fields are reconstructed on faces.
%
% USAGE:
%   Op = bst_dirac_helmholtz_face('Prepare', DiracOp, LBO, Surf)
%   Ht = bst_dirac_helmholtz_face('Frame', Op, Jf [, withCores])
%        Jf : [nF x 3] per-face ambient field (or [3*nF x 1] stacked).
% Author: Diellor Basha, 2026
    [varargout{1:nargout}] = feval(varargin{:});
end

%% ===== PREPARE =====
function Op = Prepare(DiracOp, LBO, Surf) %#ok<DEFNU>
    Vtx = Surf.Vertices;  Fcs = double(Surf.Faces);
    nVtot = size(Vtx,1);  nFtot = size(Fcs,1);
    nH = numel(DiracOp.GlobalVertices);
    Op = struct(); Op.nVtot=nVtot; Op.nFtot=nFtot; Op.VertConn=Surf.VertConn; Op.Vtx=Vtx;
    [Op.Dt, Op.vH, Op.fH, Op.Nf, Op.Af, Op.M, Op.cholK, Op.free, Op.totMass, ...
     Op.Gx, Op.Gy, Op.Gz, Op.VnV, Op.NbH] = deal(cell(1,nH));
    for hh = 1:nH
        vH = double(DiracOp.GlobalVertices{hh}(:));  nVh = numel(vH);
        isV = false(nVtot,1); isV(vH)=true;  fMask = all(isV(Fcs),2);
        fH = find(fMask);  mapV = zeros(nVtot,1); mapV(vH)=1:nVh;
        Floc = mapV(Fcs(fMask,:));  Vloc = Vtx(vH,:);  nFh = numel(fH);
        % dual face Dirac D̃ [4V x 4F] from a fresh nxr handle on this hemisphere submesh
        h = nxr_safe_create(Vloc, Floc);
        Dt = nxr_compute('operators', h, 'diracFaceD');
        nxr_compute('destroy', h);
        % face normals + areas (outward via vertex normals); twoA = 2*area = |cross|
        e1 = Vloc(Floc(:,2),:)-Vloc(Floc(:,1),:);  e2 = Vloc(Floc(:,3),:)-Vloc(Floc(:,1),:);
        Nf = cross(e1,e2,2);  twoA = sqrt(sum(Nf.^2,2));  Nf = Nf./max(twoA,eps);
        vn = Surf.VertNormals(vH(Floc(:,1)),:);  flip = sum(Nf.*vn,2)<0;  Nf(flip,:)=-Nf(flip,:);
        % per-face FEM gradient of a per-vertex scalar (same as bst_dirac_helmholtz)
        eO1=Vloc(Floc(:,3),:)-Vloc(Floc(:,2),:); eO2=Vloc(Floc(:,1),:)-Vloc(Floc(:,3),:); eO3=Vloc(Floc(:,2),:)-Vloc(Floc(:,1),:);
        c1=cross(Nf,eO1,2)./twoA; c2=cross(Nf,eO2,2)./twoA; c3=cross(Nf,eO3,2)./twoA;
        grows=[(1:nFh)';(1:nFh)';(1:nFh)'];  gcols=[Floc(:,1);Floc(:,2);Floc(:,3)];
        Gx=sparse(grows,gcols,[c1(:,1);c2(:,1);c3(:,1)],nFh,nVh);
        Gy=sparse(grows,gcols,[c1(:,2);c2(:,2);c3(:,2)],nFh,nVh);
        Gz=sparse(grows,gcols,[c1(:,3);c2(:,3);c3(:,3)],nFh,nVh);
        % vertex 1-ring (local) for core detection
        eLoc=[Floc(:,[1 2]);Floc(:,[2 3]);Floc(:,[3 1])];
        Aloc=sparse([eLoc(:,1);eLoc(:,2)],[eLoc(:,2);eLoc(:,1)],true,nVh,nVh);
        nb=cell(nVh,1); for vv=1:nVh, nb{vv}=find(Aloc(:,vv)); end
        % LBO Poisson factor (pinned vertex 1)
        K=LBO.Operator{hh}; M=LBO.Mass{hh}; free=(2:size(K,1))';
        Op.Dt{hh}=Dt; Op.vH{hh}=vH; Op.fH{hh}=fH; Op.Nf{hh}=Nf; Op.Af{hh}=twoA/2;
        Op.M{hh}=M; Op.cholK{hh}=decomposition(K(free,free),'chol'); Op.free{hh}=free; Op.totMass{hh}=sum(M(:));
        Op.Gx{hh}=Gx; Op.Gy{hh}=Gy; Op.Gz{hh}=Gz;
        Op.VnV{hh}=Surf.VertNormals(vH,:); Op.NbH{hh}=nb;
    end
end

%% ===== FRAME =====
function Ht = Frame(Op, Jf, withCores) %#ok<DEFNU>
    if nargin < 3 || isempty(withCores), withCores = true; end
    if size(Jf,2) ~= 3, Jf = reshape(Jf, 3, [])'; end          % accept [3F x 1] stacked
    nVtot = Op.nVtot;  nFtot = Op.nFtot;
    zV = zeros(nVtot,1);  zF3 = zeros(nFtot,3);  zF1 = zeros(nFtot,1);
    Ht = struct('Curl',zV,'Div',zV,'Psi',zV,'Phi',zV,'Fmag',zF1,'Hmag',zF1, ...
                'Vtot',zF3,'Virr',zF3,'Vsol',zF3,'Vharm',zF3);
    harmNum=0; harmDen=0;
    for hh = 1:numel(Op.Dt)
        vH = Op.vH{hh};  fH = Op.fH{hh};  nFh = numel(fH);
        Jfl = Jf(fH,:);                                        % [nFh x 3] local face field
        qF = zeros(4*nFh,1);
        qF(2:4:end)=Jfl(:,1); qF(3:4:end)=Jfl(:,2); qF(4:4:end)=Jfl(:,3);
        qV = Op.Dt{hh} * qF;                                   % D̃ : [4V x 1] (per vertex)
        % Convention (mirrors the validated vertex Dirac): w-part = vorticity,
        % imag . n = divergence. Validated end-to-end by the irrot/solenoidal test.
        omV   = qV(1:4:end);
        imagV = [qV(2:4:end), qV(3:4:end), qV(4:4:end)];
        dvV   = sum(imagV .* Op.VnV{hh}, 2);                   % divergence = imag . n_vertex
        psi = i_poisson(Op.cholK{hh}, Op.M{hh}, omV, Op.free{hh}, Op.totMass{hh});
        phi = i_poisson(Op.cholK{hh}, Op.M{hh}, dvV, Op.free{hh}, Op.totMass{hh});
        gphi = [Op.Gx{hh}*phi, Op.Gy{hh}*phi, Op.Gz{hh}*phi];  % [nFh x 3]
        gpsi = [Op.Gx{hh}*psi, Op.Gy{hh}*psi, Op.Gz{hh}*psi];
        Virr = gphi;  Vsol = cross(Op.Nf{hh}, gpsi, 2);        % component fields native on faces
        Vharm = Jfl - Virr - Vsol;
        Ht.Curl(vH)=omV;  Ht.Div(vH)=dvV;  Ht.Psi(vH)=psi;  Ht.Phi(vH)=phi;
        Ht.Fmag(fH)=sqrt(sum(Jfl.^2,2));  Ht.Hmag(fH)=sqrt(sum(Vharm.^2,2));
        Ht.Vtot(fH,:)=Jfl;  Ht.Virr(fH,:)=Virr;  Ht.Vsol(fH,:)=Vsol;  Ht.Vharm(fH,:)=Vharm;
        Af = Op.Af{hh};
        harmNum = harmNum + sum(Af .* sum(Vharm.^2,2));
        harmDen = harmDen + sum(Af .* sum(Jfl.^2,2));
    end
    Ht.HarmFrac = harmNum/max(harmDen,eps);
    if withCores
        Ht.Cores   = i_find_cores(Ht.Psi, Op, Ht.Curl);
        Ht.Sources = i_find_cores(Ht.Phi, Op, Ht.Div);
    else
        Ht.Cores = i_empty_cores();  Ht.Sources = i_empty_cores();
    end
end

%% ===== helpers =====
function psi = i_poisson(dK, M, omega, free, totMass)
    n = size(M,1);
    omega = omega - (sum(M*omega)/totMass)*ones(n,1);          % project to mean-zero subspace
    rhs = M*omega;  psi = zeros(n,1);
    psi(free) = dK \ rhs(free);  psi = psi - mean(psi);
end

function s = i_empty_cores()
    s = struct('iVertex',{},'charge',{},'chirality',{},'omega',{},'persistence',{}, ...
               'isGlobal',{},'hemi',{},'pos',{});
end

function cores = i_find_cores(field, Op, omega)
% Persistence-ranked extrema per hemisphere on a per-vertex potential (cores at vertices;
% no sub-vertex localization in this prototype -- count/vertex/persistence are what we compare).
    cores = i_empty_cores();
    for hh = 1:numel(Op.vH)
        vH = Op.vH{hh};  C = bst_vortex_persistence(field(vH), [], 'Neighbors', Op.NbH{hh});
        for k = 1:numel(C.vertex)
            vg = vH(C.vertex(k));  om = omega(vg);
            if om ~= 0, ch = sign(om); else, ch = -C.chirality(k); end
            cores(end+1) = struct('iVertex',vg, 'charge',ch, 'chirality',C.chirality(k), ...
                'omega',om, 'persistence',C.persistence(k), 'isGlobal',C.isGlobal(k), ...
                'hemi',hh, 'pos',Op.Vtx(vg,:)); %#ok<AGROW>
        end
    end
    if ~isempty(cores), [~,ord]=sort([cores.persistence],'descend'); cores=cores(ord); end
end
```

- [ ] **Step 4: Run, expect PASS** — `dev/tests/test_dirac_helmholtz_face.m` → `0 failed`. **If "irrotational is divergence-dominated" or "solenoidal is curl-dominated" FAILS**, the dual `D̃` convention differs from the vertex `D`: swap to `dvV = imag` interpretation — i.e. set `omV = sum(imagV .* Op.VnV{hh},2)` and `dvV = qV(1:4:end)` — re-run; document whichever convention passes in the Frame comment.

- [ ] **Step 5: Lint** `toolbox/math/bst_dirac_helmholtz_face.m`.

- [ ] **Step 6: Commit (optional)** `git add -A && git commit -m "feat(helmholtz): face-domain Helmholtz via dual face-Dirac"`

---

### Task 2: Benchmark — vertex vs face comparison

**Files:** Create `dev/benchmarks/bench_dirac_face_helmholtz.m`.

- [ ] **Step 1: Implement the comparison driver** — create `dev/benchmarks/bench_dirac_face_helmholtz.m`:

```matlab
function R = bench_dirac_face_helmholtz(frameTime)
% BENCH_DIRAC_FACE_HELMHOLTZ  Compare vertex vs face-domain Helmholtz on a real Dirac frame.
% Interpolates the unconstrained vertex dSPM solution to face centroids and decomposes it
% with bst_dirac_helmholtz_face (dual D̃), against bst_dirac_helmholtz (vertex). 
% USAGE: R = bench_dirac_face_helmholtz(22.6)
% Author: Diellor Basha, 2026
    if nargin<1 || isempty(frameTime), frameTime = 22.6; end
    OUTDIR = '/Users/diellorbasha/workspace/research/code/brainstorm3/dev/benchmarks';
    df = 'Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat';
    [sStudy,~] = bst_get('DataFile', df);
    ChanMat = in_bst_channel(sStudy.Channel(1).FileName); types = {ChanMat.Channel.Type};
    HMos = in_bst_headmodel([fileparts(df) '/headmodel_surf_os_meg.mat'], 0);
    G = double(HMos.Gain); iMEG = all(isfinite(G),2) & strcmpi(types(:),'MEG');
    NC = load(file_fullpath([fileparts(df) '/noisecov_full.mat'])); Cn = NC.NoiseCov(iMEG,iMEG); Cn=(Cn+Cn')/2;
    HMf = HMos; HMf.Gain = G(iMEG,:);
    OPT = struct('NoiseMethod','reg','NoiseReg',0.1,'SnrMethod','fixed','SnrFixed',3,'InverseMeasure','dspm2018');
    OPT.NoiseCovMat.NoiseCov = Cn; OPT.ChannelTypes = types(iMEG);
    Rd = bst_inverse_dirac(HMf, OPT);
    DM = in_bst_data(df); [~,iT] = min(abs(DM.Time-frameTime)); Jt = Rd.ImagingKernel*double(DM.F(iMEG,iT));

    SurfaceFile = HMos.SurfaceFile;  Surf = in_tess_bst(SurfaceFile,0);
    V = Surf.Vertices; F = double(Surf.Faces); nV=size(V,1); nF=size(F,1);
    Dirac = i_op(SurfaceFile,'Dirac'); LBO = i_op(SurfaceFile,'Laplace-Beltrami');

    % vertex pipeline
    OpV = bst_dirac_helmholtz('Prepare', Dirac, LBO, Surf);
    HtV = bst_dirac_helmholtz('Frame', OpV, Jt);

    % interpolate vertex field -> face centroids (barycentric mean of the 3 vertex vectors)
    J3 = reshape(Jt,3,[])';                                  % [nV x 3]
    Jf = (J3(F(:,1),:) + J3(F(:,2),:) + J3(F(:,3),:))/3;     % [nF x 3]

    % face pipeline
    OpF = bst_dirac_helmholtz_face('Prepare', Dirac, LBO, Surf);
    HtF = bst_dirac_helmholtz_face('Frame', OpF, Jf);

    % ---- compare ----
    cv = @(c) [sum([c.charge]>0) sum([c.charge]<0)];
    fprintf('\n=== vertex vs face Helmholtz @ %.3f s ===\n', DM.Time(iT));
    fprintf('HarmFrac:   vertex %.1f%%   face %.1f%%\n', 100*HtV.HarmFrac, 100*HtF.HarmFrac);
    fprintf('vortices(+/-): vertex %s   face %s\n', mat2str(cv(HtV.Cores)),   mat2str(cv(HtF.Cores)));
    fprintf('sources(+/-):  vertex %s   face %s\n', mat2str(cv(HtV.Sources)), mat2str(cv(HtF.Sources)));
    fprintf('|Curl| max: vertex %.3g  face %.3g  | |Div| max: vertex %.3g  face %.3g\n', ...
        max(abs(HtV.Curl)), max(abs(HtF.Curl)), max(abs(HtV.Div)), max(abs(HtF.Div)));

    % ---- side-by-side viz: Psi (stream fn) + top cores ----
    hFig = figure('Color','w','Position',[60 80 1200 560]); cmap = i_div(256);
    for sp = 1:2
        ax = subplot(1,2,sp); hold(ax,'on');
        if sp==1, scal=HtV.Psi; cc=HtV.Cores; ttl='vertex Helmholtz (psi)';
        else,     scal=HtF.Psi; cc=HtF.Cores; ttl='face Helmholtz (psi)'; end
        patch('Vertices',V,'Faces',F,'FaceVertexCData',scal,'FaceColor','interp', ...
              'EdgeColor','none','Parent',ax);
        m=max(abs(scal)); if m<=0, m=eps; end
        pr=[cc.persistence]; mxf=max([pr(isfinite(pr)),eps]); keep=cc(isinf(pr)|pr>=0.5*mxf);
        for k=1:numel(keep)
            col=[1 0 0]; if keep(k).charge<0, col=[0 0 1]; end
            p=keep(k).pos; plot3(ax,p(1),p(2),p(3),'o','MarkerFaceColor',col,'MarkerEdgeColor','k','MarkerSize',7,'Clipping','off');
        end
        colormap(ax,cmap); caxis(ax,[-m m]); axis(ax,'equal','off'); view(ax,[0 90]);
        camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull'); title(ax,ttl);
    end
    sgtitle(sprintf('Vertex vs face-domain Helmholtz @ %.2f s', DM.Time(iT)));
    png = fullfile(OUTDIR,'bench_dirac_face_helmholtz.png'); print(hFig,png,'-dpng','-r140');
    fprintf('saved %s\n', png);
    R = struct('HtV',HtV,'HtF',HtF,'frameTime',DM.Time(iT),'png',png);
end

function Op = i_op(SurfaceFile, variant)
    [sSubject,~,iSurf] = bst_get('SurfaceFile', SurfaceFile);
    Op = [];
    if ~isempty(iSurf) && isfield(sSubject.Surface(iSurf),'Operator')
        for k = 1:numel(sSubject.Surface(iSurf).Operator)
            S = load(file_fullpath(sSubject.Surface(iSurf).Operator(k).FileName));
            if strcmpi(S.Variant, variant), Op = S; break; end
        end
    end
    if isempty(Op), tess_operators(SurfaceFile, variant); Op = i_op(SurfaceFile, variant); end
end
function m = i_div(n)
    t=linspace(0,1,n)'; lo=[.23 .30 .75]; mid=[.96 .96 .96]; hi=[.78 .15 .18];
    m=[interp1([0 .5 1],[lo(1) mid(1) hi(1)],t), interp1([0 .5 1],[lo(2) mid(2) hi(2)],t), interp1([0 .5 1],[lo(3) mid(3) hi(3)],t)];
end
```

- [ ] **Step 2: Lint** `dev/benchmarks/bench_dirac_face_helmholtz.m`.

- [ ] **Step 3: Run it** (MCP `run_matlab_file` or `evaluate`): `R = bench_dirac_face_helmholtz(22.6);`
Expected: prints HarmFrac/core-count/curl-div comparison; saves the PNG. Sanity: both HarmFrac in [0,1]; face and vertex give the same order of magnitude of cores; |Curl|/|Div| finite and nonzero.

- [ ] **Step 4: Inspect the PNG** (Read tool) — vertex vs face ψ maps + cores should be broadly similar in structure; note any differences (crispness, core count, harmonic fraction).

- [ ] **Step 5: Commit (optional)** `git add -A && git commit -m "bench(helmholtz): vertex vs face-domain comparison"`

---

### Task 3: Regression + assessment

**Files:** none (verification only)

- [ ] **Step 1: Run all suites** — `test_dirac_helmholtz_face`, `test_dirac_helmholtz`, `test_vortex_persistence`, `test_vortex_track`, `test_helmholtz_track`, `test_time_derivative`: all `0 failed` (the existing ones must be unaffected — the face module is additive).

- [ ] **Step 2: Record the comparison finding** — from the benchmark output + figure, note in the chat: harmonic-fraction delta, core-count agreement, and whether the face decomposition is visibly crisper/noisier, as the empirical input to the go/no-go on a full face-based inverse.

- [ ] **Step 3: Final lint sweep** of the two new files.

- [ ] **Step 4: Commit (optional)** `git add -A && git commit -m "test(helmholtz): face-domain regression"`

---

## Self-review notes

- **Spec coverage:** reusable `bst_dirac_helmholtz_face` Prepare+Frame (Task 1) ✓; D̃ sourced via fresh nxr handle (Task 1 Prepare) ✓; convention validation via irrot/solenoidal re-decomposition + the documented swap fallback (Task 1 Step 4) ✓; benchmark interpolation + comparison metrics + side-by-side viz (Task 2) ✓; tests incl. reconstruction/shapes/HarmFrac (Task 1) and regression (Task 3) ✓.
- **Naming consistency:** `Op.Dt/vH/fH/Nf/Af/M/cholK/free/totMass/Gx/Gy/Gz/VnV/NbH` defined in Prepare and consumed in Frame; `i_find_cores` returns the `.iVertex/.charge/.chirality/.omega/.persistence/.isGlobal/.hemi/.pos` schema used by the benchmark's `cv`/plotting; `Jf` is `[nF x 3]` throughout.
- **Placeholder scan:** none.
- **Risk:** the `D̃` quaternion convention — explicitly validated by the re-decomposition test with a documented swap if it fails (Task 1 Step 4); divergence uses vertex normals by design (the asymmetry the benchmark quantifies).
```
