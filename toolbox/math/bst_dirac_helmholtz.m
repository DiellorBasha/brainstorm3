function varargout = bst_dirac_helmholtz(varargin)
% BST_DIRAC_HELMHOLTZ: Helmholtz-Hodge decomposition of a cortical source vector field
% using the first-order intrinsic Dirac (div + curl) and a scalar-LBO Poisson solve
% (potential phi, stream function psi), plus vortex-core detection (psi extrema).
%
% The expensive part (the cotan Poisson solve) is factorized ONCE via Prepare and reused
% per frame via Frame -- so the interactive view decomposes only the displayed time step,
% and a batch process can loop Frame over the series with a single factorization.
%
% USAGE:
%   Op = bst_dirac_helmholtz('Prepare', DiracOp, LBO, Surf)
%       Builds + caches the per-hemisphere operator pieces and the Cholesky factor of the
%       pinned cotan stiffness. DiracOp/LBO are the loaded operator nodes, Surf the loaded
%       tessellation (.Vertices, .Faces, .VertConn, .VertNormals).
%
%   Ht = bst_dirac_helmholtz('Frame', Op, Jt)
%       Decomposes a single frame Jt [3nV x 1] into the Hodge components. Returns per-vertex
%       [nV x 1] scalars .Curl (vorticity) .Div .Psi .Phi .Fmag .Hmag; component vector
%       fields [nV x 3] .Vtot .Virr (grad phi) .Vsol (curl psi) .Vharm (residual); the
%       harmonic energy fraction .HarmFrac; and singular points .Cores (psi extrema) /
%       .Sources (phi extrema), each a struct array (iVertex, charge, omega).
%
%   H = bst_dirac_helmholtz(DiracOp, LBO, Surf, J)
%       Whole-series convenience (Prepare once, Frame per column). J [3nV x nT]; returns the
%       scalar fields H.Curl/Div/Psi/Phi/Fmag [nV x nT] and H.Cores (1 x nT cell). (Batch
%       primitive; the per-frame component vector fields are produced by Frame.)
%
%   psi   = bst_dirac_helmholtz('PoissonSolve', K, M, omega)   % K psi = M omega, mean-zero
%   cores = bst_dirac_helmholtz('FindCoresOp', field, Op, omega)   % per-hemisphere, persistence-ranked
%   cores = bst_dirac_helmholtz('FindCores', field, VertConn, omega)  % legacy full-mesh entry
%       Cores struct array, sorted by persistence desc: .iVertex .charge .omega
%       .persistence .isGlobal .birth .death .pos (1x3 sub-vertex xyz; NaN for the
%       legacy entry, which lacks geometry).
%
% Math (per hemisphere h): embed J as a pure-imaginary quaternion; q = D_int * psi gives, per
% face, the VORTICITY (the w-part) and -- via the imaginary part dotted with the face normal
% -- the DIVERGENCE (validated against pure gradient/rotational fields; NOT the reverse).
% Area-weighted face->vertex averaging; Poisson K psi = M omega (stream psi from vorticity),
% K phi = M div (potential phi from divergence); component fields grad(phi), n x grad(psi),
% and the harmonic residual J - grad(phi) - n x grad(psi). Cores = local extrema of psi over
% the 1-ring (a max or a min is a vortex center; handedness = sign(omega) at the core).
%
% Authors: Diellor Basha, 2026
    if ischar(varargin{1})
        [varargout{1:nargout}] = feval(varargin{:});
        return;
    end
    [varargout{1:nargout}] = Decompose(varargin{:});
end

%% ===== PREPARE: build + cache the static operator pieces and the Cholesky factor =====
function Op = Prepare(DiracOp, LBO, Surf) %#ok<DEFNU>
    Vtx = Surf.Vertices;  Fcs = double(Surf.Faces);
    nVtot = size(Vtx,1);
    nH = numel(DiracOp.FirstOrder.Intrinsic);
    Op = struct();
    Op.nVtot    = nVtot;
    Op.VertConn = Surf.VertConn;
    [Op.D, Op.vH, Op.Nf, Op.Wfv, Op.M, Op.cholK, Op.free, Op.totMass, ...
     Op.Gx, Op.Gy, Op.Gz, Op.NbH, Op.VtxH, Op.VnH] = deal(cell(1,nH));
    Op.Vtx = Vtx;
    for hh = 1:nH
        D  = DiracOp.FirstOrder.Intrinsic{hh};           % [4F x 4V]
        vH = double(DiracOp.GlobalVertices{hh}(:));
        nVh = numel(vH);
        % reconstruct the hemisphere's local faces in tess_operators' ordering
        isV = false(nVtot,1); isV(vH) = true;
        fMask = all(isV(Fcs), 2);
        mapV  = zeros(nVtot,1); mapV(vH) = 1:nVh;
        Floc  = mapV(Fcs(fMask, :));   Vloc = Vtx(vH, :);
        % per-hemisphere 1-ring (LOCAL indexing) + geometry for core detection
        eLoc = [Floc(:,[1 2]); Floc(:,[2 3]); Floc(:,[3 1])];
        Aloc = sparse([eLoc(:,1);eLoc(:,2)], [eLoc(:,2);eLoc(:,1)], true, nVh, nVh);
        nbLoc = cell(nVh,1);
        for vv = 1:nVh, nbLoc{vv} = find(Aloc(:,vv)); end
        Op.NbH{hh}  = nbLoc;
        Op.VtxH{hh} = Vloc;
        Op.VnH{hh}  = Surf.VertNormals(vH, :);
        % face normals (oriented outward via the surface vertex normals)
        e1 = Vloc(Floc(:,2),:) - Vloc(Floc(:,1),:);
        e2 = Vloc(Floc(:,3),:) - Vloc(Floc(:,1),:);
        Nf = cross(e1, e2, 2);  Af = sqrt(sum(Nf.^2,2));   % Af = 2*area (un-normalized)
        Nf = Nf ./ max(Af, eps);
        vn = Surf.VertNormals(vH(Floc(:,1)), :);         % a per-face reference normal
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
        grows = [(1:nFh)';(1:nFh)';(1:nFh)'];  gcols = [Floc(:,1);Floc(:,2);Floc(:,3)];
        Gx = sparse(grows, gcols, [c1(:,1);c2(:,1);c3(:,1)], nFh, nVh);
        Gy = sparse(grows, gcols, [c1(:,2);c2(:,2);c3(:,2)], nFh, nVh);
        Gz = sparse(grows, gcols, [c1(:,3);c2(:,3);c3(:,3)], nFh, nVh);
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

%% ===== FRAME: decompose a single time frame using the cached factor =====
function Ht = Frame(Op, Jt, withCores) %#ok<DEFNU>
    if nargin < 3 || isempty(withCores), withCores = true; end   % skip core detection when false
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
        q = Op.D{hh} * psiQ;                              % [4F x 1]
        % Dirac convention (validated against pure gradient/rotational fields): the w-part
        % is the VORTICITY and the imaginary part . n is the DIVERGENCE (not the reverse).
        omF   = q(1:4:end);                               % vorticity = w-part (per face)
        imagF = [q(2:4:end), q(3:4:end), q(4:4:end)];     % imaginary part (per face)
        divF  = sum(imagF .* Op.Nf{hh}, 2);               % divergence = imag . n_face
        omV = Op.Wfv{hh} * omF;   dvV = Op.Wfv{hh} * divF;
        psi = i_poisson(Op.cholK{hh}, Op.M{hh}, omV, Op.free{hh}, Op.totMass{hh});  % stream from vorticity
        phi = i_poisson(Op.cholK{hh}, Op.M{hh}, dvV, Op.free{hh}, Op.totMass{hh});  % potential from divergence
        % component vector fields: grad(phi) (irrotational), skew-grad(psi) (solenoidal)
        gphi = [Op.Gx{hh}*phi, Op.Gy{hh}*phi, Op.Gz{hh}*phi];   % [nF x 3]
        gpsi = [Op.Gx{hh}*psi, Op.Gy{hh}*psi, Op.Gz{hh}*psi];
        skew = cross(Op.Nf{hh}, gpsi, 2);                       % n x grad(psi)
        Virr = Op.Wfv{hh} * gphi;                               % face->vertex [nVh x 3]
        Vsol = Op.Wfv{hh} * skew;
        Jv   = [Jx Jy Jz];
        Vharm = Jv - Virr - Vsol;                               % exact residual
        Ht.Curl(vH)=omV;  Ht.Div(vH)=dvV;  Ht.Psi(vH)=psi;  Ht.Phi(vH)=phi;
        Ht.Fmag(vH)=sqrt(Jx.^2+Jy.^2+Jz.^2);  Ht.Hmag(vH)=sqrt(sum(Vharm.^2,2));
        Ht.Vtot(vH,:)=Jv;  Ht.Virr(vH,:)=Virr;  Ht.Vsol(vH,:)=Vsol;  Ht.Vharm(vH,:)=Vharm;
        av = full(sum(Op.M{hh},2));                            % lumped vertex mass
        harmNum = harmNum + sum(av .* sum(Vharm.^2,2));
        harmDen = harmDen + sum(av .* sum(Jv.^2,2));
    end
    Ht.HarmFrac = harmNum / max(harmDen, eps);
    if withCores
        Ht.Cores    = FindCoresOp(Ht.Psi, Op, Ht.Curl);      % vortex cores (sign = vorticity)
        Ht.Sources  = FindCoresOp(Ht.Phi, Op, Ht.Div);       % sources/sinks (sign = divergence)
    else
        Ht.Cores    = i_empty_cores();                       % detection skipped (GUI markers/track off)
        Ht.Sources  = i_empty_cores();
    end
end

%% ===== whole-series convenience (Prepare once, Frame per column) =====
function H = Decompose(DiracOp, LBO, Surf, J)
    Op = Prepare(DiracOp, LBO, Surf);
    nVtot = Op.nVtot;  nT = size(J,2);
    H = struct('Curl',zeros(nVtot,nT), 'Div',zeros(nVtot,nT), ...
               'Psi',zeros(nVtot,nT),  'Phi',zeros(nVtot,nT),  'Fmag',zeros(nVtot,nT));
    H.Cores = cell(1, nT);
    H.Sources = cell(1, nT);
    for t = 1:nT
        Ht = Frame(Op, J(:,t));
        H.Curl(:,t)=Ht.Curl; H.Div(:,t)=Ht.Div; H.Psi(:,t)=Ht.Psi;
        H.Phi(:,t)=Ht.Phi;   H.Fmag(:,t)=Ht.Fmag;  H.Cores{t}=Ht.Cores;  H.Sources{t}=Ht.Sources;
    end
end

%% ===== Poisson solves =====
function psi = PoissonSolve(K, M, omega) %#ok<DEFNU>
    % Convenience: factor then solve (Frame uses the cached factor instead).
    free = (2:size(K,1))';
    dK = decomposition(K(free,free), 'chol');
    psi = i_poisson(dK, M, omega, free, sum(M(:)));
end

function psi = i_poisson(dK, M, omega, free, totMass)
    n = size(M,1);
    omega = omega - (sum(M*omega) / totMass) * ones(n,1);  % project to the mean-zero subspace
    rhs = M * omega;
    psi = zeros(n,1);
    psi(free) = dK \ rhs(free);                            % pinned solve (vertex 1 fixed)
    psi = psi - mean(psi);
end

%% ===== core detection: persistence-ranked extrema, per hemisphere =====
function cores = FindCoresOp(field, Op, omega) %#ok<DEFNU>
% Per-hemisphere persistence detection + sub-vertex localization. Returns a struct
% array (sorted by persistence desc) carrying both the legacy fields (iVertex,
% charge, omega) and the new ones (persistence, isGlobal, birth, death, pos).
    cores = i_empty_cores();
    for hh = 1:numel(Op.vH)
        vH = Op.vH{hh};  nb = Op.NbH{hh};  Vloc = Op.VtxH{hh};  Vn = Op.VnH{hh};
        fl = field(vH);
        C  = bst_vortex_persistence(fl, [], 'Neighbors', nb);
        for k = 1:numel(C.vertex)
            vloc = C.vertex(k);  vg = vH(vloc);
            cores(end+1) = i_make_core(vg, omega(vg), C.chirality(k), C.persistence(k), ...
                C.isGlobal(k), C.birth(k), C.death(k), ...
                i_subvertex(vloc, fl, nb, Vloc, Vn), hh); %#ok<AGROW>
        end
    end
    if ~isempty(cores)
        [~, ord] = sort([cores.persistence], 'descend');
        cores = cores(ord);
    end
end

% Legacy entry (field + VertConn): runs on the full mesh; pos = NaN (no geometry).
function cores = FindCores(field, VertConn, omega) %#ok<DEFNU>
    nV = numel(field);
    [ii, jj] = find(VertConn);
    nb = accumarray(ii, jj, [nV 1], @(x){x}, {zeros(0,1)});
    C  = bst_vortex_persistence(field, [], 'Neighbors', nb);
    cores = i_empty_cores();
    for k = 1:numel(C.vertex)
        v = C.vertex(k);
        cores(end+1) = i_make_core(v, omega(v), C.chirality(k), C.persistence(k), ...
            C.isGlobal(k), C.birth(k), C.death(k), nan(1,3), 1); %#ok<AGROW>
    end
end

function s = i_empty_cores()
    s = struct('iVertex',{},'charge',{},'chirality',{},'omega',{},'persistence',{}, ...
               'isGlobal',{},'birth',{},'death',{},'pos',{},'hemi',{});
end

function s = i_make_core(vg, om, chirality, persistence, isGlobal, birth, death, pos, hemi)
    % charge = vorticity sign (display continuity); chirality = topological sweep sign.
    if om ~= 0, ch = sign(om); else, ch = -chirality; end   % preserve legacy fallback
    s = struct('iVertex',vg, 'charge',ch, 'chirality',chirality, 'omega',om, ...
               'persistence',persistence, 'isGlobal',logical(isGlobal), ...
               'birth',birth, 'death',death, 'pos',pos, 'hemi',hemi);
end

% Sub-vertex localization: quadratic fit of FIELD over the 1-ring in a tangent chart.
function p = i_subvertex(vloc, field, nb, Vloc, Vn)
    v0 = Vloc(vloc,:);  p = v0;
    ns = nb{vloc};
    if numel(ns) < 5, return; end
    n = Vn(vloc,:);  n = n / max(norm(n), eps);
    e0 = [1 0 0]; if abs(e0*n') > 0.9, e0 = [0 1 0]; end
    t1 = e0 - (e0*n')*n; t1 = t1/norm(t1);  t2 = cross(n, t1);
    off = Vloc(ns,:) - v0;
    u = off*t1';  w = off*t2';  g = field(ns) - field(vloc);
    A = [u, w, 0.5*u.^2, u.*w, 0.5*w.^2];
    if rcond(A'*A) < 1e-10, return; end
    c = A \ g;                       % [b; c; d; e; f]
    H = [c(3) c(4); c(4) c(5)];
    if rcond(H) < 1e-8, return; end
    uw = -H \ [c(1); c(2)];
    if norm(uw) > max(sqrt(u.^2 + w.^2)), return; end   % reject runaway -> keep vertex
    p = v0 + uw(1)*t1 + uw(2)*t2;
end
