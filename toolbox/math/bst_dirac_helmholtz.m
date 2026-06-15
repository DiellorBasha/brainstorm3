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
%       Decomposes a single frame Jt [3nV x 1] -> struct with per-vertex [nV x 1] fields
%       .Curl .Div .Psi .Phi .Fmag and .Cores (struct array (iVertex, charge, omega)).
%
%   H = bst_dirac_helmholtz(DiracOp, LBO, Surf, J)
%       Whole-series convenience (Prepare once, Frame per column). J [3nV x nT]; returns
%       H.Curl/Div/Psi/Phi/Fmag [nV x nT] and H.Cores (1 x nT cell). (Batch primitive.)
%
%   psi   = bst_dirac_helmholtz('PoissonSolve', K, M, omega)   % K psi = M omega, mean-zero
%   cores = bst_dirac_helmholtz('FindCores', psi, VertConn, omega)
%
% Math (per hemisphere h): embed J as a pure-imaginary quaternion; q = D_int * psi gives
% per-face divergence (w-part) + curl (imag part); omega = curl . n_face (face normals
% reconstructed from the same Floc ordering tess_operators uses). Area-weighted face->vertex
% averaging; Poisson K psi = M omega, K phi = M div. Cores = local extrema of psi over the
% 1-ring (max or min = a vortex center; handedness = sign(omega) at the core).
%
% Authors: Diellor Basha, 2026
    if ischar(varargin{1})
        [varargout{1:nargout}] = feval(varargin{:});
        return;
    end
    [varargout{1:nargout}] = Decompose(varargin{:});
end

%% ===== PREPARE: build + cache the per-face gradient, the consistent stiffness factor =====
% Consistent discrete Hodge: the gradient G, its OWN adjoint divergence, and the stiffness
% A = G' diag(2*area) G are one operator family -- so the decomposition residual is genuine
% harmonic + boundary, not operator mismatch. (The Dirac node is used only to partition the
% hemispheres via GlobalVertices; div/curl come from G's adjoint, not the Dirac.)
function Op = Prepare(DiracOp, LBO, Surf) %#ok<DEFNU>
    Vtx = Surf.Vertices;  Fcs = double(Surf.Faces);
    nVtot = size(Vtx,1);
    nH = numel(DiracOp.GlobalVertices);
    Op = struct();
    Op.nVtot    = nVtot;
    Op.VertConn = Surf.VertConn;
    [Op.vH, Op.Floc, Op.Nf, Op.Af, Op.Wfv, Op.Gx, Op.Gy, Op.Gz, ...
     Op.cholA, Op.free, Op.av] = deal(cell(1,nH));
    for hh = 1:nH
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
        % consistent stiffness A = G' diag(Af) G + its pinned (vertex 1) Cholesky factor
        A = Gx'*(Af.*Gx) + Gy'*(Af.*Gy) + Gz'*(Af.*Gz);
        A = (A + A')/2;                                  % kill round-off asymmetry for chol
        free = (2:nVh)';
        Op.vH{hh}=vH; Op.Floc{hh}=Floc; Op.Nf{hh}=Nf; Op.Af{hh}=Af; Op.Wfv{hh}=Wfv;
        Op.Gx{hh}=Gx; Op.Gy{hh}=Gy; Op.Gz{hh}=Gz;
        Op.cholA{hh} = decomposition(A(free,free), 'chol');
        Op.free{hh}  = free;
        Op.av{hh}    = full(sum(LBO.Mass{hh}, 2));        % lumped vertex mass (energy weights)
    end
end

%% ===== FRAME: discrete Hodge projection of one frame (gradient + its adjoint divergence) =====
function Ht = Frame(Op, Jt) %#ok<DEFNU>
    nVtot = Op.nVtot;
    z1 = zeros(nVtot,1);  z3 = zeros(nVtot,3);
    Ht = struct('Curl',z1,'Div',z1,'Psi',z1,'Phi',z1,'Fmag',z1,'Hmag',z1, ...
                'Vtot',z3,'Virr',z3,'Vsol',z3,'Vharm',z3);
    harmNum = 0;  harmDen = 0;
    for hh = 1:numel(Op.vH)
        vH = Op.vH{hh};  Floc = Op.Floc{hh};  Nf = Op.Nf{hh};  Af = Op.Af{hh};  av = Op.av{hh};
        Gx = Op.Gx{hh};  Gy = Op.Gy{hh};  Gz = Op.Gz{hh};
        Jv = [Jt(3*(vH-1)+1), Jt(3*(vH-1)+2), Jt(3*(vH-1)+3)];          % [nVh x 3] vertex field
        Jf = (Jv(Floc(:,1),:) + Jv(Floc(:,2),:) + Jv(Floc(:,3),:))/3;   % [nFh x 3] face field
        % weak divergence and weak curl (gradient's own adjoint, area-weighted)
        rhsI = Gx'*(Af.*Jf(:,1)) + Gy'*(Af.*Jf(:,2)) + Gz'*(Af.*Jf(:,3));
        nxJ  = cross(Nf, Jf, 2);
        rhsS = -(Gx'*(Af.*nxJ(:,1)) + Gy'*(Af.*nxJ(:,2)) + Gz'*(Af.*nxJ(:,3)));
        phi = i_solveA(Op.cholA{hh}, rhsI, Op.free{hh}, numel(vH));     % potential   (irrotational)
        psi = i_solveA(Op.cholA{hh}, rhsS, Op.free{hh}, numel(vH));     % stream func. (solenoidal)
        Virr = Op.Wfv{hh} * [Gx*phi, Gy*phi, Gz*phi];                   % grad(phi) -> vertex
        Vsol = Op.Wfv{hh} * cross(Nf, [Gx*psi, Gy*psi, Gz*psi], 2);     % n x grad(psi) -> vertex
        Vharm = Jv - Virr - Vsol;                                       % exact residual
        omV = rhsS ./ max(av, eps);   dvV = rhsI ./ max(av, eps);       % vorticity / divergence (marker signs)
        Ht.Curl(vH)=omV;  Ht.Div(vH)=dvV;  Ht.Psi(vH)=psi;  Ht.Phi(vH)=phi;
        Ht.Fmag(vH)=sqrt(sum(Jv.^2,2));  Ht.Hmag(vH)=sqrt(sum(Vharm.^2,2));
        Ht.Vtot(vH,:)=Jv;  Ht.Virr(vH,:)=Virr;  Ht.Vsol(vH,:)=Vsol;  Ht.Vharm(vH,:)=Vharm;
        harmNum = harmNum + sum(av .* sum(Vharm.^2,2));
        harmDen = harmDen + sum(av .* sum(Jv.^2,2));
    end
    Ht.HarmFrac = harmNum / max(harmDen, eps);
    Ht.Cores    = FindCores(Ht.Psi, Op.VertConn, Ht.Curl);    % vortex cores (sign = vorticity)
    Ht.Sources  = FindCores(Ht.Phi, Op.VertConn, Ht.Div);     % sources/sinks (sign = divergence)
end

%% ===== pinned solve of the consistent stiffness (vertex 1 fixed, mean-zero gauge) =====
function x = i_solveA(cholA, rhs, free, n)
    x = zeros(n,1);
    x(free) = cholA \ rhs(free);
    x = x - mean(x);
end

%% ===== whole-series convenience (Prepare once, Frame per column) =====
function H = Decompose(DiracOp, LBO, Surf, J)
    Op = Prepare(DiracOp, LBO, Surf);
    nVtot = Op.nVtot;  nT = size(J,2);
    H = struct('Curl',zeros(nVtot,nT), 'Div',zeros(nVtot,nT), ...
               'Psi',zeros(nVtot,nT),  'Phi',zeros(nVtot,nT),  'Fmag',zeros(nVtot,nT));
    H.Cores = cell(1, nT);
    for t = 1:nT
        Ht = Frame(Op, J(:,t));
        H.Curl(:,t)=Ht.Curl; H.Div(:,t)=Ht.Div; H.Psi(:,t)=Ht.Psi;
        H.Phi(:,t)=Ht.Phi;   H.Fmag(:,t)=Ht.Fmag;  H.Cores{t}=Ht.Cores;
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

%% ===== core detection: local extrema of psi over the 1-ring =====
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
        if omega(v) ~= 0
            ch = sign(omega(v));                                % handedness from vorticity
        else
            ch = -double(isMax(v)) + double(isMin(v));          % fallback: max=-1, min=+1
        end
        cores(end+1) = struct('iVertex', v, 'charge', ch, 'omega', omega(v)); %#ok<AGROW>
    end
end
