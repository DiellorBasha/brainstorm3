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
% 1-ring (max or min = a vortex center; handedness = sign(omega) at the core).
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
        % LBO pieces for this hemisphere
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
        if omega(v) ~= 0
            ch = sign(omega(v));                                % handedness from vorticity
        else
            ch = -double(isMax(v)) + double(isMin(v));          % fallback: max=-1, min=+1
        end
        cores(end+1) = struct('iVertex', v, 'charge', ch, 'omega', omega(v)); %#ok<AGROW>
    end
end
