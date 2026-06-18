function varargout = bst_dirac_helmholtz_face(varargin)
% BST_DIRAC_HELMHOLTZ_FACE  Face-NATIVE Helmholtz/Hodge decomposition of a PER-FACE 3D
% field, the exact dual of bst_dirac_helmholtz. The field lives on face centroids; the
% scalar potentials (stream psi, potential phi) ALSO live on faces, solved by a Poisson
% on the face Laplacian K̃ (nxr lapFace) with a barycentric dual-mesh gradient G̃ (nxr
% gradFace) for reconstruction -- so K̃ = G̃ᵀ W_F G̃ by construction and the decomposition
% is self-consistent (HarmFrac->0 for a field in the range of G̃, genus-0 harmonic = const).
%
% Math (per hemisphere): with G̃ = gradFace [3F×F], SkewG = n_f × G̃ (skew/stream gradient),
% W_F = face-area mass, the Hodge decomposition is the COUPLED variational projection of J
% onto range([G̃ SkewG]) -- jointly minimise ||J - G̃ phi - SkewG psi||_W_F. The discrete
% skew-gradient is NOT exactly divergence-free (G̃ᵀ W_F SkewG ≠ 0), so an INDEPENDENT
% phi/psi solve would leak ~few% between channels; the coupled normal system
%   [K  C ] [phi]   [G̃ᵀ W_F J ]            K = G̃ᵀ W_F G̃ (= SkewGᵀ W_F SkewG),
%   [Cᵀ K ] [psi] = [SkewGᵀ W_F J],        C = G̃ᵀ W_F SkewG,
% removes the leakage so Vharm is the genuine (genus-0 → ~0) harmonic residual.
% Components: Virr = G̃ phi (irrotational), Vsol = SkewG psi (solenoidal), Vharm = J-Virr-Vsol.
%
% USAGE:
%   Op = bst_dirac_helmholtz_face('Prepare', DiracOp, LBO, Surf)
%   Ht = bst_dirac_helmholtz_face('Frame', Op, Jf [, withCores])
%        Jf : [nF x 3] per-face ambient field (or [3*nF x 1] stacked).
% Returns per-FACE scalars .Curl (vorticity) .Div .Psi .Phi [nF x 1]; component vector
% fields [nF x 3] .Vtot .Virr .Vsol .Vharm; .HarmFrac; cores/sources (psi/phi extrema on
% the dual face-adjacency, positions at face centroids). DiracOp/LBO are accepted for API
% parity with the vertex pipeline but are unused (the operators come from nxr gradFace/lapFace).
% Author: Diellor Basha, 2026
    [varargout{1:nargout}] = feval(varargin{:});
end

%% ===== PREPARE =====
function Op = Prepare(DiracOp, LBO, Surf) %#ok<INUSD,DEFNU>
    Vtx = Surf.Vertices;  Fcs = double(Surf.Faces);
    nVtot = size(Vtx,1);  nFtot = size(Fcs,1);
    nH = numel(DiracOp.GlobalVertices);
    Op = struct(); Op.nVtot=nVtot; Op.nFtot=nFtot; Op.Vtx=Vtx; Op.Fcs=Fcs;
    [Op.G, Op.SkewG, Op.WF, Op.cholA, Op.freeA, Op.nFh, Op.fH, Op.vH, Op.Nf, Op.Af, Op.Cf, Op.NbF] = deal(cell(1,nH));
    for hh = 1:nH
        vH = double(DiracOp.GlobalVertices{hh}(:));  nVh = numel(vH);
        isV = false(nVtot,1); isV(vH)=true;  fMask = all(isV(Fcs),2);
        fH = find(fMask);  mapV = zeros(nVtot,1); mapV(vH)=1:nVh;
        Floc = mapV(Fcs(fMask,:));  Vloc = Vtx(vH,:);  nFh = numel(fH);
        % nxr face-native gradient on this closed hemisphere submesh
        h = nxr_safe_create(Vloc, Floc);
        G = nxr_compute('operators', h, 'gradFace');          % [3F x F]
        nxr_compute('destroy', h);
        % per-face geometry (normals outward via vertex normals, areas, centroids)
        e1 = Vloc(Floc(:,2),:)-Vloc(Floc(:,1),:);  e2 = Vloc(Floc(:,3),:)-Vloc(Floc(:,1),:);
        Nf = cross(e1,e2,2);  twoA = sqrt(sum(Nf.^2,2));  Nf = Nf./max(twoA,eps);  Af = twoA/2;
        vn = Surf.VertNormals(vH(Floc(:,1)),:);  flip = sum(Nf.*vn,2)<0;  Nf(flip,:)=-Nf(flip,:);
        Cf = (Vloc(Floc(:,1),:)+Vloc(Floc(:,2),:)+Vloc(Floc(:,3),:))/3;     % local centroids
        % SkewG = n_f x G : per-face block rotation Rn (cross-product matrix) times G
        SkewG = i_block_rotation(Nf, nFh) * G;
        % coupled variational Hodge normal matrix A = M' W_F M, M = [G SkewG] [3F x 2F]
        WFv = repelem(Af,3);  M = [G, SkewG];
        A = M' * (WFv .* M);  A = (A + A')/2;                 % [2F x 2F] sym PSD
        % nullspace = constants in each block: pin phi(1) and psi(1) -> SPD
        freeA = setdiff((1:2*nFh)', [1; nFh+1]);
        Op.G{hh}=G; Op.SkewG{hh}=SkewG; Op.WF{hh}=Af; Op.nFh{hh}=nFh;
        Op.cholA{hh}=decomposition(A(freeA,freeA),'chol');  Op.freeA{hh}=freeA;
        Op.fH{hh}=fH; Op.vH{hh}=vH; Op.Nf{hh}=Nf; Op.Af{hh}=Af; Op.Cf{hh}=Cf;
        Op.NbF{hh}=i_dual_adjacency(Floc, nFh);               % faces sharing an edge
    end
end

%% ===== FRAME =====
function Ht = Frame(Op, Jf, withCores) %#ok<DEFNU>
    if nargin < 3 || isempty(withCores), withCores = true; end
    if size(Jf,2) ~= 3, Jf = reshape(Jf, 3, [])'; end          % accept [3F x 1] stacked
    nFtot = Op.nFtot;
    zF1 = zeros(nFtot,1);  zF3 = zeros(nFtot,3);
    Ht = struct('Curl',zF1,'Div',zF1,'Psi',zF1,'Phi',zF1,'Fmag',zF1,'Hmag',zF1, ...
                'Vtot',zF3,'Virr',zF3,'Vsol',zF3,'Vharm',zF3);
    harmNum=0; harmDen=0;
    for hh = 1:numel(Op.G)
        fH = Op.fH{hh};  Af = Op.Af{hh};  nFh = Op.nFh{hh};
        Jl = Jf(fH,:);  Jcol = reshape(Jl', [], 1);            % [3F x 1] (x,y,z per face)
        WF = repelem(Af,3);                                    % [3F] area weights
        divS  = Op.G{hh}'    * (WF .* Jcol);                   % G' W_F J   -> [F] divergence source
        curlS = Op.SkewG{hh}'* (WF .* Jcol);                   % SkewG' W_F J -> [F] curl source
        % coupled variational Hodge solve A[phi;psi]=[divS;curlS], A pinned at phi(1),psi(1)
        rhs = [divS; curlS];  sol = zeros(2*nFh,1);
        sol(Op.freeA{hh}) = Op.cholA{hh} \ rhs(Op.freeA{hh});
        phi = sol(1:nFh);       phi = phi - mean(phi);
        psi = sol(nFh+1:2*nFh); psi = psi - mean(psi);
        Virr = reshape(Op.G{hh}    *phi, 3, [])';              % [F x 3] irrotational
        Vsol = reshape(Op.SkewG{hh}*psi, 3, [])';              % [F x 3] solenoidal
        Vharm = Jl - Virr - Vsol;
        Ht.Curl(fH)=curlS;  Ht.Div(fH)=divS;  Ht.Psi(fH)=psi;  Ht.Phi(fH)=phi;
        Ht.Fmag(fH)=sqrt(sum(Jl.^2,2));  Ht.Hmag(fH)=sqrt(sum(Vharm.^2,2));
        Ht.Vtot(fH,:)=Jl;  Ht.Virr(fH,:)=Virr;  Ht.Vsol(fH,:)=Vsol;  Ht.Vharm(fH,:)=Vharm;
        harmNum = harmNum + sum(Af .* sum(Vharm.^2,2));
        harmDen = harmDen + sum(Af .* sum(Jl.^2,2));
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
function Rn = i_block_rotation(Nf, nFh)
% Block-diagonal [3F x 3F] sparse: per face the cross-product matrix [n]_x so that
% (Rn * v) restricted to face f equals n_f x v_f. Rows/cols 3k-2:3k per face k.
    nx=Nf(:,1); ny=Nf(:,2); nz=Nf(:,3);
    r=(1:nFh)';  b0=3*(r-1);
    I=[b0+1; b0+1; b0+2; b0+2; b0+3; b0+3];
    J=[b0+2; b0+3; b0+1; b0+3; b0+1; b0+2];
    S=[-nz;  ny;   nz;  -nx;  -ny;   nx];
    Rn = sparse(I, J, S, 3*nFh, 3*nFh);
end

function NbF = i_dual_adjacency(Floc, nFh)
% Dual face-adjacency: faces sharing an edge. Build edge->faces incidence, then
% link the (<=2) faces on each edge. Returns a cell{nFh} of neighbor face lists.
    E = [Floc(:,[1 2]); Floc(:,[2 3]); Floc(:,[3 1])];
    E = sort(E,2);  fId = repmat((1:nFh)',3,1);
    [~,~,ic] = unique(E,'rows');
    A = sparse(ic, fId, true, max(ic), nFh);             % edge x face incidence
    Adj = (A'*A) > 0;  Adj(1:nFh+1:end) = false;          % face-face via shared edge
    NbF = cell(nFh,1);
    for f=1:nFh, NbF{f}=find(Adj(:,f)); end
end

function s = i_empty_cores()
    s = struct('iVertex',{},'iFace',{},'charge',{},'chirality',{},'omega',{}, ...
               'persistence',{},'isGlobal',{},'hemi',{},'pos',{});
end

function cores = i_find_cores(field, Op, omega)
% Persistence-ranked extrema of a per-FACE potential over the dual face-adjacency
% (a max or min of psi is a vortex core); charge = sign of the vorticity/divergence
% at that face; position = face centroid.
    cores = i_empty_cores();
    for hh = 1:numel(Op.fH)
        fH = Op.fH{hh};  Cf = Op.Cf{hh};
        C = bst_vortex_persistence(field(fH), [], 'Neighbors', Op.NbF{hh});
        for k = 1:numel(C.vertex)
            fl = C.vertex(k);  fg = fH(fl);  om = omega(fg);
            if om ~= 0, ch = sign(om); else, ch = -C.chirality(k); end
            cores(end+1) = struct('iVertex',fg, 'iFace',fg, 'charge',ch, ...
                'chirality',C.chirality(k), 'omega',om, 'persistence',C.persistence(k), ...
                'isGlobal',C.isGlobal(k), 'hemi',hh, 'pos',Cf(fl,:)); %#ok<AGROW>
        end
    end
    if ~isempty(cores), [~,ord]=sort([cores.persistence],'descend'); cores=cores(ord); end
end
