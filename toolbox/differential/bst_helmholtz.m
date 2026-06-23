function varargout = bst_helmholtz(varargin)
% BST_HELMHOLTZ: Helmholtz-Hodge decomposition of a cortical source vector field over the
% VERTEX domain (per-vertex field; scalar potentials on vertices via the cotan LBO) or the
% FACE domain (per-face field; potentials on faces via the nxr face Laplacian / coupled
% variational Hodge). One domain-dispatching orchestrator (the dual of the two former
% bst_dirac_helmholtz / _face files), routed by 'Domain'.
%
% SINGLE SOURCE OF TRUTH FOR GEOMETRY: per-face normals/areas/centroids and per-vertex
% positions/normals come from the canonical manifold_ node (Embedded group, schemaVersion>=2).
% Operators (Dirac D, LBO K/M, gradFace) come from the operator_ node. Face/vertex
% CONNECTIVITY (face->vertex index triples) comes from Surf.Faces (topology, not geometry).
% I/O-free: the caller resolves and passes the loaded nodes (tess_manifold / operator nodes /
% in_tess_bst); this function never touches the DB or disk.
%
% The expensive Cholesky factorization is built ONCE in Prepare and reused per Frame, so the
% interactive view decomposes only the displayed step and a batch process loops Frame.
%
% USAGE:
%   Op = bst_helmholtz('Prepare', OperatorNode, ManifoldMat, Surf, 'Domain','vertex')
%       OperatorNode = {DiracNode, LBONode}  (vertex): DiracNode.FirstOrder.Intrinsic{hh} [4F x 4V],
%                      LBONode.Operator{hh} (cotan K), LBONode.Mass{hh} (M).
%   Op = bst_helmholtz('Prepare', FaceOperatorNode, ManifoldMat, Surf, 'Domain','face')
%       FaceOperatorNode.FaceAux{hh}.GradFace [3F x F] (Hodge-Face / Dirac-Face operator).
%       ManifoldMat : db_template('manifoldmat'); Surf : tessellation (Surf.Faces/.VertConn used).
%
%   Ht = bst_helmholtz('Frame', Op, J [, withCores])
%       Decomposes a single frame; routes on Op.Domain. J is [3nV x 1] (vertex) or
%       [nF x 3] / [3nF x 1] (face). Returns per-domain scalars .Curl .Div .Psi .Phi
%       .Fmag .Hmag; component vector fields .Vtot .Virr .Vsol .Vharm; .HarmFrac; and
%       .Cores / .Sources (psi/phi extrema) struct arrays.
%
%   H  = bst_helmholtz('Decompose', OperatorNode, ManifoldMat, Surf, J)   % vertex whole-series
%   cores = bst_helmholtz('FindCoresOp', field, Op, omega)                % vertex, persistence-ranked
%   cores = bst_helmholtz('FindCores', field, VertConn, omega)            % vertex legacy full-mesh
%
% Math (per hemisphere): VERTEX -- embed J as a pure-imaginary quaternion; q = D*psiQ gives
% per face the VORTICITY (w-part) and the DIVERGENCE (imag . face-normal); two scalar-LBO
% Poisson solves recover the stream psi and potential phi; components grad(phi), n x grad(psi),
% harmonic residual. FACE -- coupled variational projection of J onto range([gradFace SkewG])
% with the face-area mass; one coupled solve recovers (phi,psi); components gradFace*phi,
% SkewG*psi, residual. Cores = local extrema of psi/phi (handedness = sign of vorticity/div).
%
% Authors: Diellor Basha, 2026
    [varargout{1:nargout}] = feval(varargin{:});
end

%% ===== PREPARE: domain dispatch =====
function Op = Prepare(OperatorNode, ManifoldMat, Surf, varargin) %#ok<DEFNU>
    Domain = 'vertex';
    for i = 1:2:numel(varargin)
        if strcmpi(varargin{i}, 'Domain'), Domain = lower(varargin{i+1}); end
    end
    switch Domain
        case 'vertex', Op = i_prepare_vertex(OperatorNode, ManifoldMat, Surf);
        case 'face',   Op = i_prepare_face(OperatorNode, ManifoldMat, Surf);
        otherwise, error('bst_helmholtz:badDomain', 'Domain must be ''vertex'' or ''face''.');
    end
    Op.Domain = Domain;
end

%% ===== FRAME: domain dispatch =====
function Ht = Frame(Op, J, withCores) %#ok<DEFNU>
    if nargin < 3 || isempty(withCores), withCores = true; end
    switch Op.Domain
        case 'vertex', Ht = i_frame_vertex(Op, J, withCores);
        case 'face',   Ht = i_frame_face(Op, J, withCores);
        otherwise, error('bst_helmholtz:badDomain', 'Op.Domain must be ''vertex'' or ''face''.');
    end
end

%% ===== whole-series convenience (vertex; Prepare once, Frame per column) =====
function H = Decompose(OperatorNode, ManifoldMat, Surf, J) %#ok<DEFNU>
    Op = Prepare(OperatorNode, ManifoldMat, Surf, 'Domain', 'vertex');
    nVtot = Op.nVtot;  nT = size(J,2);
    H = struct('Curl',zeros(nVtot,nT), 'Div',zeros(nVtot,nT), ...
               'Psi',zeros(nVtot,nT),  'Phi',zeros(nVtot,nT),  'Fmag',zeros(nVtot,nT));
    H.Cores = cell(1, nT);
    H.Sources = cell(1, nT);
    for t = 1:nT
        Ht = i_frame_vertex(Op, J(:,t), true);
        H.Curl(:,t)=Ht.Curl; H.Div(:,t)=Ht.Div; H.Psi(:,t)=Ht.Psi;
        H.Phi(:,t)=Ht.Phi;   H.Fmag(:,t)=Ht.Fmag;  H.Cores{t}=Ht.Cores;  H.Sources{t}=Ht.Sources;
    end
end

%% ============================================================================
%% VERTEX DOMAIN
%% ============================================================================
function Op = i_prepare_vertex(OperatorNode, ManifoldMat, Surf)
    DiracOp = OperatorNode{1};  LBO = OperatorNode{2};
    Fcs   = double(Surf.Faces);
    nVtot = size(Surf.Vertices, 1);
    nH    = numel(DiracOp.FirstOrder.Intrinsic);
    Op = struct();
    Op.nVtot    = nVtot;
    Op.VertConn = Surf.VertConn;
    [Op.D, Op.vH, Op.Nf, Op.Wfv, Op.M, Op.cholK, Op.free, Op.totMass, ...
     Op.Gx, Op.Gy, Op.Gz, Op.NbH, Op.VtxH, Op.VnH] = deal(cell(1,nH));
    for hh = 1:nH
        D   = DiracOp.FirstOrder.Intrinsic{hh};            % [4F x 4V]
        E   = ManifoldMat.Embedded(hh);
        vH  = double(E.GlobalVertices(:));
        gf  = double(E.GlobalFaces(:));
        nVh = numel(vH);  nFh = numel(gf);
        % local face->vertex connectivity (topology from Surf.Faces, reindexed local)
        mapV  = zeros(nVtot,1);  mapV(vH) = 1:nVh;
        Floc  = mapV(Fcs(gf, :));                          % [nFh x 3] local indices
        Vloc  = E.vertex.position;                         % [nVh x 3] canonical positions (manifold)
        % canonical per-face geometry (manifold): unit normal + TRUE area (no cross/flip, no 2A/Af hazard)
        Af    = E.face.area;                               % [nFh x 1] true area
        % orient the (consistently-wound) manifold normals OUTWARD so the divergence sign
        % (div = imag . n) and source/sink charge follow the physical convention
        Nf    = i_orient_outward(E.face.normal, E.face.centroid, Af);   % [nFh x 3]
        % per-hemisphere 1-ring (LOCAL) for core detection
        eLoc = [Floc(:,[1 2]); Floc(:,[2 3]); Floc(:,[3 1])];
        Aloc = sparse([eLoc(:,1);eLoc(:,2)], [eLoc(:,2);eLoc(:,1)], true, nVh, nVh);
        nbLoc = cell(nVh,1);
        for vv = 1:nVh, nbLoc{vv} = find(Aloc(:,vv)); end
        Op.NbH{hh}  = nbLoc;
        Op.VtxH{hh} = Vloc;
        Op.VnH{hh}  = E.vertex.normal;                     % canonical vertex normals (manifold)
        % face->vertex area-weighted incidence [nVh x nFh] (row-normalized: area scale cancels)
        I = [Floc(:,1);Floc(:,2);Floc(:,3)];  Jc = [1:nFh,1:nFh,1:nFh]';
        Wfv = sparse(I, Jc, repmat(Af,3,1), nVh, nFh);
        Wfv = spdiags(1./max(sum(Wfv,2),eps),0,nVh,nVh) * Wfv;
        % per-face FEM gradient of a per-vertex scalar: grad f|_face = sum_i f_i (n x e_i^opp)/(2A)
        eO1 = Vloc(Floc(:,3),:) - Vloc(Floc(:,2),:);
        eO2 = Vloc(Floc(:,1),:) - Vloc(Floc(:,3),:);
        eO3 = Vloc(Floc(:,2),:) - Vloc(Floc(:,1),:);
        twoA = 2*Af;                                        % denominator = 2*Area (manifold true area)
        c1 = cross(Nf, eO1, 2) ./ twoA;   c2 = cross(Nf, eO2, 2) ./ twoA;   c3 = cross(Nf, eO3, 2) ./ twoA;
        grows = [(1:nFh)';(1:nFh)';(1:nFh)'];  gcols = [Floc(:,1);Floc(:,2);Floc(:,3)];
        Gx = sparse(grows, gcols, [c1(:,1);c2(:,1);c3(:,1)], nFh, nVh);
        Gy = sparse(grows, gcols, [c1(:,2);c2(:,2);c3(:,2)], nFh, nVh);
        Gz = sparse(grows, gcols, [c1(:,3);c2(:,3);c3(:,3)], nFh, nVh);
        % LBO pieces + cached Cholesky of the pinned (vertex 1) cotan stiffness via tess_cholesky
        M = LBO.Mass{hh};
        Op.D{hh}=D; Op.vH{hh}=vH; Op.Nf{hh}=Nf; Op.Wfv{hh}=Wfv; Op.M{hh}=M;
        Op.Gx{hh}=Gx; Op.Gy{hh}=Gy; Op.Gz{hh}=Gz;
        Op.cholK{hh}  = tess_cholesky(LBO, hh, 1);    % pin vertex 1; pure getter (I/O-free)
        Op.free{hh}   = Op.cholK{hh}.free;
        Op.totMass{hh}= sum(M(:));
    end
end

function Ht = i_frame_vertex(Op, Jt, withCores)
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
        omF   = q(1:4:end);                               % vorticity = w-part (per face)
        imagF = [q(2:4:end), q(3:4:end), q(4:4:end)];     % imaginary part (per face)
        divF  = sum(imagF .* Op.Nf{hh}, 2);               % divergence = imag . n_face
        omV = Op.Wfv{hh} * omF;   dvV = Op.Wfv{hh} * divF;
        psi = i_poisson(Op.cholK{hh}, Op.M{hh}, omV, Op.free{hh}, Op.totMass{hh});  % stream from vorticity
        phi = i_poisson(Op.cholK{hh}, Op.M{hh}, dvV, Op.free{hh}, Op.totMass{hh});  % potential from divergence
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
        Ht.Cores    = i_empty_cores_vertex();
        Ht.Sources  = i_empty_cores_vertex();
    end
end

%% ===== Poisson solve (vertex): mean-zero project -> cached pinned solve -> recenter =====
function psi = i_poisson(dK, M, omega, free, totMass) %#ok<INUSD>
    n = size(M,1);
    omega = omega - (sum(M*omega) / totMass) * ones(n,1);   % project to mean-zero
    x = tess_cholesky('solve', dK, M*omega);                % shared permuted Cholesky solve
    psi = x - mean(x);                                      % recenter
end

%% ===== vertex core detection: persistence-ranked extrema, per hemisphere =====
function cores = FindCoresOp(field, Op, omega) %#ok<DEFNU>
    cores = i_empty_cores_vertex();
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
    cores = i_empty_cores_vertex();
    for k = 1:numel(C.vertex)
        v = C.vertex(k);
        cores(end+1) = i_make_core(v, omega(v), C.chirality(k), C.persistence(k), ...
            C.isGlobal(k), C.birth(k), C.death(k), nan(1,3), 1); %#ok<AGROW>
    end
end

function s = i_empty_cores_vertex()
    s = struct('iVertex',{},'charge',{},'chirality',{},'omega',{},'persistence',{}, ...
               'isGlobal',{},'birth',{},'death',{},'pos',{},'hemi',{});
end

function s = i_make_core(vg, om, chirality, persistence, isGlobal, birth, death, pos, hemi)
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

%% ============================================================================
%% FACE DOMAIN
%% ============================================================================
function Op = i_prepare_face(OperatorNode, ManifoldMat, Surf)
    nFtot = size(Surf.Faces, 1);
    nH = numel(ManifoldMat.Embedded);
    Op = struct(); Op.nFtot = nFtot;
    [Op.G, Op.SkewG, Op.WF, Op.cholA, Op.freeA, Op.nFh, Op.fH, Op.vH, Op.Nf, Op.Af, Op.Cf, Op.NbF] = deal(cell(1,nH));
    for hh = 1:nH
        E   = ManifoldMat.Embedded(hh);
        vH  = double(E.GlobalVertices(:));
        fH  = double(E.GlobalFaces(:));
        nFh = numel(fH);
        % gradFace from the operator node (NOT a fresh nxr build); geometry from the manifold
        G   = OperatorNode.FaceAux{hh}.GradFace;            % [3F x F]
        Af  = E.face.area;                                 % true area
        Cf  = E.face.centroid;                             % barycentric centroid
        Nf  = i_orient_outward(E.face.normal, Cf, Af);     % outward (physical divergence sign)
        % local face->vertex connectivity (for the dual adjacency only)
        nVh = numel(vH);  mapV = zeros(max([vH;1]),1); mapV(vH) = 1:nVh;
        Floc = mapV(double(Surf.Faces(fH,:)));
        % SkewG = n_f x G : per-face block rotation times G
        SkewG = i_block_rotation(Nf, nFh) * G;
        % coupled variational Hodge normal matrix A = M' W_F M, M = [G SkewG]
        WFv = repelem(Af,3);  Mm = [G, SkewG];
        A = Mm' * (WFv .* Mm);  A = (A + A')/2;             % [2F x 2F] sym PSD
        freeA = setdiff((1:2*nFh)', [1; nFh+1]);           % pin phi(1), psi(1)
        Op.G{hh}=G; Op.SkewG{hh}=SkewG; Op.WF{hh}=Af; Op.nFh{hh}=nFh;
        Op.cholA{hh}=decomposition(A(freeA,freeA),'chol');  Op.freeA{hh}=freeA;
        Op.fH{hh}=fH; Op.vH{hh}=vH; Op.Nf{hh}=Nf; Op.Af{hh}=Af; Op.Cf{hh}=Cf;
        Op.NbF{hh}=i_dual_adjacency(Floc, nFh);
    end
end

function Ht = i_frame_face(Op, Jf, withCores)
    if size(Jf,2) ~= 3, Jf = reshape(Jf, 3, [])'; end       % accept [3F x 1] stacked
    nFtot = Op.nFtot;
    zF1 = zeros(nFtot,1);  zF3 = zeros(nFtot,3);
    Ht = struct('Curl',zF1,'Div',zF1,'Psi',zF1,'Phi',zF1,'Fmag',zF1,'Hmag',zF1, ...
                'Vtot',zF3,'Virr',zF3,'Vsol',zF3,'Vharm',zF3);
    harmNum=0; harmDen=0;
    for hh = 1:numel(Op.G)
        fH = Op.fH{hh};  Af = Op.Af{hh};  nFh = Op.nFh{hh};
        Jl = Jf(fH,:);  Jcol = reshape(Jl', [], 1);          % [3F x 1] (x,y,z per face)
        WF = repelem(Af,3);                                  % [3F] area weights
        divS  = Op.G{hh}'    * (WF .* Jcol);                 % G' W_F J   -> [F] divergence source
        curlS = Op.SkewG{hh}'* (WF .* Jcol);                 % SkewG' W_F J -> [F] curl source
        sol = i_pinned_solve(Op.cholA{hh}, [divS; curlS], Op.freeA{hh});   % coupled pinned solve + recenter
        phi = sol(1:nFh);       phi = phi - mean(phi);
        psi = sol(nFh+1:2*nFh); psi = psi - mean(psi);
        Virr = reshape(Op.G{hh}    *phi, 3, [])';            % [F x 3] irrotational
        Vsol = reshape(Op.SkewG{hh}*psi, 3, [])';            % [F x 3] solenoidal
        Vharm = Jl - Virr - Vsol;
        Ht.Curl(fH)=curlS;  Ht.Div(fH)=divS;  Ht.Psi(fH)=psi;  Ht.Phi(fH)=phi;
        Ht.Fmag(fH)=sqrt(sum(Jl.^2,2));  Ht.Hmag(fH)=sqrt(sum(Vharm.^2,2));
        Ht.Vtot(fH,:)=Jl;  Ht.Virr(fH,:)=Virr;  Ht.Vsol(fH,:)=Vsol;  Ht.Vharm(fH,:)=Vharm;
        harmNum = harmNum + sum(Af .* sum(Vharm.^2,2));
        harmDen = harmDen + sum(Af .* sum(Jl.^2,2));
    end
    Ht.HarmFrac = harmNum/max(harmDen,eps);
    if withCores
        Ht.Cores   = i_find_cores_face(Ht.Psi, Op, Ht.Curl);
        Ht.Sources = i_find_cores_face(Ht.Phi, Op, Ht.Div);
    else
        Ht.Cores = i_empty_cores_face();  Ht.Sources = i_empty_cores_face();
    end
end

function Rn = i_block_rotation(Nf, nFh)
% Block-diagonal [3F x 3F] sparse: per face the cross-product matrix [n]_x so that
% (Rn * v) restricted to face f equals n_f x v_f.
    nx=Nf(:,1); ny=Nf(:,2); nz=Nf(:,3);
    r=(1:nFh)';  b0=3*(r-1);
    I=[b0+1; b0+1; b0+2; b0+2; b0+3; b0+3];
    J=[b0+2; b0+3; b0+1; b0+3; b0+1; b0+2];
    S=[-nz;  ny;   nz;  -nx;  -ny;   nx];
    Rn = sparse(I, J, S, 3*nFh, 3*nFh);
end

function NbF = i_dual_adjacency(Floc, nFh)
% Dual face-adjacency: faces sharing an edge. Returns a cell{nFh} of neighbor face lists.
    E = [Floc(:,[1 2]); Floc(:,[2 3]); Floc(:,[3 1])];
    E = sort(E,2);  fId = repmat((1:nFh)',3,1);
    [~,~,ic] = unique(E,'rows');
    A = sparse(ic, fId, true, max(ic), nFh);
    Adj = (A'*A) > 0;  Adj(1:nFh+1:end) = false;
    NbF = cell(nFh,1);
    for f=1:nFh, NbF{f}=find(Adj(:,f)); end
end

function s = i_empty_cores_face()
    s = struct('iVertex',{},'iFace',{},'charge',{},'chirality',{},'omega',{}, ...
               'persistence',{},'isGlobal',{},'hemi',{},'pos',{});
end

function cores = i_find_cores_face(field, Op, omega)
% Persistence-ranked extrema of a per-FACE potential over the dual face-adjacency;
% charge = sign of vorticity/divergence at that face; position = face centroid.
    cores = i_empty_cores_face();
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

%% ===== orient consistently-wound face normals outward (both domains) =====
function Nf = i_orient_outward(Nf, Cf, Af)
% Globally flip the per-face normals to point OUTWARD. The manifold (nxr/geometry-central)
% normals are consistently oriented by mesh winding but the global sign is gauge-dependent;
% for a closed surface the divergence theorem gives sum_f Af * n_f . (c_f - c0) = 3*Volume,
% positive iff n_f is outward. One robust global sign per hemisphere -- no per-face VertNormals
% flip (so no VertNormals noise), restoring the physical divergence/source-sink convention.
    c0 = mean(Cf, 1);
    if sum(Af .* sum(Nf .* (Cf - c0), 2)) < 0
        Nf = -Nf;
    end
end

%% ===== shared cached-factor solver (both domains) =====
function x = i_pinned_solve(dChol, rhs, free)
% Solve a pinned SPD system with a cached Cholesky factor: x(free) = dChol\rhs(free),
% pinned entries 0, then recenter to mean-zero. Shared by the vertex dual-Poisson and the
% face coupled variational solve.
    x = zeros(size(rhs));
    x(free) = dChol \ rhs(free);
    x = x - mean(x);
end
