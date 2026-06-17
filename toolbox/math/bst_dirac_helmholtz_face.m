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
        % INTRINSIC dual face Dirac D̃_int [4V x 4F] (centroid immersion root) from a fresh
        % nxr handle on this hemisphere submesh -- consistent with the cotan/intrinsic pipeline
        h = nxr_safe_create(Vloc, Floc);
        Dt = nxr_compute('operators', h, 'diracFaceIntrinsicD');
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
