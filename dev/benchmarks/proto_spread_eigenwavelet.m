function R = proto_spread_eigenwavelet(SurfaceFile, A, Opts)
% PROTO_SPREAD_EIGENWAVELET: Layer-1 heat-eigenwavelet dynamic analysis of a longitudinal surface
% field series A = [nVert x nT]. The heat eigenwavelet is the LBO diffusion scale-space (the
% linearized reaction-diffusion core), so its bands are diffusion scales (coarse -> fine). Tracks:
%   - per-band energy E(m,t)         -> which scales carry/grow the signal
%   - scale centroid scaleC(t)       -> focal (fine) <-> diffuse (coarse) progression
%   - total energy etot(t)           -> global accumulation (validate vs Centiloid)
%   - epicenter epi(t)               -> peak of the coarse-band energy (the amyloid center)
%
% USAGE: R = proto_spread_eigenwavelet(SurfaceFile, A, Opts)   Opts: .Nf=6 .nModes=300
%
% Author: Diellor Basha, 2026 (prototype)
    if nargin<3, Opts=struct(); end
    if ~isfield(Opts,'Nf'), Opts.Nf=6; end
    if ~isfield(Opts,'nModes'), Opts.nModes=300; end
    Eig=tess_eigen(SurfaceFile,'Laplace-Beltrami','nModes',Opts.nModes);
    Op =tess_operators(SurfaceFile,'Laplace-Beltrami');
    lmax=max(cellfun(@(L)max(L),Eig.Lambda)); lmin=min(cellfun(@(L)min(L(L>1e-10)),Eig.Lambda));
    frame=bst_eigenwavelet('Design','heat',Opts.Nf,[lmin lmax]);
    sW=in_tess_bst(SurfaceFile); V=sW.Vertices; nV=size(V,1); nT=size(A,2);

    % decompose each timepoint -> per-vertex per-band coefficients
    a1=A(:,1); a1(~isfinite(a1))=0;
    W1=bst_eigenwavelet('Analysis', a1, Eig, Op, frame); nB=size(W1,3);
    Wv=zeros(nV,nB,nT); E=zeros(nB,nT);
    for t=1:nT
        a=A(:,t); a(~isfinite(a))=0;
        Wt=bst_eigenwavelet('Analysis', a, Eig, Op, frame);
        for m=1:nB, Wv(:,m,t)=Wt(:,1,m); E(m,t)=sum(Wt(:,1,m).^2); end
    end
    etot=sum(E,1);
    scaleC=sum((1:nB)'.*E,1)./max(sum(E,1),eps);          % 1=coarsest .. nB=finest
    coarseB=1:max(1,round(nB/3));                          % coarse bands = the diffuse center
    epi=zeros(1,nT);
    for t=1:nT, Ev=sum(Wv(:,coarseB,t).^2,2); [~,epi(t)]=max(Ev); end
    R=struct('E',E,'etot',etot,'scaleC',scaleC,'epi',epi,'nB',nB,'V',V,'Wv',Wv,'frame',frame);
end
