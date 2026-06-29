function R = proto_fiber_connectome(trkFile, SurfaceFile, Opts)
% PROTO_FIBER_CONNECTOME: build a FIBER-RESOLUTION (vertex-level) structural connectome operator on a
% cortical surface from a tractography .trk, plus its low-frequency eigenbasis. Unlike the ROI-binned
% connectome (which averages every fiber under a parcel together), this keeps the fiber endpoints at
% vertex resolution and bridges the sulcal-disconnection problem (raw endpoints leave ~46% of vertices
% isolated) with a continuous-kernel mesh smoothing W = S^p Wraw S^p' (S = row-normalized mesh
% averaging). The eigenbasis supports network heat kernels exp(-t*L) (low-pass) and spectral graph
% wavelets g(s*lambda) (band-pass, e.g. Mexican hat s*lambda*exp(-s*lambda)) seeded at any region -
% the proper foundation for sampling the connectome by scale (long/medium/short-range) rather than
% ROI-binning, and for differential operators (gradient/divergence/Hodge) over the network.
%
% Pipeline: parse .trk endpoints -> MNI (DSI ICBM152 2mm affine) -> subject SCS -> nearest surface
% vertex -> raw vertex connectome -> mesh-smooth -> largest connected component -> symmetric normalized
% Laplacian -> low modes via the normalized-adjacency trick (eigs 'largestreal', lambda_L = 1 -
% lambda_N; no factorization, ~5s for 250 modes on 15k vertices).
%
% USAGE: R = proto_fiber_connectome('.../whole_brain.trk', icbm152_cortex_mid_low, ...
%               struct('nModes',250,'smoothHops',3,'seedRegion','entorhinal L'))
%   R: .Wsm (smoothed connectome) .keep (largest-CC vertices) .Phi [nk x nModes] .lam [nModes]
%      .V .Faces .seedDemoPng
%
% Validated (29Jun2026): whole_brain.trk (105221 streamlines) on ICBM152 cortex_mid_low (15002 vtx) ->
% smoothing connects 99% (was 46% isolated), 250 modes lambda in [0 0.925], Fiedler 0.13. Heat low-pass
% spreads entorhinal to its bilateral ventral-temporal targets and coarsens; Mexican-hat band-pass
% isolates long-range (coarse) vs local (fine) network with center-surround.
%
% Author: Diellor Basha, 2026 (prototype; ICBM152 template tractography)
    if nargin<3, Opts=struct(); end
    if ~isfield(Opts,'nModes'),     Opts.nModes=250;  end
    if ~isfield(Opts,'smoothHops'), Opts.smoothHops=3; end
    if ~isfield(Opts,'edgeThr'),    Opts.edgeThr=1e-4; end
    here=bst_fileparts(mfilename('fullpath'));
    sSubj=bst_get('SurfaceFile', SurfaceFile); sMri=in_mri_bst(sSubj.Anatomy(sSubj.iAnatomy).FileName);
    sIc=in_tess_bst(SurfaceFile); V=sIc.Vertices; Fa=sIc.Faces; nV=size(V,1);

    % --- fiber endpoints -> nearest vertices ---
    EP=local_read_trk_endpoints(trkFile);
    toMNI=@(P)[-P(:,1)+79.5, -P(:,2)+81.5, P(:,3)-72];      % DSI ICBM152 2mm voxel(*2=mm) -> MNI mm
    scs1=cs_convert(sMri,'mni','scs',toMNI(EP(:,1:3))/1000);
    scs2=cs_convert(sMri,'mni','scs',toMNI(EP(:,4:6))/1000);
    v1=dsearchn(V,scs1); v2=dsearchn(V,scs2);
    Wraw=sparse([v1;v2],[v2;v1],1,nV,nV); Wraw=Wraw-spdiags(diag(Wraw),0,nV,nV);

    % --- continuous-kernel mesh smoothing (bridge sulcal disconnection) ---
    ii=[Fa(:,1);Fa(:,2);Fa(:,3)]; jj=[Fa(:,2);Fa(:,3);Fa(:,1)];
    A=double(sparse([ii;jj],[jj;ii],1,nV,nV)>0);
    Sm=A+speye(nV); Sm=spdiags(1./sum(Sm,2),0,nV,nV)*Sm; Sp=Sm^Opts.smoothHops;
    Wsm=Sp*Wraw*Sp'; Wsm=(Wsm+Wsm')/2; Wsm=Wsm-spdiags(diag(Wsm),0,nV,nV);
    Wsm(Wsm<Opts.edgeThr*max(Wsm(:)))=0;

    % --- largest connected component + symmetric normalized Laplacian ---
    cc=conncomp(graph(Wsm>0)); tb=tabulate(cc); [~,big]=max(tb(:,2)); keep=find(cc==big); nk=numel(keep);
    Wk=Wsm(keep,keep); d=full(sum(Wk,2)); Dm12=spdiags(1./sqrt(d),0,nk,nk); Nrm=Dm12*Wk*Dm12; Nrm=(Nrm+Nrm')/2;
    o2.maxit=400; [Phi,Ln]=eigs(Nrm,Opts.nModes,'largestreal',o2);   % largest of N = lowest of L
    lam=1-real(diag(Ln)); [lam,o]=sort(lam); Phi=real(Phi(:,o));

    R=struct('Wsm',Wsm,'keep',keep,'Phi',Phi,'lam',lam,'V',V,'Faces',Fa,'seedDemoPng','');
    fprintf('fiber connectome: %d vtx (largest CC of %d), %d modes, lambda in [%.3f %.3f]\n', nk, nV, Opts.nModes, lam(1), lam(end));
    if isfield(Opts,'seedRegion') && ~isempty(Opts.seedRegion)
        R.seedDemoPng=local_seed_demo(sIc, keep, Phi, lam, Opts.seedRegion, here);
        fprintf('seed demo -> %s\n', R.seedDemoPng);
    end
end

function EP=local_read_trk_endpoints(trk)
    fid=fopen(trk,'r','l'); fseek(fid,36,'bof'); nSc=fread(fid,1,'int16');
    fseek(fid,238,'bof'); nPr=fread(fid,1,'int16'); fseek(fid,988,'bof'); nT=fread(fid,1,'int32');
    fseek(fid,1000,'bof'); EP=zeros(nT,6);
    for i=1:nT
        m=fread(fid,1,'int32'); pts=reshape(fread(fid,(3+nSc)*m,'float32'),3+nSc,m)';
        if nPr>0, fread(fid,nPr,'float32'); end
        EP(i,:)=[pts(1,1:3) pts(end,1:3)];
    end
    fclose(fid);
end

function png=local_seed_demo(sIc, keep, Phi, lam, seedRegion, here)
    V=sIc.Vertices; Fa=sIc.Faces; nV=size(V,1); nk=numel(keep);
    ai=find(strcmp({sIc.Atlas.Name},'Desikan-Killiany'),1); Sc=sIc.Atlas(ai).Scouts;
    ent=Sc(find(strcmp({Sc.Label},seedRegion),1)).Vertices;
    f2k=zeros(nV,1); f2k(keep)=1:nk; seedK=f2k(ent); seedK=seedK(seedK>0);
    delta=zeros(nk,1); delta(seedK)=1/numel(seedK); phiSeed=Phi'*delta;
    cmapH=hot(64); cmapD=[ [linspace(0,1,32)';ones(32,1)],[linspace(0,1,32)';linspace(1,0,32)'],[ones(32,1);linspace(1,0,32)'] ];
    ts=[2 8 30 100]; ss=[20 5 2 1.2];                          % heat times / wavelet scales (coarse->fine)
    f=figure('Visible','off','Position',[10 10 1550 720]);
    for i=1:4
        ht=Phi*(exp(-ts(i)*lam).*phiSeed); uH=zeros(nV,1); uH(keep)=ht; cn=ht; cn(seedK)=0;
        ax=subplot(2,4,i); patch('Faces',Fa,'Vertices',V,'FaceVertexCData',uH,'FaceColor','interp','EdgeColor','none');
        clim([0 max(cn)]); colormap(ax,cmapH); view([0 -90]); axis equal off vis3d; camlight headlight; lighting gouraud; material dull;
        title(sprintf('HEAT  t=%g',ts(i)),'FontSize',10);
        g=ss(i)*lam.*exp(-ss(i)*lam); ps=Phi*(g.*phiSeed); uM=zeros(nV,1); uM(keep)=ps; mn=ps; mn(seedK)=0; m=max(abs(mn));
        ax=subplot(2,4,4+i); patch('Faces',Fa,'Vertices',V,'FaceVertexCData',uM,'FaceColor','interp','EdgeColor','none');
        clim([-m m]); colormap(ax,cmapD); view([0 -90]); axis equal off vis3d; camlight headlight; lighting gouraud; material dull;
        title(sprintf('MEXHAT  s=%g (band~%.2f)',ss(i),1/ss(i)),'FontSize',10);
    end
    png=fullfile(here,'fiber_connectome_entorhinal.png'); print(f,png,'-dpng','-r110'); close(f);
end
