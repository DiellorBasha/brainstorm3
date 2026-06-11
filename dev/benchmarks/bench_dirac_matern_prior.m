function R = bench_dirac_matern_prior(base)
% BENCH_DIRAC_MATERN_PRIOR  EXPERIMENTAL: does a curvature-aware Matern/smoothness
% source prior on the Dirac modes beat the standard inverses?
%
% A source prior is a window h(lambda) on the Dirac eigenvalues = the assumed
% spatial power spectrum of cortical activity (R = diag(h(lambda)) in the mode
% basis). MNE: h=const (white). LORETA: h=1/lambda^2 (smoothness). Matern (SPDE):
% h=(lambda+kappa^2)^(-alpha), measured by the CURVATURE-AWARE Dirac operator.
%
% Two regimes (a smoothness prior is EXPECTED to help extended, not focal, sources):
%   TEST 1  FOCAL sources    -> peak localization error  (focal is the wrong regime
%                               for smoothness; included for completeness)
%   TEST 2  EXTENDED patches  -> patch recovery (Dice overlap, centroid, area) --
%                               smoothness's home turf
% Methods: vertex MNE / dSPM / sLORETA, and Dirac amplitude Matern + Dirac
% Matern+dSPM (prior AND normalization, the fair competitor to sLORETA).
% Sanity: Dirac-white prior reproduces vertex MNE.
%
% USAGE:  R = bench_dirac_matern_prior
%
% Authors: Diellor Basha, 2026

    if nargin<1 || isempty(base), base='Subject01/S01_AEF_20131218_01_notch/'; end
    OUTDIR='/Users/diellorbasha/workspace/research/code/brainstorm3/dev/benchmarks';
    SNR=3; NoiseReg=0.1; fprintf('\n=== EXPERIMENTAL: Dirac Matern/smoothness prior vs MNE/dSPM/sLORETA ===\n');

    % ---- setup ----
    HMos=in_bst_headmodel([base 'headmodel_surf_os_meg.mat'],0);
    ChanMat=in_bst_channel([base 'channel_ctf_acc1.mat']); types={ChanMat.Channel.Type};
    G=double(HMos.Gain); iMEG=all(isfinite(G),2)&strcmpi(types(:),'MEG'); G=G(iMEG,:); nCh=size(G,1);
    Srf=in_tess_bst(HMos.SurfaceFile); Vtx=Srf.Vertices; Nrm=Srf.VertNormals; nV=size(Vtx,1);
    VertConn=Srf.VertConn;
    NC=load(file_fullpath([base 'noisecov_full.mat'])); Cn=NC.NoiseCov(iMEG,iMEG); Cn=(Cn+Cn')/2;
    HMf=HMos; HMf.Gain=G;
    CompHM=bst_dirac(HMf,'nModes',400,'Tau',0.5); Gm=double(CompHM.Gain); lam=double(CompHM.Eigenvalues);

    % ---- whitener (shared) ----
    OPTd=struct('NoiseMethod','reg','NoiseReg',NoiseReg,'SnrMethod','fixed','SnrFixed',SNR,'InverseMeasure','amplitude');
    OPTd.NoiseCovMat.NoiseCov=Cn; OPTd.ChannelTypes=types(iMEG);
    Rdi=bst_inverse_dirac(HMf,OPTd); W=Rdi.Whitener;

    % ---- vertex gold-standard kernels ----
    HMm=struct('Gain',G,'GridLoc',Vtx,'GridOrient',Nrm,'GridAtlas',[],'HeadModelType','surface', ...
               'SurfaceFile',HMos.SurfaceFile,'MEGMethod','os_meg','EEGMethod','','ECOGMethod','','SEEGMethod','');
    OPTv=struct('InverseMethod','minnorm','SourceOrient',{{'free'}},'DataTypes',{{'MEG'}}, ...
                'NoiseMethod','reg','NoiseReg',NoiseReg,'SnrMethod','fixed','SnrFixed',SNR,'UseDepth',0,'ComputeKernel',1,'DisplayMessages',0);
    OPTv.NoiseCovMat.NoiseCov=Cn; OPTv.ChannelTypes=types(iMEG);
    K.vMNE  = vkern(HMm,OPTv,'amplitude');
    K.vdSPM = vkern(HMm,OPTv,'dspm2018');
    K.vsLOR = vkern(HMm,OPTv,'sloreta');

    % ---- Dirac prior kernels ----
    med=median(lam(lam>0)); kappa2=med; alpha=2;             % moderate Matern smoothness
    hWhite=@(L) ones(size(L)); hMat=@(L)(L+kappa2).^(-alpha);
    K.dWhite   = dkern(CompHM,W,Gm,lam,hWhite(lam),SNR,'amplitude');   % sanity (=MNE)
    K.dMatAmp  = dkern(CompHM,W,Gm,lam,hMat(lam),  SNR,'amplitude');   % Matern, amplitude
    K.dMatDSPM = dkern(CompHM,W,Gm,lam,hMat(lam),  SNR,'dspm2018');    % Matern + dSPM normalization

    mags=@(J)sqrt(sum(reshape(J,3,[]).^2,1))';
    vt0=round(0.4*nV); d0=G(:,(vt0-1)*3+(1:3))*Nrm(vt0,:)';
    R.sanity_corr=corr(mags(K.dWhite*d0),mags(K.vMNE*d0));
    fprintf('SANITY: Dirac-white vs vertex-MNE map corr = %.3f (basis-truncation gap from 1)\n', R.sanity_corr);

    methods={'vMNE','vdSPM','vsLOR','dMatAmp','dMatDSPM'};
    labels ={'MNE','dSPM','sLORETA','Matern(amp)','Matern+dSPM'};
    [Uc,Dc]=eig(Cn); dc=max(real(diag(Dc)),0);

    % =================== TEST 1: FOCAL sources (peak error) ===================
    vfoc=round(linspace(300,nV-300,60)); nM=numel(methods);
    errPk=zeros(numel(vfoc),nM);
    for i=1:numel(vfoc)
        vt=vfoc(i); d=G(:,(vt-1)*3+(1:3))*Nrm(vt,:)';
        d=d+pnoise(Uc,dc,nCh,vt,norm(d),SNR);
        for j=1:nM, [~,pk]=max(mags(K.(methods{j})*d)); errPk(i,j)=norm(Vtx(pk,:)-Vtx(vt,:))*1e3; end
    end
    R.errPk=errPk;
    fprintf('\nTEST 1 FOCAL (%d srcs): peak error mm (median [IQR])\n', numel(vfoc));
    for j=1:nM, fprintf('   %-12s %5.1f [%4.1f %5.1f]\n', labels{j}, median(errPk(:,j)), prctile(errPk(:,j),25), prctile(errPk(:,j),75)); end

    % =================== TEST 2: EXTENDED patches (Dice/centroid/area) ===================
    seeds=round(linspace(400,nV-400,30)); nP=numel(seeds);
    dice=zeros(nP,nM); cerr=zeros(nP,nM); area=zeros(nP,nM);
    for i=1:nP
        patch=grow_patch(VertConn, seeds(i), 5);            % ~1.5 cm^2 geodesic patch
        Jtrue=zeros(3*nV,1);
        Jtrue((patch-1)*3+1)=Nrm(patch,1); Jtrue((patch-1)*3+2)=Nrm(patch,2); Jtrue((patch-1)*3+3)=Nrm(patch,3);
        d=G*Jtrue; d=d+pnoise(Uc,dc,nCh,seeds(i),norm(d),SNR);
        trueC=mean(Vtx(patch,:),1); trueSet=false(nV,1); trueSet(patch)=true;
        for j=1:nM
            m=mags(K.(methods{j})*d); estSet=m>0.5*max(m);
            dice(i,j)=2*sum(estSet&trueSet)/(sum(estSet)+sum(trueSet)+eps);
            wc=sum(Vtx.*m,1)/sum(m); cerr(i,j)=norm(wc-trueC)*1e3;
            area(i,j)=sum(estSet)/numel(patch);
        end
    end
    R.methods=methods; R.labels=labels; R.dice=dice; R.cerr=cerr; R.area=area;
    fprintf('\nTEST 2 EXTENDED (%d patches, ~%d verts each): Dice | centroid err mm | area ratio\n', nP, numel(grow_patch(VertConn,seeds(1),5)));
    for j=1:nM
        fprintf('   %-12s Dice=%.2f [%.2f %.2f] | cErr=%4.1f | area=%4.1f\n', labels{j}, ...
            median(dice(:,j)), prctile(dice(:,j),25), prctile(dice(:,j),75), median(cerr(:,j)), median(area(:,j)));
    end
    [~,bD]=max(median(dice,1)); fprintf('   --> best patch Dice: %s\n', labels{bD});

    % ---- figure ----
    f=figure('Color','w','Position',[40 40 1500 560],'Visible','off');
    subplot(1,3,1); boxplot(errPk,'Labels',labels); ylabel('peak error (mm)');
    title('TEST 1: focal localization'); grid on; set(gca,'XTickLabelRotation',25);
    subplot(1,3,2); boxplot(dice,'Labels',labels); ylabel('patch Dice overlap');
    title('TEST 2: extended-patch recovery'); grid on; set(gca,'XTickLabelRotation',25); ylim([0 1]);
    subplot(1,3,3); boxplot(cerr,'Labels',labels); ylabel('patch centroid error (mm)');
    title('TEST 2: extended centroid'); grid on; set(gca,'XTickLabelRotation',25);
    sgtitle('EXPERIMENTAL: Matern/smoothness prior -- focal vs extended sources','FontWeight','bold');
    print(f,[OUTDIR '/bench_dirac_matern_prior.png'],'-dpng','-r110'); close(f);
    fprintf('Saved %s/bench_dirac_matern_prior.png\n', OUTDIR);
end

% ---- vertex kernel via bst_inverse_linear_2018 ----
function Kv = vkern(HMm, OPTbase, measure)
    OPT=OPTbase; OPT.InverseMeasure=measure; Rr=bst_inverse_linear_2018(HMm,OPT); Kv=Rr.ImagingKernel;
end

% ---- Dirac weighted-prior kernel (measure: amplitude/dspm2018) ----
function Kd = dkern(CompHM, W, Gm, lam, hvec, SNR, measure) %#ok<INUSL>
    Dh=sqrt(max(hvec(:),0)); Dh=Dh/max(Dh);
    Gbar=W*(Gm.*Dh.'); [U,Ss,V]=svd(Gbar,'econ'); s=diag(Ss); s2=s.^2;
    Lam=SNR^2/mean(s2); g=(Lam*s)./(Lam*s2+1);
    KcW = Dh.*(V*(g.*(U')));                                   % [2K x nCh] whitened-data
    KvtxW = bst_dirac(CompHM,'Reconstruct', KcW.').';          % [3nVert x nCh] whitened-data
    switch lower(measure)
        case 'amplitude', Knorm=KvtxW;
        case 'dspm2018'
            rowN2=sum(KvtxW.^2,2); Sn=sqrt(sum(reshape(rowN2,3,[]),1))';
            Sn3=reshape(repmat(Sn',3,1),[],1); Sn3(Sn3<eps)=eps; Knorm=KvtxW./Sn3;
    end
    Kd = Knorm * W;                                            % raw-data kernel
end

% ---- grow a geodesic patch by BFS on the vertex adjacency ----
function patch = grow_patch(VertConn, seed, nRings)
    patch=seed; frontier=seed;
    for r=1:nRings
        nbrs=find(any(VertConn(frontier,:),1)); frontier=setdiff(nbrs,patch); patch=union(patch,frontier);
    end
    patch=patch(:);
end

% ---- reproducible pseudo-noise at a target SNR ----
function e = pnoise(Uc, dc, nCh, seed, sig, SNR)
    e=Uc*(sqrt(dc).*cos((1:nCh)'*seed)); e=e/norm(e)*sig/SNR;
end
