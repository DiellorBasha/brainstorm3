function R = bench_dirac_matern_prior(base)
% BENCH_DIRAC_MATERN_PRIOR  EXPERIMENTAL: does a curvature-aware Matern/smoothness
% source prior on the Dirac modes beat the standard inverses?
%
% A source prior is a window h(lambda) on the Dirac eigenvalues = the assumed
% spatial power spectrum of cortical activity (R = diag(h(lambda)) in the mode
% basis). Standard MNE uses h=const (white). LORETA uses h=1/lambda^2 (smoothness).
% Matern uses h=(lambda+kappa^2)^(-alpha) (tunable smoothness/correlation length),
% measured by the CURVATURE-AWARE Dirac operator.
%
% Validation protocol (experimental method must first reproduce the gold standard):
%   (A) sanity: Dirac-white prior reproduces vertex MNE (the established match)
%   (B) simulate ground-truth focal sources (+noise); measure peak localization
%       error and spatial spread for each method, with statistics over sources.
% Methods compared:
%   vertex (bst_inverse_linear_2018): MNE, dSPM, sLORETA
%   Dirac (this script, amplitude):   white(=MNE), LORETA(1/lambda^2), Matern
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
    K.vMNE   = vkern(HMm,OPTv,'amplitude');
    K.vdSPM  = vkern(HMm,OPTv,'dspm2018');
    K.vsLOR  = vkern(HMm,OPTv,'sloreta');

    % ---- Dirac prior kernels (amplitude, different h(lambda)) ----
    med=median(lam(lam>0)); eps0=1e-3*max(lam);
    hWhite =@(L) ones(size(L));
    hLOR   =@(L) 1./(L+eps0).^2;                 % classic LORETA smoothness prior
    hMat   =@(L,k2,a) (L+k2).^(-a);              % Matern (SPDE) spectral prior
    K.dWhite = dkern(CompHM,W,Gm,lam,hWhite(lam),SNR);
    K.dLOR   = dkern(CompHM,W,Gm,lam,hLOR(lam),SNR);
    kappa2=[0.3 1 3]*med; alpha=2; matKeys={}; matLabels={};
    for q=1:numel(kappa2)
        kk=sprintf('dMat%d',q); matKeys{end+1}=kk; %#ok<AGROW>
        matLabels{end+1}=sprintf('Matern \\kappa^2=%.2g',kappa2(q)); %#ok<AGROW>
        K.(kk)=dkern(CompHM,W,Gm,lam,hMat(lam,kappa2(q),alpha),SNR);
    end

    % ---- (A) sanity: Dirac-white reproduces vertex MNE ----
    mags=@(J)sqrt(sum(reshape(J,3,[]).^2,1))';
    vt0=round(0.4*nV); d0=G(:,(vt0-1)*3+(1:3))*Nrm(vt0,:)';
    R.sanity_corr = corr(mags(K.dWhite*d0), mags(K.vMNE*d0));
    fprintf('(A) sanity: Dirac-white vs vertex-MNE source-map corr = %.3f (should be ~1)\n', R.sanity_corr);

    % ---- (B) simulate ground-truth sources, measure localization ----
    methods = [{'vMNE','vdSPM','vsLOR','dWhite','dLOR'}, matKeys];
    labels  = [{'MNE','dSPM','sLORETA','Dirac-white','Dirac-LORETA'}, matLabels];
    verts=round(linspace(300,nV-300,60)); nM=numel(methods);
    errPk=zeros(numel(verts),nM); spread=zeros(numel(verts),nM);
    [Uc,Dc]=eig(Cn); dc=max(real(diag(Dc)),0);
    for i=1:numel(verts)
        vt=verts(i); o=Nrm(vt,:)'; d=G(:,(vt-1)*3+(1:3))*o;
        e=Uc*(sqrt(dc).*cos((1:nCh)'*vt)); e=e/norm(e)*norm(d)/SNR; d=d+e;   % reproducible pseudo-noise
        dist=sqrt(sum((Vtx-Vtx(vt,:)).^2,2))*1e3;
        for j=1:nM
            m=mags(K.(methods{j})*d); [~,pk]=max(m);
            errPk(i,j)=norm(Vtx(pk,:)-Vtx(vt,:))*1e3;
            spread(i,j)=sum(m.*dist)/sum(m);                 % magnitude-weighted spatial spread
        end
    end
    R.methods=methods; R.labels=labels; R.errPk=errPk; R.spread=spread;
    fprintf('\n(B) %d simulated sources (SNR=%g):  peak error (mm)   |  spread (mm)\n', numel(verts), SNR);
    for j=1:nM
        fprintf('   %-18s  med=%5.1f IQR[%4.1f %5.1f]  | %5.1f\n', labels{j}, ...
            median(errPk(:,j)), prctile(errPk(:,j),25), prctile(errPk(:,j),75), median(spread(:,j)));
    end
    [~,best]=min(median(errPk,1)); fprintf('   --> best median peak error: %s\n', labels{best});

    % ---- figure ----
    f=figure('Color','w','Position',[40 40 1400 560],'Visible','off');
    subplot(1,2,1); boxplot(errPk,'Labels',labels); ylabel('peak localization error (mm)');
    title(sprintf('Localization (%d sources, SNR=%g)',numel(verts),SNR)); grid on; set(gca,'XTickLabelRotation',30);
    subplot(1,2,2);
    bar([median(errPk,1); median(spread,1)]'); set(gca,'XTickLabel',labels,'XTickLabelRotation',30);
    ylabel('mm'); legend('median peak err','median spread'); grid on; title('Summary: error & spread by method');
    sgtitle('EXPERIMENTAL: curvature-aware Matern/smoothness prior vs MNE/dSPM/sLORETA','FontWeight','bold');
    print(f,[OUTDIR '/bench_dirac_matern_prior.png'],'-dpng','-r110'); close(f);
    fprintf('Saved %s/bench_dirac_matern_prior.png\n', OUTDIR);
end

% ---- vertex kernel via bst_inverse_linear_2018 ----
function Kv = vkern(HMm, OPTbase, measure)
    OPT=OPTbase; OPT.InverseMeasure=measure; Rr=bst_inverse_linear_2018(HMm,OPT); Kv=Rr.ImagingKernel;
end

% ---- Dirac weighted-prior amplitude kernel: R = diag(h(lambda)) ----
function Kd = dkern(CompHM, W, Gm, lam, hvec, SNR)
    Dh = sqrt(max(hvec(:),0));                    % [2K x 1] prior sqrt
    Dh = Dh / max(Dh);                            % scale-free (overall scale absorbed by SNR->Lambda)
    Gbar = W * (Gm .* Dh.');                       % whitened, prior-weighted mode-forward [nCh x 2K]
    [U,Ss,V]=svd(Gbar,'econ'); s=diag(Ss); s2=s.^2;
    Lam = SNR^2/mean(s2); g=(Lam*s)./(Lam*s2+1);
    Kc = Dh .* (V*(g.*(U'*W)));                    % mode-coeff kernel [2K x nCh] (c_hat = R^{1/2} x_hat)
    Kd = bst_dirac(CompHM,'Reconstruct', Kc.').';  % vertex kernel [3nVert x nCh]
end
