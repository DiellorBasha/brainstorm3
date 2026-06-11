function R = demo_mne_dirac_filters(base)
% DEMO_MNE_DIRAC_FILTERS  Show that every step of the (Dirac) MNE inverse is a
% SPECTRAL FILTER, and visualize each one on its natural axis:
%   - noise regularization / whitening  -> filter on the NOISE-COVARIANCE eigenvalues d
%   - leadfield SVD                      -> defines the OBSERVABILITY axis (sing. values s)
%   - SNR -> Lambda (Tikhonov/Wiener)    -> low-pass window on s
%   - amplitude min-norm gain            -> band-pass on s
%   - dSPM / sLORETA                     -> source-side resolution reshaping
% and links the observability axis to the Dirac GEOMETRY axis (eigenvalues lambda).
%
% The point: MNE = a fixed cascade of spectral windows. Seeing the shapes/formulae
% exposes the knobs (and where a custom g(lambda)/g(s) window could surpass MNE).
%
% USAGE:  R = demo_mne_dirac_filters
%
% Authors: Diellor Basha, 2026

    if nargin < 1 || isempty(base), base = 'Subject01/S01_AEF_20131218_01_notch/'; end
    OUTDIR = '/Users/diellorbasha/workspace/research/code/brainstorm3/dev/benchmarks';
    NoiseReg = 0.1; SNR0 = 3;

    % ---- inputs ----
    HMos=in_bst_headmodel([base 'headmodel_surf_os_meg.mat'],0);
    ChanMat=in_bst_channel([base 'channel_ctf_acc1.mat']); types={ChanMat.Channel.Type};
    G=double(HMos.Gain); iMEG=all(isfinite(G),2)&strcmpi(types(:),'MEG'); G=G(iMEG,:); nCh=size(G,1);
    Srf=in_tess_bst(HMos.SurfaceFile); Vtx=Srf.Vertices; Nrm=Srf.VertNormals; nV=size(Vtx,1);
    NC=load(file_fullpath([base 'noisecov_full.mat'])); Cn=NC.NoiseCov(iMEG,iMEG); Cn=(Cn+Cn')/2;
    HMf=HMos; HMf.Gain=G;
    CompHM=bst_dirac(HMf,'nModes',400,'Tau',0.5); Gm=double(CompHM.Gain); lamMode=double(CompHM.Eigenvalues);

    % ================= FILTER 1: noise regularization + whitening =================
    [Un,Dc]=eig(Cn); d=real(diag(Dc)); [d,ix]=sort(d,'descend'); Un=Un(:,ix);
    ridge = mean(d)*NoiseReg;                       % Hamalainen ridge
    alphas = [0 0.01 0.1 1];                        % regularization knob family
    Wn = Un*diag(1./sqrt(d + ridge))*Un';           % the actual whitener C_reg^{-1/2}
    R.F1.formula = 'w(d) = (d + alpha*mean(d))^{-1/2},  alpha=NoiseReg';
    R.F1.ridge = ridge; R.F1.cond = max(d)/min(d); R.F1.plateau = 1/sqrt(ridge);

    % ================= FILTER 2/3: observability SVD + SNR->Lambda =================
    Gt = Wn*Gm;                                     % whitened mode-forward
    [U,Ssvd,V] = svd(Gt,'econ'); s=diag(Ssvd); s2=s.^2; %#ok<ASGLU>
    dkMode = sqrt(sum(Gt.^2,1))';                    % per-Dirac-mode observability (links to lambda)
    Lambda = SNR0^2/mean(s2);
    fWiener = @(ss,Lam) (Lam*ss.^2)./(Lam*ss.^2+1);  % filter factor (on the solution)
    gAmp    = @(ss,Lam) (Lam*ss)./(Lam*ss.^2+1);     % gain on whitened data (amplitude)
    R.F2.formula = 's_i = svd(C_reg^{-1/2} * Gm) ;  Lambda = SNR^2 / mean(s^2)';
    R.F3.formula = 'Wiener filter factor f_i = s_i^2 / (s_i^2 + 1/Lambda)';
    R.Fgain.formula = 'amplitude gain g_i = Lambda*s_i / (Lambda*s_i^2 + 1) = f_i / s_i ; peak at s=1/sqrt(Lambda)';
    R.Lambda=Lambda; R.s=s; R.lamMode=lamMode; R.dkMode=dkMode;

    % ================= effect on coefficients (real data @ M100) =================
    DataMat=in_bst_data(local_find(base,'deviant_average')); F=double(DataMat.F(iMEG,:)); Time=DataMat.Time;
    win=Time>=0.06 & Time<=0.16; gfp=sqrt(sum((Wn*F).^2,1)); [~,tpk]=max(gfp.*win);
    dM = Wn*F(:,tpk); a = U'*dM;                      % observable data components a_i = u_i'*d~
    cNaive = a./s;                                    % unregularized inverse coeff (noise blows up)
    cReg   = gAmp(s,Lambda).*a;                       % regularized (amplitude) coeff
    R.coeff.naive=cNaive; R.coeff.reg=cReg; R.M100ms=Time(tpk)*1e3;

    % ================= FILTER 5: measure reshaping (PSF radial profile) =================
    OPT=struct('NoiseMethod','reg','NoiseReg',NoiseReg,'SnrMethod','fixed','SnrFixed',SNR0);
    OPT.NoiseCovMat.NoiseCov=Cn; OPT.ChannelTypes=types(iMEG);
    Kk=struct();
    for m={'amplitude','dspm2018','sloreta'}, OPT.InverseMeasure=m{1}; Rr=bst_inverse_dirac(HMf,OPT); Kk.(strrep(m{1},'2018',''))=Rr.ImagingKernel; end
    vt=round(0.4*nV); dsrc=G(:,(vt-1)*3+(1:3))*Nrm(vt,:)';
    dist=sqrt(sum((Vtx-Vtx(vt,:)).^2,2))*1e3; edges=0:5:80; ctr=(edges(1:end-1)+edges(2:end))/2;
    psf=struct(); meas={'amplitude','dspm','sloreta'};
    for j=1:3
        mg=sqrt(sum(reshape(Kk.(meas{j})*dsrc,3,nV).^2,1))'; mg=mg/max(mg);
        prof=arrayfun(@(b) mean(mg(dist>=edges(b)&dist<edges(b+1))), 1:numel(ctr));
        psf.(meas{j})=prof;
    end
    R.psf=psf; R.psfCtr=ctr;

    % ================= report =================
    fprintf('\n==== MNE steps as spectral filters ====\n');
    fprintf('F1 whiten+reg : %s | cond(C)=%.1e | ridge=%.2e -> gain plateau %.2e\n', R.F1.formula, R.F1.cond, ridge, R.F1.plateau);
    fprintf('F2 observ.    : %s | Lambda=%.3e\n', R.F2.formula, Lambda);
    fprintf('F3 Wiener     : %s | half-power at s^2=1/Lambda (kept DOF=%d)\n', R.F3.formula, sum(s2>1/Lambda));
    fprintf('Gain amplitude: %s | peak s=%.3g\n', R.Fgain.formula, 1/sqrt(Lambda));
    fprintf('F5 measures   : amplitude / dSPM (1/sqrt(trace K_v K_v'')) / sLORETA (R_v^{-1/2})\n');
    fprintf('coeff effect  : naive max|a/s|=%.2e vs regularized max|g*a|=%.2e (noise tamed)\n', max(abs(cNaive)), max(abs(cReg)));

    % ================= figure =================
    f=figure('Color','w','Position',[40 40 1500 900],'Visible','off');
    % A: noise regularization filter (gain family over the alpha knob)
    subplot(2,3,1);
    dd=logspace(log10(min(d)),log10(max(d)),300); cols1=lines(numel(alphas));
    for q=1:numel(alphas)
        loglog(dd, 1./sqrt(dd+alphas(q)*mean(d)),'LineWidth',2,'Color',cols1(q,:)); hold on;
    end
    loglog(d, 1./sqrt(d+ridge),'k.','MarkerSize',6);     % actual noise eigenvalues (alpha=0.1)
    xlabel('noise eigenvalue d'); ylabel('whitening gain 1/\surd(d+\alpha\langle d\rangle)'); grid on;
    legend([arrayfun(@(a)sprintf('\\alpha=%g',a),alphas,'uni',0), {'actual d'}],'Location','southwest');
    title('F1: noise reg = floor on whitening gain');
    % B: observability vs Dirac lambda
    subplot(2,3,2);
    loglog(lamMode+eps, dkMode,'.','Color',[.3 .3 .6]); grid on;
    xlabel('Dirac eigenvalue \lambda (geometry axis)'); ylabel('per-mode observability ||C^{-1/2}G\Psi_k||');
    title('B: geometry axis \lambda \rightarrow observability (loose low-pass)');
    % C: Wiener filter factor vs s, several SNR
    subplot(2,3,3);
    ss=logspace(log10(min(s)),log10(max(s)),300);
    cols=lines(4); snrs=[3 10 30 100];
    for q=1:4, Lam=snrs(q)^2/mean(s2); semilogx(ss,fWiener(ss,Lam),'LineWidth',2,'Color',cols(q,:)); hold on; end
    semilogx(ss, double(ss.^2>1/Lambda),'k--','LineWidth',1);   % TSVD limit
    xlabel('whitened singular value s (observability)'); ylabel('filter factor f(s)'); grid on; ylim([0 1.05]);
    legend([arrayfun(@(x)sprintf('SNR %d',x),snrs,'uni',0), {'TSVD'}],'Location','northwest');
    title('F3: Tikhonov low-pass  f=s^2/(s^2+1/\Lambda)');
    % D: amplitude gain on data (band-pass)
    subplot(2,3,4);
    semilogx(ss, gAmp(ss,Lambda)/max(gAmp(ss,Lambda)),'LineWidth',2,'Color',[.2 .5 .2]); hold on;
    xline(1/sqrt(Lambda),'k:'); xlabel('s'); ylabel('gain g(s) (norm.)'); grid on;
    title('Amplitude gain g=\Lambda s/(\Lambda s^2+1): BAND-PASS');
    % E: effect on coefficients
    subplot(2,3,5);
    semilogy(abs(cNaive),'r','LineWidth',1.2); hold on; semilogy(abs(cReg),'b','LineWidth',2);
    xlabel('observable direction index'); ylabel('|coefficient|'); grid on;
    legend('naive a/s (noise blows up)','regularized g\cdota','Location','northwest');
    title(sprintf('E: effect on coeffs (real data @ %.0f ms)', R.M100ms));
    % F: measure PSF radial profiles
    subplot(2,3,6);
    plot(ctr,psf.amplitude,'LineWidth',2); hold on; plot(ctr,psf.dspm,'LineWidth',2); plot(ctr,psf.sloreta,'LineWidth',2);
    xlabel('distance from source (mm)'); ylabel('PSF (norm.)'); grid on;
    legend('amplitude','dSPM','sLORETA'); title('F5: measures reshape the resolution (PSF)');
    sgtitle('MNE \leftrightarrow Dirac: every inverse step is a spectral filter','FontWeight','bold');
    print(f,[OUTDIR '/demo_mne_dirac_filters.png'],'-dpng','-r110'); close(f);
    fprintf('Saved %s/demo_mne_dirac_filters.png\n', OUTDIR);
end

function fn = local_find(base, tag)
    sStudy=bst_get('StudyWithCondition', base(1:end-1)); fn='';
    for i=1:numel(sStudy.Data), if contains(sStudy.Data(i).FileName,tag), fn=sStudy.Data(i).FileName; return; end; end
    fn=sStudy.Data(1).FileName;
end
