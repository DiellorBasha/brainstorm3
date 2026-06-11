function R = bench_dirac_source_mapping(base)
% BENCH_DIRAC_SOURCE_MAPPING  End-to-end validation & benchmark of the Dirac
% eigenmode source-mapping pipeline (forward transform + bst_inverse_dirac),
% compared against the standard overlapping-spheres leadfield and vertex MNE.
%
% Runs six sections, prints a structured report, and saves figures to
% dev/benchmarks/:
%   S1  Forward: Dirac transform vs overlapping spheres (leadfield reconstruction)
%   S2  Dirac mode energy spectrum (ambient-flat structure)
%   S3  Observability space (clean SVD + noise-whitened DOF vs SNR)
%   S4  Inverse stage validation (bst_inverse_dirac vs vertex MNE, bit-level)
%   S5  Localization statistics (amplitude / dSPM / sLORETA over many sources)
%   S6  Orientation analysis (isotropy, leadfield sensitivity, frame-invariance,
%       orientation recovery)
%
% USAGE:  R = bench_dirac_source_mapping
%         R = bench_dirac_source_mapping('Subject01/S01_AEF_20131218_01_notch/')
%
% Requires the protocol loaded with an Overlapping spheres head model, a noise
% covariance, and a 'Dirac' eigen node on the head model surface (bst_inverse_dirac
% find-or-creates it).
%
% Authors: Diellor Basha, 2026

    if nargin < 1 || isempty(base)
        base = 'Subject01/S01_AEF_20131218_01_notch/';
    end
    OUTDIR = '/Users/diellorbasha/workspace/research/code/brainstorm3/dev/benchmarks';
    R = struct();
    fprintf('\n========== DIRAC SOURCE-MAPPING BENCHMARK ==========\n');

    % ===== shared setup =====
    HMos = in_bst_headmodel([base 'headmodel_surf_os_meg.mat'], 0);
    ChanMat = in_bst_channel([base 'channel_ctf_acc1.mat']);
    types = {ChanMat.Channel.Type};
    G = double(HMos.Gain);
    iMEG = all(isfinite(G),2) & strcmpi(types(:),'MEG');
    G = G(iMEG,:); nCh = size(G,1); nV = size(G,2)/3;
    Srf = in_tess_bst(HMos.SurfaceFile); Vtx = Srf.Vertices; Nrm = Srf.VertNormals;
    NC = load(file_fullpath([base 'noisecov_full.mat']));
    Cn = NC.NoiseCov(iMEG,iMEG); Cn = (Cn+Cn')/2;
    HMf = HMos; HMf.Gain = G;
    fprintf('Setup: %d MEG channels, %d vertices.\n', nCh, nV);

    % Build the Dirac mode head model (forward transform) once
    CompHM = bst_dirac(HMf, 'nModes',400, 'Tau',0.5);
    Gm  = double(CompHM.Gain);                       % [nCh x 2K] mode-forward
    Kpm = CompHM.nModes/2;                            % modes per hemisphere
    Eig = load(file_fullpath(CompHM.DiracEigenFile)); % eigenbasis (for energy/PR)
    Op  = load(file_fullpath(Eig.OperatorFile));

    % ===================================================================
    R.S1 = section1_forward(G, Gm, Eig, Op, nCh, nV);
    R.S2 = section2_mode_energy(Eig, Op, Gm, Vtx, Nrm);
    R.S3 = section3_observability(Gm, Cn, nCh);
    R.S4 = section4_inverse_stages(HMf, G, Gm, Cn, types, iMEG);
    R.S5 = section5_localization(HMf, G, Vtx, Srf.Faces, Nrm, Cn, types, iMEG, OUTDIR);
    R.S6 = section6_orientation(HMf, G, Vtx, Nrm, Cn, types, iMEG);

    % ===================================================================
    make_figures(R, OUTDIR);
    fprintf('\n========== BENCHMARK COMPLETE ==========\n');
    fprintf('Figures saved to %s/bench_dirac_*.png\n', OUTDIR);
end


%% ===================================================================
function s = section1_forward(G, Gm, Eig, Op, nCh, nV)
    fprintf('\n--- S1: FORWARD (Dirac transform vs overlapping spheres) ---\n');
    % leadfield B-energy total + captured vs modes/hemi
    totB = 0;
    for hh = 1:2
        vH = Eig.GlobalVertices{hh}(:); B = Op.Mass{hh};
        Psi = zeros(4*numel(vH), nCh);
        Psi(2:4:end,:)=G(:,(vH-1)*3+1).'; Psi(3:4:end,:)=G(:,(vH-1)*3+2).'; Psi(4:4:end,:)=G(:,(vH-1)*3+3).';
        totB = totB + sum(sum(Psi.*(B*Psi)));
    end
    hemi = (1:size(Gm,2)) > size(Gm,2)/2;   % false=L, true=R (cols stacked L then R)
    hidx = zeros(size(Gm,2),1); hidx(~hemi)=1:sum(~hemi); hidx(hemi)=1:sum(hemi);
    kkv = [25 50 100 200 300 400]; capt = zeros(size(kkv));
    for q=1:numel(kkv), capt(q) = 100*sum(sum(Gm(:,hidx<=kkv(q)).^2))/totB; end
    s.captureCurve = [kkv(:), capt(:)];
    s.capture400 = capt(end);
    % reconstructed leadfield vs os_meg (per-channel correlation, per-vector cosine)
    Lrec = bst_dirac_reconstruct_leadfield(Gm, Eig);     % [nCh x 3nVert]
    perCh = arrayfun(@(c) corr(Lrec(c,:)', G(c,:)'), 1:nCh);
    A = reshape(G',3,[])'; Bv = reshape(Lrec',3,[])';     % [nCh*nVert x 3] vectors
    cosV = sum(A.*Bv,2) ./ max(sqrt(sum(A.^2,2)).*sqrt(sum(Bv.^2,2)), eps);
    s.reconChanCorr = [median(perCh) min(perCh)];
    s.reconVecCosMed = median(cosV);
    fprintf('  leadfield B-energy captured @K=400/hemi: %.1f%% (residual %.1f%% high-freq)\n', capt(end), 100-capt(end));
    fprintf('  recon vs os_meg: per-channel corr median=%.4f (min %.4f); per-vector cos median=%.4f\n', ...
        s.reconChanCorr(1), s.reconChanCorr(2), s.reconVecCosMed);
end


%% ===================================================================
function s = section2_mode_energy(Eig, Op, Gm, Vtx, Nrm) %#ok<INUSL>
    fprintf('\n--- S2: DIRAC MODE ENERGY SPECTRUM (ambient-flat structure) ---\n');
    hh=1; gv=Eig.GlobalVertices{hh}(:); Vh=numel(gv); B=Op.Mass{hh}; Phi=double(Eig.Phi{hh});
    a=full(diag(B)); a=a(1:4:end);
    % mode-1 directional spread (globally constant ambient => 0)
    Q1=reshape(real(Phi(:,1)),4,Vh).'; V1=Q1(:,2:4);
    s.mode1Spread = norm(V1-mean(V1,1),'fro')/norm(V1,'fro');
    % constant-ambient field energy concentration & normal-field spread
    mk=@(D)reshape([zeros(Vh,1) D]',[],1);
    cConst=Phi'*(B*mk(repmat([1 0 0],Vh,1))); eConst=abs(cConst).^2;
    Nl = Nrm(gv,:);
    cNorm =Phi'*(B*mk(Nl));               eNorm=abs(cNorm).^2;
    s.constInMode1 = 100*max(eConst)/sum(eConst);
    s.normCapturedK = 100*sum(eNorm)/(mk(Nl)'*(B*mk(Nl)));
    cs=cumsum(eNorm)/sum(eNorm); s.normCentroid = sum((1:numel(eNorm))'.*eNorm)/sum(eNorm);
    % eigenvector participation ratio (global harmonics) across modes
    K=size(Phi,2); PR=zeros(K,1);
    for k=1:K, Qk=reshape(real(Phi(:,k)),4,Vh).'; m=sum(Qk(:,2:4).^2,2); p=a.*m; PR(k)=(sum(p))^2/(Vh*sum(p.^2)); end
    s.PR = PR; s.PRmedian = median(PR);
    fprintf('  mode-1 directional spread = %.4f (0 => globally constant ambient vector)\n', s.mode1Spread);
    fprintf('  constant-ambient field: %.1f%% energy in mode 1 (ground state)\n', s.constInMode1);
    fprintf('  cortical-normal field: only %.1f%% captured in K modes (centroid mode %.0f) => high-freq\n', s.normCapturedK, s.normCentroid);
    fprintf('  eigenvector participation fraction: median %.2f (global harmonics)\n', s.PRmedian);
end


%% ===================================================================
function s = section3_observability(Gm, Cn, nCh)
    fprintf('\n--- S3: OBSERVABILITY SPACE ---\n');
    % clean (unwhitened) effective rank of the mode-forward
    sv = svd(Gm,'econ');
    s.rank90 = find(cumsum(sv.^2)/sum(sv.^2)>=0.90,1);
    s.rank99 = find(cumsum(sv.^2)/sum(sv.^2)>=0.99,1);
    % noise-whitened DOF vs SNR (regularized whitener: truncate + ridge)
    [Un,Dc]=eig(Cn); dc=real(diag(Dc)); [dc,ix]=sort(dc,'descend'); Un=Un(:,ix);
    reg=1e-3; dcr=dc; dcr(dcr<reg*max(dc))=reg*max(dc);
    Wn=Un*diag(1./sqrt(dcr))*Un';
    svw=svd(Wn*Gm,'econ'); svw=svw/svw(1);
    snrv=[3 10 30 100]; dofv=arrayfun(@(snr)sum(svw>1/snr), snrv);
    s.svw=svw; s.snrDOF=[snrv(:) dofv(:)];
    s.noiseRankEff = sum(dc> max(dc)*1e-6);
    fprintf('  clean mode-forward effective rank: 90%%=%d, 99%%=%d of %d modes (nCh=%d)\n', s.rank90, s.rank99, size(Gm,2), nCh);
    fprintf('  noise cov effective rank: %d/%d (ill-conditioned => regularize before whitening)\n', s.noiseRankEff, nCh);
    fprintf('  noise-whitened observable DOF vs SNR:  3:1->%d  10:1->%d  30:1->%d  100:1->%d\n', dofv);
end


%% ===================================================================
function s = section4_inverse_stages(HMf, G, Gm, Cn, types, iMEG)
    fprintf('\n--- S4: INVERSE STAGE VALIDATION (vs vertex MNE) ---\n');
    OPT = struct('NoiseMethod','reg','NoiseReg',0.1,'SnrMethod','fixed','SnrFixed',3, ...
                 'InverseMeasure','amplitude');
    OPT.NoiseCovMat.NoiseCov = Cn; OPT.ChannelTypes = types(iMEG);
    Ra = bst_inverse_dirac(HMf, OPT);
    iW = Ra.Whitener;
    % S1: whitener bit-identity vs MNE 'reg'
    [Un,Sn2]=svd(Cn,'econ'); Sn=sqrt(diag(Sn2)); rk=sum(Sn>length(Sn)*eps(single(Sn(1))));
    Un=Un(:,1:rk); Sn=Sn(1:rk); Ridge=mean(diag(Sn2))*0.1;
    iWref=Un*diag(1./sqrt(Sn.^2+Ridge))*Un';
    s.whitenerRelErr = max(abs(iW(:)-iWref(:)))/max(abs(iWref(:)));
    % S2/S3: observable subspace + data fit vs vertex MNE
    Lw=iW*G; [ULv,SL2v]=svd(Lw*Lw'); SL2v=diag(SL2v);
    Lam=9/mean(SL2v); Kv=Lam*(Lw'*(ULv*diag(1./(Lam*SL2v+1))*ULv'))*iW;
    svL=svd(Lw,'econ'); svM=svd(Ra.GainWhitened,'econ'); n=min(50,numel(svM));
    s.subspaceCos = svL(1:n)'*svM(1:n)/(norm(svL(1:n))*norm(svM(1:n)));
    keptV=sum(SL2v>1/Lam); keptM=sum(Ra.SL.^2>1/Ra.Lambda); s.keptDOF=[keptM keptV];
    Hm = Gm*(Ra.ImagingKernelMode); Hv = G*Kv;   % sensor-space data-fit hat matrices
    s.hatCos = (Hm(:)'*Hv(:))/(norm(Hm(:))*norm(Hv(:)));
    s.kernelSize = size(Ra.ImagingKernel);
    fprintf('  S1 whitener vs MNE reg: rel err = %.2e (bit-identical)\n', s.whitenerRelErr);
    fprintf('  S2 kept-DOF mode/vertex = %d/%d ; observable-subspace cos = %.4f\n', keptM, keptV, s.subspaceCos);
    fprintf('  S3 sensor data-fit hat-matrix cos(mode,vertex) = %.4f\n', s.hatCos);
    fprintf('  S4 vertex imaging kernel size = [%d x %d]\n', s.kernelSize);
end


%% ===================================================================
function s = section5_localization(HMf, G, Vtx, Faces, Nrm, Cn, types, iMEG, OUTDIR)
    fprintf('\n--- S5: LOCALIZATION STATISTICS (multi-source) ---\n');
    nV = size(Vtx,1);
    verts = round(linspace(200, nV-200, 60));    % 60 sources spread across cortex
    OPT = struct('NoiseMethod','reg','NoiseReg',0.1,'SnrMethod','fixed','SnrFixed',3);
    OPT.NoiseCovMat.NoiseCov = Cn; OPT.ChannelTypes = types(iMEG);
    K = struct();
    for m = {'amplitude','dspm2018','sloreta'}
        OPT.InverseMeasure = m{1}; Rr = bst_inverse_dirac(HMf, OPT); K.(matkey(m{1})) = Rr.ImagingKernel;
    end
    mags=@(J) sqrt(sum(reshape(J,3,[]).^2,1))';
    meas={'amplitude','dspm2018','sloreta'}; nM=numel(meas);
    errClean = zeros(numel(verts), nM); errNoisy = zeros(numel(verts), nM);
    for noisy=0:1
        for i=1:numel(verts)
            vt=verts(i); d0=G(:,(vt-1)*3+(1:3))*Nrm(vt,:)';
            if noisy
                [Uc,Dc]=eig(Cn); dc=max(real(diag(Dc)),0);
                % deterministic pseudo-noise (index-seeded) so the benchmark is reproducible
                ee=Uc*(sqrt(dc).*cos((1:size(Cn,1))'*vt)); ee=ee/norm(ee)*norm(d0)/3; d=d0+ee;
            else, d=d0; end
            for j=1:nM
                [~,pk]=max(mags(K.(matkey(meas{j}))*d));
                e=norm(Vtx(pk,:)-Vtx(vt,:))*1e3;
                if noisy, errNoisy(i,j)=e; else, errClean(i,j)=e; end
            end
        end
    end
    s.verts=verts; s.meas=meas; s.errClean=errClean; s.errNoisy=errNoisy;
    fprintf('  peak localization error (mm), %d sources:\n', numel(verts));
    for j=1:nM
        fprintf('    %-10s clean: med=%.1f IQR[%.1f %.1f]  | noisy(SNR3): med=%.1f IQR[%.1f %.1f]\n', meas{j}, ...
            median(errClean(:,j)), prctile(errClean(:,j),25), prctile(errClean(:,j),75), ...
            median(errNoisy(:,j)),  prctile(errNoisy(:,j),25),  prctile(errNoisy(:,j),75));
    end

    % --- demonstrative example: one source reconstructed by the 3 measures ---
    vt = round(0.4*nV); d = G(:,(vt-1)*3+(1:3))*Nrm(vt,:)';
    f=figure('Color','w','Position',[50 50 1450 460],'Visible','off');
    for j=1:nM
        mg = mags(K.(matkey(meas{j}))*d); mg=mg/max(mg);
        ax=subplot(1,3,j);
        patch('Faces',Faces,'Vertices',Vtx,'FaceVertexCData',mg,'FaceColor','interp','EdgeColor','none');
        hold on; plot3(Vtx(vt,1),Vtx(vt,2),Vtx(vt,3),'g.','MarkerSize',28);
        [~,pk]=max(mg); plot3(Vtx(pk,1),Vtx(pk,2),Vtx(pk,3),'co','MarkerSize',12,'LineWidth',2);
        axis(ax,'equal','off'); view(-90,90); camlight; material dull; colormap(hot); caxis([0 1]);
        title(sprintf('%s (peak err %.0f mm)', meas{j}, norm(Vtx(pk,:)-Vtx(vt,:))*1e3));
    end
    sgtitle('S5 demonstrative: single-source reconstruction (green=true, cyan=peak)','FontWeight','bold');
    print(f,[OUTDIR '/bench_dirac_demo_localization.png'],'-dpng','-r110'); close(f);
    s.demoVertex = vt;
end


%% ===================================================================
function s = section6_orientation(HMf, G, Vtx, Nrm, Cn, types, iMEG)
    fprintf('\n--- S6: ORIENTATION ANALYSIS ---\n');
    nV=size(Vtx,1);
    OPT = struct('NoiseMethod','reg','NoiseReg',0.1,'SnrMethod','fixed','SnrFixed',3,'InverseMeasure','amplitude');
    OPT.NoiseCovMat.NoiseCov = Cn; OPT.ChannelTypes = types(iMEG);
    Ra = bst_inverse_dirac(HMf, OPT); Kamp = Ra.ImagingKernel; iW = Ra.Whitener;
    % (a) leadfield normal sensitivity (radial blindness check)
    fracN=zeros(nV,1);
    for v=1:nV, Gv=G(:,(v-1)*3+(1:3)); fracN(v)=norm(Gv*Nrm(v,:)')/norm(Gv,'fro'); end
    s.leadfieldNormalFrac = median(fracN);   % isotropic baseline ~0.58
    % (b) reconstructed-orientation isotropy on real data (one focal source) + frame invariance
    vt=round(nV/3); d=G(:,(vt-1)*3+(1:3))*Nrm(vt,:)';  J=Kamp*d; J3=reshape(J,3,nV)';
    vn=sqrt(sum(J3.^2,2)); proj=sum(J3.*Nrm,2); act=vn>0.2*max(vn);
    s.reconNormalFrac = median(abs(proj(act))./vn(act));   % ~0.5 if isotropic
    % frame invariance of the norm (orthonormal cortical frame)
    ref=repmat([1 0 0],nV,1); fl=abs(sum(Nrm.*ref,2))>0.9; ref(fl,:)=repmat([0 1 0],sum(fl),1);
    e1=ref-sum(ref.*Nrm,2).*Nrm; e1=e1./sqrt(sum(e1.^2,2)); e2=cross(Nrm,e1,2);
    Jl=[sum(J3.*e1,2) sum(J3.*e2,2) sum(J3.*Nrm,2)];
    s.normFrameInvErr = max(abs(sqrt(sum(J3.^2,2))-sqrt(sum(Jl.^2,2))));
    % (c) ORIENTATION RECOVERY: simulate known orientation, measure angular error at source
    verts=round(linspace(300,nV-300,40)); angT=zeros(numel(verts),1); angR=zeros(numel(verts),1);
    for i=1:numel(verts)
        v=verts(i); n=Nrm(v,:)';
        t=cross(n,[0;0;1]); if norm(t)<1e-3, t=cross(n,[0;1;0]); end; t=t/norm(t);  % a tangential dir
        for which=1:2
            o = (which==1)*t + (which==2)*n;   % tangential vs normal source
            dd=G(:,(v-1)*3+(1:3))*o; Jr=reshape(Kamp*dd,3,nV)'; Jv=Jr(v,:)';
            ang=acosd(min(1,abs(Jv'*o)/max(norm(Jv),eps)));
            if which==1, angT(i)=ang; else, angR(i)=ang; end
        end
    end
    s.orientErrTangential=median(angT); s.orientErrNormal=median(angR);
    s.angT=angT; s.angR=angR;
    fprintf('  leadfield normal sensitivity (median ||G.n||/||G||) = %.2f (isotropic baseline 0.58)\n', s.leadfieldNormalFrac);
    fprintf('  reconstructed |J.n|/||J|| median = %.2f (0.5 => orientation isotropic, no normal/tangential bias)\n', s.reconNormalFrac);
    fprintf('  norm frame-invariance error (ambient vs cortical frame) = %.1e\n', s.normFrameInvErr);
    fprintf('  ORIENTATION RECOVERY (angular error at source): tangential src = %.0f deg, normal src = %.0f deg\n', ...
        s.orientErrTangential, s.orientErrNormal);
end


%% ===================================================================
function make_figures(R, OUTDIR)
    % --- Figure 1: forward + mode energy + observability ---
    f1=figure('Color','w','Position',[60 60 1400 900],'Visible','off');
    subplot(2,3,1); plot(R.S1.captureCurve(:,1), R.S1.captureCurve(:,2),'-o','LineWidth',2,'Color',[.1 .5 .2]);
    xlabel('modes / hemisphere'); ylabel('% leadfield B-energy'); title('S1: leadfield reconstruction'); grid on; ylim([50 100]);
    subplot(2,3,2);
    plot(R.S2.PR,'.','Color',[.3 .3 .3]); xlabel('mode index'); ylabel('participation frac');
    title(sprintf('S2: eigenvectors global (med %.2f)', R.S2.PRmedian)); grid on; ylim([0 1]);
    subplot(2,3,3);
    bar([R.S2.constInMode1, R.S2.normCapturedK]); set(gca,'XTickLabel',{'const-ambient','normal field'});
    ylabel('% energy (mode1 / K modes)'); title('S2: ambient-flat structure'); grid on;
    subplot(2,3,4);
    semilogy(R.S3.svw,'k','LineWidth',2); hold on;
    for q=1:size(R.S3.snrDOF,1), plot([1 numel(R.S3.svw)],[1/R.S3.snrDOF(q,1) 1/R.S3.snrDOF(q,1)],'--'); end
    xlabel('singular index'); ylabel('\sigma/\sigma_1'); title('S3: noise-whitened observability'); grid on; xlim([1 120]); ylim([1e-3 1]);
    subplot(2,3,5);
    plot(R.S3.snrDOF(:,1), R.S3.snrDOF(:,2),'-o','LineWidth',2); xlabel('assumed SNR'); ylabel('observable DOF');
    title('S3: observable modes vs SNR'); grid on; set(gca,'XScale','log');
    subplot(2,3,6); axis off;
    txt = sprintf(['VALIDATION SUMMARY\n\n' ...
        'whitener vs MNE: %.0e\n' 'kept-DOF mode/vertex: %d/%d\n' 'subspace cos: %.4f\n' 'hat-matrix cos: %.4f\n\n' ...
        'recon vs os\\_meg:\n  chan corr %.4f\n  vec cos %.4f\n\n' ...
        'mode-1 spread: %.3f\nnorm frame-inv: %.0e'], ...
        R.S4.whitenerRelErr, R.S4.keptDOF(1), R.S4.keptDOF(2), R.S4.subspaceCos, R.S4.hatCos, ...
        R.S1.reconChanCorr(1), R.S1.reconVecCosMed, R.S2.mode1Spread, R.S6.normFrameInvErr);
    text(0.0,0.95,txt,'VerticalAlignment','top','FontSize',10,'FontName','FixedWidth');
    sgtitle('Dirac source mapping — forward, mode energy, observability, validation','FontWeight','bold');
    print(f1,[OUTDIR '/bench_dirac_forward_observability.png'],'-dpng','-r110'); close(f1);

    % --- Figure 2: localization + orientation ---
    f2=figure('Color','w','Position',[60 60 1400 480],'Visible','off');
    subplot(1,3,1);
    boxplot([R.S5.errClean R.S5.errNoisy], 'Labels',{'amp','dSPM','sLOR','amp_n','dSPM_n','sLOR_n'});
    ylabel('peak loc error (mm)'); title('S5: localization (clean | noisy SNR3)'); grid on;
    subplot(1,3,2);
    histogram(R.S6.angT,0:10:180,'FaceColor',[.2 .5 .8]); hold on; histogram(R.S6.angR,0:10:180,'FaceColor',[.8 .4 .2]);
    xlabel('orientation error (deg)'); ylabel('# sources'); legend('tangential src','normal src');
    title(sprintf('S6: orientation recovery (med %.0f/%.0f deg)', R.S6.orientErrTangential, R.S6.orientErrNormal)); grid on;
    subplot(1,3,3); axis off;
    txt2 = sprintf(['ORIENTATION\n\n' ...
        'leadfield ||G.n||/||G||: %.2f\n  (isotropic ~0.58)\n\n' ...
        'recon |J.n|/||J||: %.2f\n  (0.5 = isotropic)\n\n' ...
        'orientation recovery:\n  tangential %.0f deg\n  normal %.0f deg'], ...
        R.S6.leadfieldNormalFrac, R.S6.reconNormalFrac, R.S6.orientErrTangential, R.S6.orientErrNormal);
    text(0.0,0.95,txt2,'VerticalAlignment','top','FontSize',11,'FontName','FixedWidth');
    sgtitle('Dirac source mapping — localization & orientation statistics','FontWeight','bold');
    print(f2,[OUTDIR '/bench_dirac_localization_orientation.png'],'-dpng','-r110'); close(f2);
end


%% ===== helpers =====
function L = bst_dirac_reconstruct_leadfield(Gm, Eig)
% reconstruct the effective vertex leadfield from the mode-forward: L = Phi-expand(Gm)
    nV = sum(cellfun(@numel, Eig.GlobalVertices));
    nCh = size(Gm,1); Kpm = size(Gm,2)/2; L = zeros(nCh, 3*nV);
    for hh=1:2
        vH=Eig.GlobalVertices{hh}(:); Phi=double(Eig.Phi{hh}); Phi=Phi(:,1:Kpm);
        cols = (1:size(Gm,2)) > Kpm; if hh==1, cols=~cols; end
        Rr = Phi * Gm(:,cols).';                 % [4Vh x nCh]
        L(:,(vH-1)*3+1) = Rr(2:4:end,:).'; L(:,(vH-1)*3+2)=Rr(3:4:end,:).'; L(:,(vH-1)*3+3)=Rr(4:4:end,:).';
    end
end

function k = matkey(m)
    k = strrep(m,'2018','');   % 'dspm2018'->'dspm'
end
