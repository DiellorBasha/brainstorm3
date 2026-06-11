function R = analyze_alpha_dirac(dataFile)
% ANALYZE_ALPHA_DIRAC  Dirac-mode scale analysis of an alpha-band data block.
%
% Computes, for the alpha block:
%   1. amplitude Dirac MODE-COEFFICIENT time series c(t) = KernelMode * M(t),
%      the scale-spectrogram E(lambda,t)=|c_k(t)|^2, and the lambda-energy
%      CENTROID(t) -- the diffusion test (centroid drops within a cycle =>
%      point->spread diffusion; flat => standing oscillation).
%   2. the Hilbert amplitude envelope: saves the SENSOR envelope as a Brainstorm
%      data file, and maps the (correct) SOURCE envelope |K*hilbert(M)| onto cortex.
%   3. the time-averaged lambda-spectrum over the high-alpha window (21.5-24 s) ->
%      the characteristic spatial SCALE of the alpha activation.
% dSPM is used only for the display map; amplitude mode coeffs for the scale axis.
%
% USAGE:  R = analyze_alpha_dirac('Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat')
% Author: Diellor Basha, 2026

    if nargin<1 || isempty(dataFile)
        dataFile='Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat';
    end
    OUTDIR='/Users/diellorbasha/workspace/research/code/brainstorm3/dev/benchmarks';
    WIN=[21.5 24];                                   % high-alpha window (s)

    % ---- study / inputs (same kernel the user has been visualizing) ----
    [sStudy,iStudy]=bst_get('DataFile', dataFile);
    ChannelFile=sStudy.Channel(1).FileName; ChanMat=in_bst_channel(ChannelFile); types={ChanMat.Channel.Type};
    HMos=in_bst_headmodel([fileparts(dataFile) '/headmodel_surf_os_meg.mat'],0);
    G=double(HMos.Gain); iMEG=all(isfinite(G),2)&strcmpi(types(:),'MEG'); G=G(iMEG,:);
    Srf=in_tess_bst(HMos.SurfaceFile); Vtx=Srf.Vertices; Faces=Srf.Faces; nV=size(Vtx,1);
    NC=load(file_fullpath([fileparts(dataFile) '/noisecov_full.mat'])); Cn=NC.NoiseCov(iMEG,iMEG); Cn=(Cn+Cn')/2;
    HMf=HMos; HMf.Gain=G;

    OPT=struct('NoiseMethod','reg','NoiseReg',0.1,'SnrMethod','fixed','SnrFixed',3);
    OPT.NoiseCovMat.NoiseCov=Cn; OPT.ChannelTypes=types(iMEG);
    OPT.InverseMeasure='amplitude'; Ra=bst_inverse_dirac(HMf,OPT);
    OPT.InverseMeasure='dspm2018';  Rd=bst_inverse_dirac(HMf,OPT);
    Kmode=Ra.ImagingKernelMode; Kdspm=Rd.ImagingKernel; lam=double(Ra.Eigenvalues);

    % ---- data ----
    DataMat=in_bst_data(dataFile); Time=DataMat.Time; M=double(DataMat.F(iMEG,:));
    inW = Time>=WIN(1)&Time<=WIN(2);
    fprintf('Block %.1f-%.1f s, %d samples (%.0f Hz), %d MEG ch. alpha window %.1f-%.1f s.\n', ...
        Time(1),Time(end),numel(Time),1/median(diff(Time)),numel(iMEG),WIN(1),WIN(2));

    % ---- (1) mode-coefficient time series + scale spectrogram + centroid ----
    c = Kmode*M;                                     % [2K x nTime] amplitude mode coeffs
    E = c.^2;                                        % mode energy
    [lamS,ord]=sort(lam,'ascend'); Es=E(ord,:);
    nz = lamS>1e-9;                                  % exclude the ~0 kernel modes for the centroid
    cent = sum(lamS(nz).*Es(nz,:),1)./max(sum(Es(nz,:),1),eps);   % lambda-energy centroid(t)
    R.lambda=lamS; R.E=Es; R.centroid=cent; R.Time=Time;

    % ---- dSPM source movie + alpha source-point time series ----
    Jd=Kdspm*M; md=squeeze(sqrt(sum(reshape(Jd,3,nV,[]).^2,1)));  % [nV x nTime] dSPM magnitude
    [~,pv]=max(max(md(:,inW),[],2)); srcTS=md(pv,:);
    R.peakVertex=pv;
    fprintf('alpha dSPM peak vertex %d (%.0f %.0f %.0f mm)\n', pv, Vtx(pv,:)*1e3);

    % ---- (2) Hilbert envelopes ----
    envSensor = abs(hilbert(DataMat.F.')).';         % [nChanAll x nTime] sensor envelope (all chans)
    Ma = hilbert(M.').';                             % analytic MEG signal
    Ja = Kdspm*Ma;                                   % complex source (dSPM)
    envSrc = squeeze(sqrt(sum(abs(reshape(Ja,3,nV,[])).^2,1)));   % source amplitude envelope [nV x nTime]
    R.envSrcPeakMap = envSrc(:, find(inW,1) + round(0.5*sum(inW)) );   % a mid-window envelope map

    % save the SENSOR envelope as a Brainstorm data file (user request)
    EnvMat=DataMat; EnvMat.F=envSensor;
    EnvMat.Comment=[DataMat.Comment ' | alpha env (|hilbert|)'];
    EnvMat=bst_history('add',EnvMat,'process','Hilbert amplitude envelope (alpha band)');
    OutFile=bst_process('GetNewFilename', bst_fileparts(file_fullpath(dataFile)), 'data_alphaEnv');
    bst_save(OutFile,EnvMat,'v6');
    try
        sNew=db_template('Data'); sNew.FileName=file_short(OutFile); sNew.Comment=EnvMat.Comment; sNew.DataType='recordings'; sNew.BadTrial=0;
        sStudy.Data(end+1)=sNew; bst_set('Study',iStudy,sStudy); panel_protocols('UpdateNode','Study',iStudy);
        fprintf('Saved + registered sensor envelope: %s\n', file_short(OutFile));
    catch e
        fprintf('Saved sensor envelope file (registration skipped: %s): %s\n', e.message, file_short(OutFile));
    end

    % ---- (3) lambda-spectrum over the alpha window -> characteristic scale ----
    spec = mean(Es(:,inW),2); spec=spec/sum(spec);   % normalized lambda-energy spectrum
    cumS=cumsum(spec); k50=find(cumS>=0.5,1); k90=find(cumS>=0.9,1);
    lamCentWin = sum(lamS(nz).*spec(nz))/sum(spec(nz));
    R.spec=spec; R.k50=k50; R.k90=k90; R.lamCentWin=lamCentWin;
    fprintf('(3) alpha-window scale: energy 50%%@mode %d, 90%%@mode %d (of %d); lambda-centroid=%.3e\n', k50,k90,numel(lam),lamCentWin);

    % ================= figures =================
    tsub = Time>=WIN(1)-0.3 & Time<=WIN(2)+0.3;
    f1=figure('Color','w','Position',[40 40 1300 760],'Visible','off');
    subplot(3,1,1);
    imagesc(Time(tsub), 1:nnz(nz), log10(Es(nz,tsub)+eps)); axis xy;
    ylabel('Dirac mode (sorted by \lambda \uparrow)'); title('(1) scale-spectrogram  log_{10}|c_k(t)|^2'); colorbar;
    subplot(3,1,2);
    plot(Time(tsub), cent(tsub),'LineWidth',1.5); xlabel('time (s)'); ylabel('\lambda centroid');
    title('\lambda-energy centroid(t) -- DIFFUSION test (drops within cycle => spread)'); grid on; xlim(WIN+[-0.3 0.3]);
    subplot(3,1,3);
    plot(Time(tsub), srcTS(tsub)/max(srcTS(tsub)),'LineWidth',1.5); xlabel('time (s)'); ylabel('dSPM @peak (norm.)');
    title(sprintf('alpha source-point time series (vertex %d)',pv)); grid on; xlim(WIN+[-0.3 0.3]);
    sgtitle('Alpha Dirac-mode scale analysis','FontWeight','bold');
    print(f1,[OUTDIR '/analyze_alpha_dirac_spectrogram.png'],'-dpng','-r110'); close(f1);

    f2=figure('Color','w','Position',[60 60 1300 460],'Visible','off');
    subplot(1,2,1);
    semilogy(lamS, spec,'.','Color',[.2 .4 .7]); hold on; xline(lamCentWin,'k:');
    xlabel('Dirac eigenvalue \lambda'); ylabel('norm. energy'); grid on;
    title(sprintf('(3) alpha-window \\lambda-spectrum (50%%@%d, 90%%@%d modes)',k50,k90));
    subplot(1,2,2);
    mp=R.envSrcPeakMap; mp=mp/max(mp);
    patch('Faces',Faces,'Vertices',Vtx,'FaceVertexCData',mp,'FaceColor','interp','EdgeColor','none');
    hold on; plot3(Vtx(pv,1),Vtx(pv,2),Vtx(pv,3),'co','MarkerSize',12,'LineWidth',2);
    axis equal off; view(180,-90); camlight; material dull; colormap(hot); clim([0 1]);
    title('source amplitude envelope |K\cdothilbert(M)| (occipital view)');
    sgtitle('Alpha scale spectrum + source envelope map','FontWeight','bold');
    print(f2,[OUTDIR '/analyze_alpha_dirac_scale.png'],'-dpng','-r110'); close(f2);
    fprintf('Saved %s/analyze_alpha_dirac_{spectrogram,scale}.png\n', OUTDIR);
end
