function R = compare_dirac_vs_mne(base, dataFile)
% COMPARE_DIRAC_VS_MNE  Head-to-head of Dirac-basis dSPM vs standard vertex dSPM
% on a real evoked recording (TutorialAuditory M100). Everything except the
% source basis is matched: same noise covariance, SNR, regularization, no depth
% weighting, unconstrained (free) orientation, dSPM measure.
%
% USAGE:  R = compare_dirac_vs_mne
%         R = compare_dirac_vs_mne(base, dataFile)
%
% Authors: Diellor Basha, 2026

    if nargin < 1 || isempty(base), base = 'Subject01/S01_AEF_20131218_01_notch/'; end
    if nargin < 2 || isempty(dataFile)
        sStudy = bst_get('StudyWithCondition', base(1:end-1));
        dataFile = '';
        for i=1:numel(sStudy.Data)
            if contains(sStudy.Data(i).FileName,'deviant_average'), dataFile = sStudy.Data(i).FileName; break; end
        end
        if isempty(dataFile), dataFile = sStudy.Data(1).FileName; end
    end
    OUTDIR = '/Users/diellorbasha/workspace/research/code/brainstorm3/dev/benchmarks';
    fprintf('\n=== Dirac dSPM vs standard dSPM on %s ===\n', dataFile);

    % ---- shared inputs ----
    HMos = in_bst_headmodel([base 'headmodel_surf_os_meg.mat'], 0);
    ChanMat = in_bst_channel([base 'channel_ctf_acc1.mat']); types = {ChanMat.Channel.Type};
    G = double(HMos.Gain); iMEG = all(isfinite(G),2) & strcmpi(types(:),'MEG'); G = G(iMEG,:);
    Srf = in_tess_bst(HMos.SurfaceFile); Vtx = Srf.Vertices; Nrm = Srf.VertNormals; Faces = Srf.Faces; nV = size(Vtx,1);
    NC = load(file_fullpath([base 'noisecov_full.mat'])); Cn = NC.NoiseCov(iMEG,iMEG); Cn=(Cn+Cn')/2;
    DataMat = in_bst_data(dataFile); F = double(DataMat.F(iMEG,:)); Time = DataMat.Time;

    % ---- Dirac unconstrained dSPM ----
    OPTd = struct('NoiseMethod','reg','NoiseReg',0.1,'SnrMethod','fixed','SnrFixed',3,'InverseMeasure','dspm2018');
    OPTd.NoiseCovMat.NoiseCov = Cn; OPTd.ChannelTypes = types(iMEG);
    Rd = bst_inverse_dirac(struct('Gain',G,'SurfaceFile',HMos.SurfaceFile,'HeadModelType','surface'), OPTd);
    Kd = Rd.ImagingKernel;

    % ---- standard vertex unconstrained dSPM (matched options, no depth weighting) ----
    HMm = struct('Gain',G,'GridLoc',Vtx,'GridOrient',Nrm,'GridAtlas',[], ...
                 'HeadModelType','surface','SurfaceFile',HMos.SurfaceFile, ...
                 'MEGMethod','os_meg','EEGMethod','','ECOGMethod','','SEEGMethod','');
    OPTm = struct('InverseMethod','minnorm','InverseMeasure','dspm2018','SourceOrient',{{'free'}}, ...
                  'DataTypes',{{'MEG'}},'NoiseMethod','reg','NoiseReg',0.1,'SnrMethod','fixed','SnrFixed',3, ...
                  'UseDepth',0,'ComputeKernel',1,'DisplayMessages',0);
    OPTm.NoiseCovMat.NoiseCov = Cn; OPTm.ChannelTypes = types(iMEG);
    Rm = bst_inverse_linear_2018(HMm, OPTm);
    Km = Rm.ImagingKernel;

    % ---- apply to data, per-vertex norm maps ----
    Jd = Kd*F; Js = Km*F;
    nd = squeeze(sqrt(sum(reshape(Jd,3,nV,[]).^2,1)));   % [nV x nT]
    ns = squeeze(sqrt(sum(reshape(Js,3,nV,[]).^2,1)));

    % ---- M100 peak (global max in 60-160 ms) ----
    win = Time>=0.06 & Time<=0.16;
    gfpD = max(nd,[],1); [~,it] = max(gfpD.*win); tpk = it; tms = Time(tpk)*1e3;
    [~,pd]=max(nd(:,tpk)); [~,ps]=max(ns(:,tpk));

    % ---- comparison metrics ----
    R.dataFile=dataFile; R.tpk_ms=tms;
    R.mapCorr = corr(nd(:,tpk), ns(:,tpk));
    R.peakDist_mm = norm(Vtx(pd,:)-Vtx(ps,:))*1e3;
    thr=0.5; dD = nd(:,tpk)>thr*max(nd(:,tpk)); dS = ns(:,tpk)>thr*max(ns(:,tpk));
    R.diceTop50 = 2*sum(dD&dS)/(sum(dD)+sum(dS));
    R.tcCorrAtPeak = corr(nd(pd,:)', ns(pd,:)');
    % whole movie spatial-temporal correlation
    R.movieCorr = corr(nd(:), ns(:));
    fprintf('M100 peak at %.0f ms\n', tms);
    fprintf('  spatial map corr (Dirac vs MNE) @peak = %.3f\n', R.mapCorr);
    fprintf('  peak-vertex distance = %.1f mm ; top-50%% Dice = %.2f\n', R.peakDist_mm, R.diceTop50);
    fprintf('  time-course corr at peak vertex = %.3f ; whole-movie corr = %.3f\n', R.tcCorrAtPeak, R.movieCorr);

    % ---- figure ----
    f=figure('Color','w','Position',[40 40 1500 760],'Visible','off');
    maps={nd(:,tpk)/max(nd(:,tpk)), ns(:,tpk)/max(ns(:,tpk))}; tit={'Dirac dSPM','standard dSPM'};
    views=[-90 90; 90 -90];   % top of left, top of right-ish (two views)
    for s=1:2
        subplot(2,3,s);
        patch('Faces',Faces,'Vertices',Vtx,'FaceVertexCData',maps{s},'FaceColor','interp','EdgeColor','none');
        axis equal off; view(-90,90); camlight; material dull; colormap(hot); clim([0 1]);
        title(sprintf('%s @ %.0f ms', tit{s}, tms));
    end
    subplot(2,3,3);
    scatter(ns(:,tpk)/max(ns(:,tpk)), nd(:,tpk)/max(nd(:,tpk)), 4, 'filled', 'MarkerFaceAlpha',0.2); hold on;
    plot([0 1],[0 1],'k--'); xlabel('standard dSPM (norm.)'); ylabel('Dirac dSPM (norm.)');
    title(sprintf('per-vertex @peak (r=%.3f)', R.mapCorr)); axis square; grid on;
    subplot(2,3,[4 5]);
    plot(Time*1e3, nd(pd,:)/max(nd(pd,:)),'LineWidth',2); hold on;
    plot(Time*1e3, ns(ps,:)/max(ns(ps,:)),'LineWidth',2);
    xline(tms,'k:'); xlabel('time (ms)'); ylabel('source mag (norm.)'); xlim([Time(1) Time(end)]*1e3);
    legend('Dirac (its peak vtx)','MNE (its peak vtx)','Location','best');
    title('Source time courses at each method''s peak vertex'); grid on;
    subplot(2,3,6); axis off;
    txt=sprintf(['MATCHED dSPM (unconstrained)\nbasis is the only difference\n\n' ...
        'M100 peak: %.0f ms\n\nspatial map corr: %.3f\npeak dist: %.1f mm\ntop-50%% Dice: %.2f\n' ...
        'peak-vtx TC corr: %.3f\nwhole-movie corr: %.3f'], ...
        tms, R.mapCorr, R.peakDist_mm, R.diceTop50, R.tcCorrAtPeak, R.movieCorr);
    text(0,0.95,txt,'VerticalAlignment','top','FontSize',11,'FontName','FixedWidth');
    sgtitle(sprintf('Dirac dSPM vs standard dSPM on real evoked data (%s)', strrep(DataMat.Comment,'_','\_')),'FontWeight','bold');
    print(f,[OUTDIR '/compare_dirac_vs_mne.png'],'-dpng','-r110'); close(f);
    fprintf('Saved %s/compare_dirac_vs_mne.png\n', OUTDIR);
end
