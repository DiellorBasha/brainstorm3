function R = bench_face_dirac_inverse(frameTime)
% BENCH_FACE_DIRAC_INVERSE  Decisive end-to-end test of the face-Dirac eigenbasis inverse.
% Compares, on a real frame: (1) face Dirac inverse (bst_inverse_dirac, face leadfield +
% 'Dirac-Face' eigenbasis), (2) vertex Dirac inverse (bst_inverse_dirac, vertex), and
% (3) plain face wMNE (no eigenbasis). The question: does the face Dirac eigenbasis --
% which captures the observable subspace less efficiently than vertex (wide-root, non-
% smooth modes) -- still localize the alpha source where the vertex/ wMNE do?
% USAGE: R = bench_face_dirac_inverse(22.6)
% Author: Diellor Basha, 2026
    if nargin<1 || isempty(frameTime), frameTime = 22.6; end
    OUTDIR = '/Users/diellorbasha/workspace/research/code/brainstorm3/dev/benchmarks';
    df = 'Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat';
    [sStudy,~] = bst_get('DataFile', df);
    ChanMat = in_bst_channel(sStudy.Channel(1).FileName); types = {ChanMat.Channel.Type};
    HMos = in_bst_headmodel([fileparts(df) '/headmodel_surf_os_meg.mat'], 0);
    G = double(HMos.Gain); iMEG = all(isfinite(G),2) & strcmpi(types(:),'MEG');
    NC = load(file_fullpath([fileparts(df) '/noisecov_full.mat'])); Cn = NC.NoiseCov(iMEG,iMEG); Cn=(Cn+Cn')/2;
    OPT = struct('NoiseMethod','reg','NoiseReg',0.1,'SnrMethod','fixed','SnrFixed',3, ...
                 'InverseMeasure','dspm2018','nModes',400,'Tau',0.5);
    OPT.NoiseCovMat.NoiseCov = Cn;  OPT.ChannelTypes = types(iMEG);
    DM = in_bst_data(df); [~,iT] = min(abs(DM.Time-frameTime)); d = double(DM.F(iMEG,iT));
    Gv = G(iMEG,:);
    SurfaceFile = HMos.SurfaceFile;  Surf = in_tess_bst(SurfaceFile,0);
    V = Surf.Vertices; F = double(Surf.Faces);

    % ---- (1) vertex Dirac inverse ----
    HMv = HMos; HMv.Gain = Gv;
    Rv = bst_inverse_dirac(HMv, OPT);  Pv = i_pow(Rv.ImagingKernel * d);     % [nV]

    % ---- (2) face Hodge-eigenbasis inverse (the working face vector basis) ----
    [Lf, FG] = bst_face_leadfield(SurfaceFile, ChanMat.Channel(iMEG), HMos.Param(iMEG), 'Mode','unconstrained');
    HMf = struct('Gain',Lf,'SurfaceFile',SurfaceFile,'HeadModelType','surface','isFaceBased',1, ...
                 'FaceBasis','hodge', 'GridLoc',FG.Centroids);
    Rf = bst_inverse_dirac(HMf, OPT);  Pf = i_pow(Rf.ImagingKernel * d);     % [nF]

    % ---- (3) plain face wMNE (no eigenbasis) ----
    Cr = Cn + 0.1*mean(diag(Cn))*eye(size(Cn));  W = inv(chol(Cr,'lower'));
    GGt = (W*Lf)*(W*Lf)';  lam = trace(GGt)/(size(Lf,1)*9);
    Jw = (W*Lf)' * ((GGt + lam*eye(size(Lf,1))) \ (W*d));  Pw = i_pow3(Jw);   % [nF]

    % ---- compare ----
    Cf = FG.Centroids;
    [~,iv]=max(Pv); [~,ifd]=max(Pf); [~,ifw]=max(Pw);
    pkV=V(iv,:); pkFd=Cf(ifd,:); pkFw=Cf(ifw,:);
    sepDirac = 1000*norm(pkV-pkFd);  sepWmne = 1000*norm(pkV-pkFw);  sepFF = 1000*norm(pkFd-pkFw);
    PvF = (Pv(F(:,1))+Pv(F(:,2))+Pv(F(:,3)))/3;
    ccDV = corr(Pf, PvF);  ccDW = corr(Pf, Pw);
    fprintf('\n=== face-Hodge inverse vs vertex-Dirac vs face-wMNE @ %.3f s ===\n', DM.Time(iT));
    fprintf('peak: vertex-Dirac %s\n      face-Hodge   %s  (%.1f mm from vertex)\n      face-wMNE    %s  (%.1f mm from vertex)\n', ...
        mat2str(pkV,3), mat2str(pkFd,3), sepDirac, mat2str(pkFw,3), sepWmne);
    fprintf('face-Hodge vs face-wMNE peak separation: %.1f mm\n', sepFF);
    fprintf('corr(face-Hodge power, vertex-Dirac@faces) = %.3f   corr(face-Hodge, face-wMNE) = %.3f\n', ccDV, ccDW);

    % ---- 3-panel power maps ----
    hFig = figure('Color','w','Position',[40 80 1500 480]);
    pan = {Pv,'interp',pkV,'vertex Dirac |J|'; Pf,'flat',pkFd,'face Hodge |J|'; Pw,'flat',pkFw,'face wMNE |J|'};
    for sp=1:3
        ax=subplot(1,3,sp); hold(ax,'on'); scal=pan{sp,1};
        patch('Vertices',V,'Faces',F,'FaceVertexCData',scal,'FaceColor',pan{sp,2},'EdgeColor','none','Parent',ax);
        m=prctile(scal,99.5); if m<=0, m=max(scal); end
        pk=pan{sp,3}; plot3(ax,pk(1),pk(2),pk(3),'o','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',9,'Clipping','off');
        colormap(ax,hot(256)); clim(ax,[0 m]); axis(ax,'equal','off'); view(ax,[0 90]);
        camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull'); title(ax,pan{sp,4});
    end
    sgtitle(sprintf('Face Hodge-eigenbasis inverse @ %.2f s | peak vs vertex %.1f mm, vs wMNE %.1f mm', DM.Time(iT), sepDirac, sepFF));
    png = fullfile(OUTDIR,'bench_face_dirac_inverse.png'); print(hFig,png,'-dpng','-r130');
    fprintf('saved %s\n', png);
    R = struct('Pv',Pv,'Pf',Pf,'Pw',Pw,'sepDirac_mm',sepDirac,'sepFF_mm',sepFF,'corrDV',ccDV,'corrDW',ccDW,'png',png);
end
function P = i_pow(J),  P = sqrt(sum(reshape(J,3,[]).^2,1))'; end
function P = i_pow3(J), P = sqrt(sum(reshape(J,3,[]).^2,1))'; end
