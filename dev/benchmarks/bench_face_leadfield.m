function R = bench_face_leadfield(frameTime)
% BENCH_FACE_LEADFIELD  Face vs vertex UNCONSTRAINED inverse on a real frame.
% Builds the full-unconstrained face leadfield (raw Sarvas at centroids, tess_manifold
% geometry) and runs a self-contained whitened minimum-norm (wMNE) on it, against the
% vertex leadfield. NOTE: bst_inverse_dirac is NOT used here -- it reconstructs in the
% stored VERTEX Dirac eigenbasis (it would silently return a vertex-sized kernel for a
% face leadfield); the face-Dirac eigenbasis inverse is the next phase. The inline wMNE
% genuinely lives in the face column space, so this validates the leadfield end-to-end.
% Compares observability and source-power localization (observable subspace ~tens of
% DOF is basis-invariant, so the maps should co-localize; the face map is finer).
% USAGE: R = bench_face_leadfield(22.6)
% Author: Diellor Basha, 2026
    if nargin<1 || isempty(frameTime), frameTime = 22.6; end
    OUTDIR = '/Users/diellorbasha/workspace/research/code/brainstorm3/dev/benchmarks';
    df = 'Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat';
    [sStudy,~] = bst_get('DataFile', df);
    ChanMat = in_bst_channel(sStudy.Channel(1).FileName); types = {ChanMat.Channel.Type};
    HMos = in_bst_headmodel([fileparts(df) '/headmodel_surf_os_meg.mat'], 0);
    G = double(HMos.Gain); iMEG = all(isfinite(G),2) & strcmpi(types(:),'MEG');
    NC = load(file_fullpath([fileparts(df) '/noisecov_full.mat'])); Cn = NC.NoiseCov(iMEG,iMEG); Cn=(Cn+Cn')/2;
    DM = in_bst_data(df); [~,iT] = min(abs(DM.Time-frameTime)); d = double(DM.F(iMEG,iT));
    Gv = G(iMEG,:);

    SurfaceFile = HMos.SurfaceFile;  Surf = in_tess_bst(SurfaceFile,0);
    V = Surf.Vertices; F = double(Surf.Faces);

    % whitener from the (regularized) noise covariance, shared by both inverses
    Cr = Cn + 0.1*mean(diag(Cn))*eye(size(Cn));  W = inv(chol(Cr,'lower'));   % W*Cr*W' = I
    dW = W*d;  SNR = 3;
    wMNE = @(Gw) i_wmne(Gw, dW, SNR);

    % ---- vertex unconstrained wMNE ----
    Jv = wMNE(W*Gv);                                                        % [3nV]
    Pv = sqrt(sum(reshape(Jv,3,[]).^2,1))';                                % per-vertex power [nV]

    % ---- face unconstrained wMNE ----
    Channel = ChanMat.Channel(iMEG);  Param = HMos.Param(iMEG);
    [Lf, FG] = bst_face_leadfield(SurfaceFile, Channel, Param, 'Mode','unconstrained');
    Jf = wMNE(W*Lf);                                                        % [3nF]
    Pf = sqrt(sum(reshape(Jf,3,[]).^2,1))';                                % per-face power [nF]
    HMv = struct('Gain',Gv);                                               % for the rank print below

    % ---- compare ----
    er = @(M) sum(svd(M) > 1e-3*max(svd(M)));
    [~,ivPk] = max(Pv); [~,ifPk] = max(Pf);
    pkV = V(ivPk,:); pkF = FG.Centroids(ifPk,:); dPk = 1000*norm(pkV-pkF);   % mm
    % project vertex power to faces (barycentric) for a correlation on a common domain
    PvF = (Pv(F(:,1))+Pv(F(:,2))+Pv(F(:,3)))/3;  cc = corr(Pf, PvF);
    fprintf('\n=== face vs vertex unconstrained inverse @ %.3f s ===\n', DM.Time(iT));
    fprintf('gain: vertex [%d x %d]  face [%d x %d]\n', size(HMv.Gain,1),size(HMv.Gain,2), size(Lf,1),size(Lf,2));
    fprintf('whitened-leadfield effective rank: vertex %d  face %d\n', er(HMv.Gain), er(Lf));
    fprintf('peak power location: vertex %s  face %s  | separation %.1f mm\n', mat2str(pkV,3), mat2str(pkF,3), dPk);
    fprintf('corr(face power, vertex power@faces) = %.3f\n', cc);

    % ---- side-by-side power maps ----
    hFig = figure('Color','w','Position',[60 80 1200 560]);
    for sp = 1:2
        ax = subplot(1,2,sp); hold(ax,'on');
        if sp==1, scal=Pv; fc='interp'; pk=pkV; ttl='vertex unconstrained |J| (dSPM)';
        else,     scal=Pf; fc='flat';   pk=pkF; ttl='face unconstrained |J| (dSPM)'; end
        patch('Vertices',V,'Faces',F,'FaceVertexCData',scal,'FaceColor',fc,'EdgeColor','none','Parent',ax);
        m = prctile(scal,99.5); if m<=0, m=max(scal); end
        plot3(ax,pk(1),pk(2),pk(3),'o','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',9,'Clipping','off');
        colormap(ax,hot(256)); clim(ax,[0 m]); axis(ax,'equal','off'); view(ax,[0 90]);
        camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull'); title(ax,ttl);
    end
    sgtitle(sprintf('Face vs vertex full-unconstrained inverse @ %.2f s (peak sep %.1f mm)', DM.Time(iT), dPk));
    png = fullfile(OUTDIR,'bench_face_leadfield.png'); print(hFig,png,'-dpng','-r140');
    fprintf('saved %s\n', png);
    R = struct('Jv',Jv,'Jf',Jf,'Pv',Pv,'Pf',Pf,'pkSep_mm',dPk,'corr',cc,'png',png);
end

function J = i_wmne(Gw, dW, SNR)
% Whitened minimum-norm (identity source cov): J = Gw' (Gw Gw' + lambda I)^-1 dW,
% lambda = trace(Gw Gw')/(nCh * SNR^2) so regularization scales with the leadfield.
    GGt = Gw * Gw';  nCh = size(Gw,1);
    lambda = trace(GGt)/(nCh * SNR^2);
    J = Gw' * ((GGt + lambda*eye(nCh)) \ dW);
end
