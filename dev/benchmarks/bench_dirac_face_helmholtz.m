function R = bench_dirac_face_helmholtz(frameTime)
% BENCH_DIRAC_FACE_HELMHOLTZ  Compare vertex vs face-NATIVE Helmholtz on a real Dirac frame.
% Interpolates the unconstrained vertex dSPM solution to face centroids and decomposes it
% with bst_helmholtz Domain=face (face-native gradFace coupled Hodge), against
% bst_helmholtz Domain=vertex. Now apples-to-apples (both self-consistent, genus-0 harmonic~0).
% USAGE: R = bench_dirac_face_helmholtz(22.6)
% Author: Diellor Basha, 2026
    if nargin<1 || isempty(frameTime), frameTime = 22.6; end
    OUTDIR = '/Users/diellorbasha/workspace/research/code/brainstorm3/dev/benchmarks';
    df = 'Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat';
    [sStudy,~] = bst_get('DataFile', df);
    ChanMat = in_bst_channel(sStudy.Channel(1).FileName); types = {ChanMat.Channel.Type};
    HMos = in_bst_headmodel([fileparts(df) '/headmodel_surf_os_meg.mat'], 0);
    G = double(HMos.Gain); iMEG = all(isfinite(G),2) & strcmpi(types(:),'MEG');
    NC = load(file_fullpath([fileparts(df) '/noisecov_full.mat'])); Cn = NC.NoiseCov(iMEG,iMEG); Cn=(Cn+Cn')/2;
    HMf = HMos; HMf.Gain = G(iMEG,:);
    OPT = struct('NoiseMethod','reg','NoiseReg',0.1,'SnrMethod','fixed','SnrFixed',3,'InverseMeasure','dspm2018');
    OPT.NoiseCovMat.NoiseCov = Cn; OPT.ChannelTypes = types(iMEG);
    Rd = bst_inverse_dirac(HMf, OPT);
    DM = in_bst_data(df); [~,iT] = min(abs(DM.Time-frameTime)); Jt = Rd.ImagingKernel*double(DM.F(iMEG,iT));

    SurfaceFile = HMos.SurfaceFile;  Surf = in_tess_bst(SurfaceFile,0);
    V = Surf.Vertices; F = double(Surf.Faces);
    Dirac = i_op(SurfaceFile,'Dirac'); LBO = i_op(SurfaceFile,'Laplace-Beltrami');
    Mani  = tess_manifold(SurfaceFile);
    FaceOp = i_op(SurfaceFile,'Hodge-Face');                 % carries FaceAux.GradFace

    % vertex pipeline
    OpV = bst_helmholtz('Prepare', {Dirac, LBO}, Mani, Surf, 'Domain','vertex');
    HtV = bst_helmholtz('Frame', OpV, Jt);

    % interpolate vertex field -> face centroids (barycentric mean of the 3 vertex vectors)
    J3 = reshape(Jt,3,[])';                                  % [nV x 3]
    Jf = (J3(F(:,1),:) + J3(F(:,2),:) + J3(F(:,3),:))/3;     % [nF x 3]

    % face pipeline (intrinsic dual)
    OpF = bst_helmholtz('Prepare', FaceOp, Mani, Surf, 'Domain','face');
    HtF = bst_helmholtz('Frame', OpF, Jf);

    % ---- compare (count only persistence-significant cores; the 3-valence face dual
    %      yields many tiny-persistence raw extrema, so gate both sides identically) ----
    GATE = 0.1;                                            % keep persistence >= 10% of per-hemi max
    cv = @(c) i_count(bst_persistence_gate(c, GATE));
    fprintf('\n=== vertex vs face Helmholtz @ %.3f s ===\n', DM.Time(iT));
    fprintf('HarmFrac:   vertex %.1f%%   face %.1f%%\n', 100*HtV.HarmFrac, 100*HtF.HarmFrac);
    fprintf('vortices(+/-) [persist>=%.0f%%]: vertex %s   face %s\n', 100*GATE, mat2str(cv(HtV.Cores)),   mat2str(cv(HtF.Cores)));
    fprintf('sources(+/-)  [persist>=%.0f%%]: vertex %s   face %s\n', 100*GATE, mat2str(cv(HtV.Sources)), mat2str(cv(HtF.Sources)));
    fprintf('|Curl| max: vertex %.3g  face %.3g  | |Div| max: vertex %.3g  face %.3g\n', ...
        max(abs(HtV.Curl)), max(abs(HtF.Curl)), max(abs(HtV.Div)), max(abs(HtF.Div)));

    % ---- side-by-side viz: Psi (stream fn) + top cores ----
    hFig = figure('Color','w','Position',[60 80 1200 560]); cmap = i_divmap(256);
    for sp = 1:2
        ax = subplot(1,2,sp); hold(ax,'on');
        if sp==1, scal=HtV.Psi; cc=HtV.Cores; ttl='vertex Helmholtz (psi)'; fc='interp';
        else,     scal=HtF.Psi; cc=HtF.Cores; ttl='face Helmholtz (psi, face-native)'; fc='flat'; end
        patch('Vertices',V,'Faces',F,'FaceVertexCData',scal,'FaceColor',fc, ...
              'EdgeColor','none','Parent',ax);
        m=prctile(abs(scal(scal~=0)),99); if isempty(m)||m<=0, m=max(abs(scal)); end  % robust clim (outliers wash out max)
        if m<=0, m=eps; end
        if ~isempty(cc)
            pr=[cc.persistence]; mxf=max([pr(isfinite(pr)),eps]); keep=cc(isinf(pr)|pr>=0.5*mxf);
            for k=1:numel(keep)
                col=[1 0 0]; if keep(k).charge<0, col=[0 0 1]; end
                p=keep(k).pos; plot3(ax,p(1),p(2),p(3),'o','MarkerFaceColor',col,'MarkerEdgeColor','k','MarkerSize',7,'Clipping','off');
            end
        end
        colormap(ax,cmap); clim(ax,[-m m]); axis(ax,'equal','off'); view(ax,[0 90]);
        camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull'); title(ax,ttl);
    end
    sgtitle(sprintf('Vertex vs face-domain (intrinsic) Helmholtz @ %.2f s', DM.Time(iT)));
    png = fullfile(OUTDIR,'bench_dirac_face_helmholtz.png'); print(hFig,png,'-dpng','-r140');
    fprintf('saved %s\n', png);
    R = struct('HtV',HtV,'HtF',HtF,'frameTime',DM.Time(iT),'png',png);
end

function Op = i_op(SurfaceFile, variant)
    [sSubject,~,iSurf] = bst_get('SurfaceFile', SurfaceFile);
    Op = [];
    if ~isempty(iSurf) && isfield(sSubject.Surface(iSurf),'Operator')
        for k = 1:numel(sSubject.Surface(iSurf).Operator)
            S = load(file_fullpath(sSubject.Surface(iSurf).Operator(k).FileName));
            if strcmpi(S.Variant, variant), Op = S; break; end
        end
    end
    if isempty(Op), tess_operators(SurfaceFile, variant); Op = i_op(SurfaceFile, variant); end
end
function c = i_count(mk)
    if isempty(mk), c = [0 0]; else, c = [sum([mk.charge]>0) sum([mk.charge]<0)]; end
end
function m = i_divmap(n)
    t=linspace(0,1,n)'; lo=[.23 .30 .75]; mid=[.96 .96 .96]; hi=[.78 .15 .18];
    m=[interp1([0 .5 1],[lo(1) mid(1) hi(1)],t), interp1([0 .5 1],[lo(2) mid(2) hi(2)],t), interp1([0 .5 1],[lo(3) mid(3) hi(3)],t)];
end
