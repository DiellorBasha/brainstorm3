function R = bench_dirac_streamribbon_real(frameTime, dataFile)
% BENCH_DIRAC_STREAMRIBBON_REAL  Streamribbon the REAL solenoidal Dirac field.
%
% Pipeline (mirrors view_helmholtz + the dirac benchmarks):
%   1. build the dSPM unconstrained kernel with bst_inverse_dirac
%   2. at the chosen frame, J = K*M  ->  ambient 3D source vectors [nV x 3]
%   3. Helmholtz/Hodge split (bst_helmholtz) -> solenoidal field Vsol + cores
%   4. seed nested rings around the strongest posterior vortex core and ribbon Vsol
%
% USAGE:  R = bench_dirac_streamribbon_real()           % frame ~22.6 s (alpha burst)
%         R = bench_dirac_streamribbon_real(22.6)
% Author: Diellor Basha, 2026

    if nargin<1 || isempty(frameTime), frameTime = 22.6; end
    if nargin<2 || isempty(dataFile)
        dataFile = 'Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat';
    end
    OUTDIR = '/Users/diellorbasha/workspace/research/code/brainstorm3/dev/benchmarks';

    % ---- inverse kernel (same recipe as analyze_alpha_dirac) -------------------
    [sStudy,~] = bst_get('DataFile', dataFile);
    ChanMat = in_bst_channel(sStudy.Channel(1).FileName); types = {ChanMat.Channel.Type};
    HMos = in_bst_headmodel([fileparts(dataFile) '/headmodel_surf_os_meg.mat'], 0);
    G = double(HMos.Gain); iMEG = all(isfinite(G),2) & strcmpi(types(:),'MEG'); G = G(iMEG,:);
    NC = load(file_fullpath([fileparts(dataFile) '/noisecov_full.mat']));
    Cn = NC.NoiseCov(iMEG,iMEG); Cn = (Cn+Cn')/2;
    HMf = HMos; HMf.Gain = G;
    OPT = struct('NoiseMethod','reg','NoiseReg',0.1,'SnrMethod','fixed','SnrFixed',3, ...
                 'InverseMeasure','dspm2018');
    OPT.NoiseCovMat.NoiseCov = Cn; OPT.ChannelTypes = types(iMEG);
    Rd = bst_inverse_dirac(HMf, OPT);  Kdspm = Rd.ImagingKernel;

    % ---- data frame ------------------------------------------------------------
    DataMat = in_bst_data(dataFile); Time = DataMat.Time; M = double(DataMat.F(iMEG,:));
    [~, iT] = min(abs(Time - frameTime));
    Jt = Kdspm * M(:,iT);                                  % [3nV x 1]

    % ---- surface + operators + Helmholtz split ---------------------------------
    SurfaceFile = HMos.SurfaceFile;
    Surf = in_tess_bst(SurfaceFile, 0); V = Surf.Vertices; F = double(Surf.Faces);
    nV = size(V,1);
    if isfield(Surf,'VertNormals') && ~isempty(Surf.VertNormals)
        VN = Surf.VertNormals;  VN = VN ./ max(sqrt(sum(VN.^2,2)),eps);
    else
        VN = vertexNormals(V, F);
    end
    Dirac = loadOp(SurfaceFile, 'Dirac');
    LBO   = loadOp(SurfaceFile, 'Laplace-Beltrami');
    Mani  = tess_manifold(SurfaceFile);
    Op    = bst_helmholtz('Prepare', {Dirac, LBO}, Mani, Surf, 'Domain','vertex');
    Ht    = bst_helmholtz('Frame', Op, Jt);
    Vsol  = Ht.Vsol;                                       % [nV x 3] solenoidal field

    % ---- pick the strongest vortex core in the posterior cortex ----------------
    cores = Ht.Cores;
    if isempty(cores)
        [~, core] = max(sqrt(sum(Vsol.^2,2)));             % fallback: max |Vsol|
        fprintf('No vortex cores found; seeding at max |Vsol| vertex %d.\n', core);
    else
        cv = [cores.iVertex];  om = abs([cores.omega]);
        post = V(cv,1) < 0;                                % SCS: x<0 = posterior
        if ~any(post), post = true(size(cv)); end
        om(~post) = -inf;
        [~, ic] = max(om);  core = cv(ic);
        fprintf('%d cores (%d posterior). Seeding strongest posterior core: vertex %d, |omega|=%.3g.\n', ...
                numel(cores), sum(post), core, abs(cores(ic).omega));
    end
    PATCH_R = 0.014;
    coreVerts = core;
    if ~isempty(cores)
        cv = [cores.iVertex];
        coreVerts = cv(sqrt(sum((V(cv,:) - V(core,:)).^2,2)) <= PATCH_R);   % cores in patch
    end
    fprintf('Frame %.3f s (idx %d).  harmonic fraction %.1f%%.  core @ [%.0f %.0f %.0f] mm\n', ...
            Time(iT), iT, 100*Ht.HarmFrac, V(core,:)*1e3);

    % ---- smoothed display surface (unwrap folds for legibility) ----------------
    VertConn = tess_vertconn(V, F);
    Vsm  = tess_smooth(V, 0.5, 80, VertConn, 1, F);
    VNsm = vertexNormals(Vsm, F);

    % ---- trace on folded geometry, display on smoothed (shared renderer) -------
    opts = struct('SIGMA',0.006, 'PATCH_R',PATCH_R, 'STEP',0.0010, 'MAXSTEP',140, ...
                  'seedRadii',linspace(0.3,1.4,8), 'CoreVerts',coreVerts, ...
                  'HW',0.0018, 'LIFT',0.0011, 'surfAlpha',0.25, ...
                  'DispV',Vsm, 'DispVN',VNsm, ...
                  'outPng',fullfile(OUTDIR,'bench_dirac_streamribbon_real.png'), ...
                  'titleStr',sprintf('Real Dirac solenoidal field @ %.2f s  —  S01 alpha vortex on SMOOTHED cortex (rings=rotation, banking=depth)', Time(iT)));
    R = streamribbon_surface(V, F, VN, Vsol, core, opts);
    R.Ht = Ht; R.frameIdx = iT; R.frameTime = Time(iT);
end


% ===========================================================================
function Op = loadOp(SurfaceFile, variant)
    [sSubject,~,iSurf] = bst_get('SurfaceFile', SurfaceFile);
    Op = [];
    if ~isempty(iSurf) && isfield(sSubject.Surface(iSurf),'Operator')
        for k = 1:numel(sSubject.Surface(iSurf).Operator)
            S = load(file_fullpath(sSubject.Surface(iSurf).Operator(k).FileName));
            if strcmpi(S.Variant, variant), Op = S; break; end
        end
    end
    if isempty(Op), tess_operators(SurfaceFile, variant); Op = loadOp(SurfaceFile, variant); end
end

function VN = vertexNormals(V, F)
    fn = cross(V(F(:,2),:)-V(F(:,1),:), V(F(:,3),:)-V(F(:,1),:));
    VN = zeros(size(V));
    for j=1:3, VN(:,j) = accumarray(F(:),[fn(:,j);fn(:,j);fn(:,j)],[size(V,1) 1]); end
    VN = VN ./ max(sqrt(sum(VN.^2,2)),eps);
end
