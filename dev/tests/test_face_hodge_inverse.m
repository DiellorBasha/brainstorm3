function test_face_hodge_inverse()
% Face Hodge vector eigenbasis ('Hodge-Face') end-to-end: smooth full-rank basis
% (scalar lapFace modes lifted to {grad, n x grad}), B-orthonormal, captures the
% observable subspace, and -- the decisive gate -- LOCALIZES (unlike Dirac-Face's 43mm).
% Author: Diellor Basha, 2026
    nFail = 0;  K = 400;
    df = 'Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat';
    [sStudy,~] = bst_get('DataFile', df);
    BaseHM = in_bst_headmodel([fileparts(df) '/headmodel_surf_os_meg.mat'], 0);
    ChanMat = in_bst_channel(sStudy.Channel(1).FileName); types = {ChanMat.Channel.Type};
    iMEG = find(strcmpi(types,'MEG'));
    SurfaceFile = BaseHM.SurfaceFile;
    Surf = in_tess_bst(SurfaceFile,0); V = Surf.Vertices; F = double(Surf.Faces);

    % --- eigenbasis: pure-imaginary, B-orthonormal ---
    Eig = tess_eigen(SurfaceFile, 'Hodge-Face', 'K', K, 'Tau', 0.5);
    Op  = tess_operators(SurfaceFile, 'Hodge-Face', 'NoSave', true);
    for hh = 1:2
        Phi = Eig.Phi{hh}; B = Op.Mass{hh};
        nFail = nFail + chk(sprintf('h%d Phi [4nF x K], w-rows zero',hh), ...
            size(Phi,2)>=K && max(abs(Phi(1:4:end,:)),[],'all') < 1e-12);
        G = Phi(:,1:K)'*B*Phi(:,1:K);
        nFail = nFail + chk(sprintf('h%d B-orthonormal',hh), norm(G-eye(K),'fro')/sqrt(K) < 1e-6);
    end

    % --- leadfield + transform: observability retained (Dirac-Face managed only ~100/178) ---
    [Lf, FG] = bst_face_leadfield(SurfaceFile, ChanMat.Channel(iMEG), BaseHM.Param(iMEG), 'Mode','unconstrained');
    HMf = struct('Gain',Lf,'SurfaceFile',SurfaceFile,'HeadModelType','surface','isFaceBased',1,'FaceBasis','hodge','GridLoc',FG.Centroids);
    Comp = bst_dirac(HMf, 'nModes', K, 'Tau', 0.5);
    er = @(M) sum(svd(M) > 1e-3*max(svd(M)));
    rL = er(Lf); rM = er(Comp.Gain);
    fprintf('  observability: face LF=%d   Hodge mode-forward=%d (%.0f%%)\n', rL, rM, 100*rM/rL);
    nFail = nFail + chk('Hodge mode-forward retains >= 80%% observability', rM >= 0.80*rL);

    % --- THE GATE: the Hodge face inverse localizes near the vertex Dirac inverse ---
    NC = load(file_fullpath([fileparts(df) '/noisecov_full.mat'])); Cn = NC.NoiseCov(iMEG,iMEG); Cn=(Cn+Cn')/2;
    OPT = struct('NoiseMethod','reg','NoiseReg',0.1,'SnrMethod','fixed','SnrFixed',3, ...
                 'InverseMeasure','dspm2018','nModes',K,'Tau',0.5);
    OPT.NoiseCovMat.NoiseCov = Cn; OPT.ChannelTypes = types(iMEG);
    DM = in_bst_data(df); [~,iT] = min(abs(DM.Time-22.6)); d = double(DM.F(iMEG,iT));
    G = double(BaseHM.Gain); ivMEG = all(isfinite(G),2)&strcmpi(types(:),'MEG');
    HMv = struct('Gain',G(ivMEG,:),'SurfaceFile',SurfaceFile,'HeadModelType','surface');
    Rv = bst_inverse_dirac(HMv, OPT); Pv = i_p(Rv.ImagingKernel*d);
    Rf = bst_inverse_dirac(HMf, OPT); Pf = i_p(Rf.ImagingKernel*d);
    [~,iv]=max(Pv); [~,iff]=max(Pf); sep = 1000*norm(V(iv,:) - FG.Centroids(iff,:));
    PvF = (Pv(F(:,1))+Pv(F(:,2))+Pv(F(:,3)))/3; cc = corr(Pf, PvF);
    fprintf('  [GATE] face-Hodge peak %.1f mm from vertex-Dirac; power corr %.3f\n', sep, cc);
    nFail = nFail + chk('Hodge inverse peak < 20 mm from vertex-Dirac', sep < 20);
    nFail = nFail + chk('Hodge inverse power corr > 0.70', cc > 0.70);

    fprintf('\n==== test_face_hodge_inverse: %d failed ====\n', nFail);
    if nFail > 0, error('test_face_hodge_inverse FAILED'); end
end
function P = i_p(J), P = sqrt(sum(reshape(J,3,[]).^2,1))'; end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
