function test_bst_inverse_dirac(HeadModelFile, NoiseCovFile, ChannelFile)
% TEST_BST_INVERSE_DIRAC  Regression test for the staged Dirac-basis MNE inverse.
%
% Verifies bst_inverse_dirac Stages 1-5 against the canonical vertex-space MNE
% (bst_inverse_linear_2018), on the live protocol's overlapping-spheres model:
%   T1  Stage 1 whitener bit-identical to MNE 'reg' formula
%   T2  whitening identity: iW*C_reg*iW' is a rank-(noiseRank) projector
%   T3  observable subspace basis-invariant (whitened mode-forward vs vertex)
%   T4  Stage 2 regularization corner matches vertex MNE (kept DOF, filter shape)
%   T5  Stage 3 data fit identical to vertex MNE (focal-source sensor refit corr)
%   T6  Stage 5 dSPM removes the depth bias (peak loc error << amplitude)
%
% USAGE:  test_bst_inverse_dirac    % uses the TutorialAuditory defaults below
%
% Requires: the protocol loaded, an Overlapping spheres head model, a noise
% covariance, and a 'Dirac' eigen node on the head model's surface.
%
% Authors: Diellor Basha, 2026

    if nargin < 1 || isempty(HeadModelFile)
        base = 'Subject01/S01_AEF_20131218_01_notch/';
        HeadModelFile = [base 'headmodel_surf_os_meg.mat'];
        NoiseCovFile  = [base 'noisecov_full.mat'];
        ChannelFile   = [base 'channel_ctf_acc1.mat'];
    end
    % (the Dirac eigenbasis is auto-found/created by bst_dirac from the surface)
    tol = struct('whit',1e-9, 'proj',1e-4, 'cos',1e-3, 'dof',3, 'fit',1e-3);
    nPass = 0; nFail = 0;

    % ---- assemble inputs (MEG data channels only) ----
    HMos = in_bst_headmodel(HeadModelFile, 0);
    ChanMat = in_bst_channel(ChannelFile);
    allTypes = {ChanMat.Channel.Type};
    G = double(HMos.Gain);
    iMEG = all(isfinite(G),2) & strcmpi(allTypes(:),'MEG');
    G = G(iMEG,:); nCh = size(G,1);
    HMf = HMos; HMf.Gain = G;
    CompHM = bst_dirac(HMf, 'nModes',400, 'Tau',0.5);
    Srf = in_tess_bst(HMos.SurfaceFile); Vtx = Srf.Vertices; Nrm = Srf.VertNormals;
    NC = load(file_fullpath(NoiseCovFile));

    OPT = struct('NoiseMethod','reg','NoiseReg',0.1,'SnrMethod','fixed','SnrFixed',3);
    OPT.NoiseCovMat.NoiseCov = NC.NoiseCov(iMEG,iMEG);
    OPT.ChannelTypes = allTypes(iMEG);

    % pass the VERTEX head model (new interface; transformed internally)
    OPT.InverseMeasure = 'amplitude';
    Ra = bst_inverse_dirac(HMf, OPT);
    OPT.InverseMeasure = 'dspm2018';
    Rd = bst_inverse_dirac(HMf, OPT);
    OPT.InverseMeasure = 'sloreta';
    Rs = bst_inverse_dirac(HMf, OPT);

    iW = Ra.Whitener;
    C  = (OPT.NoiseCovMat.NoiseCov + OPT.NoiseCovMat.NoiseCov')/2;

    % ---- T1: whitener bit-identical to MNE 'reg' (single modality) ----
    [Un,Sn2] = svd(C,'econ'); Sn = sqrt(diag(Sn2));
    rk = sum(Sn > length(Sn)*eps(single(Sn(1))));
    Un = Un(:,1:rk); Sn = Sn(1:rk); Ridge = mean(diag(Sn2))*OPT.NoiseReg;
    iW_ref = Un*diag(1./sqrt(Sn.^2+Ridge))*Un';
    [nPass,nFail] = check('T1 whitener == MNE reg', ...
        max(abs(iW(:)-iW_ref(:)))/max(abs(iW_ref(:))) < tol.whit, nPass, nFail);

    % ---- T2: whitening identity (projector of rank = noiseRank) ----
    P = iW*(C + Ridge*eye(nCh))*iW';
    [nPass,nFail] = check('T2 iW*C_reg*iW'' is projector', ...
        (abs(max(eig(P))-1) < tol.proj) && (abs(trace(P)-Ra.NoiseRankKept) < tol.dof), nPass, nFail);

    % ---- T3: observable subspace basis-invariant (mode vs vertex) ----
    Lw = iW*G; svL = svd(Lw,'econ'); svM = svd(Ra.GainWhitened,'econ');
    n = min(50,numel(svM)); csp = svL(1:n)'*svM(1:n)/(norm(svL(1:n))*norm(svM(1:n)));
    [nPass,nFail] = check('T3 observable subspace cos>0.999', (1-csp) < tol.cos, nPass, nFail);

    % ---- T4: Stage 2 corner matches vertex MNE ----
    [~,SL2v] = svd(Lw*Lw'); SL2v = diag(SL2v);
    Lam_v = OPT.SnrFixed^2/mean(SL2v);
    keptV = sum(SL2v > 1/Lam_v);  keptM = sum(Ra.SL.^2 > 1/Ra.Lambda);
    [nPass,nFail] = check('T4 kept-DOF mode~vertex', abs(keptM-keptV) <= tol.dof, nPass, nFail);

    % ---- T5: Stage 3 data fit identical to vertex MNE (focal source) ----
    [ULv,SL2v2] = svd(Lw*Lw'); SL2v2 = diag(SL2v2);
    Kv = Lam_v*(Lw'*(ULv*diag(1./(Lam_v*SL2v2+1))*ULv'))*iW;     % vertex MNE kernel
    Gm = double(CompHM.Gain);
    vt = 5000; d = G(:,(vt-1)*3+(1:3))*Nrm(vt,:)';
    dfit_m = Gm*(Ra.ImagingKernelMode*d); dfit_v = G*(Kv*d);
    [nPass,nFail] = check('T5 mode vs vertex sensor-refit corr==1', ...
        (1-corr(dfit_m,dfit_v)) < tol.fit, nPass, nFail);

    % ---- T6: dSPM removes depth bias (peak loc error) ----
    mags = @(J) sqrt(sum(reshape(J,3,[]).^2,1))';
    verts = [5000 1200 8000 15000 18000]; eA = []; eD = []; eS = [];
    for v = verts
        dd = G(:,(v-1)*3+(1:3))*Nrm(v,:)';
        [~,pa] = max(mags(Ra.ImagingKernel*dd)); [~,pd] = max(mags(Rd.ImagingKernel*dd));
        [~,ps] = max(mags(Rs.ImagingKernel*dd));
        eA(end+1) = norm(Vtx(pa,:)-Vtx(v,:))*1e3; %#ok<AGROW>
        eD(end+1) = norm(Vtx(pd,:)-Vtx(v,:))*1e3; %#ok<AGROW>
        eS(end+1) = norm(Vtx(ps,:)-Vtx(v,:))*1e3; %#ok<AGROW>
    end
    [nPass,nFail] = check(sprintf('T6 dSPM<<amplitude (%.0f vs %.0f mm)',median(eD),median(eA)), ...
        (median(eD) < median(eA)) && (median(eD) < 15), nPass, nFail);
    [nPass,nFail] = check(sprintf('T7 sLORETA localizes (%.0f mm median)',median(eS)), ...
        (median(eS) < median(eA)) && (median(eS) < 15), nPass, nFail);

    fprintf('\n==== test_bst_inverse_dirac: %d passed, %d failed ====\n', nPass, nFail);
    if nFail > 0, error('test_bst_inverse_dirac: %d test(s) FAILED.', nFail); end
end

% -----------------------------------------------------------------------------
function [nPass,nFail] = check(name, cond, nPass, nFail)
    if cond
        fprintf('  PASS  %s\n', name); nPass = nPass + 1;
    else
        fprintf('  FAIL  %s\n', name); nFail = nFail + 1;
    end
end
