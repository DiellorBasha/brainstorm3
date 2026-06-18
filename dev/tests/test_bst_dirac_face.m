function test_bst_dirac_face()
% Phase 3: bst_dirac face branch (Transform / Reconstruct). Validated by an INDEPENDENT
% re-derivation from the 'Dirac-Face' eigen + operator nodes (exact, not a round-trip --
% reconstruct drops the quaternion w-part, so T/R are not mutual inverses).
% Author: Diellor Basha, 2026
    nFail = 0;  K = 200;  Tau = 0.5;
    df = 'Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat';
    [sStudy,~] = bst_get('DataFile', df);
    BaseHM = in_bst_headmodel([fileparts(df) '/headmodel_surf_os_meg.mat'], 0);
    ChanMat = in_bst_channel(sStudy.Channel(1).FileName);
    iMEG = find(strcmpi({ChanMat.Channel.Type},'MEG'));
    SurfaceFile = BaseHM.SurfaceFile;
    [Lf, FG] = bst_face_leadfield(SurfaceFile, ChanMat.Channel(iMEG), BaseHM.Param(iMEG), 'Mode','unconstrained');
    nCh = numel(iMEG);  nF = size(Lf,2)/3;
    HMf = struct('Gain',Lf, 'SurfaceFile',SurfaceFile, 'HeadModelType','surface', ...
                 'isFaceBased',1, 'FaceBasis','dirac', 'GridLoc',FG.Centroids);  % validate the Dirac-Face machinery

    Comp = bst_dirac(HMf, 'nModes', K, 'Tau', Tau);
    nFail = nFail + chk('Transform Gain [nCh x 2K]', isequal(size(Comp.Gain),[nCh 2*K]));
    nFail = nFail + chk('isFaceBased flag set', isfield(Comp,'isFaceBased') && Comp.isFaceBased==1);
    nFail = nFail + chk('HemiGlobalFaces present', isfield(Comp,'HemiGlobalFaces') && numel(Comp.HemiGlobalFaces)==2);

    % --- independent re-derivation of the Transform from the nodes ---
    Eig = tess_eigen(SurfaceFile, 'Dirac-Face', 'K', K, 'Tau', Tau);
    Op  = tess_operators(SurfaceFile, 'Dirac-Face', 'Tau', Tau, 'NoSave', true);
    Lref = zeros(nCh, 2*K);
    for hh = 1:2
        fH = Eig.GlobalFaces{hh}(:);  nFh = numel(fH);
        Phi = double(Eig.Phi{hh}(:,1:K));  B = Op.Mass{hh};
        Psi = zeros(4*nFh, nCh);
        Psi(2:4:end,:) = Lf(:, (fH-1)*3+1).';
        Psi(3:4:end,:) = Lf(:, (fH-1)*3+2).';
        Psi(4:4:end,:) = Lf(:, (fH-1)*3+3).';
        Lref(:, (hh-1)*K + (1:K)) = Psi' * (B * Phi);
    end
    relT = norm(Comp.Gain - Lref,'fro') / max(norm(Lref,'fro'),eps);
    fprintf('  Transform vs independent re-derivation: rel err = %.2e\n', relT);
    nFail = nFail + chk('Transform matches independent re-derivation', relT < 1e-10);

    % --- Reconstruct: independent expansion imag(Phi*c) scattered to faces ---
    rng(3); C = randn(2, 2*K);                          % 2 coefficient rows
    J = bst_dirac(Comp, 'Reconstruct', C);
    nFail = nFail + chk('Reconstruct [m x 3nF]', isequal(size(J),[2 3*nF]));
    Jref = zeros(2, 3*nF);
    for hh = 1:2
        fH = Eig.GlobalFaces{hh}(:);  Phi = double(Eig.Phi{hh}(:,1:K));
        cols = ((hh-1)*K + (1:K));
        R = Phi * C(:, cols).';                          % [4nFh x 2]
        Jref(:, (fH-1)*3+1) = R(2:4:end,:).';
        Jref(:, (fH-1)*3+2) = R(3:4:end,:).';
        Jref(:, (fH-1)*3+3) = R(4:4:end,:).';
    end
    relR = norm(J - Jref,'fro') / max(norm(Jref,'fro'),eps);
    fprintf('  Reconstruct vs independent re-derivation: rel err = %.2e\n', relR);
    nFail = nFail + chk('Reconstruct matches independent re-derivation', relR < 1e-10);

    % --- observability preservation (REPORTED finding, not gated here): how much of the
    %     observable subspace the K modes retain. The face Dirac modes are not smooth-ordered
    %     (wide-root spectrum), so they capture observability far less efficiently than the
    %     vertex Dirac -- quantified end-to-end in the Phase-4 benchmark. Gate only on the
    %     implementation being sane (rank positive, <= ceiling). ---
    er = @(M) sum(svd(M) > 1e-3*max(svd(M)));
    rL = er(Lf);  rM = er(Comp.Gain);
    fprintf('  observability: face leadfield rank=%d   mode-forward (2K=%d) rank=%d  (%.0f%% retained @K=%d)\n', ...
        rL, 2*K, rM, 100*rM/rL, K);
    nFail = nFail + chk('mode-forward rank in (0, ceiling]', rM > 0 && rM <= rL + 1);

    fprintf('\n==== test_bst_dirac_face: %d failed ====\n', nFail);
    if nFail > 0, error('test_bst_dirac_face FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
