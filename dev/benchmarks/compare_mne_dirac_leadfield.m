function compare_mne_dirac_leadfield()
% COMPARE_MNE_DIRAC_LEADFIELD: Vector-by-vector comparison of the per-vertex
% leadfield for two example MEG sensors, standard OS-MEG vs Dirac eigenmode.
%
% Standalone analysis (plain .mat loads; no Brainstorm DB needed). Reuses the
% pure helper bst_dirac_eigenmode_field for the Dirac reconstruction.
%
% For each sensor it builds:
%   gStd(v,:)  = raw unconstrained leadfield 3-vector at vertex v   [nVert x 3]
%   gDir(v,:)  = Dirac eigenmode leadfield reconstructed to vertex v [nVert x 3]
% and reports, per vertex: magnitude ratio |gDir|/|gStd| and the angle between
% the two vectors; plus overall energy capture and correlation. Saves a PNG per
% sensor with the distributions.
%
% Authors: Diellor Basha, 2026

    repo     = '/Users/diellorbasha/workspace/research/code/brainstorm3';
    addpath(fullfile(repo,'toolbox','forward'));   % bst_dirac_eigenmode_field
    protocol = '/Users/diellorbasha/workspace/library/datasets/brainstorm_db/tmp_aggregate/TutorialAuditory';
    studyDir = fullfile(protocol,'data','Subject01','S01_AEF_20131218_01_notch');
    outDir   = fullfile(repo,'dev','benchmarks');

    HMstd = load(fullfile(studyDir,'headmodel_surf_os_meg.mat'));            % .Gain [nCh x 3nVert]
    HMdir = load(fullfile(studyDir,'headmodel_dirac_eigenmode_260610_1108.mat')); % .Gain [nCh x 2K] + maps
    S     = load(fullfile(protocol,'anat', HMstd.SurfaceFile));             % .Vertices, .DiracEigen
    nVert = size(S.Vertices,1);
    fprintf('nVert=%d | std Gain %s | dirac Gain %s | nModes=%d\n', ...
        nVert, mat2str(size(HMstd.Gain)), mat2str(size(HMdir.Gain)), HMdir.nModes);

    % --- resolve sensor channel indices by name ---
    chd = dir(fullfile(studyDir,'channel_*.mat'));
    Ch  = load(fullfile(chd(1).folder, chd(1).name));
    names = {Ch.Channel.Name};
    want = {'MLC11','MZP01'};   % MLC11 = the 1st MEG channel (user wrote "MCL11"); MZP01 = 274th MEG
    for k = 1:numel(want)
        idx(k) = find(strcmpi(names, want{k}), 1); %#ok<AGROW>
        assert(~isempty(idx(k)), 'channel %s not found', want{k});
    end
    fprintf('MCL11 -> row %d | MZP01 -> row %d\n', idx(1), idx(2));

    tol = 1e-15;
    for k = 1:numel(want)
        s = idx(k);
        gStd = reshape(HMstd.Gain(s,:), 3, nVert).';                         % [nVert x 3]
        Jrow = bst_dirac_eigenmode_field(S.DiracEigen, HMdir.Gain(s,:), HMdir);
        gDir = reshape(Jrow, 3, nVert).';                                    % [nVert x 3]

        mStd = sqrt(sum(gStd.^2,2));  mDir = sqrt(sum(gDir.^2,2));
        good = mStd > tol & isfinite(mStd) & all(isfinite(gDir),2);          % usable vertices

        % per-vertex angle (deg) and magnitude ratio
        cosang = sum(gStd(good,:).*gDir(good,:),2) ./ max(mStd(good).*mDir(good), tol);
        ang    = acosd(min(1,max(-1,cosang)));
        ratio  = mDir(good) ./ mStd(good);

        % overall metrics (restricted to usable vertices, all 3 components)
        a = gStd(good,:); b = gDir(good,:);
        capture = 1 - norm(b-a,'fro')/norm(a,'fro');
        cc      = corr(a(:), b(:));

        fprintf('\n==== sensor %s (row %d), %d usable vertices ====\n', want{k}, s, nnz(good));
        fprintf('  energy captured by Dirac basis : %.1f%%   (corr %.3f)\n', 100*capture, cc);
        fprintf('  magnitude ratio |Dir|/|Std|    : median %.2f, IQR [%.2f, %.2f]\n', ...
            median(ratio), prctile(ratio,25), prctile(ratio,75));
        fprintf('  vector angle (deg)             : median %.1f, 90th pct %.1f, frac<30deg %.0f%%\n', ...
            median(ang), prctile(ang,90), 100*mean(ang<30));

        % --- figure: distributions ---
        f = figure('Visible','off','Position',[100 100 1100 320]);
        subplot(1,3,1);
        loglog(mStd(good), mDir(good), '.', 'MarkerSize', 3); hold on;
        lim = [min(mStd(good)) max(mStd(good))]; loglog(lim, lim, 'r-');
        xlabel('|Std| leadfield'); ylabel('|Dirac| leadfield'); axis tight; title('per-vertex magnitude');
        subplot(1,3,2);
        histogram(ang, 0:5:180); xlabel('angle Std vs Dirac (deg)'); ylabel('# vertices'); title('direction change');
        subplot(1,3,3);
        histogram(ratio(ratio<3), 0:0.05:3); xlabel('|Dirac|/|Std|'); ylabel('# vertices'); title('magnitude ratio');
        sgtitle(sprintf('%s : OS-MEG (Std) vs Dirac eigenmode leadfield (K=%d)', want{k}, HMdir.nModes));
        png = fullfile(outDir, sprintf('compare_leadfield_%s.png', want{k}));
        saveas(f, png); close(f);
        fprintf('  saved %s\n', png);
    end
end
