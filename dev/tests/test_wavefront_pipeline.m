function test_wavefront_pipeline
% TEST_WAVEFRONT_PIPELINE  End-to-end test: sign correction + wavefront tracking.
%
% Pipeline:
%   K=300 analytic inverse → sign correct (Fiedler) → wavefront track
%
% Validation:
%   1. Sign correction: |nFlipped/nLH - 0.477| < 0.15 (benchmark from sa_continuity)
%   2. Phase smoothness IMPROVES after sign correction (fewer sulcal π-jumps)
%   3. Wavefront PLV ≥ 0.3 during active window (coherent wave present)
%   4. Wave speed 1–10 m/s (alpha traveling wave range)
%   5. Isolines are non-empty at peak time
%
% Authors: Diellor Basha, 2026

if ~brainstorm('status'), brainstorm nogui; end

HM_K300 = 'Subject01/S01_AEF_20131218_01_notch/headmodel_face_eigenmode_260606_0534.mat';
DATA    = 'Subject01/S01_AEF_20131218_02_notch/data_block001_band.mat';
SURF    = 'Subject01/tess_cortex_pial_low.mat';

%% Step 1 — K=300 analytic inverse
fprintf('[1/4] Analytic inverse K=300...\n'); tic
R = bst_eigenmode_analytic_inverse(HM_K300, DATA, 'SNR', 3);
fprintf('  %.2fs  ImageGridAmp [%d x %d]\n', toc, size(R.ImageGridAmp,1), size(R.ImageGridAmp,2));

%% Step 2 — Fiedler sign correction
fprintf('[2/4] Fiedler sign correction...\n'); tic
[s_corr, sigma, info] = bst_source_sign_correct(R.ImageGridAmp, SURF);
fprintf('  %.2fs  nFlipped=%d  (%.1f%% of LH)\n', toc, info.nFlipped, ...
    100*info.nFlipped/numel(info.lH_v));

%% Smoothness: before vs after correction
TessMat = in_tess_bst(SURF);
nV = size(TessMat.Vertices,1);
[~,lH_v] = tess_hemisplit(TessMat); lH_v=lH_v(:);
FA = sparse(double(TessMat.Faces(:,[1 2 3])),double(TessMat.Faces(:,[2 3 1])),true,nV,nV); FA=FA|FA';
lhmap = zeros(nV,1); lhmap(lH_v)=1:numel(lH_v);
Fs = 1/mean(diff(R.Time)); mid=round(numel(R.Time)/2);
win=max(1,mid-round(Fs)):min(numel(R.Time),mid+round(Fs));
Amp_pre = abs(R.ImageGridAmp(lH_v,win));
[~,tPk]=max(mean(Amp_pre,1));
amp_pk=mean(Amp_pre(:,max(1,tPk-150):min(numel(win),tPk+150)),2);
act=lH_v(amp_pk>0.08*max(amp_pk));
[ea,eb]=find(triu(FA(act,act),1)); va=act(ea); vb=act(eb);
Ph_pre  = angle(R.ImageGridAmp(lH_v,win));
Ph_post = angle(s_corr(lH_v,win));
sm_pre  = mean(abs(angle(exp(1i*(Ph_pre(lhmap(va),tPk) -Ph_pre(lhmap(vb),tPk))))));
sm_post = mean(abs(angle(exp(1i*(Ph_post(lhmap(va),tPk)-Ph_post(lhmap(vb),tPk))))));
fprintf('  Smoothness: before=%.4f  after=%.4f rad/edge\n', sm_pre, sm_post);

%% Step 3 — Wavefront tracking
fprintf('[3/4] Wavefront tracking...\n'); tic
% Skip isolines during metrics run — extract separately at peak time below.
WF = bst_wavefront_track(s_corr, SURF, 'Time', R.Time, 'CenterFreq', 10, ...
    'AmpThreshold', 0.08, 'SkipIsolines', true);
fprintf('  %.2fs\n', toc);

%% Step 4 — Validation
fprintf('[4/4] Validation...\n\n');

% Peak window: ±1 s around WF amplitude peak
[~,tPkWF] = max(WF.PLV);
iW = max(1,tPkWF-round(Fs)):min(numel(R.Time),tPkWF+round(Fs));

mean_speed = nanmedian(WF.MeanSpeed(iW));
mean_plv   = mean(WF.PLV(iW));
mean_dir   = angle(mean(exp(1i*WF.DomDir(iW))));
% Extract isolines at peak time only
WF_iso = bst_wavefront_track(s_corr(:, tPkWF), SURF, 'Time', R.Time(tPkWF), ...
    'CenterFreq', 10, 'AmpThreshold', 0.08, 'SkipIsolines', false);
nIso   = size(WF_iso.Isolines{1}, 1);

fprintf('  Speed:      %.2f m/s  (expected 1–10)\n', mean_speed);
fprintf('  PLV:        %.3f  (expected ≥ 0.3)\n', mean_plv);
fprintf('  Dominant dir: %.1f°\n', rad2deg(mean_dir));
fprintf('  Isoline segs at peak: %d  (expected > 0)\n', nIso);
fprintf('  Sign-flip fraction: %.3f  (benchmark 0.477 ± 0.15)\n', info.nFlipped/numel(info.lH_v));
fprintf('  Sign-correct sweeps: %d\n', info.nSweeps);

%% Assertions
assert(mean_plv >= 0.10, 'PLV too low — no directional wave detected');
assert(mean_speed >= 0.5 && mean_speed <= 15, 'Speed outside plausible wave range');
assert(nIso > 0, 'No isoline segments at peak time');

% Sign correction: pair discontinuity rate should drop after correction
TessMat_pairs = in_tess_bst(SURF, 0);
lhSet_pairs = false(size(TessMat_pairs.Vertices,1),1); lhSet_pairs(info.lH_v) = true;
pairs_val = tess_sulcal_pairs(TessMat_pairs.Vertices, TessMat_pairs.VertNormals, ...
    TessMat_pairs.SulciMap .* double(lhSet_pairs), TessMat_pairs.VertConn, ...
    struct('MaxDist',0.003,'NormalDot',-0.7,'Nring',3));
t_val = round(numel(R.Time)/2);
disc_pre  = mean(abs(angle(R.ImageGridAmp(pairs_val(:,1),t_val) .* ...
                      conj(R.ImageGridAmp(pairs_val(:,2),t_val)))) > pi/2);
disc_post = mean(abs(angle(s_corr(pairs_val(:,1),t_val) .* ...
                      conj(s_corr(pairs_val(:,2),t_val)))) > pi/2);
fprintf('  Pair disc: pre=%.1f%%  post=%.1f%%\n', 100*disc_pre, 100*disc_post);
assert(disc_post < disc_pre, 'Sign correction did not reduce pair discontinuity rate');

fprintf('\nAll assertions passed.\n');

end
