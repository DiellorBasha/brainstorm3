function test_phase_recovery_v3
% TEST_PHASE_RECOVERY_V3  Phase recovery using the Eigen-MNE kernel.
%
% This test uses the eigenmode-parameterised inverse (Eigen-MNE) built into
% Brainstorm.  The source model is:
%
%   s(x,t) = Σ_k θ_k(t) · ψ_k(x)       [smooth scalar field]
%
% θ_k(t) are estimated from sensors via the composed eigenmode leadfield.
% The ImagingKernel already encodes the full Φ · M_eigenspace · L_tilde^-1
% chain — applying it to sensor data returns a smooth signed scalar at
% every vertex directly.
%
% This is the correct parameterisation for wave detection:
%   • Signed scalar — Hilbert gives meaningful phase at ω₀
%   • Smooth by construction — per-vertex sign flips impossible
%   • No post-processing needed — the correct basis was used from the start
%
% Compares: Eigen-MNE  vs  standard constrained MNE  vs  ground truth.
%
% Tests
% -----
%   1. Sign error census  — fraction of vertices with wrong sign vs constrained MNE
%   2. Phase error        — mean |Φ_rec − Φ_true| at active vertices
%   3. Spatial smoothness — mean |ΔΦ| per mesh edge
%   4. Wave direction     — gradient correlation with planted wave
%
% Same synthetic wave as v2: traveling alpha (10 Hz) across gyral LH patch,
% propagation along −Y, speed 3 m/s, SNR 10 dB.

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);
if ~brainstorm('status'), brainstorm nogui; end

%% ── Load kernels and surface ──────────────────────────────────────────────

sStudies = bst_get('ProtocolStudies');
s6       = [sStudies.Study];  s6 = s6(6);
rc       = {s6.Result.Comment};
rf       = {s6.Result.FileName};

iEig   = find(strcmp(rc, 'Eigen-MNE: MEG 2018'), 1);
iConst = find(strcmp(rc, 'MN: MEG(Constr) 2018'), 1);
assert(~isempty(iEig),   'Eigen-MNE kernel not found in Study 6.');
assert(~isempty(iConst), 'Constrained MNE kernel not found in Study 6.');

RE = in_bst_results(rf{iEig},   0);   % Eigen-MNE
RC = in_bst_results(rf{iConst}, 0);   % Constrained MNE

K_eig   = RE.ImagingKernel;   % [nV x nCh]
K_const = RC.ImagingKernel;   % [nV x nCh]
nV      = size(K_eig, 1);

fprintf('Eigen-MNE   kernel: %s   nModes=%d\n', mat2str(size(K_eig)),   RE.nModes);
fprintf('Constrained kernel: %s\n', mat2str(size(K_const)));

TessMat = in_tess_bst(RE.SurfaceFile);
n_hat   = TessMat.VertNormals;
Vtx_mm  = TessMat.Vertices * 1000;

% Channel mapping — kernels use GoodChannel subset
hm   = load(file_fullpath('Subject01/S01_AEF_20131218_01_notch/headmodel_surf_os_meg.mat'),'Gain');
K_f  = load(file_fullpath(rf{iEig}), 'GoodChannel');
Gain = hm.Gain(K_f.GoodChannel, :);      % [nCh x 3nV]

NC = load(file_fullpath('Subject01/S01_AEF_20131218_01_notch/noisecov_full.mat'),'NoiseCov');
noise_std = sqrt(max(diag(NC.NoiseCov(K_f.GoodChannel, K_f.GoodChannel)), 0));

%% ── Synthetic traveling wave ──────────────────────────────────────────────

Fs     = 600;
t      = (0:Fs-1) / Fs;
f0     = 10;
om     = 2*pi*f0;
A_src  = 1e-9;
SNR_dB = 10;
e_prop = [0 -1 0];
v_ms   = 3.0;
k_rad  = om / (v_ms * 1000);

[~, lH_v] = tess_hemisplit(TessMat); lH_v = lH_v(:);
inL        = false(nV,1); inL(lH_v) = true;
gyral_lh   = lH_v(TessMat.SulciMap(lH_v) == 0);

phi_true   = k_rad * (Vtx_mm(gyral_lh, :) * e_prop');   % [nPatch x 1]
nT         = length(t);

% Build synthetic source matrix [3nV x nT] — constrained (along n̂)
J_synth = zeros(3*nV, nT);
for ki = 1:numel(gyral_lh)
    v = gyral_lh(ki);
    J_synth((3*v-2):(3*v), :) = n_hat(v,:)' * (A_src * cos(om*t + phi_true(ki)));
end

% Forward + noise
sensors  = Gain * J_synth;
sig_rms  = sqrt(mean(sensors(:).^2));
nse_rms  = sig_rms / 10^(SNR_dB/20);
sensors  = sensors + nse_rms * (noise_std./max(noise_std,eps)) .* randn(size(sensors));

%% ── Apply kernels ─────────────────────────────────────────────────────────

s_eig   = K_eig   * sensors;   % [nV x nT]  Eigen-MNE
s_const = K_const * sensors;   % [nV x nT]  Constrained MNE

%% ── Analytic signal (Hilbert) ─────────────────────────────────────────────

A_eig   = abs  (s_eig   + 1i*imag(hilbert(s_eig'  )')); % [nV x nT]
Ph_eig  = angle(s_eig   + 1i*imag(hilbert(s_eig'  )'));
A_con   = abs  (s_const + 1i*imag(hilbert(s_const' )')); % [nV x nT]
Ph_con  = angle(s_const + 1i*imag(hilbert(s_const' )'));

%% ── Amplitude mask and peak time ─────────────────────────────────────────

amp_eig_mean = mean(A_eig(gyral_lh,:), 2);
amp_con_mean = mean(A_con(gyral_lh,:), 2);
thresh_eig   = 0.10 * max(amp_eig_mean);
thresh_con   = 0.10 * max(amp_con_mean);

act_eig = gyral_lh(amp_eig_mean > thresh_eig);
act_con = gyral_lh(amp_con_mean > thresh_con);

[~, t_pk_e] = max(mean(A_eig(act_eig,:), 1));
[~, t_pk_c] = max(mean(A_con(act_con,:), 1));

fprintf('\nActive vertices (10%% threshold):\n');
fprintf('  Eigen-MNE:   %d / %d gyral\n', numel(act_eig), numel(gyral_lh));
fprintf('  Constrained: %d / %d gyral\n', numel(act_con), numel(gyral_lh));

%% ── TEST 1 — Sign error census ────────────────────────────────────────────
fprintf('\n─── TEST 1: Sign error census ───\n');

act_both    = intersect(act_eig, act_con);
phi_map     = zeros(nV,1); phi_map(gyral_lh) = phi_true;
s_true_pk_e = A_src * cos(om*t(t_pk_e) + phi_map(act_both));
s_true_pk_c = A_src * cos(om*t(t_pk_c) + phi_map(act_both));

se_eig  = sign(s_eig  (act_both, t_pk_e)) == sign(s_true_pk_e);
se_con  = sign(s_const(act_both, t_pk_c)) == sign(s_true_pk_c);

fprintf('  Vertices evaluated (both methods active): %d\n', numel(act_both));
fprintf('  Eigen-MNE   correct sign: %.1f%%  (%d errors)\n', ...
    100*mean(se_eig),  sum(~se_eig));
fprintf('  Constrained correct sign: %.1f%%  (%d errors)\n', ...
    100*mean(se_con),  sum(~se_con));

%% ── TEST 2 — Phase error ──────────────────────────────────────────────────
fprintf('\n─── TEST 2: Phase error at active vertices ───\n');

Ph_true_e = om*t(t_pk_e) + phi_map(act_eig);
Ph_true_c = om*t(t_pk_c) + phi_map(act_con);

err_eig = abs(angle(exp(1i*(Ph_eig(act_eig, t_pk_e) - Ph_true_e))));
err_con = abs(angle(exp(1i*(Ph_con(act_con, t_pk_c) - Ph_true_c))));

fprintf('  Eigen-MNE:   mean=%.4f rad  median=%.4f rad\n', mean(err_eig), median(err_eig));
fprintf('  Constrained: mean=%.4f rad  median=%.4f rad\n', mean(err_con), median(err_con));

% Breakdown: sign-error subset vs sign-correct subset (Eigen-MNE)
act_both_eig = intersect(act_eig, act_both);
se_eig_local = se_eig(ismember(act_both, act_both_eig));
err_eig_local = err_eig(ismember(act_eig, act_both_eig));
if any(~se_eig_local)
    fprintf('  Eigen-MNE at SIGN-ERROR vertices: mean=%.4f rad\n', ...
        mean(err_eig_local(~se_eig_local)));
    fprintf('  Eigen-MNE at SIGN-OK   vertices:  mean=%.4f rad\n', ...
        mean(err_eig_local(se_eig_local)));
end

%% ── TEST 3 — Spatial smoothness ──────────────────────────────────────────
fprintf('\n─── TEST 3: Spatial smoothness (mean |ΔΦ| per mesh edge) ───\n');

VC = TessMat.VertConn;
[ei,ej] = find(triu(VC(act_eig, act_eig) > 0, 1));
ga = act_eig(ei); gb = act_eig(ej);

dphi_true = abs(angle(exp(1i*(Ph_true_e(ei) - Ph_true_e(ej)))));
dphi_eig  = abs(angle(exp(1i*(Ph_eig(ga,t_pk_e) - Ph_eig(gb,t_pk_e)))));

% For constrained MNE use active vertices subset
[ei_c,ej_c] = find(triu(VC(act_con,act_con)>0,1));
ga_c=act_con(ei_c); gb_c=act_con(ej_c);
dphi_con = abs(angle(exp(1i*(Ph_con(ga_c,t_pk_c)-Ph_con(gb_c,t_pk_c)))));

fprintf('  Ground truth:  mean edge ΔΦ = %.4f rad\n', mean(dphi_true));
fprintf('  Eigen-MNE:     mean edge ΔΦ = %.4f rad\n', mean(dphi_eig));
fprintf('  Constrained:   mean edge ΔΦ = %.4f rad\n', mean(dphi_con));

%% ── TEST 4 — Wave direction recovery ─────────────────────────────────────
fprintf('\n─── TEST 4: Wave direction recovery (gradient correlation) ───\n');

edge_v         = Vtx_mm(gb,:) - Vtx_mm(ga,:);
true_grad      = k_rad * (edge_v * e_prop');
rec_grad_eig   = Ph_eig(gb,t_pk_e) - Ph_eig(ga,t_pk_e);

edge_v_c       = Vtx_mm(gb_c,:) - Vtx_mm(ga_c,:);
true_grad_c    = k_rad * (edge_v_c * e_prop');
rec_grad_con   = Ph_con(gb_c,t_pk_c) - Ph_con(ga_c,t_pk_c);

corr_eig = corr(true_grad, rec_grad_eig);
corr_con = corr(true_grad_c, rec_grad_con);

fprintf('  Eigen-MNE:   r = %.4f\n', corr_eig);
fprintf('  Constrained: r = %.4f\n', corr_con);

fprintf('\n=======================================================\n');
fprintf(' TEST_PHASE_RECOVERY_V3 COMPLETE\n');
fprintf('=======================================================\n');

end
