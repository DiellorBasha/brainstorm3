function test_phase_recovery_v2
% TEST_PHASE_RECOVERY_V2  Signed normal component + eigenmode reconstruction.
%
% Corrected pipeline (v1 used the norm, which was wrong):
%
%   J(x,t)  [3D unconstrained MNE]
%       ↓  project onto outward normal n̂(x)
%   s(x,t)  [signed scalar, oscillates through zero]
%       ↓  Hilbert per vertex
%   A(x,t), Φ(x,t)  [instantaneous amplitude and phase at ω₀]
%
% The norm is NOT used — it rectifies the signal to 2ω and destroys phase.
%
% Sign ambiguity source: the MNE inverse can assign wrong sign at individual
% vertices, especially near sulcal fundi (ill-conditioned).  This is an
% INVERSE PROBLEM artefact, not a geometric one — FreeSurfer normals are
% consistently outward so the geometric sign is well-defined.
%
% Eigenmode hypothesis (to be tested here):
%   Reconstructing s via scalar LBO eigenmodes forces spatial smoothness.
%   Isolated MNE sign errors — which require an abrupt local sign flip —
%   cannot survive the smooth-mode projection, which averages over nearby
%   vertices.  Multi-vertex sign errors persist but isolated ones are
%   suppressed.
%
% Tests
% -----
%   1. Sign error census — how many vertices get a wrong sign from the MNE?
%      Do the errors cluster (sulcal) or scatter (random)?
%   2. Phase error: direct signed s vs eigenmode reconstruction
%      at ALL vertices; focus on the sign-error subset
%   3. Spatial smoothness: phase gradient magnitude (lower = smoother)
%   4. Wave direction: gradient correlation with planted wave
%
% Synthetic source: smooth traveling wave across a gyral patch of the LH,
% s(x,t) = A(x) · cos(ω₀t + k · (x·ê))  where ê is the propagation direction.

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);
addpath(fullfile(repoRoot,'dev','benchmarks','sign_ambiguity'));
if ~brainstorm('status'), brainstorm nogui; end

%% ── Setup ─────────────────────────────────────────────────────────────────

TessMat = in_tess_bst('Subject01/tess_cortex_mid_low.mat');
nV      = size(TessMat.Vertices,1);
n_hat   = TessMat.VertNormals;     % outward normals
Vtx_mm  = TessMat.Vertices * 1000;

hm = load(file_fullpath('Subject01/S01_AEF_20131218_01_notch/headmodel_surf_os_meg.mat'),'Gain');
K  = load(file_fullpath('Subject01/S01_AEF_20131218_01_notch/results_MN_MEG_KERNEL_260605_0111.mat'), ...
          'ImagingKernel','GoodChannel');
Gain = hm.Gain(K.GoodChannel,:);
W    = K.ImagingKernel;

NC = load(file_fullpath('Subject01/S01_AEF_20131218_01_notch/noisecov_full.mat'),'NoiseCov');
noise_std = sqrt(max(diag(NC.NoiseCov(K.GoodChannel,K.GoodChannel)),0));

Fs   = 600;
t    = (0:Fs-1)/Fs;
f0   = 10;
om   = 2*pi*f0;

%% ── Scalar LBO eigenmodes (left hemisphere) ───────────────────────────────

clear mex
h_lbo = nxr_compute('create', TessMat.Vertices, double(TessMat.Faces));
dec   = nxr_compute('assembleDECOperators', h_lbo);
L_sc  = dec.d0' * dec.hodge1 * dec.d0;
M_sc  = dec.hodge0;
nxr_compute('destroy', h_lbo);

[~,lH_v] = tess_hemisplit(TessMat); lH_v = lH_v(:);
inL = false(nV,1); inL(lH_v)=true;

nModes = 20;
[Phi, Lam] = eigs(L_sc(lH_v,lH_v), M_sc(lH_v,lH_v), nModes, 'smallestabs');
lam = real(diag(Lam)); [lam,so]=sort(lam,'ascend'); Phi=Phi(:,so);
Phi_nm  = Phi(:,2:end);   % skip DC (λ=0)
M_diag  = full(diag(M_sc(lH_v,lH_v)));

fprintf('LBO modes: %d  (λ range [%.0f … %.0f])\n\n', nModes-1, lam(2), lam(end));

%% ── Synthetic traveling wave ──────────────────────────────────────────────
% Active vertices: gyral LH patch (SulciMap==0) for clean ground truth.
% Propagation: anterior-posterior (−Y direction).

gyral_lh = lH_v(TessMat.SulciMap(lH_v)==0);
e_prop   = [0 -1 0];           % propagation direction
v_ms     = 3.0;                 % wave speed m/s
k_rad    = om / (v_ms*1000);    % wavenumber rad/mm
A_src    = 1e-9;                % 1 nAm amplitude
SNR_dB   = 10;

phi_true = k_rad * (Vtx_mm(gyral_lh,:) * e_prop');   % [nPatch x 1]

% Synthetic source matrix [3nV x nT]
nT = length(t);
J_synth = zeros(3*nV, nT);
for ki = 1:numel(gyral_lh)
    v  = gyral_lh(ki);
    J_synth((3*v-2):(3*v), :) = n_hat(v,:)' * (A_src * cos(om*t + phi_true(ki)));
end

% Forward model + noise
sensors = Gain * J_synth;
sig_rms = sqrt(mean(sensors(:).^2));
nse_rms = sig_rms / 10^(SNR_dB/20);
sensors = sensors + nse_rms * (noise_std./max(noise_std,eps)) .* randn(size(sensors));

% Inverse: J_rec [nV x 3 x nT], signed normal component s [nV x nT]
J_rec3 = reshape(W * sensors, 3, nV, []);
J_rec  = permute(J_rec3, [2 1 3]);                      % [nV x 3 x nT]
s_rec  = squeeze(sum(J_rec .* n_hat, 2));               % [nV x nT]  signed

%% ── Analytic signal of s (direct approach) ────────────────────────────────
% Hilbert per vertex along time axis.
A_dir  = abs(s_rec + 1i*imag(hilbert(s_rec')')); % [nV x nT]
Ph_dir = angle(s_rec + 1i*imag(hilbert(s_rec')'));

%% ── Eigenmode reconstruction approach ────────────────────────────────────
% 1. Project s_rec onto LBO modes → θ_k(t) [nModes-1 x nT]
% 2. Reconstruct: s_eig = Φ_nm * θ  [nLH x nT]
% 3. Hilbert → A_eig, Ph_eig

theta    = (Phi_nm .* M_diag)' * s_rec(lH_v,:);    % [nModes-1 x nT]
s_eig    = Phi_nm * theta;                           % [nLH x nT]  smooth
A_eig_lh = abs(s_eig + 1i*imag(hilbert(s_eig')')); % [nLH x nT]
Ph_eig_lh= angle(s_eig + 1i*imag(hilbert(s_eig')'));

%% ── Amplitude threshold: active vertices ─────────────────────────────────
amp_mean  = mean(A_dir(gyral_lh,:), 2);
amp_th    = 0.10 * max(amp_mean);
active    = gyral_lh(amp_mean > amp_th);   % global indices
active_lh = find(ismember(lH_v, active));  % local LH indices

fprintf('Active vertices (>10%% peak, gyral LH): %d / %d\n', numel(active), numel(gyral_lh));

%% ── TEST 1 — Sign error census ────────────────────────────────────────────
fprintf('\n─── TEST 1: MNE sign error census ───\n');

% Ground truth sign at t_peak (moment of maximum source amplitude)
[~, t_pk] = max(mean(A_dir(active,:),1));
s_true_pk = A_src * cos(om*t(t_pk) + phi_true(ismember(gyral_lh, active)));
s_rec_pk  = s_rec(active, t_pk);

sign_correct = sign(s_rec_pk) == sign(s_true_pk);
frac_correct = mean(sign_correct);
fprintf('  Sign correct at active vertices: %.1f%%\n', 100*frac_correct);
fprintf('  Sign errors: %d / %d vertices\n', sum(~sign_correct), numel(active));

% Do errors cluster at sulcal-adjacent vertices?
is_sulcal_adj = false(numel(active),1);
VC = TessMat.VertConn;
for ki=1:numel(active)
    nb = find(VC(active(ki),:)>0);
    is_sulcal_adj(ki) = any(TessMat.SulciMap(nb)>0);
end
fprintf('  Sign errors at sulcal-adjacent vertices: %.1f%% of errors\n', ...
    100*mean(~sign_correct & is_sulcal_adj) / max(mean(~sign_correct),eps));

%% ── TEST 2 — Phase error: direct vs eigenmode ─────────────────────────────
fprintf('\n─── TEST 2: Phase error at active vertices ───\n');

Ph_true  = om*t(t_pk) + phi_true(ismember(gyral_lh,active));  % [nActive x 1]
ph_dir_pk  = Ph_dir(active, t_pk);
ph_eig_pk  = Ph_eig_lh(active_lh, t_pk);

err_dir = abs(angle(exp(1i*(ph_dir_pk  - Ph_true))));
err_eig = abs(angle(exp(1i*(ph_eig_pk  - Ph_true))));

fprintf('  ALL active vertices:\n');
fprintf('    Direct signed s:     mean=%.4f rad  median=%.4f rad\n', mean(err_dir), median(err_dir));
fprintf('    Eigenmode recon:     mean=%.4f rad  median=%.4f rad\n', mean(err_eig), median(err_eig));

% Focus on sign-error vertices
if any(~sign_correct)
    ae_dir_bad = err_dir(~sign_correct);
    ae_eig_bad = err_eig(~sign_correct);
    fprintf('  SIGN-ERROR vertices (%d):\n', sum(~sign_correct));
    fprintf('    Direct signed s:     mean=%.4f rad  (expect ≈π=%.4f)\n', mean(ae_dir_bad), pi);
    fprintf('    Eigenmode recon:     mean=%.4f rad\n', mean(ae_eig_bad));
end

%% ── TEST 3 — Spatial smoothness of the phase field ───────────────────────
fprintf('\n─── TEST 3: Spatial smoothness (mean |Δφ| per mesh edge) ───\n');

[ei,ej] = find(triu(VC(active,active)>0,1));
ga=active(ei); gb=active(ej);
% Build local-LH index map:  global vertex v → position in lH_v vector
lh_map = zeros(nV,1); lh_map(lH_v) = 1:numel(lH_v);

dphi_dir  = abs(angle(exp(1i*(Ph_dir(ga,t_pk)     - Ph_dir(gb,t_pk)))));
dphi_eig  = abs(angle(exp(1i*(Ph_eig_lh(lh_map(ga),t_pk) - ...
                               Ph_eig_lh(lh_map(gb),t_pk)))));
act_map   = zeros(nV,1); act_map(active) = 1:numel(active);
dphi_true = abs(angle(exp(1i*(Ph_true(act_map(ga)) - Ph_true(act_map(gb))))));

fprintf('  Ground truth:        mean edge Δφ = %.4f rad\n', mean(dphi_true));
fprintf('  Direct signed s:     mean edge Δφ = %.4f rad\n', mean(dphi_dir));
fprintf('  Eigenmode recon:     mean edge Δφ = %.4f rad\n', mean(dphi_eig));

%% ── TEST 4 — Wave direction recovery ─────────────────────────────────────
fprintf('\n─── TEST 4: Wave direction recovery ───\n');

edge_v = Vtx_mm(gb,:) - Vtx_mm(ga,:);
true_grad_edge = k_rad * (edge_v * e_prop');   % true phase diff per edge

corr_dir  = corr(true_grad_edge, Ph_dir(gb,t_pk)           - Ph_dir(ga,t_pk));
corr_eig  = corr(true_grad_edge, Ph_eig_lh(lh_map(gb),t_pk) - ...
                                  Ph_eig_lh(lh_map(ga),t_pk));

fprintf('  Gradient correlation with true wave direction:\n');
fprintf('    Direct signed s:   r = %.4f\n', corr_dir);
fprintf('    Eigenmode recon:   r = %.4f\n', corr_eig);

fprintf('\n=======================================================\n');
fprintf(' TEST_PHASE_RECOVERY_V2 COMPLETE\n');
fprintf('=======================================================\n');

end
