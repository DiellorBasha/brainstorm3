function test_phase_recovery
% TEST_PHASE_RECOVERY  Diagnostic tests for temporal phase recovery.
%
% Validates whether the pipeline
%   synthetic source  →  forward model  →  [+ noise]  →  MN inverse  →  phase
% faithfully recovers temporal phase on the cortical surface.
%
% Tests
% -----
%   1. Baseline  — single gyral vertex, noise-free then at SNR 20/10/0 dB
%   2. Sign flip — two opposing sulcal wall vertices, same planted phase:
%                  signed J·n̂ should give ~π error; |J| should give ~0 error
%   3. Delay     — two vertices with known phase delay Δφ; measure recovered Δφ
%   4. Wave map  — smooth spatial phase gradient across a gyral patch;
%                  compare recovered gradient direction and magnitude
%
% Requires: Brainstorm running (nogui), TutorialAuditory protocol loaded,
%           sa_sulcal_walls in path (dev/benchmarks/sign_ambiguity/).
%
% Run: test_phase_recovery()

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);
addpath(fullfile(repoRoot,'dev','benchmarks','sign_ambiguity'));
if ~brainstorm('status'), brainstorm nogui; end

%% ── Common setup ──────────────────────────────────────────────────────────

TessMat = in_tess_bst('Subject01/tess_cortex_mid_low.mat');
nV      = size(TessMat.Vertices,1);
n_hat   = TessMat.VertNormals;
Vtx_mm  = TessMat.Vertices * 1000;

% Headmodel and kernel (272-channel MEG subset)
hm   = load(file_fullpath('Subject01/S01_AEF_20131218_01_notch/headmodel_surf_os_meg.mat'),'Gain');
K    = load(file_fullpath('Subject01/S01_AEF_20131218_01_notch/results_MN_MEG_KERNEL_260605_0111.mat'), ...
           'ImagingKernel','GoodChannel');
Gain = hm.Gain(K.GoodChannel,:);   % [272 x 3nV]
W    = K.ImagingKernel;             % [3nV x 272]

% Noise covariance (MEG channels only, diagonal for simplicity)
NC   = load(file_fullpath('Subject01/S01_AEF_20131218_01_notch/noisecov_full.mat'),'NoiseCov');
noise_std = sqrt(max(diag(NC.NoiseCov(K.GoodChannel, K.GoodChannel)), 0));  % [272 x 1]

% Signal parameters
Fs     = 600;
t      = (0:Fs-1) / Fs;       % 1 s at 600 Hz
omega  = 2*pi*10;              % 10 Hz alpha

% Sulcal wall pairs (LH, anti-aligned normals)
[~,lH_v] = tess_hemisplit(TessMat); lH_v = lH_v(:);
SM_lh    = TessMat.SulciMap;  SM_lh(~ismember((1:nV)',lH_v)) = 0;
sw_pairs = sa_sulcal_walls(TessMat.Vertices, n_hat, SM_lh, TessMat.VertConn, ...
           struct('MaxDist',0.003,'NormalDot',-0.7,'Nring',3));

% High-sensitivity gyral vertex (best round-trip)
gain_norm = sqrt(sum(reshape(sum(Gain.^2,1),3,nV),1));
[~,iGyral] = max(gain_norm);

fprintf('=======================================================\n');
fprintf(' TEST_PHASE_RECOVERY\n');
fprintf(' Gyral test vertex : %d  [%.0f %.0f %.0f] mm\n', iGyral, Vtx_mm(iGyral,:));
fprintf(' Sulcal pairs      : %d\n', size(sw_pairs,1));
fprintf(' Signal            : 10 Hz, 1 s @ %d Hz\n', Fs);
fprintf('=======================================================\n\n');

%% ── Helper: forward→inverse round-trip for a temporal source vector ──────
%  J_patch  : [nVpatch x nTime] amplitudes (each vertex oscillates along n̂)
%  verts    : [nVpatch x 1]     global vertex indices
%  snr_db   : scalar            sensor SNR (Inf = noise-free)
%  Returns recovered J [nV x 3 x nTime] and its normal component J_n [nV x nTime]

    function [J_n_rec, J_norm_rec] = forward_inverse(J_patch, verts, snr_db)
        nT = size(J_patch, 2);
        % Build 3D source matrix [3nV x nT]
        J3 = zeros(3*nV, nT);
        for vi = 1:numel(verts)
            v = verts(vi);
            J3((3*v-2):(3*v), :) = n_hat(v,:)' * J_patch(vi,:);
        end
        % Forward
        sensors = Gain * J3;                    % [272 x nT]
        % Add noise
        if isfinite(snr_db)
            sig_rms = sqrt(mean(sensors(:).^2));
            nse_rms = sig_rms / 10^(snr_db/20);
            sensors = sensors + nse_rms * (noise_std ./ max(noise_std,eps)) .* randn(size(sensors));
        end
        % Inverse
        J_rec3 = reshape(W * sensors, 3, nV, []);   % [3 x nV x nT]
        J_rec  = permute(J_rec3, [2 1 3]);           % [nV x 3 x nT]
        % Normal component and norm
        J_n_rec   = squeeze(sum(J_rec .* n_hat, 2));  % [nV x nT]
        J_norm_rec = squeeze(sqrt(sum(J_rec.^2, 2))); % [nV x nT]
    end

% Phase from signed analytic signal (J_n)
phase_signed = @(x) angle(x + 1i*imag(hilbert(x')')); % [nV x nT] → [nV x nT]

% Phase from norm squared (sign-free, at 2ω) → divide by 2
phase_norm   = @(x) angle(abs(x).^2 + 1i*imag(hilbert((abs(x).^2)')')) / 2;

%% ═══════════════════════════════════════════════════════════════════════════
%% TEST 1 — Baseline phase recovery at a single gyral vertex
%% ═══════════════════════════════════════════════════════════════════════════
fprintf('─── TEST 1: Baseline phase recovery (gyral vertex) ───\n');

phi_true = pi/3;   % planted phase
A_src    = 1e-9;   % 1 nAm
J_patch1 = A_src * cos(omega*t + phi_true);   % [1 x nT]

snr_levels = [Inf 20 10 0];
for snr_db = snr_levels
    [J_n1, ~] = forward_inverse(J_patch1, iGyral, snr_db);
    % Phase at target and at peak-amplitude vertex
    ph_n1 = phase_signed(J_n1);
    [~,iPk1] = max(mean(abs(J_n1).^2, 2));
    ph_at_target = mean(angle(exp(1i*(ph_n1(iGyral,:) - (omega*t + phi_true)))));
    ph_at_peak   = mean(angle(exp(1i*(ph_n1(iPk1,:)   - (omega*t + phi_true)))));
    dist_peak    = norm(Vtx_mm(iPk1,:) - Vtx_mm(iGyral,:));
    lbl = sprintf('%+.0f dB', snr_db); if isinf(snr_db), lbl='noise-free'; end
    fprintf('  SNR=%-10s  |phase_err(target)|=%.4f rad  peak_dist=%.1f mm  |phase_err(peak)|=%.4f rad\n', ...
        lbl, abs(ph_at_target), dist_peak, abs(ph_at_peak));
end
fprintf('\n');

%% ═══════════════════════════════════════════════════════════════════════════
%% TEST 2 — Sign flip at opposing sulcal wall vertices
%% ═══════════════════════════════════════════════════════════════════════════
fprintf('─── TEST 2: Sign flip at sulcal wall pairs ───\n');
fprintf('  Ground truth: same phase (Δφ = 0) planted at both walls.\n');
fprintf('  Prediction: signed J·n̂ → Δφ ≈ π;  norm |J| → Δφ ≈ 0.\n\n');

phi_true2 = pi/4;
A_src2    = 1e-9;

n_test = min(5, size(sw_pairs,1));
dph_signed = zeros(n_test,1);
dph_norm   = zeros(n_test,1);

for k = 1:n_test
    vi = sw_pairs(k,1);  vj = sw_pairs(k,2);
    J_patch2 = A_src2 * [cos(omega*t + phi_true2); cos(omega*t + phi_true2)];
    [J_n2, Jnorm2] = forward_inverse(J_patch2, [vi;vj], Inf);   % noise-free

    ph_n2    = phase_signed(J_n2);
    ph_nm2   = phase_norm(Jnorm2);

    % Phase difference at the two planted vertices
    dph_signed(k) = mean(abs(angle(exp(1i*(ph_n2(vi,:)  - ph_n2(vj,:))))));
    dph_norm(k)   = mean(abs(angle(exp(1i*(ph_nm2(vi,:) - ph_nm2(vj,:))))));
end

fprintf('  Pair   n̂(i)·n̂(j)   Δφ signed    Δφ norm\n');
for k = 1:n_test
    nd = sum(n_hat(sw_pairs(k,1),:).*n_hat(sw_pairs(k,2),:));
    fprintf('  %2d     %+.3f       %.4f rad    %.4f rad\n', k, nd, dph_signed(k), dph_norm(k));
end
fprintf('\n  Mean   |Δφ signed| = %.4f rad  (%.1f° from π: %.2f rad)\n', ...
    mean(dph_signed), rad2deg(abs(mean(dph_signed)-pi)), abs(mean(dph_signed)-pi));
fprintf('  Mean   |Δφ norm|   = %.4f rad\n\n', mean(dph_norm));

%% ═══════════════════════════════════════════════════════════════════════════
%% TEST 3 — Phase delay recovery between two vertices
%% ═══════════════════════════════════════════════════════════════════════════
fprintf('─── TEST 3: Phase delay recovery ───\n');

delays_rad = [0.2 0.5 1.0 pi/2];   % planted delays in radians
A_src3     = 1e-9;

% Three vertex pairs: (a) two nearby gyral, (b) two distant gyral, (c) sulcal
[~, iGyral2] = maxk(gain_norm, 10);
v_pairs3 = [iGyral2(1) iGyral2(2);    % nearby gyral
            iGyral2(1) iGyral2(10);   % distant gyral
            sw_pairs(1,1) sw_pairs(1,2)];  % sulcal wall pair
pair_labels = {'nearby gyral','distant gyral','sulcal pair'};

for pp = 1:3
    vi = v_pairs3(pp,1);  vj = v_pairs3(pp,2);
    dist_mm = norm(Vtx_mm(vi,:)-Vtx_mm(vj,:));
    fprintf('  %s  (dist=%.0f mm):\n', pair_labels{pp}, dist_mm);
    fprintf('    Planted Δφ  |  Recovered (signed)  |  Recovered (norm)  |  Err (signed)  |  Err (norm)\n');
    for dph_true = delays_rad
        J3 = A_src3 * [cos(omega*t); cos(omega*t + dph_true)];
        [J_n3, Jnorm3] = forward_inverse(J3, [vi;vj], 20);  % SNR=20dB

        ph_n3  = phase_signed(J_n3);
        ph_nm3 = phase_norm(Jnorm3);

        dph_rec_s = mean(angle(exp(1i*(ph_n3(vi,:) - ph_n3(vj,:)))));
        dph_rec_n = mean(angle(exp(1i*(ph_nm3(vi,:)- ph_nm3(vj,:)))));
        fprintf('    %+.2f rad    |  %+.2f rad              |  %+.2f rad            |  %.4f rad      |  %.4f rad\n', ...
            dph_true, dph_rec_s, dph_rec_n, abs(dph_rec_s-dph_true), abs(dph_rec_n-dph_true));
    end
    fprintf('\n');
end

%% ═══════════════════════════════════════════════════════════════════════════
%% TEST 4 — Smooth spatial phase gradient recovery
%% ═══════════════════════════════════════════════════════════════════════════
fprintf('─── TEST 4: Spatial phase gradient recovery ───\n');

% Define wave: propagate along posterior-anterior direction (−Y)
e_prop = [0 -1 0];   % propagation direction in mm space
v_ms   = 3.0;        % wave speed m/s = 3 mm/ms
k_rad  = omega / (v_ms * 1000);  % wavenumber rad/mm

% Active patch: LH vertices within gyral crowns (SulciMap == 0)
gyral_lh = lH_v(TessMat.SulciMap(lH_v) == 0);
phi_spatial = k_rad * (Vtx_mm(gyral_lh,:) * e_prop');  % [nPatch x 1]
A_src4 = 1e-9;
J_patch4 = A_src4 * cos(omega*t + phi_spatial);   % [nPatch x nT]

[J_n4, Jnorm4] = forward_inverse(J_patch4, gyral_lh, 10);  % SNR=10dB

% Recover phase at active vertices
ph_n4  = phase_signed(J_n4);
ph_nm4 = phase_norm(Jnorm4);

% Phase error at planted vertices (SNR≥threshold)
amp4 = mean(abs(J_n4(gyral_lh,:)).^2, 2);
thresh4 = 0.10 * max(amp4);
active4 = gyral_lh(amp4 > thresh4);

if ~isempty(active4)
    ph_true4  = k_rad * (Vtx_mm(active4,:) * e_prop') + omega*t(1);
    ph_rec_s4 = ph_n4(active4, 1);
    ph_rec_n4 = ph_nm4(active4, 1);
    err_s = mean(abs(angle(exp(1i*(ph_rec_s4 - ph_true4)))));
    err_n = mean(abs(angle(exp(1i*(ph_rec_n4 - ph_true4)))));
    fprintf('  Active vertices above 10%% threshold: %d / %d\n', numel(active4), numel(gyral_lh));
    fprintf('  Mean phase error (signed J·n̂): %.4f rad\n', err_s);
    fprintf('  Mean phase error (norm |J|):   %.4f rad\n', err_n);

    % Gradient direction test: compute mean gradient of recovered vs true phase
    % Use edge-based finite differences on the mesh
    VC = TessMat.VertConn;
    [ei,ej] = find(triu(VC(active4,active4)>0,1));
    gi=active4(ei); gj=active4(ej);
    edge_v = Vtx_mm(gj,:)-Vtx_mm(gi,:);
    edge_u = edge_v ./ max(sqrt(sum(edge_v.^2,2)),eps);

    % True gradient direction: uniform, = e_prop projected to edges
    dph_true_edge  = k_rad * (edge_v * e_prop');          % true phase diff per edge
    dph_rec_s_edge = angle(exp(1i*(ph_rec_s4(ei)-ph_rec_s4(ej))));
    dph_rec_n_edge = angle(exp(1i*(ph_rec_n4(ei)-ph_rec_n4(ej))));

    % Correlation of recovered gradient with true gradient
    corr_s = corr(dph_true_edge, dph_rec_s_edge);
    corr_n = corr(dph_true_edge, dph_rec_n_edge);
    fprintf('  Gradient correlation, signed: %.4f\n', corr_s);
    fprintf('  Gradient correlation, norm:   %.4f\n', corr_n);
else
    fprintf('  No vertices above amplitude threshold at SNR=10dB.\n');
end

fprintf('\n=======================================================\n');
fprintf(' ALL TESTS COMPLETE\n');
fprintf('=======================================================\n');

end
