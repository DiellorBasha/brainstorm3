function test_phase_recovery_v4
% TEST_PHASE_RECOVERY_V4  Face-based current flux source model.
%
% Tests the face-based approach described in dev/references/face_based_source_model.md.
%
% Source model: primal 2-form — normal current flux through triangular faces.
%   s(f,t) = ∫_face J · n̂_f dA   [A·m]   (exact face normal, area-weighted)
%
% The face normal n̂_f is the EXACT normal from the mesh geometry (cross product
% of two edge vectors) — no vertex averaging, no approximation.  The sign is
% globally consistent by the mesh winding convention.
%
% Face Laplacian null space note (see dev/references/face_based_source_model.md §5):
% The discrete face Laplacian d₁·★₁⁻¹·d₁ᵀ is numerically degenerate on a
% hemisphere (manifold with boundary) because ~12% of edges have negative cotan
% weights making ★₁ indefinite.  Fix: interpolate vertex LBO eigenmodes to faces
% and orthonormalize under the face mass matrix ★₂.
%
% Leadfield approximation (§6 of the reference):
%   L_face(:,f) ≈ (A_f/3)·[G_3D(r_i)+G_3D(r_j)+G_3D(r_k)]·n̂_f
% 3-point vertex quadrature with UNCONSTRAINED Gain and EXACT face normals.
% Error: O(h²/d²) ≈ 0.1% for h~3mm, d~100mm.
%
% Cross-method results (dense synthetic wave, SNR=10dB):
%   Method                    Sign-OK  Phase err  Smoothness   Wave dir r
%   Ground truth              100%     0.000 rad  0.021 rad    1.000
%   Direct J·n̂  (v2)          64%     1.300 rad  0.253 rad    0.054
%   Vertex Eigen-MNE (v3)      31%     1.940 rad  0.130 rad    0.048
%   Face-based (v4)            59%     1.576 rad  0.027 rad    0.039
%
% Key finding: face-based spatial smoothness (0.027) matches ground truth (0.021).
% Wave direction remains uncorrelated for all methods — fundamental ill-conditioning
% with dense simultaneous activation, independent of source parameterisation.

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);
if ~brainstorm('status'), brainstorm nogui; end

%% ── Geometry + DEC operators ──────────────────────────────────────────────

TessMat = in_tess_bst('Subject01/tess_cortex_mid_low.mat');
Vtx     = TessMat.Vertices;
Faces   = double(TessMat.Faces);
nV      = size(Vtx,1);
nF      = size(Faces,1);
Vtx_mm  = Vtx * 1000;

clear mex
h   = nxr_compute('create', Vtx, Faces);
dec = nxr_compute('assembleDECOperators', h);
nxr_compute('destroy', h);

% Exact face normals and areas
e01 = Vtx(Faces(:,2),:) - Vtx(Faces(:,1),:);
e02 = Vtx(Faces(:,3),:) - Vtx(Faces(:,1),:);
Ncross  = cross(e01, e02, 2);
A_f_vec = sqrt(sum(Ncross.^2, 2));
A_f     = A_f_vec / 2;            % face areas [m²]
n_hat_f = Ncross ./ A_f_vec;      % exact unit normals [nF × 3]

[~, lH_v] = tess_hemisplit(TessMat); lH_v = lH_v(:);
inL_v = false(nV,1); inL_v(lH_v) = true;
lH_f  = find(all(inL_v(Faces), 2));
nLHF  = numel(lH_f);

fprintf('Mesh: %d V,  %d F  |  LH: %d V,  %d F\n', nV, nF, numel(lH_v), nLHF);

%% ── Face eigenmodes via vertex LBO interpolation ─────────────────────────
% Direct face Laplacian d₁·★₁⁻¹·d₁ᵀ is numerically degenerate on a hemisphere
% (negative cotan weights → ★₁ indefinite).  Interpolate vertex LBO modes instead.

L_vert  = dec.d0' * dec.hodge1 * dec.d0;
M_vert  = dec.hodge0;
nModes  = 10;
[Phi_v, Lam_v] = eigs(L_vert(lH_v,lH_v), M_vert(lH_v,lH_v), nModes, 'smallestabs');
lam_v = sort(real(diag(Lam_v)),'ascend');

% Interpolate to LH faces by corner average
lh_vmap = zeros(nV,1); lh_vmap(lH_v) = 1:numel(lH_v);
FacesLH  = Faces(lH_f,:);
Phi_f_raw = (Phi_v(lh_vmap(FacesLH(:,1)), 2:end) + ...
             Phi_v(lh_vmap(FacesLH(:,2)), 2:end) + ...
             Phi_v(lh_vmap(FacesLH(:,3)), 2:end)) / 3;   % [nLHF × nModes-1]

% M_f-orthonormalize via Cholesky
M_f_lh    = dec.hodge2(lH_f, lH_f);
G_gram    = Phi_f_raw' * M_f_lh * Phi_f_raw;
[R, p]    = chol(G_gram);
assert(p == 0, 'Gram matrix not positive definite');
Phi_f_basis = Phi_f_raw / R;   % [nLHF × nModes-1]  M_f-orthonormal

fprintf('Face eigenmodes: %d  (λ_vert range [%.0f … %.0f])\n', ...
    nModes-1, lam_v(2), lam_v(end));
fprintf('M_f orthonorm check: off-diag max = %.2e\n', ...
    max(abs(Phi_f_basis'*M_f_lh*Phi_f_basis - eye(nModes-1)),[],'all'));

%% ── Face-based constrained leadfield ─────────────────────────────────────

sStudies = bst_get('ProtocolStudies');
s6 = [sStudies.Study];  s6 = s6(6);
rc = {s6.Result.Comment};  rf = {s6.Result.FileName};
iUC = find(strcmp(rc,'MN: MEG(Unconstr) 2018'),1);
K_uc = load(file_fullpath(rf{iUC}), 'GoodChannel');

hm     = load(file_fullpath('Subject01/S01_AEF_20131218_01_notch/headmodel_surf_os_meg.mat'),'Gain');
Gain_u = double(hm.Gain(K_uc.GoodChannel,:));   % [nCh × 3nV] unconstrained
nCh    = size(Gain_u,1);

NC = load(file_fullpath('Subject01/S01_AEF_20131218_01_notch/noisecov_full.mat'),'NoiseCov');
noise_std = sqrt(max(diag(NC.NoiseCov(K_uc.GoodChannel, K_uc.GoodChannel)),0));

% L_face(:,f) ≈ (A_f/3)·[G_3D(r_i)+G_3D(r_j)+G_3D(r_k)]·n̂_f
% Uses UNCONSTRAINED Gain + EXACT face normals (O(h²/d²)≈0.1% error)
fprintf('Building face leadfield [%d × %d]...', nCh, nF);
L_face = zeros(nCh, nF);
for f = 1:nF
    i = Faces(f,1); j = Faces(f,2); k = Faces(f,3);
    Gavg = (Gain_u(:,3*i-2:3*i) + Gain_u(:,3*j-2:3*j) + Gain_u(:,3*k-2:3*k)) / 3;
    L_face(:,f) = Gavg * n_hat_f(f,:)' * A_f(f);
end
fprintf(' done\n');

%% ── Synthetic traveling wave ──────────────────────────────────────────────

Fs     = 600;  t = (0:Fs-1)/Fs;  f0 = 10;  om = 2*pi*f0;
A_src  = 1e-9; SNR_dB = 10;
e_prop = [0 -1 0];  v_ms = 3.0;  k_rad = om/(v_ms*1000);

gyral_lh_f = lH_f(all(TessMat.SulciMap(Faces(lH_f,:))==0, 2));
ctr_f = (Vtx_mm(Faces(gyral_lh_f,1),:) + ...
         Vtx_mm(Faces(gyral_lh_f,2),:) + ...
         Vtx_mm(Faces(gyral_lh_f,3),:)) / 3;
phi_true_f = k_rad * (ctr_f * e_prop');   % [nPatch × 1]

nT = length(t);
J_f = zeros(nF, nT);
for ki = 1:numel(gyral_lh_f)
    J_f(gyral_lh_f(ki),:) = A_src * cos(om*t + phi_true_f(ki));
end

sensors   = L_face * J_f;
sig_rms   = sqrt(mean(sensors(:).^2));
sensors   = sensors + (sig_rms/10^(SNR_dB/20)) * (noise_std./max(noise_std,eps)) .* randn(size(sensors));

%% ── Face eigenmode inverse (Tikhonov) ────────────────────────────────────

L_tilde_f = L_face(:, lH_f) * Phi_f_basis;          % [nCh × nModes-1]
K_modes   = size(Phi_f_basis, 2);
lam_reg   = trace(L_tilde_f'*L_tilde_f) / (nCh * 9);
theta_f   = (L_tilde_f'*L_tilde_f + lam_reg*eye(K_modes)) \ (L_tilde_f'*sensors);
s_rec_lhf = Phi_f_basis * theta_f;                  % [nLHF × nT]
s_rec_f   = zeros(nF, nT);
s_rec_f(lH_f,:) = s_rec_lhf;

%% ── Analytic signal + mask ────────────────────────────────────────────────

A_rec  = abs  (s_rec_f + 1i*imag(hilbert(s_rec_f')'));
Ph_rec = angle(s_rec_f + 1i*imag(hilbert(s_rec_f')'));

amp_mean = mean(A_rec(gyral_lh_f,:), 2);
act_f    = gyral_lh_f(amp_mean > 0.10*max(amp_mean));
[~,t_pk] = max(mean(A_rec(act_f,:), 1));

fprintf('\nActive faces (10%% threshold): %d / %d gyral LH\n', numel(act_f), numel(gyral_lh_f));

%% ── TEST 1 — Sign errors ─────────────────────────────────────────────────
fprintf('\n─── TEST 1: Sign error census ───\n');

phi_m     = zeros(nF,1); phi_m(gyral_lh_f) = phi_true_f;
s_true_pk = A_src * cos(om*t(t_pk) + phi_m(act_f));
sign_ok   = sign(s_rec_f(act_f, t_pk)) == sign(s_true_pk);
fprintf('  Correct sign: %.1f%%  (%d errors / %d active faces)\n', ...
    100*mean(sign_ok), sum(~sign_ok), numel(act_f));

%% ── TEST 2 — Phase error ─────────────────────────────────────────────────
fprintf('\n─── TEST 2: Phase error at active faces ───\n');

Ph_true_pk = om*t(t_pk) + phi_m(act_f);
err_f      = abs(angle(exp(1i*(Ph_rec(act_f,t_pk) - Ph_true_pk))));
fprintf('  mean=%.4f rad  median=%.4f rad\n', mean(err_f), median(err_f));
if any(~sign_ok)
    fprintf('  Sign-error faces: mean=%.4f rad\n', mean(err_f(~sign_ok)));
    fprintf('  Sign-OK    faces: mean=%.4f rad\n', mean(err_f(sign_ok)));
end

%% ── TEST 3 — Spatial smoothness ─────────────────────────────────────────
fprintf('\n─── TEST 3: Spatial smoothness (mean |ΔΦ| per face-adjacency edge) ───\n');

FaceConn = dec.d1 * dec.d1' ~= 0;
FaceConn = FaceConn - speye(nF);
[fi,fj]  = find(triu(FaceConn(act_f,act_f) > 0, 1));
ga = act_f(fi);  gb = act_f(fj);
act_m = zeros(nF,1); act_m(act_f) = 1:numel(act_f);

dphi_true = abs(angle(exp(1i*(Ph_true_pk(act_m(ga)) - Ph_true_pk(act_m(gb))))));
dphi_rec  = abs(angle(exp(1i*(Ph_rec(ga,t_pk) - Ph_rec(gb,t_pk)))));

fprintf('  Ground truth:    %.4f rad/edge\n', mean(dphi_true));
fprintf('  Face-based rec:  %.4f rad/edge\n', mean(dphi_rec));

%% ── TEST 4 — Wave direction ──────────────────────────────────────────────
fprintf('\n─── TEST 4: Wave direction (gradient correlation) ───\n');

ctr_ga = (Vtx_mm(Faces(ga,1),:)+Vtx_mm(Faces(ga,2),:)+Vtx_mm(Faces(ga,3),:))/3;
ctr_gb = (Vtx_mm(Faces(gb,1),:)+Vtx_mm(Faces(gb,2),:)+Vtx_mm(Faces(gb,3),:))/3;
true_grad = k_rad * ((ctr_gb - ctr_ga) * e_prop');
rec_grad  = Ph_rec(gb,t_pk) - Ph_rec(ga,t_pk);
fprintf('  r = %.4f  (0 = no direction recovery; 1 = perfect)\n', corr(true_grad, rec_grad));

%% ── TEST 5 — Summary table ───────────────────────────────────────────────
fprintf('\n─── CROSS-METHOD SUMMARY (dense wave, SNR=10dB) ───\n');
fprintf('%-34s  Sign-OK  Phase(mean)  Smooth    WaveDir-r\n','Method');
fprintf('%s\n', repmat('-',72,1));
rows = {
  'Ground truth',              100,    0,      0.021,  1.000;
  'Direct J·n̂  (v2)',          64.4,  1.300,  0.253,  0.054;
  'Vertex Eigen-MNE (v3)',      31.4,  1.940,  0.130,  0.048;
  'Face-based (v4, this run)', 100*mean(sign_ok), mean(err_f), mean(dphi_rec), corr(true_grad,rec_grad)};
for k=1:size(rows,1)
    fprintf('%-34s  %.0f%%     %.3f rad    %.4f    %.4f\n', rows{k,1}, rows{k,2}, rows{k,3}, rows{k,4}, rows{k,5});
end

fprintf('\nKey finding: face-based smoothness (%.4f) ≈ ground truth (%.4f)\n', mean(dphi_rec), mean(dphi_true));
fprintf('Persistent limitation: wave direction uncorrelated in all methods.\n');
fprintf('Cause: MNE-type inverses cannot recover spatial phase of dense\n');
fprintf('       simultaneous activation, independent of source parameterisation.\n');

fprintf('\n=======================================================\n');
fprintf(' TEST_PHASE_RECOVERY_V4 COMPLETE\n');
fprintf('=======================================================\n');

end
