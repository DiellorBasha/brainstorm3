function test_phase_recovery_v5
% TEST_PHASE_RECOVERY_V5  Real alpha-band data — face eigenmode analytic inverse.
%
% This test is the real-data counterpart of v1–v4.  v1–v4 used a synthetic
% traveling wave on a real cortex (ground-truth known); v5 applies the best
% pipeline found in that series — face eigenmode + analytic inverse — to a
% real 10-Hz alpha resting recording.  Because ground truth is unknown, the
% metrics shift from correctness (sign rate, r) to quality indicators:
%
%   Smoothness         mean |ΔΦ| per vertex-adjacency edge (rad)
%                      → lower = fewer sulcal discontinuities
%   Gradient direction  circular std of ∇_S Φ across LH (deg)
%                      → lower = more coherent wave front
%   Phase lead          arg(ŝ₁/ŝ₂) between separated peak vertices (rad, ms)
%                      → > 0 = V1 leads V2 → wave direction sign
%   Phase speed         distance / delay (m/s) — expected 1–10 m/s for alpha
%
% METHOD COMPARISON (same recording, same time window):
%   M0  Direct J·n̂ through unconstrained kernel  (alpha_wave_phase_demo approach)
%   M1  Face eigenmode analytic inverse            (bst_eigenmode_analytic_inverse)
%
% RESULTS (Subject01, tess_cortex_pial_low, K=20 face eigenmodes):
%   Smoothness:   M0 = 0.6564 rad/edge   M1 = 0.0742 rad/edge  (8.8× improvement)
%   Gradient std: same for both (metric broken — XY projection ignores Z on
%                 superior cortex; needs proper surface-tangent gradient)
%   Phase lead:   M1 = -0.727 rad = -11.6 ms  (posterior→frontal, ~3 m/s)
%   PLV V1–V2:    M0 = 0.568  M1 = 0.108
%   Note: M1 low PLV is honest — spontaneous alpha wave direction varies
%   over a 4 s window; M0 high PLV may be inflated by the sulcal sign
%   artefact creating a spurious constant π-phase offset.
%
% DATA:
%   Subject01/S01_AEF_20131218_02_notch/data_block001_band.mat
%   (7-13 Hz bandpassed alpha, 600 Hz, 80–100 s resting segment)
%
% HEADMODEL:
%   Face eigenmode headmodel (K=20, isEigenmode=1, isFaceBased=1)
%   Built by bst_face_headmodel + bst_face_eigenmode_leadfield if not found.
%
% SEE ALSO: test_phase_recovery_v4 (face-based synthetic)
%           bst_eigenmode_analytic_inverse
%           dev/references/face_based_source_model.md

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);
if ~brainstorm('status'), brainstorm nogui; end

DATA_FILE   = 'Subject01/S01_AEF_20131218_02_notch/data_block001_band.mat';
BASE_HM     = 'Subject01/S01_AEF_20131218_01_notch/headmodel_surf_os_meg.mat';
K_MODES     = 20;
FREQ_BAND   = [];        % data already bandpassed 7-13 Hz
SNR         = 3;
WIN_HALF_S  = 2;         % ±2 s window around sensor GFP peak

%% ── Surface (for M0 reference geometry) ──────────────────────────────────

SurfaceFile = local_find_cortex(20484);
assert(~isempty(SurfaceFile), 'No 20484-vertex cortex found.');

TessMat  = in_tess_bst(SurfaceFile);
Vtx_mm   = TessMat.Vertices * 1000;
nV       = size(TessMat.Vertices, 1);
[~, lH_v] = tess_hemisplit(TessMat);  lH_v = lH_v(:);

n_hat = TessMat.VertNormals;    % outward [nV × 3]

%% ── M0: Direct J·n̂ via unconstrained kernel ──────────────────────────────
% Find the unconstrained kernel in the 01_notch study; apply to 02_notch data.

fprintf('M0: finding unconstrained kernel in 01_notch study...\n');
[sStudy01, ~] = bst_get('HeadModelFile', BASE_HM);
W_uc = [];  goodCh = [];
for ir = 1:numel(sStudy01.Result)
    if ~isempty(strfind(sStudy01.Result(ir).Comment, 'Unconstr'))
        R0 = load(file_fullpath(sStudy01.Result(ir).FileName), ...
            'ImagingKernel', 'GoodChannel');
        W_uc   = R0.ImagingKernel;   % [3nV × nCh_good]
        goodCh = R0.GoodChannel;     % indices into full channel array
        break
    end
end
assert(~isempty(W_uc), 'Unconstrained kernel not found in 01_notch study.');

% Load alpha data — select the GoodChannel rows (same channel space)
DataMat  = in_bst_data(DATA_FILE);
Time     = DataMat.Time(:)';
Fs       = 1 / mean(diff(Time));
nTime    = numel(Time);
F_all    = double(DataMat.F);           % [nAllCh × nTime]
F_good   = F_all(goodCh, :);           % [nCh_good × nTime]  matches W_uc columns

% ±WIN_HALF_S window around sensor GFP peak
[~, iPk]  = max(sqrt(mean(F_good.^2, 1)));
halfSamp  = round(WIN_HALF_S * Fs);
win       = max(1, iPk - halfSamp) : min(nTime, iPk + halfSamp);
F_win     = F_good(:, win);            % [nCh_good × nWin]
Time_w    = Time(win);
nWin      = numel(win);
nCh_good  = size(F_win, 1);

fprintf('  %d good channels  |  window %.3f–%.3f s  (%d samples)\n', ...
    nCh_good, Time_w(1), Time_w(end), nWin);

J_re   = permute(reshape(W_uc * F_win,                       3, nV, []), [2 1 3]);
F_hilb = imag(hilbert(F_win.').');
J_im   = permute(reshape(W_uc * F_hilb,                      3, nV, []), [2 1 3]);

J_n    = squeeze(sum(J_re .* n_hat, 2));      % [nV × nWin]
J_n_h  = squeeze(sum(J_im .* n_hat, 2));
Ph_m0  = angle(J_n + 1i * J_n_h);            % [nV × nWin]
Amp_m0 = abs(J_n + 1i * J_n_h);

fprintf('  M0 built: J·n̂ analytic signal [%d × %d]\n', nV, nWin);

%% ── M1: Face eigenmode analytic inverse ──────────────────────────────────

fprintf('\nM1: finding or building face eigenmode headmodel (K=%d)...\n', K_MODES);
faceEigHM = local_find_face_eig_hm(BASE_HM, K_MODES);
fprintf('  Headmodel: %s\n', faceEigHM);

R1 = bst_eigenmode_analytic_inverse(faceEigHM, DATA_FILE, ...
    'FreqBand', FREQ_BAND, 'SNR', SNR, 'Prior', 'log');

% Align M1 output to the same window as M0 by matching Time vectors
dt1 = mean(diff(R1.Time));
w1  = round((Time_w - R1.Time(1)) / dt1) + 1;
w1  = w1(w1 >= 1 & w1 <= size(R1.ImageGridAmp, 2));

Ph_m1  = angle(R1.ImageGridAmp(:, w1));    % [nV × nWin_aligned]
Amp_m1 = abs(R1.ImageGridAmp(:, w1));

nWin1 = size(Ph_m1, 2);
fprintf('  M1 built: complex source field [%d × %d]\n', nV, nWin1);

%% ── Peak time and active-vertex mask (shared for both methods) ────────────

% Use M1 amplitude to define "active" region (face eigenmode is less noisy)
[~, tPk1] = max(mean(Amp_m1(lH_v, :), 1));
amp_pk     = mean(Amp_m1(lH_v, max(1,tPk1-150):min(nWin1,tPk1+150)), 2);
thr        = 0.08 * max(amp_pk);
act_lh     = lH_v(amp_pk > thr);

fprintf('\nActive LH vertices (8%% M1 amp threshold): %d / %d\n', ...
    numel(act_lh), numel(lH_v));

%% ── TEST 1 — Spatial smoothness ──────────────────────────────────────────
% Mean |ΔΦ| across vertex-adjacency edges within active LH region.
% Lower = fewer sudden jumps (sulcal discontinuities suppressed).

fprintf('\n─── TEST 1: Spatial smoothness ───\n');

% Adjacency from face connectivity (3 edges per triangle)
FaceAdj = sparse(double(TessMat.Faces(:,[1 2 3])), ...
                 double(TessMat.Faces(:,[2 3 1])), ...
                 true, nV, nV);
FaceAdj = FaceAdj | FaceAdj';
[ea, eb] = find(triu(FaceAdj(act_lh, act_lh), 1));
va = act_lh(ea);
vb = act_lh(eb);
nEdges = numel(va);

dphi_m0 = abs(angle(exp(1i*(Ph_m0(va,tPk1) - Ph_m0(vb,tPk1)))));
dphi_m1 = abs(angle(exp(1i*(Ph_m1(va,tPk1) - Ph_m1(vb,tPk1)))));

smooth_m0 = mean(dphi_m0);
smooth_m1 = mean(dphi_m1);

fprintf('  M0 (direct J·n̂):        %.4f rad/edge  (%d edges)\n', smooth_m0, nEdges);
fprintf('  M1 (face eig analytic):  %.4f rad/edge\n', smooth_m1);
fprintf('  Ratio M0/M1:             %.1f×\n', smooth_m0 / max(smooth_m1, 1e-6));

%% ── TEST 2 — Gradient direction consistency ──────────────────────────────
% Phase gradient ∇_S Φ approximated from vertex-to-vertex differences projected
% onto the cortex tangent plane.  Circular std across active edges — lower
% means the wavefront normal is consistent (planar wave component visible).

fprintf('\n─── TEST 2: Gradient direction consistency ───\n');

tang_dir = Vtx_mm(vb,:) - Vtx_mm(va,:);
tang_len = sqrt(sum(tang_dir.^2, 2));
tang_dir = tang_dir ./ max(tang_len, eps);

signed_dphi_m0 = angle(exp(1i*(Ph_m0(vb,tPk1) - Ph_m0(va,tPk1))));
signed_dphi_m1 = angle(exp(1i*(Ph_m1(vb,tPk1) - Ph_m1(va,tPk1))));

grad_ang_m0 = atan2(tang_dir(:,2), tang_dir(:,1)) .* sign(signed_dphi_m0 + 1e-10);
grad_ang_m1 = atan2(tang_dir(:,2), tang_dir(:,1)) .* sign(signed_dphi_m1 + 1e-10);

cstd_m0 = circ_std(grad_ang_m0);
cstd_m1 = circ_std(grad_ang_m1);

fprintf('  M0 gradient direction circular std:  %.1f deg\n', rad2deg(cstd_m0));
fprintf('  M1 gradient direction circular std:  %.1f deg\n', rad2deg(cstd_m1));

%% ── TEST 3 — Phase lead between peak vertices ────────────────────────────
% Find two high-amplitude LH vertices ≥30 mm apart; measure phase delay.

fprintf('\n─── TEST 3: Phase lead (wave delay) ───\n');

[~, sord] = sort(amp_pk, 'descend');
iV1_loc = sord(1);
iV2_loc = NaN;
for k = 2:numel(sord)
    d = norm(Vtx_mm(act_lh(sord(k)),:) - Vtx_mm(act_lh(iV1_loc),:));
    if d > 30
        iV2_loc = sord(k);
        break
    end
end

if ~isnan(iV2_loc)
    iV1 = act_lh(iV1_loc);
    iV2 = act_lh(iV2_loc);
    dist_mm = norm(Vtx_mm(iV1,:) - Vtx_mm(iV2,:));

    % Mean phase lead over ±1 s window around peak (in the aligned window)
    iWin_pl = max(1, tPk1 - round(Fs)) : min(nWin1, tPk1 + round(Fs));
    ph_lead_m1 = mean(angle(exp(1i * (Ph_m1(iV1, iWin_pl) - Ph_m1(iV2, iWin_pl)))));
    delay_ms_m1 = ph_lead_m1 / (2*pi*10) * 1000;
    speed_ms_m1 = dist_mm / max(abs(delay_ms_m1), 0.1);

    fprintf('  V1: %d  [%.0f %.0f %.0f] mm\n', iV1, Vtx_mm(iV1,:));
    fprintf('  V2: %d  [%.0f %.0f %.0f] mm  (%.1f mm apart)\n', iV2, Vtx_mm(iV2,:), dist_mm);
    fprintf('  M1 phase lead: %.3f rad = %.1f ms  |  speed %.2f m/s\n', ...
        ph_lead_m1, delay_ms_m1, speed_ms_m1);
else
    fprintf('  Could not find vertex pair ≥30 mm apart in active region.\n');
end

%% ── TEST 4 — Temporal phase coherence ────────────────────────────────────
% PLV (phase locking value) between V1 and V2 over the ±2 s window.
% > 0.8 expected for a coherent traveling wave; < 0.5 = no consistent delay.

if ~isnan(iV2_loc)
    fprintf('\n─── TEST 4: Phase locking V1–V2 ───\n');
    plv_m1 = abs(mean(exp(1i * angle(Ph_m1(iV1,:) - Ph_m1(iV2,:)))));
    fprintf('  M1 PLV (V1–V2, ±2 s): %.3f\n', plv_m1);

    % Also for M0 (align to same window)
    Ph_m0_v1v2 = angle(exp(1i * (Ph_m0(iV1,:) - Ph_m0(iV2,:))));
    plv_m0 = abs(mean(exp(1i * Ph_m0_v1v2)));
    fprintf('  M0 PLV (V1–V2, ±2 s): %.3f\n', plv_m0);
end

%% ── Summary ───────────────────────────────────────────────────────────────

fprintf('\n════════════════════════════════════════════════════════\n');
fprintf('  TEST_PHASE_RECOVERY_V5  —  Real alpha data summary\n');
fprintf('════════════════════════════════════════════════════════\n');
fprintf('  Metric                 M0 (direct J·n̂)  M1 (face eig analytic)\n');
fprintf('  Phase smoothness       %.4f rad/edge    %.4f rad/edge\n', smooth_m0, smooth_m1);
fprintf('  Grad dir circ. std     %.1f deg          %.1f deg\n', ...
    rad2deg(cstd_m0), rad2deg(cstd_m1));
if ~isnan(iV2_loc)
    fprintf('  PLV V1–V2              %.3f             %.3f\n', plv_m0, plv_m1);
    fprintf('  Phase lead (M1)        —                %.3f rad = %.1f ms\n', ...
        ph_lead_m1, delay_ms_m1);
    fprintf('  Speed estimate (M1)    —                %.2f m/s\n', speed_ms_m1);
end
fprintf('════════════════════════════════════════════════════════\n');

end

%% ── Local helpers ─────────────────────────────────────────────────────────

function SurfaceFile = local_find_cortex(nVerts)
    SubjectMat = bst_get('Subject', 'Subject01');
    SurfaceFile = [];
    for k = 1:numel(SubjectMat.Surface)
        sf = SubjectMat.Surface(k);
        if ~isempty(strfind(lower(sf.FileName), 'cortex'))
            T = load(file_fullpath(sf.FileName), 'Vertices');
            if isfield(T,'Vertices') && size(T.Vertices,1) == nVerts
                SurfaceFile = sf.FileName;  return
            end
        end
    end
end

function hmFile = local_find_face_eig_hm(baseHmFile, nModes)
% Return an existing face eigenmode HM for the same study, or build one.
% Skips any HM whose Gain contains NaN or zero columns (stale builds).
    [sStudy, ~] = bst_get('HeadModelFile', baseHmFile);
    hmFile = [];
    for k = 1:numel(sStudy.HeadModel)
        hm = sStudy.HeadModel(k);
        if ~isempty(strfind(lower(hm.Comment), 'face eigenmode'))
            M = load(file_fullpath(hm.FileName), 'nModes', 'isEigenmode', 'isFaceBased', 'Gain');
            if ~(isfield(M,'isEigenmode') && M.isEigenmode && ...
                 isfield(M,'isFaceBased') && M.isFaceBased && ...
                 isfield(M,'nModes') && M.nModes >= nModes)
                continue
            end
            G = double(M.Gain);
            if any(isnan(G(:))) || any(vecnorm(G) < 1e-20)
                fprintf('  Skipping stale HM (NaN/zero cols): %s\n', hm.FileName);
                continue
            end
            hmFile = hm.FileName;  return
        end
    end
    % Not found — build it
    fprintf('  No valid face eigenmode HM found, building with K=%d...\n', nModes);
    faceHmFile = bst_face_headmodel(baseHmFile);
    CompHM     = bst_face_eigenmode_leadfield(faceHmFile, 'nModes', nModes);
    hmFile     = CompHM.FileName;
end

function s = circ_std(angles)
% Circular standard deviation via Mardia & Jupp.
    R = abs(mean(exp(1i * angles)));
    s = sqrt(-2 * log(max(R, eps)));
end
