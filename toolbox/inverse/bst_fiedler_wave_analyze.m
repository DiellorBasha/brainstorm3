function Results = bst_fiedler_wave_analyze(Z_joint, phi_F, FaceIndices, SurfaceFile, WTInfo, varargin)
% BST_FIEDLER_WAVE_ANALYZE  Wave analysis in the joint Fiedler coordinate frame.
%
% USAGE:
%   R = bst_fiedler_wave_analyze(Z_joint, phi_F, FaceIndices, SurfaceFile, WTInfo)
%
% DESCRIPTION:
%   Extracts wave properties from the joint Fiedler field Z_joint(x,t) where
%   arg(Z_joint) ≈ ω₀t − k·φ_F(x) in Fiedler coordinates.
%
%   Source localization: the source face has MINIMUM temporal variance of
%   arg(Z_joint(x,t)) — the phase is stationary at the origin of the wave.
%
%   Wave speed from slope: at fixed burst-onset time t₀, arg(Z_joint) vs
%   φ_F(x) has slope −k = −ω₀/v. Linear regression gives v directly.
%
%   PLV: amplitude-weighted circular mean of exp(i·arg(Z_joint)) over faces
%   at each time step — directional coherence in Fiedler coordinates.
%
% INPUTS:
%   Z_joint     [nLHF × nTime]  complex — from bst_fiedler_joint_field
%   phi_F       [nLHF × 1]     real    — Fiedler longitude arg(u₁(x))
%   FaceIndices [nLHF × 1]     global face indices
%   SurfaceFile Brainstorm surface path
%   WTInfo      struct from bst_sensor_cwt (.Fs, .f_scales)
%               OR struct with just .Fs if from Hilbert pipeline
%
% OPTIONS (name-value):
%   'CenterFreq'      Hz for speed = ω₀/k (default: 10)
%   'AmpThreshold'    fraction of peak amplitude for active face mask (default: 0.08)
%   'AmpThresholdTime' temporal gate fraction (default: 0.20)
%   'Verbose'         logical (default: true)
%
% OUTPUT  Results struct:
%   .SourceFace      global face index — minimum temporal phase variance
%   .SourcePos       [1 × 3] face centroid [m]
%   .v_fiedler       phase velocity from Fiedler slope [m/s]
%   .v_fiedler_std   uncertainty on slope estimate [m/s]
%   .k_fiedler       wavenumber in Fiedler coordinates [1/rad]
%   .PLV             [1 × nTime] directional coherence in Fiedler frame
%   .DomPhase        [1 × nTime] dominant phase in Fiedler frame
%   .MeanAmp         [1 × nTime] mean |Z_joint| across active faces
%   .ValidTime       [1 × nTime] logical temporal gate
%   .ActiveMask      [nLHF × 1] logical spatial gate
%   .phi_F_range     [min max] Fiedler longitude range [rad]
%
% Authors: Diellor Basha, 2026

%% ── Parse options ────────────────────────────────────────────────────────
CenterFreq      = 10;
AmpThresh       = 0.08;
AmpThreshTime   = 0.20;
Verbose         = true;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'centerfreq',       CenterFreq    = varargin{k+1};
        case 'ampthreshold',     AmpThresh     = varargin{k+1};
        case 'ampthresholdtime', AmpThreshTime = varargin{k+1};
        case 'verbose',          Verbose       = logical(varargin{k+1});
    end
end

nLHF  = size(Z_joint, 1);
nTime = size(Z_joint, 2);
omega0 = 2*pi * CenterFreq;

%% ── Load face centroids ───────────────────────────────────────────────────
TessMat = in_tess_bst(SurfaceFile, 0);
Vtx   = TessMat.Vertices;
Faces = double(TessMat.Faces);
ctr_f = (Vtx(Faces(:,1),:)+Vtx(Faces(:,2),:)+Vtx(Faces(:,3),:))/3;
ctr_lh = ctr_f(FaceIndices, :);   % [nLHF × 3] m

%% ── Amplitude masks ──────────────────────────────────────────────────────
amp_env  = abs(Z_joint);                       % [nLHF × nTime]
ampMean  = mean(amp_env, 2);                   % [nLHF × 1]
activeF  = ampMean >= AmpThresh * max(ampMean);
amp_t    = mean(amp_env, 1);                   % [1 × nTime]
validTime = amp_t >= AmpThreshTime * max(amp_t);

%% ── Phase field ──────────────────────────────────────────────────────────
phi_joint = angle(Z_joint);   % [nLHF × nTime]  arg(Z_joint) = wave phase in Fiedler coords

%% ── Source: minimum temporal phase variance ──────────────────────────────
% Source face = where the phase is most temporally stable during high-amp periods
phi_var = var(phi_joint(:, validTime), [], 2);   % [nLHF × 1]
phi_var(~activeF) = Inf;
[~, src_local] = min(phi_var);
src_global = FaceIndices(src_local);
src_pos    = ctr_f(src_global, :);   % [1 × 3] m

%% ── Wave speed from Fiedler slope at burst peak ──────────────────────────
% At burst peak: arg(Z_joint(x,t₀)) ≈ const − k·φ_F(x)
% Linear regression: phi_joint ~ a + b·phi_F  →  k = −b  →  v = ω₀/k
[~, t_peak] = max(amp_t);
act_idx = find(activeF & phi_var < prctile(phi_var(activeF), 75));
phi_joint_peak = phi_joint(act_idx, t_peak);
phi_F_act      = phi_F(act_idx);

% Circular regression: minimize circular variance of (phi_joint - b*phi_F)
b_vals = linspace(-60, 60, 2401);   % range of candidate k values
cv = zeros(size(b_vals));
for bi = 1:numel(b_vals)
    residual = phi_joint_peak - b_vals(bi) * phi_F_act;
    cv(bi) = 1 - abs(mean(exp(1i * residual)));   % circular variance
end
[~, best_bi] = min(cv);
k_fiedler = -b_vals(best_bi);   % slope dφ_joint/dφ_F [dimensionless]

% Physical scaling: φ_F is dimensionless (Fiedler longitude in radians).
% The wave phase φ_joint = ω₀t − k_phys·d(x) where d is physical distance.
% slope b = dφ_joint/dφ_F, and k_phys = |b| · |∇φ_F| [rad/m].
% Estimate |∇φ_F| from active-region centroids: regress φ_F on position,
% the gradient magnitude is the physical longitude-to-distance scale.
ctr_act = ctr_lh(act_idx, :);                  % [nAct × 3] m
pF      = phi_F_act;
% Robust gradient: |∇φ_F| ≈ range(φ_F) / spatial extent along principal axis
% Use PCA of centroids to find the dominant spatial axis, project φ_F onto it
ctr_c = ctr_act - mean(ctr_act, 1);
[~, ~, Vpca] = svd(ctr_c, 'econ');
proj  = ctr_c * Vpca(:,1);                     % [nAct × 1] position along principal axis [m]
% Linear fit φ_F vs proj → slope = |∇φ_F| along principal axis [rad/m]
grad_phiF = abs((proj(:)' * (pF - mean(pF))) / max(proj(:)'*proj(:), eps));
grad_phiF = max(grad_phiF, eps);

k_phys = abs(k_fiedler) * grad_phiF;           % [rad/m]
v_fiedler = omega0 / max(k_phys, eps);          % [m/s]

% Bootstrap uncertainty over the slope estimate
nBoot = 20; v_boot = zeros(nBoot,1);
for bi_boot = 1:nBoot
    idx_b = randsample(numel(act_idx), round(0.7*numel(act_idx)), false);
    ph_b  = phi_joint_peak(idx_b);
    pF_b  = phi_F_act(idx_b);
    cv_b  = arrayfun(@(b) 1-abs(mean(exp(1i*(ph_b - b*pF_b)))), b_vals);
    [~,bb] = min(cv_b);
    v_boot(bi_boot) = omega0 / max(abs(b_vals(bb))*grad_phiF, eps);
end
v_fiedler_std = std(v_boot);

%% ── PLV and dominant phase over time ─────────────────────────────────────
PLV      = nan(1, nTime);
DomPhase = nan(1, nTime);
for t = 1:nTime
    if ~validTime(t) || ~any(activeF), continue; end
    A_t = amp_env(activeF, t);
    ph_t = phi_joint(activeF, t);
    wSum = sum(A_t);
    if wSum < eps, continue; end
    r          = sum(A_t .* exp(1i*ph_t)) / wSum;
    PLV(t)     = abs(r);
    DomPhase(t) = angle(r);
end

%% ── Pack output ──────────────────────────────────────────────────────────
Results.SourceFace     = src_global;
Results.SourcePos      = src_pos;
Results.v_fiedler      = v_fiedler;
Results.v_fiedler_std  = v_fiedler_std;
Results.k_fiedler      = k_fiedler;
Results.PLV            = PLV;
Results.DomPhase       = DomPhase;
Results.MeanAmp        = amp_t;
Results.ValidTime      = validTime;
Results.ActiveMask     = activeF;
Results.phi_F_range    = [min(phi_F), max(phi_F)];

if Verbose
    fprintf('bst_fiedler_wave_analyze:\n');
    fprintf('  Source: face %d  pos=[%.1f %.1f %.1f] mm\n', ...
        src_global, src_pos*1000);
    fprintf('  Wave speed: %.2f ± %.2f m/s  (k=%.2f rad⁻¹)\n', ...
        v_fiedler, v_fiedler_std, k_fiedler);
    fprintf('  PLV: max=%.3f  mean=%.3f\n', ...
        max(PLV,[],'omitnan'), mean(PLV,'omitnan'));
end
end
