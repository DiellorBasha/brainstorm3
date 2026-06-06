function Results = bst_eigenmode_poisson_sharpen(Inv, varargin)
% BST_EIGENMODE_POISSON_SHARPEN  Laplacian / Poisson deconvolution of an
%   eigenmode inverse estimate to recover a sharper source map.
%
% USAGE:
%   R = bst_eigenmode_poisson_sharpen(InvResults)
%   R = bst_eigenmode_poisson_sharpen(InvResults, 'Alpha', 0.1)
%   R = bst_eigenmode_poisson_sharpen(InvResults, 'Mode', 'laplacian')
%
% MOTIVATION:
%   The eigenmode inverse recovers smoothed coefficients
%       θ_k(t) ≈ T_k · ρ_k^true,   T_k = σ_k² / (σ_k² + α λ_k^β)
%   which act as a spatial low-pass filter — the reconstructed field is a
%   spatially blurred version of the true source.
%
%   In eigenmode space the surface Laplacian acts as:
%       −Δ_S ψ_k = λ_k ψ_k
%   so the Laplacian of the estimate is:
%       −Δ_S ŝ(x,t) = Σ_k λ_k θ_k(t) ψ_k(x)
%   which re-weights coefficients by λ_k — boosting high spatial frequencies.
%
%   The full regularised Poisson deconvolution solves:
%       (−Δ_S + α) ρ = −Δ_S ŝ
%   whose eigenmode solution is:
%       ρ_k(t) = λ_k θ_k(t) / (λ_k + α)
%
%   For α → 0  this gives pure Laplacian sharpening  (noise-amplifying).
%   For α → ∞  this recovers ŝ / α  (no sharpening, just scaling).
%   The optimal α is near the noise floor of the high-k coefficients.
%
% INPUT:
%   Inv   Results struct from bst_eigenmode_analytic_inverse or
%         bst_eigenmode_oscillator_inverse (must contain .theta,
%         .Eigenvalues, .nModes, .SurfaceFile and either .ImagingKernel
%         or the face-based .PhiVertices via HM).
%
% OPTIONS (name-value):
%   'Alpha'   Regularisation parameter (default: 'auto' — estimated from
%             the noise floor of the top-3 high-k mode coefficients)
%   'Mode'    'poisson'   ρ_k = λ_k θ_k / (λ_k + α)   (default)
%             'laplacian' ρ_k = λ_k θ_k               (no regularisation;
%                                                       useful for diagnostics)
%   'HMFile'  Headmodel file to reload Phi (needed when Inv comes from
%             bst_eigenmode_analytic_inverse and PhiVertices are not cached
%             in the Inv struct).  Optional — inferred from Inv.SurfaceFile
%             if omitted.
%
% OUTPUT:
%   Results  Struct extending Inv with additional fields:
%     .RhoGridAmp       [nVert x nTime]  sharpened source field (complex if
%                                        Inv.ImageGridAmp was complex)
%     .RhoAmplitude     [nVert x nTime]  |ρ|  amplitude envelope
%     .RhoPhase         [nVert x nTime]  arg(ρ) instantaneous phase [rad]
%     .rho_k            [K x nTime]      sharpened eigenmode coefficients
%     .Alpha            scalar           regularisation α used
%     .SharpenMode      string           'poisson' | 'laplacian'
%
% EIGENMODE INTERPRETATION:
%   The sharpened field ρ(x,t) is the Poisson-deconvolved estimate of the
%   focal source distribution.  For a traveling wave, the wavefront in
%   ρ is geometrically sharper than in ŝ — the phase gradient ∇_S arg(ρ)
%   gives a more precise wave velocity estimate.
%
%   Sharpening is meaningful only within the observable subspace: modes
%   where θ_k is dominated by noise will still be noisy after sharpening.
%   The regularisation α prevents runaway amplification of those modes.
%
% FACE-BASED NOTE:
%   If Inv was produced from a face-based headmodel (isFaceBased=1), the
%   PhiVertices matrix from the HM is used for reconstruction — exactly as
%   in bst_eigenmode_analytic_inverse.
%
% SEE ALSO: bst_eigenmode_analytic_inverse, bst_eigenmode_oscillator_inverse
%
% Authors: Diellor Basha, 2026

%% ── Parse options ────────────────────────────────────────────────────────
Alpha    = 'auto';
Mode     = 'poisson';
HMFile   = '';
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'alpha',  Alpha  = varargin{k+1};
        case 'mode',   Mode   = lower(varargin{k+1});
        case 'hmfile', HMFile = varargin{k+1};
    end
end

%% ── Validate inputs ──────────────────────────────────────────────────────
required = {'theta', 'Eigenvalues', 'nModes'};
for r = required
    if ~isfield(Inv, r{1})
        error('bst_eigenmode_poisson_sharpen:MissingField', ...
            'Inv must have field .%s (from bst_eigenmode_analytic_inverse).', r{1});
    end
end

theta_k  = Inv.theta;       % [K x nTime]  complex (analytic) or real (Kalman)
lambdas  = double(Inv.Eigenvalues(:));   % [K x 1]
K        = Inv.nModes;
nTime    = size(theta_k, 2);

fprintf('Poisson sharpen: K=%d modes, mode=%s\n', K, Mode);

%% ── Estimate noise floor and α ──────────────────────────────────────────
% Auto-α: use the mean power of the top-3 highest eigenvalue (highest-k)
% coefficients as the noise floor.  These modes are the most strongly
% regularised by the inverse and are dominated by noise after recovery.
if ischar(Alpha) && strcmpi(Alpha, 'auto')
    % α must be on the scale of the eigenvalues for w_k = λ_k/(λ_k+α) to vary.
    % Strategy: estimate SNR per mode (signal power / noise floor), then set α
    % so that modes at the SNR=1 boundary get w_k = 0.5 (half-power point).
    % Noise floor: temporal variance of the LOWEST-amplitude mode (most regularised).
    [~, kord] = sort(mean(abs(theta_k).^2, 2), 'ascend');
    noise_floor = mean(abs(theta_k(kord(1:min(3,K)), :)).^2, 'all');
    snr_k = mean(abs(theta_k).^2, 2) / max(noise_floor, eps);  % [K x 1]
    % Knee: find smallest λ_k where SNR > 1 (reliably observed)
    reliable = snr_k >= 1;
    if any(reliable)
        lambda_knee = min(lambdas(reliable));
    else
        lambda_knee = min(lambdas);
    end
    % Set α = λ_knee so that w(λ_knee) = 0.5 and w(λ_max) approaches 1
    Alpha = lambda_knee;
    fprintf('  Auto α = %.3e  (λ at SNR≥1 knee; w_range=[%.2f %.2f])\n', ...
        Alpha, min(lambdas)/(min(lambdas)+Alpha), max(lambdas)/(max(lambdas)+Alpha));
else
    fprintf('  Manual α = %.3e\n', Alpha);
end

%% ── Poisson deconvolution in eigenmode space ─────────────────────────────
% ρ_k(t) = λ_k θ_k(t) / (λ_k + α)   [Poisson]
% ρ_k(t) = λ_k θ_k(t)                [pure Laplacian]

switch Mode
    case 'poisson'
        w_k   = lambdas ./ (lambdas + Alpha);   % [K x 1]  sharpening weights
    case 'laplacian'
        w_k   = lambdas;                         % [K x 1]  un-regularised
    otherwise
        error('bst_eigenmode_poisson_sharpen:UnknownMode', ...
            'Mode must be ''poisson'' or ''laplacian''.');
end

rho_k = w_k .* theta_k;   % [K x nTime]

fprintf('  θ_k power before: %.3e\n', mean(abs(theta_k(:)).^2));
fprintf('  ρ_k power after:  %.3e  (ratio %.1f×)\n', ...
    mean(abs(rho_k(:)).^2), mean(abs(rho_k(:)).^2) / max(mean(abs(theta_k(:)).^2), eps));

%% ── Reconstruct sharpened source field ───────────────────────────────────
% Load Phi_sel [nV x K] — try cached field first (fastest), then HMFile,
% then reload from SurfaceFile eigenmodes.
Phi_sel = [];
if isfield(Inv,'PhiVertices') && ~isempty(Inv.PhiVertices)
    % bst_eigenmode_analytic_inverse caches Phi(:,idx) here
    Phi_sel = double(Inv.PhiVertices);   % [nV x K]  already column-selected

elseif ~isempty(HMFile)
    HM = in_bst_headmodel(HMFile, 0);
    if isfield(HM,'isFaceBased') && HM.isFaceBased && isfield(HM,'PhiVertices')
        Phi_sel = double(HM.PhiVertices);   % [nV x K]
    else
        Eig = in_tess_eigenmodes(HM.SurfaceFile);
        Phi = double(Eig.Vectors);
        idx = double(HM.ModeIndices(1:K));
        Phi_sel = Phi(:, idx);
    end

elseif isfield(Inv,'SurfaceFile') && ~isempty(Inv.SurfaceFile)
    Eig = in_tess_eigenmodes(Inv.SurfaceFile);
    Phi = double(Eig.Vectors);
    % Match eigenvalues to global indices
    EigAll = double(Eig.Values(:));
    idx = zeros(K,1);
    for ki = 1:K
        [~, idx(ki)] = min(abs(EigAll - lambdas(ki)));
    end
    Phi_sel = Phi(:, idx);
end

if isempty(Phi_sel)
    error('bst_eigenmode_poisson_sharpen:NoPhi', ...
        'Cannot load Phi; pass HMFile explicitly or use bst_eigenmode_analytic_inverse (which caches PhiVertices).');
end

rho_field = Phi_sel * rho_k;   % [nV x nTime]  sharpened (complex or real)

%% ── Pack results ─────────────────────────────────────────────────────────
Results = Inv;
Results.RhoGridAmp   = rho_field;
Results.RhoAmplitude = abs(rho_field);
Results.RhoPhase     = angle(rho_field);
Results.rho_k        = rho_k;
Results.Alpha        = Alpha;
Results.SharpenMode  = Mode;

fprintf('  Done: RhoGridAmp [%d x %d]  |α|_mean=%.3e\n', ...
    size(rho_field,1), size(rho_field,2), Alpha);
end
