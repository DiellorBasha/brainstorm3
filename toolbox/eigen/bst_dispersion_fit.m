function [v_est, s_alpha, DispInfo] = bst_dispersion_fit(S, lambdas, f_ax, WTInfo, varargin)
% BST_DISPERSION_FIT  Estimate wave dispersion from (λ,ω) joint spectrum.
%
% USAGE:
%   [v_est, s_alpha, DispInfo] = bst_dispersion_fit(S, lambdas, f_ax, WTInfo)
%   [v_est, s_alpha, DispInfo] = bst_dispersion_fit(S, lambdas, f_ax, WTInfo, 'TargetFreq', 10)
%
% DESCRIPTION:
%   Fits the dispersion relation ω = v·√λ to the ridge of the (λ,ω) power
%   spectrum S from bst_lambda_omega_spectrum.  Returns the estimated phase
%   velocity v and the CWT scale index closest to the target frequency.
%
% INPUTS:
%   S        [K × nFreq]   power spectrum from bst_lambda_omega_spectrum
%   lambdas  [K × 1]       LBO eigenvalues [1/m²]
%   f_ax     [1 × nFreq]   temporal frequency axis [Hz]
%   WTInfo   struct from bst_sensor_cwt (.f_scales [nScales × 1])
%
% OPTIONS (name-value):
%   'TargetFreq'  Hz for selecting dominant scale (default: 10)
%   'Model'       'linear' | 'powerlaw'  dispersion model (default: 'linear')
%   'MinR2'       minimum R² to trust the fit (default: 0.3)
%   'Verbose'     logical (default: true)
%
% OUTPUTS:
%   v_est     scalar  estimated phase velocity [m/s] at TargetFreq
%   s_alpha   integer index into WTInfo.f_scales for the selected CWT scale
%   DispInfo  struct:
%     .v_linear    phase velocity from linear fit ω = v·√λ [m/s]
%     .v_at_target phase velocity interpolated at TargetFreq
%     .R2          goodness of linear dispersion fit
%     .k_ridge     [nFreq × 1] √λ at power ridge per frequency
%     .f_ridge     f_ax trimmed to ridge support
%     .s_alpha     CWT scale index
%     .f_alpha     actual CWT center frequency at s_alpha [Hz]
%     .model       dispersion model used
%
% Authors: Diellor Basha, 2026

%% ── Parse options ────────────────────────────────────────────────────────
TargetFreq    = 10;
Model         = 'linear';
MinR2         = 0.3;
VelocityRange = [0.5 15];   % m/s — same default as bst_lambda_omega_spectrum
Verbose       = true;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'targetfreq',    TargetFreq    = varargin{k+1};
        case 'model',         Model         = lower(varargin{k+1});
        case 'minr2',         MinR2         = varargin{k+1};
        case 'velocityrange', VelocityRange = varargin{k+1};
        case 'verbose',       Verbose       = varargin{k+1};
    end
end

K    = size(S, 1);
k_ax = sqrt(lambdas(:));   % [K × 1] spatial frequencies [1/m]

%% ── Speed-constrained ridge extraction ───────────────────────────────────
try
    S_sm = imgaussfilt(double(S), 1.5);
catch
    g    = fspecial('gaussian', [5 5], 1.5);
    S_sm = conv2(double(S), g, 'same');
end

nFreq   = numel(f_ax);
ridge_k_idx = ones(1, nFreq);
for fi = 1:nFreq
    k_lo = 2*pi * f_ax(fi) / max(VelocityRange(2), eps);
    k_hi = 2*pi * f_ax(fi) / max(VelocityRange(1), eps);
    mask = k_ax >= k_lo & k_ax <= k_hi;
    if ~any(mask), continue; end
    S_col = S_sm(:, fi);  S_col(~mask) = -Inf;
    [~, ridge_k_idx(fi)] = max(S_col);
end
f_ax    = f_ax(:)';                 % force row [1 × nFreq]
k_ridge = k_ax(ridge_k_idx);
k_ridge = k_ridge(:)';              % force row [1 × nFreq]

valid = ridge_k_idx > 1 & ridge_k_idx < K & f_ax > 0 & isfinite(k_ridge);

%% ── Fit dispersion relation ──────────────────────────────────────────────
f_v = f_ax(valid); k_v = k_ridge(valid);
omega_v = 2*pi * f_v;     % angular frequency [rad/s]

switch Model
    case 'linear'
        % ω = v · √λ  →  f_Hz = (v/2π) · k  →  v = 2π · f / k
        % Weighted OLS: v = (k' * omega) / (k' * k)
        v_linear = (k_v(:)' * omega_v(:)) / max(k_v(:)' * k_v(:), eps);
        f_fit    = v_linear / (2*pi) * k_ridge;
        residual = (f_ax(valid) - f_fit(valid)).^2;
        R2       = max(0, 1 - sum(residual) / max(sum((f_ax(valid) - mean(f_ax(valid))).^2), eps));
        v_est    = v_linear;

    case 'powerlaw'
        % ω = v · λ^α  →  log(ω) = log(v) + α·log(√λ)
        lk = log(max(k_v, eps)); lo = log(max(omega_v, eps));
        A  = [ones(numel(lk),1), lk(:)];
        c  = A \ lo(:);
        v_linear = exp(c(1));   alpha_exp = c(2);
        f_fit    = exp(c(1)) / (2*pi) * k_ridge.^c(2);
        residual = (f_ax(valid) - f_fit(valid)).^2;
        R2       = max(0, 1 - sum(residual) / max(sum((f_ax(valid) - mean(f_ax(valid))).^2), eps));
        v_est    = exp(c(1)) / (2*pi) * (2*pi*TargetFreq / exp(c(1)))^((c(2)-1)/c(2));

    otherwise
        error('bst_dispersion_fit:UnknownModel', 'Model must be ''linear'' or ''powerlaw''.');
end

% Phase velocity at target frequency (from linear relation)
k_at_target   = 2*pi * TargetFreq / max(v_linear, eps);   % k = ω/v
v_at_target   = 2*pi * TargetFreq / max(k_at_target, eps);

R2 = double(R2(1));  % guarantee scalar, no NaN from empty-valid edge case
if isnan(R2), R2 = 0; end
if R2 < MinR2 && Verbose
    warning('bst_dispersion_fit:WeakFit', ...
        'R²=%.2f < MinR2=%.2f — weak dispersion signature.', R2, MinR2);
end

%% ── Select CWT scale nearest TargetFreq ─────────────────────────────────
f_scales = WTInfo.f_scales;
[~, s_alpha] = min(abs(f_scales - TargetFreq));

%% ── Pack output ──────────────────────────────────────────────────────────
DispInfo.v_linear    = v_linear;
DispInfo.v_at_target = v_at_target;
DispInfo.R2          = R2;
DispInfo.k_ridge     = k_ridge;
DispInfo.f_ridge     = f_ax;
DispInfo.s_alpha     = s_alpha;
DispInfo.f_alpha     = f_scales(s_alpha);
DispInfo.model       = Model;
v_est = v_at_target;

if Verbose
    fprintf('bst_dispersion_fit: model=%s  v=%.2f m/s  R²=%.3f\n', Model, v_linear, R2);
    fprintf('  Selected CWT scale: s=%d  f=%.2f Hz (target %.1f Hz)\n', ...
        s_alpha, f_scales(s_alpha), TargetFreq);
end
end
