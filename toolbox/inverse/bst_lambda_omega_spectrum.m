function [S, f_ax, k_ax, lambda_ax, peak_info] = bst_lambda_omega_spectrum(Theta, lambdas, WTInfo, varargin)
% BST_LAMBDA_OMEGA_SPECTRUM  Joint eigenvalue-frequency (lambda,omega) power
%   spectrum from CWT eigenmode coefficients — reveals the wave dispersion.
%
% USAGE:
%   [S, f_ax, k_ax, lambda_ax, peak_info] = bst_lambda_omega_spectrum(Theta, lambdas, WTInfo)
%   [S, f_ax, k_ax, lambda_ax, peak_info] = bst_lambda_omega_spectrum(..., 'FreqRange', [4 30])
%
% INPUTS:
%   Theta    [K x nScales x nTime] complex (from bst_eigenmode_cwt_inverse)
%            OR [K x nTime] complex (analytic inverse, no CWT dimension)
%   lambdas  [K x 1] LBO eigenvalues [1/m^2]
%   WTInfo   struct from bst_sensor_cwt (.Fs, .f_scales, .nScales)
%            OR scalar Fs [Hz] when Theta is [K x nTime]
%
% OPTIONS (name-value):
%   'FreqRange'     [f_lo f_hi] Hz (default: [1 30])
%   'Normalize'     'none'|'mode'|'freq' (default: 'mode')
%   'SmoothSigma'   Gaussian smoothing sigma before peak detection (default: 1.5)
%   'VelocityRange' [v_lo v_hi] m/s for speed-constrained ridge (default: [0.5 15])
%                   At each frequency f, only modes where 2πf/√λ ∈ VelocityRange
%                   are candidates for the ridge. Prevents noise-dominated
%                   high-λ modes from hijacking the peak.
%   'Verbose'       logical (default: true)
%
% OUTPUTS:
%   S          [K x nFreq] real — power spectrum S(k,omega)
%   f_ax       [1 x nFreq] Hz;  k_ax [K x 1] = sqrt(lambdas) [1/m]
%   lambda_ax  [K x 1] = lambdas;  peak_info struct: .ridge_k, .ridge_kax,
%              .f_peak, .k_peak, .v_estimate, .v_linear, .ridge_R2
%
% SEE ALSO: bst_sensor_cwt, bst_eigenmode_analytic_inverse
%
% Authors: Diellor Basha, 2026

%% -- Parse options --
FreqRange     = [1 30];
Normalize     = 'mode';
SmoothSigma   = 1.5;
VelocityRange = [0.5 15];   % m/s: physically plausible wave speeds for brain
Verbose       = true;
for ki = 1:2:numel(varargin)
    switch lower(varargin{ki})
        case 'freqrange',     FreqRange     = varargin{ki+1};
        case 'normalize',     Normalize     = lower(varargin{ki+1});
        case 'smoothsigma',   SmoothSigma   = varargin{ki+1};
        case 'velocityrange', VelocityRange = varargin{ki+1};
        case 'verbose',       Verbose       = logical(varargin{ki+1});
    end
end

%% -- Resolve Fs and Theta shape --
if isstruct(WTInfo)
    Fs = WTInfo.Fs;
else
    Fs = double(WTInfo);
end
K      = size(Theta, 1);
hasCWT = (ndims(Theta) == 3);

%% -- Step 1: per-mode temporal power spectrum --
if hasCWT
    nScales = size(Theta, 2);
    nTime   = size(Theta, 3);
    nFFT    = 2^nextpow2(nTime);
    S       = zeros(K, nFFT/2);
    for s = 1:nScales
        Th_s = squeeze(Theta(:, s, :));          % [K x nTime]
        Sf   = abs(fft(double(Th_s), nFFT, 2)).^2;
        S    = S + Sf(:, 1:nFFT/2);
    end
    S = S / nScales;
else
    nTime = size(Theta, 2);
    nFFT  = 2^nextpow2(nTime);
    Sf    = abs(fft(double(Theta), nFFT, 2)).^2;
    S     = Sf(:, 1:nFFT/2);
end

%% -- Step 2: frequency axis and range trim --
f_ax   = (0:nFFT/2-1) * (Fs / nFFT);          % [1 x nFFT/2]
f_mask = f_ax >= FreqRange(1) & f_ax <= FreqRange(2);
S      = S(:, f_mask);
f_ax   = f_ax(f_mask);

%% -- Step 3: normalize --
switch Normalize
    case 'mode'
        S = S ./ max(max(S, [], 2), eps);
    case 'freq'
        S = S ./ max(max(S, [], 1), eps);
    case 'none'
    otherwise
        error('bst_lambda_omega_spectrum:BadNormalize', ...
            'Normalize must be ''none'', ''mode'', or ''freq''.');
end

%% -- Step 4: spatial frequency axis --
lambdas   = double(lambdas(:));                % [K x 1]
k_ax      = sqrt(lambdas);                     % [K x 1] in 1/m
lambda_ax = lambdas;

%% -- Step 5: smooth and extract ridge --
if license('test', 'image_toolbox') && ~isempty(which('imgaussfilt'))
    S_smooth = imgaussfilt(S, SmoothSigma);
else
    hw = ceil(3 * SmoothSigma);
    gv = exp(-(-hw:hw).^2 / (2*SmoothSigma^2));  gv = gv / sum(gv);
    S_smooth = conv2(gv, gv', S, 'same');
end

% Speed-constrained ridge: at each frequency f, only consider modes where
% v = 2πf/√λ ∈ VelocityRange — prevents high-λ noise modes from dominating.
% For frequency f_j and spatial frequency k = √λ_i:  v = 2π*f_j / k_i
% The mode mask is: k_i ∈ [2π*f_j/v_max,  2π*f_j/v_min]
nFreq   = numel(f_ax);
ridge_k   = ones(nFreq, 1);    % [nFreq x 1] local mode index
ridge_kax = nan(nFreq, 1);
for fi = 1:nFreq
    k_lo = 2*pi * f_ax(fi) / max(VelocityRange(2), eps);   % k for v_max
    k_hi = 2*pi * f_ax(fi) / max(VelocityRange(1), eps);   % k for v_min
    speed_mask = (k_ax >= k_lo) & (k_ax <= k_hi);
    if ~any(speed_mask)
        continue;
    end
    S_col = S_smooth(:, fi);
    S_col(~speed_mask) = -Inf;
    [~, ridge_k(fi)] = max(S_col);
    ridge_kax(fi)    = k_ax(ridge_k(fi));
end

% Global peak within speed-constrained band
valid_fi = ~isnan(ridge_kax);
if any(valid_fi)
    [~, i_f]   = max(max(S_smooth .* (k_ax >= 2*pi*f_ax/(max(VelocityRange)+eps) & ...
                                       k_ax <= 2*pi*f_ax/(max(VelocityRange(1),eps)+eps)), [], 1));
    % Simpler: find peak in speed-masked spectrum
    S_masked = S_smooth;
    for fi2 = 1:nFreq
        k_lo2 = 2*pi*f_ax(fi2)/max(VelocityRange(2),eps);
        k_hi2 = 2*pi*f_ax(fi2)/max(VelocityRange(1),eps);
        S_masked(~(k_ax>=k_lo2 & k_ax<=k_hi2), fi2) = 0;
    end
    [~, i_f] = max(max(S_masked, [], 1));
    f_peak   = f_ax(i_f);
    [~, i_k] = max(S_masked(:, i_f));
    k_peak   = k_ax(i_k);
else
    [~, i_f] = max(max(S_smooth, [], 1));
    f_peak   = f_ax(i_f);
    [~, i_k] = max(S_smooth(:, i_f));
    k_peak   = k_ax(i_k);
end
v_estimate = 2*pi * f_peak / max(k_peak, eps);

% Linear dispersion fit ω = v·√λ on the constrained ridge
A     = ridge_kax(valid_fi);
b     = f_ax(valid_fi)';
if numel(A) >= 2
    v_lin = (A(:)' * b(:)) / max(A(:)' * A(:), eps);   % f = (v/2π)·k → v/2π
else
    v_lin = v_estimate / (2*pi);
end
f_fit  = v_lin * ridge_kax(valid_fi);
SS_res = sum((f_ax(valid_fi)' - f_fit(:)).^2);
SS_tot = sum((f_ax(valid_fi)' - mean(f_ax(valid_fi)')).^2);
R2     = double(max(0, 1 - SS_res / max(SS_tot, eps)));
if isnan(R2), R2 = 0; end

%% -- Step 6: pack peak_info --
peak_info.ridge_k    = ridge_k;
peak_info.ridge_kax  = ridge_kax;
peak_info.f_peak     = f_peak;
peak_info.k_peak     = k_peak;
peak_info.v_estimate = v_estimate;
peak_info.v_linear   = v_lin;
peak_info.ridge_R2   = R2;

if Verbose
    fprintf('lambda-omega spectrum: K=%d modes, f=[%g..%g] Hz\n', K, FreqRange(1), FreqRange(2));
    fprintf('  Peak: f=%.2g Hz, sqrt(lam)=%.3g m^-1 -> v_est=%.3g m/s\n', f_peak, k_peak, v_estimate);
    fprintf('  Ridge linear fit (w=v*sqrt(lam)): v=%.3g m/s  R2=%.2f\n', v_lin * 2*pi, R2);
end
end
