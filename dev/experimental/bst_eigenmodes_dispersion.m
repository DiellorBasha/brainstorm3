function Out = bst_eigenmodes_dispersion(P, lambdas, Freqs, varargin)
% BST_EIGENMODES_DISPERSION: Wave-vs-diffusion discrimination on a (lambda,omega) spectrum.
%
% USAGE:  Out = bst_eigenmodes_dispersion(P, lambdas, Freqs)
%         Out = bst_eigenmodes_dispersion(P, lambdas, Freqs, 'MinPowerFrac', 0)
%
% DESCRIPTION:
%     Discriminates the dynamical regime behind a stationary eigenmode
%     (lambda, omega) power spectrum by reducing each mode's spectrum to two
%     features and testing how they scale with the spatial frequency lambda:
%
%       Wave (dispersion curve omega = c*sqrt(lambda)): the PEAK frequency
%         f*(k) scales as sqrt(lambda_k). Weighted through-origin fit
%         f* = a*sqrt(lambda) gives speed c = 2*pi*a.
%       Diffusion (wedge; Lorentzian half-width ~ lambda): the spectral
%         BANDWIDTH w(k) scales as lambda_k. Weighted through-origin fit
%         w = b*lambda gives diffusivity alpha = 2*pi*b (proportional).
%
%     Regime = the model with the higher weighted R^2.
%     c is in m/s when lambdas are metric LBO eigenvalues (Brainstorm surfaces
%     are in metres); the regime decision and R^2 are unit-independent.
%
% INPUTS:
%     P       : [K x nFreq] non-negative power per eigenmode per frequency.
%     lambdas : [K x 1] eigenvalues (rad^2/m^2), aligned with rows of P.
%     Freqs   : [1 x nFreq] frequencies (Hz).
%
% OPTIONS (name-value):
%     'MinPowerFrac' : drop modes whose total power is below this fraction of
%                      the max per-mode total power (default 0 = keep all;
%                      power-weighting already down-weights weak modes).
%
% OUTPUTS:
%     Out.PeakFreq  [K x 1], Out.Bandwidth [K x 1], Out.Weights [K x 1]
%     Out.c (m/s), Out.alpha, Out.R2wave, Out.R2diff
%     Out.Regime ('wave'|'diffusion'), Out.Margin (|R2wave - R2diff|)
%
% SEE ALSO: bst_eigenmodes_wavelet, process_eigenmodes_dispersion

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Diellor Basha, 2026

%% ===== PARSE OPTIONS =====
MinPowerFrac = 0;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'minpowerfrac', MinPowerFrac = varargin{i+1};
    end
end

P       = max(double(P), 0);        % power is non-negative
lambdas = double(lambdas(:));        % [K x 1]
Freqs   = double(Freqs(:)');         % [1 x nFreq]

%% ===== PER-MODE FEATURES =====
p = sum(P, 2);                        % [K x 1] total power per mode (weight)
psafe = p; psafe(psafe == 0) = eps;
[~, ipk] = max(P, [], 2);            % peak-frequency index per mode
PeakFreq = Freqs(ipk)';             % [K x 1]
fbar = (P * Freqs') ./ psafe;        % power-weighted mean frequency [K x 1]
fvar = (P * (Freqs'.^2)) ./ psafe - fbar.^2;
Bandwidth = sqrt(max(fvar, 0));     % power-weighted std of frequency [K x 1]

%% ===== WEIGHTS (valid modes only) =====
valid = (lambdas > 0) & isfinite(PeakFreq) & isfinite(Bandwidth) & (p > MinPowerFrac * max(p));
w = p; w(~valid) = 0;

sl = sqrt(max(lambdas, 0));         % [K x 1]

%% ===== WAVE FIT (f* = a*sqrt(lambda)) =====
a = sum(w .* sl .* PeakFreq) / max(sum(w .* lambdas), eps);
c = 2*pi*a;
R2wave = weighted_r2(sl, PeakFreq, w);

%% ===== DIFFUSION FIT (w = b*lambda) =====
b = sum(w .* lambdas .* Bandwidth) / max(sum(w .* lambdas.^2), eps);
alpha = 2*pi*b;
R2diff = weighted_r2(lambdas, Bandwidth, w);

%% ===== REGIME =====
if R2wave >= R2diff
    Regime = 'wave';
else
    Regime = 'diffusion';
end

Out = struct('PeakFreq', PeakFreq, 'Bandwidth', Bandwidth, 'Weights', w, ...
             'c', c, 'alpha', alpha, 'R2wave', R2wave, 'R2diff', R2diff, ...
             'Regime', Regime, 'Margin', abs(R2wave - R2diff));
end


%% ===== WEIGHTED SQUARED PEARSON CORRELATION =====
function r2 = weighted_r2(x, y, w)
    sw = sum(w);
    if sw <= 0
        r2 = 0; return;
    end
    mx = sum(w .* x) / sw;
    my = sum(w .* y) / sw;
    cxy = sum(w .* (x - mx) .* (y - my));
    cxx = sum(w .* (x - mx).^2);
    cyy = sum(w .* (y - my).^2);
    denom = cxx * cyy;
    if denom <= 0
        r2 = 0;
    else
        r2 = cxy^2 / denom;
    end
end
