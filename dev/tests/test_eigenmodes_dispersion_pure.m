function test_eigenmodes_dispersion_pure
% Verify wave-vs-diffusion discrimination + parameter recovery on synthetic spectra.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% ---- Synthetic WAVE: peak frequency f* = c0*sqrt(lambda)/(2*pi), fixed bandwidth ----
K = 30; nF = 200; Freqs = linspace(1, 50, nF);
sl = linspace(0.3, 5, K)';          % sqrt(lambda) per mode
lambdas = sl.^2;
c0 = 30;                             % m/s
fpk = c0 * sl / (2*pi);             % peak freq (Hz), ~1.4..23.9
P = zeros(K, nF);
for k = 1:K
    P(k,:) = exp(-((Freqs - fpk(k)).^2) / (2 * 2^2));   % Gaussian peak, sigma 2 Hz
end
Out = bst_eigenmodes_dispersion(P, lambdas, Freqs);
assert(strcmp(Out.Regime, 'wave'), 'Expected wave regime, got %s.', Out.Regime);
assert(Out.R2wave > Out.R2diff, 'R2wave should exceed R2diff for a wave.');
assert(abs(Out.c - c0)/c0 < 0.15, 'Wave speed should recover c0 (got %.2f vs %.2f).', Out.c, c0);

% ---- Synthetic DIFFUSION: per-mode Lorentzian half-width gamma = alpha0*lambda ----
% NOTE: the Lorentzians peak at f=0, so with Freqs starting >0 the per-mode peak
% frequency is ~constant (R2wave -> 0); the diffusion win is driven by bandwidth
% scaling with lambda, which is the intended discriminator.
K2 = 30; nF2 = 200; Freqs2 = linspace(1, 50, nF2);
lambdas2 = linspace(0.5, 30, K2)';
alpha0 = 0.05;
P2 = zeros(K2, nF2);
for k = 1:K2
    g = alpha0 * lambdas2(k);                       % half-width (Hz)
    P2(k,:) = 1 ./ (Freqs2.^2 + g^2);               % Lorentzian (peaks at f=0)
end
Out2 = bst_eigenmodes_dispersion(P2, lambdas2, Freqs2);
assert(strcmp(Out2.Regime, 'diffusion'), 'Expected diffusion regime, got %s.', Out2.Regime);
assert(Out2.R2diff > Out2.R2wave, 'R2diff should exceed R2wave for diffusion.');

% ---- Shape / field checks ----
assert(isequal(size(Out.PeakFreq), [K 1]) && isequal(size(Out.Bandwidth), [K 1]), 'feature shapes');
assert(isfinite(Out.c) && isfinite(Out.alpha), 'c and alpha must be finite.');


% ---- MinPowerFrac drops low-power modes (use a spectrum with one weak mode) ----
Pv = ones(5, nF);
Pv(1,:) = 0.01 * Pv(1,:);          % mode 1 has ~1% of the others' total power
lamv = ((1:5)').^2;
Outv = bst_eigenmodes_dispersion(Pv, lamv, Freqs, 'MinPowerFrac', 0.5);
assert(Outv.Weights(1) == 0, 'Weakest mode should be dropped by MinPowerFrac.');
assert(all(Outv.Weights(2:5) > 0), 'Modes above the power floor should be kept.');

fprintf('ALL TESTS PASSED: test_eigenmodes_dispersion_pure\n');
end
