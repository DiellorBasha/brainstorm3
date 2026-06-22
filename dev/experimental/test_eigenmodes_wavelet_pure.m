function test_eigenmodes_wavelet_pure
% Verify the complex Morlet wavelet tensor: shape, complexity, amplitude
% localization, phase rate, and the default frequency grid.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sfreq = 200; T = 2; nT = sfreq*T; t = (0:nT-1)/sfreq;
f0 = 10; K = 8; kSig = 3;
Coeffs = zeros(K, nT);
Coeffs(kSig,:) = cos(2*pi*f0*t);          % one mode carries a 10 Hz sinusoid
Freqs = 2:2:40;                            % grid includes 10 Hz

[W, Fout] = bst_eigenmodes_wavelet(Coeffs, sfreq, Freqs);

assert(isequal(size(W), [K, nT, numel(Freqs)]), 'W must be [K x nTime x nFreq].');
assert(~isreal(W), 'W must be complex.');
assert(isequal(Fout(:)', Freqs), 'Freqs must pass through unchanged.');

% Amplitude localizes to (kSig, f0).
A = reshape(mean(abs(W), 2), size(W,1), []);   % [K x nFreq], time-averaged amplitude (shape-safe for K=1)
[~, imax] = max(A(:)); [kmax, fmax] = ind2sub(size(A), imax);
assert(kmax == kSig, 'Peak mode wrong (got %d, expected %d).', kmax, kSig);
[~, if0] = min(abs(Freqs - f0));
assert(fmax == if0, 'Peak frequency wrong (got bin %d).', fmax);
others = setdiff(1:K, kSig);
assert(max(A(others, if0)) < 0.1 * A(kSig, if0), 'Signal leaked to other modes.');

% Phase advances at ~2*pi*f0/sfreq per sample in the central (edge-free) region.
wsig = squeeze(W(kSig, :, if0));
ph = unwrap(angle(wsig(:)'));
mid = round(nT*0.4):round(nT*0.6);
dph = mean(diff(ph(mid)));
assert(abs(dph - 2*pi*f0/sfreq) < 0.2*(2*pi*f0/sfreq), 'Phase rate wrong (got %.4f).', dph);

% Empty Freqs -> default 40-frequency log grid.
[W2, F2] = bst_eigenmodes_wavelet(Coeffs, sfreq, []);
assert(numel(F2) == 40, 'Default grid must have 40 frequencies.');
assert(isequal(size(W2), [K, nT, 40]), 'Default tensor must be [K x nTime x 40].');

% Single-mode case (K=1) preserves the [1 x nTime x nFreq] complex shape.
[W1, ~] = bst_eigenmodes_wavelet(cos(2*pi*f0*t), sfreq, Freqs);
assert(isequal(size(W1), [1, nT, numel(Freqs)]), 'K=1 tensor must be [1 x nTime x nFreq].');
assert(~isreal(W1), 'K=1 tensor must be complex.');

fprintf('ALL TESTS PASSED: test_eigenmodes_wavelet_pure\n');
end
