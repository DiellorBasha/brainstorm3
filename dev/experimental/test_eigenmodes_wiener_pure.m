function test_eigenmodes_wiener_pure
% Verify FFT-domain application of a per-mode magnitude gain G(k,f).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sfreq = 200; T = 4;
t = 0:(1/sfreq):(T - 1/sfreq);
nTime = numel(t);
K  = 3;
fA = 10; fB = 40;
Coeffs    = repmat(cos(2*pi*fA*t) + cos(2*pi*fB*t), K, 1);   % [K x nTime]
GainFreqs = 0:1:(sfreq/2);                                   % [1 x nGF], ascending
nGF       = numel(GainFreqs);

% Power at frequency f via in-phase/quadrature projection.
powAt = @(x, f) (mean(x .* cos(2*pi*f*t)))^2 + (mean(x .* sin(2*pi*f*t)))^2;

% --- Identity: Gain == 1 returns the input ---
Yid = bst_eigenmodes_wiener(Coeffs, sfreq, ones(K, nGF), GainFreqs);
assert(isequal(size(Yid), size(Coeffs)), 'Output size must match input.');
assert(isreal(Yid), 'Output must be real.');
assert(max(abs(Yid(:) - Coeffs(:))) < 1e-6, 'Identity gain must return the input.');

% --- Null: Gain == 0 returns ~0 ---
Y0 = bst_eigenmodes_wiener(Coeffs, sfreq, zeros(K, nGF), GainFreqs);
assert(max(abs(Y0(:))) < 1e-6, 'Zero gain must return ~0.');

% --- Selectivity: pass fA (gain 1 in 5..15 Hz), kill fB ---
Gsel = zeros(K, nGF);
Gsel(:, GainFreqs >= 5 & GainFreqs <= 15) = 1;
Ysel = bst_eigenmodes_wiener(Coeffs, sfreq, Gsel, GainFreqs);
pA_in  = powAt(Coeffs(1,:), fA); pB_in  = powAt(Coeffs(1,:), fB);
pA_out = powAt(Ysel(1,:),  fA);  pB_out = powAt(Ysel(1,:),  fB);
assert(pA_out > 0.5  * pA_in, 'Passband tone fA should be largely preserved.');
assert(pB_out < 0.05 * pB_in, 'Stopband tone fB should be largely removed.');

% --- Zero-phase: a single passband tone keeps its phase (no quadrature, no inversion) ---
xin   = cos(2*pi*fA*t);
Gpass = double(GainFreqs >= 5 & GainFreqs <= 15);   % [1 x nGF]
yph   = bst_eigenmodes_wiener(xin, sfreq, Gpass, GainFreqs);
c_cos = mean(yph .* cos(2*pi*fA*t));
c_sin = mean(yph .* sin(2*pi*fA*t));
assert(abs(c_sin) < 0.05 * abs(c_cos), 'Zero-phase: quadrature component must be negligible.');
assert(c_cos > 0, 'Zero-phase: in-phase component must stay positive (no inversion).');

% --- Errors ---
threwRow = false;
try, bst_eigenmodes_wiener(Coeffs, sfreq, ones(K+1, nGF), GainFreqs); catch, threwRow = true; end
assert(threwRow, 'Row mismatch between Gain and Coeffs should error.');
threwLen = false;
try, bst_eigenmodes_wiener(Coeffs, sfreq, ones(K, nGF), GainFreqs(1:end-1)); catch, threwLen = true; end
assert(threwLen, 'GainFreqs length not matching Gain columns should error.');
threwOrd = false;
try, bst_eigenmodes_wiener(Coeffs, sfreq, ones(K, nGF), fliplr(GainFreqs)); catch, threwOrd = true; end
assert(threwOrd, 'Non-ascending GainFreqs should error.');

fprintf('ALL TESTS PASSED: test_eigenmodes_wiener_pure\n');
end
