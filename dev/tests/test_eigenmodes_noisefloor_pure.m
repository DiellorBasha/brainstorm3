function test_eigenmodes_noisefloor_pure
% Verify SNR / power-subtraction / Wiener-gain / reliable-mode-cutoff math.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

rng(3);
K = 30; nF = 40;
S = zeros(K, nF); S(1:10,:) = rand(10, nF) + 0.5;   % signal only in low modes
N = 0.2*rand(K, nF) + 0.1;                            % positive noise floor
Pdata = S + N;

Out = bst_eigenmodes_noisefloor(Pdata, N);            % defaults: Alpha=1, Floor=0, SnrThresh=1

% Power subtraction recovers the signal power exactly (Alpha=1, Floor=0).
assert(max(abs(Out.CleanPSD(:) - S(:))) < 1e-12, 'CleanPSD should equal S.');
% SNR = (S+N)/N.
SNRexp = Pdata ./ N;
assert(max(abs(Out.SNR(:) - SNRexp(:))) < 1e-12, 'SNR mismatch.');
% Wiener gain in [0,1], and exactly 0 where Pdata <= N (modes 11..30 have S=0 -> equal).
assert(all(Out.Gain(:) >= 0 & Out.Gain(:) <= 1), 'Gain must be in [0,1].');
below = (Pdata <= N);
assert(all(Out.Gain(below) == 0), 'Gain must be 0 where Pdata<=N.');

% Floor clamp respected with over-subtraction.
Out2 = bst_eigenmodes_noisefloor(Pdata, N, 'Alpha', 2, 'Floor', 0.1);
assert(all(Out2.CleanPSD(:) >= 0.1*N(:) - 1e-12), 'CleanPSD must respect the spectral floor.');

% Reliable-mode cutoff on a descending-SNR ramp at frequency 1.
Pramp = N;
Pramp(:,1) = N(:,1) .* (3 - (0:K-1)'*0.1);            % SNR(k,1) = 3 - 0.1*(k-1)
Outr = bst_eigenmodes_noisefloor(Pramp, N, 'SnrThresh', 1);
% SNR>=1  <=>  3 - 0.1*(k-1) >= 1  <=>  k <= 21
assert(Outr.Kstar(1) == 21, 'Kstar ramp wrong (got %d, expected 21).', Outr.Kstar(1));


% Reliable-mode cutoff is 0 when no mode meets the threshold (SNR=1 everywhere < 2).
OutZ = bst_eigenmodes_noisefloor(N, N, 'SnrThresh', 2);
assert(all(OutZ.Kstar == 0), 'Kstar must be 0 when no mode meets the threshold.');

fprintf('ALL TESTS PASSED: test_eigenmodes_noisefloor_pure\n');
end
