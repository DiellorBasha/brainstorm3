function test_eigenmodes_filter_gain_pure
% Verify each per-mode transfer-function gain.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

lam = ((0:9)').^2;          % lambda_0=0, strictly increasing
K = numel(lam);

% lowpass: keep modes 1..4
hLP = bst_eigenmodes_filter_gain(lam, 'lowpass', 'CutoffMode', 4);
assert(isequal(hLP, [ones(4,1); zeros(K-4,1)]), 'lowpass mask wrong.');
% highpass: keep modes 4..end
hHP = bst_eigenmodes_filter_gain(lam, 'highpass', 'CutoffMode', 4);
assert(isequal(hHP, [zeros(3,1); ones(K-3,1)]), 'highpass mask wrong.');
% bandpass: keep modes 3..6
hBP = bst_eigenmodes_filter_gain(lam, 'bandpass', 'ModeRange', [3 6]);
expBP = zeros(K,1); expBP(3:6) = 1;
assert(isequal(hBP, expBP), 'bandpass mask wrong.');
% heat: t->0 is ~identity; large t suppresses high lambda monotonically
hH0 = bst_eigenmodes_filter_gain(lam, 'heat', 'DiffusionTime', 1e-12);
assert(max(abs(hH0 - 1)) < 1e-6, 'heat t->0 should be ~1.');
hH1 = bst_eigenmodes_filter_gain(lam, 'heat', 'DiffusionTime', 1);
assert(abs(hH1(1)-1) < 1e-12 && hH1(end) < 1e-6 && all(diff(hH1) <= 0), 'heat should decay with lambda.');
% inverse_heat: clamped at MaxGain
hIH = bst_eigenmodes_filter_gain(lam, 'inverse_heat', 'DiffusionTime', 1, 'MaxGain', 5);
assert(all(hIH <= 5 + 1e-9) && any(abs(hIH - 5) < 1e-9), 'inverse_heat should clamp at MaxGain.');
% tikhonov: 1 at lambda=0, decreasing, in (0,1]
hT = bst_eigenmodes_filter_gain(lam, 'tikhonov', 'RegBeta', 1);
assert(abs(hT(1)-1) < 1e-12 && all(hT > 0 & hT <= 1) && all(diff(hT) <= 0), 'tikhonov shape wrong.');
% custom: passthrough of the supplied handle
hC = bst_eigenmodes_filter_gain(lam, 'custom', 'TransferFn', @(l) 2*ones(size(l)));
assert(all(hC == 2), 'custom passthrough wrong.');
% unknown type errors
threw = false;
try, bst_eigenmodes_filter_gain(lam, 'bogus'); catch, threw = true; end
assert(threw, 'unknown filter type should error.');


% Guard error paths.
threwH = false;  try, bst_eigenmodes_filter_gain(lam, 'heat', 'DiffusionTime', 0);          catch, threwH = true;  end
assert(threwH, 'heat with DiffusionTime<=0 should error.');
threwIH = false; try, bst_eigenmodes_filter_gain(lam, 'inverse_heat', 'DiffusionTime', 0);  catch, threwIH = true; end
assert(threwIH, 'inverse_heat with DiffusionTime<=0 should error.');
threwT = false;  try, bst_eigenmodes_filter_gain(lam, 'tikhonov', 'RegBeta', -1);           catch, threwT = true;  end
assert(threwT, 'tikhonov with RegBeta<0 should error.');
threwC = false;  try, bst_eigenmodes_filter_gain(lam, 'custom', 'TransferFn', @(l) [1;2]);  catch, threwC = true;  end
assert(threwC, 'custom with wrong-length output should error.');

fprintf('ALL TESTS PASSED: test_eigenmodes_filter_gain_pure\n');
end
