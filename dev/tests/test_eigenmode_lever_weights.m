function test_eigenmode_lever_weights
% Pure-math tests for panel_eigenmodes('BuildWeights', ...).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(fullfile(repoRoot, 'toolbox', 'gui'));

K = 200;

% single -> delta at center
w = panel_eigenmodes('BuildWeights', 'single', 1, K, 42, K);
assert(isequal(size(w), [1 K]), 'weights must be 1xK row');
assert(w(42) == 1 && sum(w) == 1, 'single must be a delta at the center');

% box -> 1 on [kLo,kHi], 0 elsewhere
w = panel_eigenmodes('BuildWeights', 'box', 30, 55, 42, K);
assert(all(w(30:55) == 1), 'box interior must be 1');
assert(sum(w) == 26, 'box must keep exactly 26 modes');
assert(w(29) == 0 && w(56) == 0, 'box must be 0 outside the band');

% tapered -> <=1 everywhere, interior 1, monotone shoulders, 0 outside
w = panel_eigenmodes('BuildWeights', 'tapered', 30, 55, 42, K);
assert(all(w >= 0 & w <= 1), 'tapered weights in [0,1]');
assert(w(29) == 0 && w(56) == 0, 'tapered 0 outside band');
mid = round((30+55)/2);
assert(abs(w(mid) - 1) < 1e-9, 'tapered interior reaches 1');
assert(w(30) <= w(31) && w(31) <= w(mid), 'tapered rising shoulder is monotone');

% gain -> Gaussian bell, peak at center, symmetric falloff, in (0,1]
w = panel_eigenmodes('BuildWeights', 'gain', 30, 55, 42, K);
assert(all(w >= 0 & w <= 1 + 1e-12), 'gain weights in [0,1]');
[~, ipk] = max(w);
assert(ipk == 42, 'gain peak at the center mode');
assert(w(42) > w(20) && w(42) > w(64), 'gain falls off away from center');

% tapered endpoints taper to 0 (open-interval pass-band)
w = panel_eigenmodes('BuildWeights', 'tapered', 30, 55, 42, K);
assert(w(30) == 0 && w(55) == 0, 'tapered tapers to 0 at the band endpoints');

% tapered degenerates to a delta when kLo==kHi (n=1 guard)
w = panel_eigenmodes('BuildWeights', 'tapered', 42, 42, 42, K);
assert(w(42) == 1 && sum(w) == 1, 'tapered n=1 must be a delta');

% tapered n=2 must not blank the band (n<=2 guard)
w = panel_eigenmodes('BuildWeights', 'tapered', 42, 43, 42, K);
assert(all(w(42:43) == 1), 'tapered n=2 must keep both modes (no silent zeros)');

% gain kLo==kHi falls back to a sigma=1 bell centered at the mode
w = panel_eigenmodes('BuildWeights', 'gain', 42, 42, 42, K);
assert(w(42) > w(41) && w(42) > w(43), 'gain kLo==kHi -> sigma=1 bell at center');

% clamping inside the builder: out-of-range band is clamped to [1,K]
w = panel_eigenmodes('BuildWeights', 'box', -5, K+100, 42, K);
assert(all(w == 1), 'box clamped to full range keeps all modes');

fprintf('ALL TESTS PASSED: test_eigenmode_lever_weights\n');
end
