function test_eigfilter_design_pure
% Pure logic for the eigenmode filter-design panel: kernel->paired-weights,
% slider-range heuristic, and the delta point-spread (impulse response).
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3');
if ~brainstorm('status'); brainstorm nogui; end
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3/toolbox/math/eigfilter');

% Synthetic 2-component eigenmodes (LH verts 1:4, RH verts 5:8), 4 paired ranks.
nV = 8; m = 4;
Vectors = zeros(nV, 2*m);
for k = 1:m
    for n = 1:4; Vectors(n, k)     = cos(pi/4*(n-0.5)*(k-1)); end
    for n = 1:4; Vectors(4+n, m+k) = cos(pi/4*(n-0.5)*(k-1)); end
end
Values   = [1;4;9;16;  1;4;9;16];          % both components share the same spectrum
Eig = struct('Vectors',Vectors,'Values',Values,'nModes',2*m, ...
    'Component',[ones(m,1);2*ones(m,1)],'CompRank',[(1:m)';(1:m)'],'nComponents',2);

% --- KernelPairedWeights: W(k) = g(lambda_paired_k) ---
[W, g] = panel_eigenmodes('KernelPairedWeights', Eig, 'heat', struct('t',0.1));
assert(isequal(size(W),[1 m]), 'W must be 1 x Kpaired.');
lamPaired = [1;4;9;16];
assert(max(abs(W(:) - exp(-0.1*lamPaired))) < 1e-12, 'paired weights must equal g(lambda_paired).');
assert(isa(g,'function_handle'), 'must return the kernel handle.');

% --- KernelSliderRange: log heuristic around the default, clamped to range ---
[smin, smax] = panel_eigenmodes('KernelSliderRange', 0.01, [0 Inf]);
assert(smin > 0 && smin < 0.01 && smax > 0.01, 'range must bracket the default.');

% --- DeltaPointSpread: heat point-spread is smooth, concentrated, finite ---
M = speye(nV);                                   % identity mass for the test
[~, gf] = panel_eigenmodes('KernelPairedWeights', Eig, 'heat', struct('t',0.05));
ps = panel_eigenmodes('DeltaPointSpread', Eig, M, gf, 3);   % impulse at vertex 3 (LH)
assert(isequal(size(ps),[nV 1]), 'point-spread must be [nVert x 1].');
assert(all(isfinite(ps)), 'point-spread must be finite.');
assert(abs(ps(3)) == max(abs(ps)), 'heat point-spread must peak at the impulse vertex.');

% --- GetDisplayColumn: synthesis = PairedGrid*W; delta = point-spread ---
% Set up panel state for the synthetic surface (no real DB needed: use a fake file).
global GlobalData;
panel_eigenmodes('ResetState', 'fake_surf.mat', m);
GlobalData.UserModes.CacheSurfaceFile = 'fake_surf.mat';
GlobalData.UserModes.CacheEig  = Eig;
GlobalData.UserModes.CacheMass = speye(nV);
GlobalData.UserModes.WeightMode = 'kernel';
GlobalData.UserModes.KernelName = 'heat';
GlobalData.UserModes.KernelParams = struct('t',0.05);
panel_eigenmodes('RecomputeWeights');   % fills Weights + KernelFn

PairedGrid = zeros(nV, m);
for k = 1:m; PairedGrid(:,k) = sum(Eig.Vectors(:, Eig.CompRank==k), 2); end

GlobalData.UserModes.DisplayMode = 'synthesis';
colS = panel_eigenmodes('GetDisplayColumn', 'fake_surf.mat', PairedGrid);
assert(isequal(size(colS),[nV 1]), 'synthesis column shape.');
assert(max(abs(colS - PairedGrid*GlobalData.UserModes.Weights(:))) < 1e-12, 'synthesis = PairedGrid*W.');

GlobalData.UserModes.DisplayMode = 'delta';
GlobalData.UserModes.DeltaVertex = 6;   % RH vertex
colD = panel_eigenmodes('GetDisplayColumn', 'fake_surf.mat', PairedGrid);
assert(isequal(size(colD),[nV 1]) && all(isfinite(colD)), 'delta column shape/finite.');
assert(abs(colD(6)) == max(abs(colD)), 'delta point-spread peaks at the impulse vertex.');

% --- ParamRange: scalar params adjustable; ideal's 2-element 'band' is not (no crash) ---
mh = bst_eigfilter_kernel('info', 'heat');
[smn, smx, okh] = panel_eigenmodes('ParamRange', mh.params.t);
assert(okh && smn > 0 && smx > smn, 'heat t must be adjustable with valid bounds.');
mi = bst_eigfilter_kernel('info', 'ideal');
[~, ~, oki] = panel_eigenmodes('ParamRange', mi.params.band);
assert(~oki, 'ideal band (non-scalar default, no range field) must be flagged non-adjustable.');

disp('ALL TESTS PASSED');
end
