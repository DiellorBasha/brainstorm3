function test_eigenmode_lever_integration
% Headless end-to-end on a synthetic mesh: ApplyToColumn must equal the
% analytic band-limited reconstruction, be a no-op when inactive or on a
% mismatched surface, and round-trip a band-limited field at full band.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(fullfile(repoRoot, 'toolbox', 'gui'));
addpath(fullfile(repoRoot, 'toolbox', 'math'));
addpath(fullfile(repoRoot, 'toolbox', 'anatomy'));
global GlobalData; %#ok<GVMIS>

% Small closed mesh; eigenmodes struct + mass matrix come back together
[V, F] = tess_sphere(642);
[Eig, ~, M] = tess_eigenmodes(V, F, 'nModes', 60, 'MassType', 'barycentric', ...
                              'RemoveDC', 1, 'Verbose', 0);
K = Eig.nModes;                       % actual count after RemoveDC

% Inject the cache directly (ApplyToColumn uses it when SurfaceFile matches)
SurfaceFile = '/synthetic/sphere.mat';
panel_eigenmodes('ResetState', SurfaceFile, K);
panel_eigenmodes('SetCache', SurfaceFile, Eig, M);

% A band-limited field (lives in span(Phi)) so full-band reconstruction is exact
rng(1);
u = Eig.Vectors * randn(K, 1);

% --- Inactive => no-op ---
panel_eigenmodes('SetActive', 0);
uF = panel_eigenmodes('ApplyToColumn', SurfaceFile, u);
assert(isequal(uF, u), 'inactive lever must be a no-op');

% --- Active, full-band box => round-trip identity (u in span(Phi)) ---
panel_eigenmodes('SetActive', 1);
panel_eigenmodes('SetWindowShape', 'box');
panel_eigenmodes('SetBand', 1, K);
uF = panel_eigenmodes('ApplyToColumn', SurfaceFile, u);
assert(max(abs(uF - u)) < 1e-6 * max(abs(u)), 'full band must round-trip a band-limited field to u');

% --- Active, low band => equals analytic reconstruction and shrinks energy ---
panel_eigenmodes('SetBand', 1, 10);
uLow = panel_eigenmodes('ApplyToColumn', SurfaceFile, u);
W = panel_eigenmodes('GetWeights');
analytic = Eig.Vectors * (W(:) .* (Eig.Vectors' * (M * u)));
assert(max(abs(uLow - analytic)) < 1e-9, 'ApplyToColumn must equal analytic band-limited reconstruction');
assert(norm(uLow) < norm(u), 'low-band reconstruction must shrink overall energy');

% --- Surface mismatch => no-op ---
uOther = panel_eigenmodes('ApplyToColumn', '/synthetic/other.mat', u);
assert(isequal(uOther, u), 'mismatched surface must be a no-op');

% IsActive query reflects active flag + surface match
panel_eigenmodes('SetActive', 1);
assert(panel_eigenmodes('IsActive', SurfaceFile), 'IsActive true when active+matching');
assert(~panel_eigenmodes('IsActive', '/synthetic/other.mat'), 'IsActive false on surface mismatch');
panel_eigenmodes('SetActive', 0);
assert(~panel_eigenmodes('IsActive', SurfaceFile), 'IsActive false when inactive');

fprintf('ALL TESTS PASSED: test_eigenmode_lever_integration\n');
end
