function test_io_eigenmodes_roundtrip
% Verify out_tess_eigenmodes -> in_tess_eigenmodes round-trips an embedded struct.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

[V, F] = tess_sphere(162);
nV = size(V, 1);

Em = struct();
Em.Vectors     = randn(nV, 15);
Em.Values      = sort(rand(15, 1));
Em.nModes      = 15;
Em.MassType    = 'barycentric';
Em.Sigma       = -1e-8;
Em.Tolerance   = 1e-10;
Em.nRemoved    = 2;
Em.ComputeTime = 1.23;
Em.Component   = [ones(8,1); 2*ones(7,1)];   % two components (non-default, to make persistence assertions real)
Em.CompRank    = [(1:8)'; (1:7)'];            % per-component rank (non-default)
Em.nComponents = 2;
[Em.Laplacian, Em.MassMatrix] = tess_laplacian(V, F, 'MassType', Em.MassType);   % whole-mesh operators (sparse)
Em.LaplacianType = 'Laplace-Beltrami';

tmpDir = tempname; mkdir(tmpDir);
cleanup = onCleanup(@() rmdir(tmpDir, 's'));

% Surface WITH eigenmodes (file must be named like a Brainstorm tess file).
SurfFile = fullfile(tmpDir, 'tess_cortex_test.mat');
bst_save(SurfFile, struct('Vertices', V, 'Faces', F, 'Comment', 'test cortex'), 'v7');
out_tess_eigenmodes(SurfFile, Em, V, F, false);

[Em2, isComputed] = in_tess_eigenmodes(SurfFile);
assert(isComputed, 'in_tess_eigenmodes reported not computed.');
assert(Em2.nModes == 15, 'nModes mismatch (%d).', Em2.nModes);
assert(isa(Em2.Vectors, 'double'), 'Vectors not converted to double.');
assert(isequal(size(Em2.Vectors), [nV, 15]), 'Vectors size mismatch.');
assert(max(abs(Em2.Values - Em.Values)) < 1e-12, 'Values changed on round-trip.');
assert(max(abs(Em2.Vectors(:) - Em.Vectors(:))) < 1e-5, 'Vectors exceeded single-precision round-trip error.');
fprintf('PASSED: round-trip preserves eigenmodes.\n');

% Surface WITHOUT eigenmodes -> isComputed false.
SurfFile2 = fullfile(tmpDir, 'tess_cortex_empty.mat');
bst_save(SurfFile2, struct('Vertices', V, 'Faces', F, 'Comment', 'empty'), 'v7');
[~, isComp2] = in_tess_eigenmodes(SurfFile2);
assert(~isComp2, 'Empty surface wrongly reported computed.');
fprintf('PASSED: missing eigenmodes report isComputed=false.\n');

% ----- New per-hemisphere metadata round-trips -----
assert(isfield(Em2, 'Component') && isfield(Em2, 'CompRank'), ...
    'Component/CompRank must persist through save/load.');
assert(isequal(Em2.Component(:), Em.Component(:)), 'Component not preserved.');
assert(isequal(Em2.CompRank(:),  Em.CompRank(:)),  'CompRank not preserved.');
assert(isequal(Em2.nComponents, Em.nComponents), 'nComponents not preserved.');
fprintf('PASSED: Component/CompRank survive round-trip.\n');

% ----- Whole-mesh operators (MassMatrix, Laplacian) round-trip -----
assert(isfield(Em2, 'MassMatrix') && isfield(Em2, 'Laplacian'), 'MassMatrix/Laplacian must persist through save/load.');
assert(issparse(Em2.MassMatrix) && isequal(size(Em2.MassMatrix), [nV, nV]), 'MassMatrix must be sparse [nV x nV].');
assert(full(max(max(abs(Em2.MassMatrix - Em.MassMatrix)))) < 1e-12, 'MassMatrix changed on round-trip.');
assert(full(max(max(abs(Em2.Laplacian  - Em.Laplacian))))  < 1e-12, 'Laplacian changed on round-trip.');
assert(strcmp(Em2.LaplacianType, 'Laplace-Beltrami'), 'LaplacianType not preserved.');
assert(strcmp(Em2.MassType, Em.MassType), 'MassType not preserved.');
fprintf('PASSED: MassMatrix/Laplacian (+types) survive round-trip.\n');

% ----- Backward compatibility: a struct WITHOUT the metadata loads with defaults -----
TessTmp = load(SurfFile);
TessTmp.Eigenmodes = rmfield(TessTmp.Eigenmodes, intersect({'Component','CompRank','nComponents','MassMatrix','Laplacian','LaplacianType'}, fieldnames(TessTmp.Eigenmodes)));
bst_save(SurfFile, TessTmp, 'v7');
[EigOld, isOld] = in_tess_eigenmodes(SurfFile);
assert(isOld, 'Old-format eigenmodes must still load.');
nOld = size(EigOld.Vectors, 2);
assert(isequal(EigOld.Component(:), ones(nOld,1)), 'Missing Component must default to a single component.');
assert(isequal(EigOld.CompRank(:),  (1:nOld)'),    'Missing CompRank must default to 1..K.');
fprintf('PASSED: backward-compat defaults for old files without Component/CompRank.\n');

fprintf('ALL TESTS PASSED: test_io_eigenmodes_roundtrip\n');
end
