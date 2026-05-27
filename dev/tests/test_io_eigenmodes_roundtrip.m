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

fprintf('ALL TESTS PASSED: test_io_eigenmodes_roundtrip\n');
end
