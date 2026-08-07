% TEST_TESS_OPERATORS: pencil assembly on the real cortex + guards.
S = load(fullfile(fileparts(mfilename('fullpath')), 'oracle_lbo_sub0002.mat'));
SurfaceFile = '/Users/diellorbasha/workspace/research/code/brainstorm3/dev/verify/phase0/bst_userdir_clean/.brainstorm/local_db/omega-tutorial-cortical-flow/anat/sub-0002/tess_cortex_pial_low.mat';
% --- happy path: parity of the assembled pencil, both hemispheres ---
[Op, Ms, vH] = tess_operators(SurfaceFile, 'Laplace-Beltrami');
for hh = 1:2
    assert(isequal(vH{hh}(:), S.vH{hh}(:)), 'hemisphere vertex order mismatch');
    assert(norm(Op{hh} - S.A{hh}, 'fro') / norm(S.A{hh}, 'fro') <= 1e-12, 'stiffness parity');
    assert(norm(Ms{hh} - S.B{hh}, 'fro') / norm(S.B{hh}, 'fro') <= 1e-12, 'mass parity');
end
% --- recipe-struct call (reproducibility invariant) ---
recipe = struct('Name', 'Laplace-Beltrami', 'Tau', []);
[Op2, Ms2] = tess_operators(SurfaceFile, recipe);
assert(isequal(Op2{1}, Op{1}) && isequal(Ms2{2}, Ms{2}), 'recipe call must reproduce pencil');
% --- guard: unknown variant ---
ok = false;
try, tess_operators(SurfaceFile, 'Dirac'); catch err, ok = strcmp(err.identifier, 'tess_operators:unknownVariant'); end
assert(ok, 'Dirac must error with unknownVariant until the variant decision');
% --- guard: missing Structures atlas ---
T = load(SurfaceFile); T2 = rmfield(T, 'Atlas');
ok = false;
try, tess_operators(T2, 'Laplace-Beltrami'); catch err, ok = strcmp(err.identifier, 'tess_operators:noHemisphereLabels'); end
assert(ok, 'missing atlas must raise noHemisphereLabels');
% --- guard: non-manifold hemisphere (duplicate a face) ---
T3 = load(SurfaceFile); T3.Faces = [T3.Faces; T3.Faces(1, [2 1 3])];
ok = false;
try, tess_operators(T3, 'Laplace-Beltrami'); catch err, ok = strcmp(err.identifier, 'tess_operators:nonManifold'); end
assert(ok, 'non-manifold input must raise nonManifold');
disp('test_tess_operators PASSED');
