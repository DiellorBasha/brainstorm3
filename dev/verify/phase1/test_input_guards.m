% TEST_INPUT_GUARDS: spec-validation guards added in the Phase 1 fix wave.
%   - NaN vertex -> invalidInput (tess_laplacian, tess_massmatrix)
%   - degenerate (zero-area, duplicated-point) triangle -> degenerateFace
%     (tess_laplacian, tess_massmatrix)
%   - tess_operators recipe with nonempty Tau on Laplace-Beltrami -> unexpectedParameter
%   - tess_operators recipe struct missing .Name -> unknownVariant

% --- NaN vertex guard: tess_laplacian ---
V = [0 0 0; 1 0 0; 0 1 0]; F = [1 2 3];
Vnan = V; Vnan(2,1) = NaN;
ok = false;
try, tess_laplacian(Vnan, F); catch err, ok = strcmp(err.identifier, 'tess_laplacian:invalidInput'); end
assert(ok, 'tess_laplacian must raise invalidInput on NaN vertex');

% --- NaN vertex guard: tess_massmatrix ---
ok = false;
try, tess_massmatrix(Vnan, F); catch err, ok = strcmp(err.identifier, 'tess_massmatrix:invalidInput'); end
assert(ok, 'tess_massmatrix must raise invalidInput on NaN vertex');

% --- Inf vertex guard: tess_laplacian (also invalidInput) ---
Vinf = V; Vinf(3,2) = Inf;
ok = false;
try, tess_laplacian(Vinf, F); catch err, ok = strcmp(err.identifier, 'tess_laplacian:invalidInput'); end
assert(ok, 'tess_laplacian must raise invalidInput on Inf vertex');

% --- degenerate (zero-area) face guard: tess_laplacian ---
% duplicated point v2==v3 -> zero-area triangle
Vdeg = [0 0 0; 1 0 0; 1 0 0];
ok = false;
try, tess_laplacian(Vdeg, F); catch err, ok = strcmp(err.identifier, 'tess_laplacian:degenerateFace'); end
assert(ok, 'tess_laplacian must raise degenerateFace on zero-area triangle');

% --- degenerate (zero-area) face guard: tess_massmatrix ---
ok = false;
try, tess_massmatrix(Vdeg, F); catch err, ok = strcmp(err.identifier, 'tess_massmatrix:degenerateFace'); end
assert(ok, 'tess_massmatrix must raise degenerateFace on zero-area triangle');

% --- tess_operators: LBO recipe with nonempty Tau -> unexpectedParameter ---
SurfaceFile = '/Users/diellorbasha/workspace/research/code/brainstorm3/dev/verify/phase0/bst_userdir_clean/.brainstorm/local_db/omega-tutorial-cortical-flow/anat/sub-0002/tess_cortex_pial_low.mat';
recipeTau = struct('Name', 'Laplace-Beltrami', 'Tau', 0.5);
ok = false;
try, tess_operators(SurfaceFile, recipeTau); catch err, ok = strcmp(err.identifier, 'tess_operators:unexpectedParameter'); end
assert(ok, 'tess_operators must raise unexpectedParameter on nonempty Tau for Laplace-Beltrami');

% --- tess_operators: recipe struct without .Name -> unknownVariant ---
recipeNoName = struct('Tau', []);
ok = false;
try, tess_operators(SurfaceFile, recipeNoName); catch err, ok = strcmp(err.identifier, 'tess_operators:unknownVariant'); end
assert(ok, 'tess_operators must raise unknownVariant when recipe struct has no .Name field');

% --- sanity: the guards must NOT fire on the real cortex happy path ---
[Op, Ms, vH] = tess_operators(SurfaceFile, 'Laplace-Beltrami'); %#ok<ASGLU>
assert(isequal(vH{1}(:), sort(vH{1}(:))), 'vH{1} must be sorted ascending');
assert(isequal(vH{2}(:), sort(vH{2}(:))), 'vH{2} must be sorted ascending');

disp('test_input_guards PASSED');
