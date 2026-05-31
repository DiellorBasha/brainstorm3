function test_laplacian_backend_select
% Verify backend selection: 'auto' tracks load state; 'nxr' errors when not
% loaded; 'matlab' always works; auto==nxr when loaded, auto==matlab when not.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

codeDir  = fileparts(repoRoot);
nxrStage = fullfile(codeDir, 'nxr-compute', 'dist', 'plugin', 'nxr-compute');
assert(isfolder(nxrStage), 'Run nxr-compute/scripts/package-plugin.sh first (missing %s).', nxrStage);
destDir = fullfile(bst_get('UserPluginsDir'), 'nxr-compute');
if isfolder(destDir), file_delete(destDir, 1, 3); end
mkdir(destDir);
copyfile(fullfile(nxrStage, '*'), destDir);

[V, F] = tess_sphere(162);

% --- With nxr UNLOADED ---
bst_plugin('Unload', 'nxr-compute');
La = tess_laplacian(V, F, 'Backend', 'auto');     % must use MATLAB
Lm = tess_laplacian(V, F, 'Backend', 'matlab');
assert(isequal(La, Lm), 'auto (unloaded) did not equal MATLAB backend.');
threw = false;
try
    tess_laplacian(V, F, 'Backend', 'nxr');        % must error when not loaded
catch ME
    threw = strcmp(ME.identifier, 'tess_laplacian:nxrNotLoaded');
end
assert(threw, 'Backend ''nxr'' did not raise nxrNotLoaded when unloaded.');
fprintf('PASSED: unloaded -> auto uses MATLAB; explicit nxr errors.\n');

% --- With nxr LOADED ---
[isOk, errMsg] = bst_plugin('Load', 'nxr-compute');
assert(isOk, 'Load failed: %s', errMsg);
cleanup = onCleanup(@() bst_plugin('Unload', 'nxr-compute'));
La2 = tess_laplacian(V, F, 'Backend', 'auto');
Ln  = tess_laplacian(V, F, 'Backend', 'nxr');
assert(isequal(La2, Ln), 'auto (loaded) did not equal nxr backend.');
fprintf('PASSED: loaded -> auto uses nxr.\n');

fprintf('ALL TESTS PASSED: test_laplacian_backend_select\n');
end
