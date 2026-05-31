function test_nxr_plugin_lifecycle
% Stage the locally-built nxr-compute plugin, load it, run a compute, unload.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));   % .../brainstorm3
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% Locate the staged plugin built by NXR/scripts/package-plugin.sh
codeDir  = fileparts(repoRoot);             % .../research/code
nxrStage = fullfile(codeDir, 'nxr-compute', 'dist', 'plugin', 'nxr-compute');
assert(isfolder(nxrStage), ...
    'Staged plugin not found at %s. Run nxr-compute/scripts/package-plugin.sh first.', nxrStage);

% Stage into the Brainstorm user plugins dir (offline "install")
destDir = fullfile(bst_get('UserPluginsDir'), 'nxr-compute');
if isfolder(destDir)
    file_delete(destDir, 1, 3);
end
mkdir(destDir);
copyfile(fullfile(nxrStage, '*'), destDir);

% Load and verify it is reported as loaded
[isOk, errMsg] = bst_plugin('Load', 'nxr-compute');
assert(isOk, 'bst_plugin Load failed: %s', errMsg);
PlugDesc = bst_plugin('GetInstalled', 'nxr-compute');
assert(~isempty(PlugDesc) && PlugDesc.isLoaded, 'Plugin not reported as loaded.');
fprintf('PASSED: nxr-compute staged and loaded.\n');

% Prove the binary actually computes (context builds K and M)
[V, F] = tess_sphere(162);
ctx = nxr.manifold.context(V, F);
assert(isfield(ctx, 'K') && size(ctx.K, 1) == size(V, 1), 'nxr stiffness K wrong size.');
assert(isfield(ctx, 'M') && size(ctx.M, 1) == size(V, 1), 'nxr mass M wrong size.');
fprintf('PASSED: nxr.manifold.context returns K/M (%dx%d).\n', size(ctx.K,1), size(ctx.K,2));

% Unload restores state
[isOkU, errMsgU] = bst_plugin('Unload', 'nxr-compute');
assert(isOkU, 'bst_plugin Unload failed: %s', errMsgU);
PlugDesc2 = bst_plugin('GetInstalled', 'nxr-compute');
assert(isempty(PlugDesc2) || ~PlugDesc2.isLoaded, 'Plugin still loaded after Unload.');
fprintf('ALL TESTS PASSED: test_nxr_plugin_lifecycle\n');
end
