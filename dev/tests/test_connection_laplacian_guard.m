function test_connection_laplacian_guard
% tess_connection_laplacian must error cleanly (nxrNotLoaded) when the nxr-compute
% plugin is not loaded, since there is no MATLAB fallback. Uses a synthetic sphere
% (the guard fires before any assembly). Restores the plugin's loaded state after.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% Record whether nxr is installed/loaded so we can restore it afterward.
PlugDesc = bst_plugin('GetInstalled', 'nxr-compute');
wasLoaded = ~isempty(PlugDesc) && isfield(PlugDesc, 'isLoaded') && PlugDesc.isLoaded;
if isempty(PlugDesc)
    fprintf('SKIP: nxr-compute not installed; cannot exercise the loaded/unloaded guard.\n');
    return;
end
% Ensure it is reloaded at the end regardless of assertion outcome.
restore = onCleanup(@() bst_plugin('Load', 'nxr-compute'));

% Unload so the guard precondition (not loaded) holds.
if wasLoaded
    bst_plugin('Unload', 'nxr-compute');
end

[V, F] = tess_sphere(162);
threw = false;
try
    tess_connection_laplacian(V, F);
catch ME
    threw = strcmp(ME.identifier, 'tess_connection_laplacian:nxrNotLoaded');
    if ~threw
        rethrow(ME);
    end
end
assert(threw, 'tess_connection_laplacian did not raise nxrNotLoaded when the plugin was unloaded.');

fprintf('PASSED: errors with nxrNotLoaded when nxr-compute is not loaded.\n');
fprintf('ALL TESTS PASSED: test_connection_laplacian_guard\n');
end
