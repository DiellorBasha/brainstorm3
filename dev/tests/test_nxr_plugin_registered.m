function test_nxr_plugin_registered
% Verify nxr-compute is registered in bst_plugin GetSupported with expected fields.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

PlugDesc = bst_plugin('GetSupported', 'nxr-compute');
assert(~isempty(PlugDesc) && strcmp(PlugDesc.Name, 'nxr-compute'), ...
    'nxr-compute not found in GetSupported.');
assert(strcmp(PlugDesc.Category, 'Anatomy'), 'Unexpected Category: %s', PlugDesc.Category);
assert(PlugDesc.AutoLoad == 0, 'nxr-compute must be AutoLoad=0 (install-on-demand).');
assert(~isempty(PlugDesc.URLinfo), 'URLinfo must be set.');
% On Apple Silicon the TestFile must point at the mac arm binary
if strcmp(bst_get('OsType'), 'mac64arm')
    assert(strcmp(PlugDesc.TestFile, 'nxr_compute.mexmaca64'), ...
        'TestFile must be nxr_compute.mexmaca64 on mac64arm, got: %s', PlugDesc.TestFile);
end
fprintf('ALL TESTS PASSED: test_nxr_plugin_registered\n');
end
