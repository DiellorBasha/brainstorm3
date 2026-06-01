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
% Cross-platform release (v0.1.0): URLzip + per-OS TestFile must be wired for
% each supported platform.
osType = bst_get('OsType');
expTestFile = struct('mac64arm', 'nxr_compute.mexmaca64', ...
                     'win64',    'nxr_compute.mexw64', ...
                     'linux64',  'nxr_compute.mexa64');
if isfield(expTestFile, osType)
    assert(~isempty(PlugDesc.URLzip), 'URLzip must be set for %s.', osType);
    assert(~isempty(strfind(PlugDesc.URLzip, 'v0.1.0')), ...
        'URLzip should reference the v0.1.0 release, got: %s', PlugDesc.URLzip);
    assert(strcmp(PlugDesc.TestFile, expTestFile.(osType)), ...
        'TestFile must be %s on %s, got: %s', expTestFile.(osType), osType, PlugDesc.TestFile);
end
fprintf('ALL TESTS PASSED: test_nxr_plugin_registered\n');
end
