function test_process_eigenmodes_options
% Verify the Compute-eigenmodes process exposes an opt-in, unchecked repair box.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sProcess = process_eigenmodes('GetDescription');
assert(isfield(sProcess.options, 'repair'), 'Missing "repair" option.');
assert(strcmp(sProcess.options.repair.Type, 'checkbox'), 'repair option must be a checkbox.');
assert(~isfield(sProcess.options, 'fixmesh'), 'Old "fixmesh" option still present.');
assert(sProcess.options.repair.Value == 0, 'Repair must default to 0 (no auto-repair).');
assert(~isempty(regexpi(sProcess.options.repair.Comment, 'risky')), ...
    'Repair comment must warn that repair is risky.');
fprintf('ALL TESTS PASSED: test_process_eigenmodes_options\n');
end
