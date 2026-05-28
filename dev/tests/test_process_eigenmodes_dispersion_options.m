function test_process_eigenmodes_dispersion_options
% Verify the dispersion process is a timefreq-in / matrix-out Sources process.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sProcess = process_eigenmodes_dispersion('GetDescription');
assert(strcmp(sProcess.SubGroup, 'Sources'), 'SubGroup must be Sources.');
assert(abs(sProcess.Index - 336.8) < 1e-9, 'Index must be 336.8.');
assert(isequal(sProcess.InputTypes, {'timefreq'}), 'InputTypes must be {timefreq}.');
assert(isequal(sProcess.OutputTypes, {'matrix'}), 'OutputTypes must be {matrix}.');
assert(isfield(sProcess.options, 'minpowerfrac'), 'Missing minpowerfrac option.');
assert(isfield(sProcess.options, 'label_info'), 'Missing label_info option.');
fprintf('ALL TESTS PASSED: test_process_eigenmodes_dispersion_options\n');
end
