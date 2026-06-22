function test_process_eigenmodes_coeffsfilter_options
% Verify the coefficient-filter process is a matrix-in / matrix-out Sources process.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sProcess = process_eigenmodes_coeffsfilter('GetDescription');
assert(strcmp(sProcess.SubGroup, 'Sources'), 'SubGroup must be Sources.');
assert(abs(sProcess.Index - 336.9) < 1e-9, 'Index must be 336.9.');
assert(isequal(sProcess.InputTypes, {'matrix'}), 'InputTypes must be {matrix}.');
assert(isequal(sProcess.OutputTypes, {'matrix'}), 'OutputTypes must be {matrix}.');
for f = {'filtertype','cutoffmode','moderange_low','moderange_high','diffusiontime','regbeta','dorecon'}
    assert(isfield(sProcess.options, f{1}), 'Missing option: %s', f{1});
end
assert(strcmp(sProcess.options.filtertype.Value, 'heat'), 'Default filtertype must be heat.');
assert(sProcess.options.dorecon.Value == 0, 'Default dorecon must be 0 (off).');
fprintf('ALL TESTS PASSED: test_process_eigenmodes_coeffsfilter_options\n');
end
