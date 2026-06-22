function test_process_eigenmodes_denoise_options
% Verify the denoise process is a 2-input Sources process with the right knobs.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sProcess = process_eigenmodes_denoise('GetDescription');
assert(strcmp(sProcess.SubGroup, 'Sources'), 'SubGroup must be Sources.');
assert(sProcess.nInputs == 2, 'Must be a 2-input process (Files A = data, Files B = empty-room).');
assert(abs(sProcess.Index - 336.6) < 1e-9, 'Index must be 336.6.');
for f = {'nmodes','noisewin','alpha','snrthresh','floorfrac'}
    assert(isfield(sProcess.options, f{1}), 'Missing option: %s', f{1});
end
fprintf('ALL TESTS PASSED: test_process_eigenmodes_denoise_options\n');
end
