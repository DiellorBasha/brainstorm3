function test_process_eigenmodes_wavelet_options
% Verify the wavelet process is a matrix-in / timefreq-out Sources process.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sProcess = process_eigenmodes_wavelet('GetDescription');
assert(strcmp(sProcess.SubGroup, 'Sources'), 'SubGroup must be Sources.');
assert(abs(sProcess.Index - 336.7) < 1e-9, 'Index must be 336.7.');
assert(isequal(sProcess.InputTypes, {'matrix'}), 'InputTypes must be {matrix}.');
assert(isequal(sProcess.OutputTypes, {'timefreq'}), 'OutputTypes must be {timefreq}.');
for f = {'flo','fhi','nfreqs','morletfc','morletfwhmtc'}
    assert(isfield(sProcess.options, f{1}), 'Missing option: %s', f{1});
end
assert(isfield(sProcess.options, 'label_info'), 'Missing label_info option.');
fprintf('ALL TESTS PASSED: test_process_eigenmodes_wavelet_options\n');
end
