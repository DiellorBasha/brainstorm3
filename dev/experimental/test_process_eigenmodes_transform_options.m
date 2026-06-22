function test_process_eigenmodes_transform_options
% Verify the transform process exposes only transform options (no regularization
% knobs), and defaults vertex reconstruction off.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sProcess = process_eigenmodes_transform('GetDescription');
assert(strcmp(sProcess.SubGroup, 'Sources'), 'SubGroup must be Sources.');
assert(isequal(sProcess.InputTypes, {'data','raw'}), 'InputTypes must be {data,raw}.');

assert(isfield(sProcess.options, 'nmodes'),  'Missing nmodes option.');
assert(isfield(sProcess.options, 'dorecon'), 'Missing dorecon option.');
assert(strcmp(sProcess.options.dorecon.Type, 'checkbox'), 'dorecon must be a checkbox.');
assert(sProcess.options.dorecon.Value == 0, 'Vertex reconstruction must default OFF.');
assert(sProcess.options.nmodes.Value{1} == 0, 'nmodes must default to 0 (auto).');

% Transform is unregularized: it must NOT expose inverse-method knobs.
assert(~isfield(sProcess.options, 'method'),     'Transform must not expose a method option.');
assert(~isfield(sProcess.options, 'prioralpha'), 'Transform must not expose a prior option.');
assert(~isfield(sProcess.options, 'snr'),        'Transform must not expose an SNR option.');

fprintf('ALL TESTS PASSED: test_process_eigenmodes_transform_options\n');
end
