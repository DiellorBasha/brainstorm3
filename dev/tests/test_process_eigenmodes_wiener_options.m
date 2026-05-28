function test_process_eigenmodes_wiener_options
% Verify the Wiener-filter process is a 2-input data->matrix Sources process.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sProcess = process_eigenmodes_wiener('GetDescription');
assert(strcmp(sProcess.Category, 'Custom'),  'Category must be Custom.');
assert(strcmp(sProcess.SubGroup, 'Sources'), 'SubGroup must be Sources.');
assert(abs(sProcess.Index - 336.95) < 1e-9,  'Index must be 336.95.');
assert(isequal(sProcess.InputTypes,  {'data','raw'}), 'InputTypes must be {data,raw}.');
assert(isequal(sProcess.OutputTypes, {'matrix'}),     'OutputTypes must be {matrix}.');
assert(sProcess.nInputs == 2, 'nInputs must be 2 (data + empty-room).');
for f = {'nmodes','noisewin','alpha','gainfloor','domirror','dorecon','savegain'}
    assert(isfield(sProcess.options, f{1}), 'Missing option: %s', f{1});
end
assert(sProcess.options.alpha.Value{1} == 1,     'Default alpha must be 1.');
assert(sProcess.options.gainfloor.Value{1} == 0, 'Default gainfloor must be 0.');
assert(sProcess.options.domirror.Value == 1,     'Default domirror must be 1 (on).');
assert(sProcess.options.dorecon.Value  == 0,     'Default dorecon must be 0 (off).');
assert(sProcess.options.savegain.Value == 0,     'Default savegain must be 0 (off).');
fprintf('ALL TESTS PASSED: test_process_eigenmodes_wiener_options\n');
end
