function test_process_import_bids_options
% Verify the BIDS importer exposes the cortex downsampling method/level options
% (default icosphere/ico5) and the GetIcoVertexCount level->count helper.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sProcess = process_import_bids('GetDescription');

% New options exist, old one preserved
assert(isfield(sProcess.options, 'downsamplemethod'), 'Missing option: downsamplemethod');
assert(isfield(sProcess.options, 'icolevel'),         'Missing option: icolevel');
assert(isfield(sProcess.options, 'nvertices'),        'nvertices option must still exist');

% Correct widget types
assert(strcmp(sProcess.options.downsamplemethod.Type, 'radio_linelabel'), 'downsamplemethod must be radio_linelabel.');
assert(strcmp(sProcess.options.icolevel.Type,         'radio_linelabel'), 'icolevel must be radio_linelabel.');

% Defaults: icosphere / ico5
assert(strcmp(sProcess.options.downsamplemethod.Value, 'icosphere'), 'Default downsamplemethod must be icosphere.');
assert(strcmp(sProcess.options.icolevel.Value,         'ico5'),      'Default icolevel must be ico5.');

% GetIcoVertexCount mappings (total cortex vertices)
assert(process_import_bids('GetIcoVertexCount', 'ico3') == 1284,  'ico3 -> 1284');
assert(process_import_bids('GetIcoVertexCount', 'ico4') == 5124,  'ico4 -> 5124');
assert(process_import_bids('GetIcoVertexCount', 'ico5') == 20484, 'ico5 -> 20484');
assert(process_import_bids('GetIcoVertexCount', 'ico6') == 81924, 'ico6 -> 81924');

% Unknown level errors
threw = false;
try, process_import_bids('GetIcoVertexCount', 'ico9'); catch, threw = true; end
assert(threw, 'Unknown ico level should error.');

fprintf('ALL TESTS PASSED: test_process_import_bids_options\n');
end
