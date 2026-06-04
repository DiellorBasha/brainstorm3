function test_process_parser_engine_skip
% ParseProcessFolder must treat a process_*.m that defines NO GetDescription
% subfunction as a shared engine / support function (driven by sibling menu
% processes, e.g. process_eigenmodes_freq) and skip it SILENTLY -- while a file
% whose GetDescription EXISTS but errors is a broken menu process and must still
% be reported as "Invalid ... function".
%
% The discrimination lives in panel_process_select('IsProcessEngineFile', file).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

%% ===== Discriminator: real files =====
% Shared engine: no GetDescription -> engine (skip).
fEngine = which('process_eigenmodes_freq');
assert(~isempty(fEngine), 'process_eigenmodes_freq must be on the path.');
assert(isequal(panel_process_select('IsProcessEngineFile', fEngine), true), ...
    'process_eigenmodes_freq (no GetDescription) must be classified as a shared engine.');

% Real menu process: has GetDescription -> NOT an engine (would be registered/warned).
fMenu = which('process_eigenmodes_psd');
assert(~isempty(fMenu), 'process_eigenmodes_psd must be on the path.');
assert(isequal(panel_process_select('IsProcessEngineFile', fMenu), false), ...
    'process_eigenmodes_psd (has GetDescription) must NOT be classified as an engine.');

%% ===== Discriminator: temp files (both branches) =====
% (a) Engine stub: only a Compute subfunction, no GetDescription -> engine.
fA = fullfile(tempdir, 'process_zz_engine_stub.m');
writeLines(fA, {
    'function varargout = process_zz_engine_stub(varargin)'
    'eval(macro_method);'
    'end'
    ''
    'function y = Compute(x) %#ok<DEFNU>'
    '    y = x;'
    'end'});
% (b) Broken menu process: GetDescription EXISTS but errors -> NOT an engine.
fB = fullfile(tempdir, 'process_zz_broken_desc.m');
writeLines(fB, {
    'function varargout = process_zz_broken_desc(varargin)'
    'eval(macro_method);'
    'end'
    ''
    'function sProcess = GetDescription() %#ok<DEFNU>'
    '    error(''boom in GetDescription'');'
    'end'});
cleanup = onCleanup(@() cellfun(@delete, {fA, fB}));

assert(isequal(panel_process_select('IsProcessEngineFile', fA), true), ...
    'A process_*.m with only a Compute subfunction must be classified as an engine.');
assert(isequal(panel_process_select('IsProcessEngineFile', fB), false), ...
    'A process_*.m whose GetDescription exists (even if it errors) must NOT be an engine.');
fprintf('PASSED: engine/menu-process discrimination (real + temp files).\n');

%% ===== End-to-end: a forced re-parse no longer flags the engine =====
out = evalc("panel_process_select('ParseProcessFolder', 1)");
assert(isempty(strfind(out, 'process_eigenmodes_freq')), ...
    'ParseProcessFolder must skip process_eigenmodes_freq silently (no "Invalid" message).');
fprintf('PASSED: forced re-parse does not flag process_eigenmodes_freq.\n');

fprintf('ALL TESTS PASSED: test_process_parser_engine_skip\n');
end


function writeLines(filePath, lines)
    fid = fopen(filePath, 'w');
    assert(fid > 0, 'Could not open %s for writing.', filePath);
    for i = 1:numel(lines)
        fprintf(fid, '%s\n', lines{i});
    end
    fclose(fid);
end
