function test_eigenmodes_freq_processes
% Validate GetDescription/FormatComment of the mode-frequency processes (no DB).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end
addpath(fullfile(repoRoot, 'toolbox', 'process', 'functions'));

% ---- PSD process ----
sPsd = process_eigenmodes_psd('GetDescription');
assert(~isempty(sPsd.Comment), 'PSD comment empty');
assert(strcmpi(sPsd.InputTypes{1}, 'results'), 'PSD input must be results');
assert(strcmpi(sPsd.OutputTypes{1}, 'timefreq'), 'PSD output must be timefreq');
assert(isfield(sPsd.options, 'nmodes') && isequal(sPsd.options.nmodes.Value{1}, 300), 'PSD nmodes default 300');
assert(isfield(sPsd.options, 'win_length'), 'PSD must expose window length');
assert(isfield(sPsd.options, 'win_overlap'), 'PSD must expose window overlap');
assert(isfield(sPsd.options, 'units'), 'PSD must expose units');
assert(ischar(process_eigenmodes_psd('FormatComment', sPsd)), 'PSD FormatComment must return char');

fprintf('ALL TESTS PASSED: test_eigenmodes_freq_processes\n');
end
