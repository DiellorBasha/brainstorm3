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

% ---- FFT process ----
sFft = process_eigenmodes_fft('GetDescription');
assert(~isempty(sFft.Comment), 'FFT comment empty');
assert(strcmpi(sFft.InputTypes{1}, 'results'), 'FFT input must be results');
assert(strcmpi(sFft.OutputTypes{1}, 'timefreq'), 'FFT output must be timefreq');
assert(isfield(sFft.options, 'nmodes') && isequal(sFft.options.nmodes.Value{1}, 300), 'FFT nmodes default 300');
assert(isfield(sFft.options, 'units'), 'FFT must expose units');
assert(~isfield(sFft.options, 'win_length'), 'FFT must not expose window length');
assert(~isfield(sFft.options, 'win_overlap'), 'FFT must not expose window overlap');
assert(~isfield(sFft.options, 'win_std'), 'FFT must not expose window std');
assert(ischar(process_eigenmodes_fft('FormatComment', sFft)), 'FFT FormatComment must return char');

fprintf('ALL TESTS PASSED: test_eigenmodes_freq_processes\n');
end
