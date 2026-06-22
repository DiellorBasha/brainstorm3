function test_eigenmodes_freq_compute
% Pure test of the engine's spectral core: shape, freqs, labels, and the
% kernel-swap == project-first parity (including bad-segment handling).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end
addpath(fullfile(repoRoot, 'toolbox', 'process', 'functions'));

rng(11);
nChan = 18; nModes = 5; sfreq = 200; nT = 4000;
F  = randn(nChan, nT);
MK = randn(nModes, nChan);            % mode kernel [nModes x nChan]
Values = (1:nModes)';

% ---- PSD: shape, frequency vector, row labels ----
[TF, Freqs, RowNames] = process_eigenmodes_freq('Compute', ...
    F, sfreq, MK, [], 'psd', 1, 50, 'mean', 'physical', Values);
assert(ndims(TF) == 3 && size(TF,1) == nModes && size(TF,2) == 1, 'PSD TF shape wrong');
assert(size(TF,3) == numel(Freqs), 'PSD freq dim mismatch');
assert(issorted(Freqs) && all(Freqs >= 0), 'freqs must be ascending and non-negative');
assert(numel(RowNames) == nModes && ischar(RowNames{1}), 'row labels wrong');
assert(all(isfinite(TF(:))) && all(TF(:) >= 0), 'PSD power must be finite and non-negative');

% ---- Parity (no bad segments): kernel-swap == project-first ----
TFref = bst_psd(MK*F, sfreq, 1, 50, [], [], 'mean', 'physical');
assert(max(abs(TF(:) - TFref(:))) < 1e-9, 'PSD kernel-swap != project-first');

% ---- Parity WITH a bad segment that drops a window ----
BadSeg = [round(1.2*sfreq); round(1.8*sfreq)];   % [start; stop] samples
[TFb] = process_eigenmodes_freq('Compute', ...
    F, sfreq, MK, BadSeg, 'psd', 1, 50, 'mean', 'physical', Values);
TFbref = bst_psd(MK*F, sfreq, 1, 50, BadSeg, [], 'mean', 'physical');
assert(max(abs(TFb(:) - TFbref(:))) < 1e-9, 'bad-segment handling differs from project-first');
assert(max(abs(TFb(:) - TF(:))) > 0, 'bad segment should change the result');

% ---- FFT: single-window transform, same parity ----
[TF2, Freqs2] = process_eigenmodes_freq('Compute', ...
    F, sfreq, MK, [], 'fft', [], [], [], 'physical', Values);
assert(size(TF2,1) == nModes && size(TF2,2) == 1 && size(TF2,3) == numel(Freqs2), 'FFT shape wrong');
TF2ref = bst_psd(MK*F, sfreq, [], 0, [], [], [], 'physical');
assert(max(abs(TF2(:) - TF2ref(:))) < 1e-9, 'FFT kernel-swap != project-first');

fprintf('ALL TESTS PASSED: test_eigenmodes_freq_compute\n');
end
