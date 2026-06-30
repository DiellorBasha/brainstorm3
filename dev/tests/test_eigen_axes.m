% test_eigen_axes - bst_eigen('Axes') assembles eigenbasis x time x temporal-frequency
% Run in a live Brainstorm session (protocol with a cortex surface).
SURF = 'sub-MTL0002/tess_cortex_pial_low.mat';
ax = bst_eigen('Axes', struct('SurfaceFile',SURF, 'Variant','Laplace-Beltrami', 'nModes',20, ...
                              'TimeWindow',[0 0.99], 'SampleRate',100));
% --- temporal/spectral axis ---
assert(ax.nT == 100,                       'nT should be 100 (1 s @ 100 Hz)');
assert(abs(ax.Fs - 100) < 1e-9,            'Fs should be 100 Hz');
assert(abs(ax.tlag(2) - 0.01) < 1e-12,     'tlag step = 1/Fs');
assert(ax.NFFT == ax.nT,                   'NFFT == nT');
assert(abs(ax.omega(2) - (ax.Fs/ax.NFFT)) < 1e-9, 'omega step = Fs/NFFT');
assert(numel(ax.omega) == ax.NFFT && numel(ax.tlag) == ax.nT, 'axis lengths');
% --- spatial axis present ---
assert(iscell(ax.Phi) && ~isempty(ax.Phi{1}) && iscell(ax.Mass) && iscell(ax.GlobalVertices), 'eigenbasis fields');
assert(numel(ax.Lambda{1}) == size(ax.Phi{1},2), 'Lambda matches Phi columns');
disp('OK');
