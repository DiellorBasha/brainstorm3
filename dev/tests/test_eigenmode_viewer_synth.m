function test_eigenmode_viewer_synth
% Pure synthesis: view_eigenmodes('SynthColumn', PairedGrid, W) == PairedGrid*W(:).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(fullfile(repoRoot, 'toolbox', 'gui'));
nV = 50; K = 8;
PairedGrid = reshape(1:(nV*K), nV, K) / 100;
W = zeros(1,K); W(3) = 1;
col = view_eigenmodes('SynthColumn', PairedGrid, W);
assert(isequal(size(col), [nV 1]), 'column is nV x 1');
assert(max(abs(col - PairedGrid(:,3))) < 1e-12, 'single -> the rank-k paired column');
W = zeros(1,K); W(2:4) = 1;
col = view_eigenmodes('SynthColumn', PairedGrid, W);
assert(max(abs(col - sum(PairedGrid(:,2:4),2))) < 1e-12, 'band -> sum of selected columns');
fprintf('ALL TESTS PASSED: test_eigenmode_viewer_synth\n');
end
