function test_benchmark_metrics_pure
% Verify the metric suite on synthetic GT vs estimate.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot, fullfile(repoRoot,'dev','benchmarks'));
if ~brainstorm('status'); brainstorm nogui; end

% 4 vertices on a line, 10 mm apart
GridLoc = [0 0 0; 10 0 0; 20 0 0; 30 0 0] / 1000;   % metres
gt = [0; 1; 0; 0];                                   % focal GT at vertex 2

% Perfect estimate
M = bst_benchmark_metrics(gt, gt, GridLoc, 2);
assert(M.LocError < 1e-9, 'perfect estimate -> zero localization error.');
assert(abs(M.Correlation - 1) < 1e-12, 'perfect estimate -> correlation 1.');
assert(abs(M.NRMSE) < 1e-9, 'perfect estimate -> NRMSE 0.');
assert(abs(M.AUC - 1) < 1e-12, 'perfect estimate -> AUC 1.');

% Estimate peaks at vertex 4 (30 mm from true vertex 2 at 10 mm) -> LE = 20 mm
est = [0; 0.1; 0; 1];
M2 = bst_benchmark_metrics(gt, est, GridLoc, 2);
assert(abs(M2.LocError - 20) < 1e-9, 'peak at vertex 4 -> 20 mm localization error.');
assert(M2.SpatialDispersion > 0, 'dispersed estimate -> positive dispersion.');

disp('ALL TESTS PASSED');
end
