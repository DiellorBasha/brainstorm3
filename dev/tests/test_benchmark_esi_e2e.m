function test_benchmark_esi_e2e
% Smoke: run a tiny sweep (1 regime x 1 SNR x 2 reps) on OMEGA and confirm the
% report struct + CSV are produced. Skips cleanly without a suitable protocol.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

outDir = fullfile(repoRoot, 'dev', 'benchmarks', 'smoke');
R = tutorial_benchmark_esi('Regimes', {'focal'}, 'SNRs', 6, 'nLoc', 2, 'nNoise', 1, ...
    'OutDir', outDir);
if isempty(R)
    disp('SKIP: no suitable OMEGA study for the benchmark.'); return;
end
assert(isfield(R,'summary') && ~isempty(R.summary), 'report must have a non-empty summary.');
assert(exist(fullfile(outDir,'summary.csv'),'file')==2, 'summary.csv must be written.');
disp('ALL TESTS PASSED');
end
