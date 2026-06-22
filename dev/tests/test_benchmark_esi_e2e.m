function test_benchmark_esi_e2e
% Smoke: tiny sweep (1 regime x 1 SNR x 2 locations x 1 noise draw) on OMEGA
% and confirm the report struct + CSV are produced. Skips cleanly without a suitable protocol.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot, fullfile(repoRoot,'dev','benchmarks'));
if ~brainstorm('status'); brainstorm nogui; end

outDir = fullfile(repoRoot, 'dev', 'benchmarks', 'smoke');
R = tutorial_benchmark_esi('Regimes', {'focal'}, 'SNRs', 6, 'nLoc', 2, 'nNoise', 1, ...
    'OutDir', outDir);
if isempty(R)
    disp('SKIP: no suitable OMEGA study for the benchmark.'); return;
end
assert(isfield(R,'summary') && ~isempty(R.summary), 'report must have a non-empty summary.');
assert(exist(fullfile(outDir,'summary.csv'),'file')==2, 'summary.csv must be written.');
% Each of the 2 realizations scores >=3 methods x 5 metrics; summary has one row per (method,metric)
methodsSeen = unique({R.summary.method});
assert(numel(methodsSeen) >= 3, 'expected at least the 3 native methods in the summary.');
assert(numel(R.summary) == numel(methodsSeen)*5, 'summary must have nMethods x 5 metrics rows.');
disp('ALL TESTS PASSED');
end
