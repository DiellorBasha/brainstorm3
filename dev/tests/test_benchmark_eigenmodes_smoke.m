function test_benchmark_eigenmodes_smoke
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot, fullfile(repoRoot,'dev','benchmarks'));
if ~brainstorm('status'); brainstorm nogui; end

out = benchmark_eigenmodes('smoke');
assert(exist(fullfile(out,'synthetic.csv'),'file')==2, 'synthetic.csv missing.');
assert(exist(fullfile(out,'stats.csv'),'file')==2, 'stats.csv missing.');
assert(exist(fullfile(out,'REPORT.md'),'file')==2, 'REPORT.md missing.');
png = dir(fullfile(out,'figures','*.png'));
assert(numel(png) >= 4, 'expect at least 4 figures (>=5 if cortex render ran).');
disp('ALL TESTS PASSED');
end
