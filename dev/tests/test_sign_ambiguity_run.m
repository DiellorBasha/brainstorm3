function test_sign_ambiguity_run
% Driver smoke test: runs end to end on the real DB, writes stats + figures.
% SKIP if no 20484V cortex (Component 1 needs at least the surface).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'dev', 'benchmarks', 'sign_ambiguity'));
if ~brainstorm('status'), brainstorm nogui; end
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute required: %s', errMsg);
bst_plugin('Load', 'nxr-compute');

out = sign_ambiguity_run();
if isempty(out)
    fprintf('SKIP: driver reported no usable cortex.\n'); return;
end
assert(exist(out.statsFile, 'file') == 2, 'stats.md not written.');
assert(~isempty(out.figures), 'No figures written.');
assert(isfield(out, 'C1') && ~isempty(out.C1), 'Component 1 stats missing.');
fprintf('PASSED: test_sign_ambiguity_run (out=%s).\n', out.outDir);
end
