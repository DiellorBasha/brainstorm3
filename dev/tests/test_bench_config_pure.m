function test_bench_config_pure
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot, fullfile(repoRoot,'dev','benchmarks'));

C = bench_config();
req = {'anatomies','methods','regimes','snr_db','k_total','nReplicates','seed','outDir','nModes_eig'};
for i=1:numel(req); assert(isfield(C,req{i}), 'config missing field %s', req{i}); end
assert(numel(C.anatomies)==2, 'expect 2 anatomies (Auditory, Neuromag).');
assert(isequal(sort(C.k_total), [600 1200 2000]), 'K-sweep must be {600,1200,2000} total.');
assert(C.nModes_eig==1000, 'eigenmodes computed at 1000 per hemisphere.');
assert(all(ismember({'focal','patch','distributed'}, C.regimes)), 'three regimes required.');

Csm = bench_config('smoke');
assert(numel(Csm.anatomies)==1, 'smoke uses 1 anatomy.');
assert(isequal(Csm.regimes, {'focal'}), 'smoke uses focal only.');
assert(numel(Csm.snr_db)==2 && Csm.nReplicates==2, 'smoke is tiny.');
disp('ALL TESTS PASSED');
end
