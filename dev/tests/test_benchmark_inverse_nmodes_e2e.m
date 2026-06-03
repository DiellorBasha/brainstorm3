function test_benchmark_inverse_nmodes_e2e
% Verify the optional nModes arg caps the eigenmode reconstruction and stays finite.
% Uses bench_fixtures to obtain a vertex-consistent (base HM, noise cov, eigenmode HM)
% set, then checks that K=600 and K=2000 give DIFFERENT, finite eig reconstructions.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot, fullfile(repoRoot,'dev','benchmarks'));
if ~brainstorm('status'); brainstorm nogui; end

if isempty(bst_get('Protocol','TutorialAuditory'))
    disp('SKIP: TutorialAuditory protocol not present.'); return;
end
C    = bench_config();
info = bench_fixtures(C.anatomies{1}, C.nModes_eig);

baseHM = in_bst_headmodel(info.baseHmFile, 1);   % [nCh x nVert]
goodMask = all(isfinite(double(baseHM.Gain)),2);
L = double(baseHM.Gain(goodMask,:));
Surf = in_tess_bst(info.surfaceFile);
if ~isfield(Surf,'VertConn') || isempty(Surf.VertConn)
    Surf.VertConn = tess_vertconn(Surf.Vertices, Surf.Faces);
end
S   = bst_benchmark_sources(struct('Vertices',Surf.Vertices,'VertConn',Surf.VertConn),'focal','nTime',5,'Seed',1);
NC  = load(file_fullpath(info.ncFile)); Cnoise = NC.NoiseCov(goodMask,goodMask);
Sim = bst_benchmark_simulate(L, S.Sources, Cnoise, 'SNR', 6, 'Seed', 1);

Est600  = bst_benchmark_inverse(Sim.F, info.baseHmFile, info.ncFile, info.chFile, goodMask, 6, 600);
Est2000 = bst_benchmark_inverse(Sim.F, info.baseHmFile, info.ncFile, info.chFile, goodMask, 6, 2000);

assert(isfield(Est600,'eig_mne_log') && all(isfinite(Est600.eig_mne_log(:))), 'eig estimate must be finite at K=600.');
assert(isfield(Est2000,'eig_mne_log') && all(isfinite(Est2000.eig_mne_log(:))), 'eig estimate must be finite at K=2000.');
assert(norm(Est600.eig_mne_log(:) - Est2000.eig_mne_log(:)) > 0, 'different K must give different eig reconstruction.');
assert(all(isfinite(Est600.wmne(:))), 'standard methods must remain finite.');
disp('ALL TESTS PASSED');
end
