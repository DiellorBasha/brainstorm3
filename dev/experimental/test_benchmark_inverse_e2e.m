function test_benchmark_inverse_e2e
% Smoke: run the comparator panel on simulated data using a real OMEGA base head
% model + eigenmode head model + noise cov; verify each method returns a finite
% [nVert x nTime] vertex estimate. Skips cleanly otherwise.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot, fullfile(repoRoot,'dev','benchmarks'));
if ~brainstorm('status'); brainstorm nogui; end

sStudies = bst_get('ProtocolStudies');
if isempty(sStudies) || ~isfield(sStudies,'Study') || isempty(sStudies.Study)
    disp('SKIP: no protocol loaded.'); return;
end
% Find a study with a base (non-eigenmode) surface head model + noise cov + cortex eigenmodes
T = [];
for iS = 1:numel(sStudies.Study)
    s = sStudies.Study(iS);
    if isempty(s.iHeadModel) || s.iHeadModel<1 || numel(s.HeadModel)<1; continue; end
    if isempty(s.NoiseCov) || isempty(s.NoiseCov(1).FileName); continue; end
    if isempty(s.Channel) || ~any(strcmpi(s.Channel.DisplayableSensorTypes,'MEG')); continue; end
    iBase = [];
    for ih = 1:numel(s.HeadModel)
        try hm = in_bst_headmodel(s.HeadModel(ih).FileName,0); catch; continue; end
        isEig = isfield(hm,'isEigenmode') && hm.isEigenmode;
        if ~isEig && strcmpi(hm.HeadModelType,'surface'); iBase = ih; break; end
    end
    if isempty(iBase); continue; end
    hmB = in_bst_headmodel(s.HeadModel(iBase).FileName,0,'SurfaceFile');
    [~, isE] = in_tess_eigenmodes(hmB.SurfaceFile);
    if isE; T = struct('iStudy',iS,'iBase',iBase); break; end
end
if isempty(T); disp('SKIP: no study with base surface HM + eigenmodes + noise cov.'); return; end

s = sStudies.Study(T.iStudy);
baseHmFile = s.HeadModel(T.iBase).FileName;
ncFile     = s.NoiseCov(1).FileName;
chFile     = bst_get('ChannelFileForStudy', s.FileName);

% Simulate a focal source through the constrained base leadfield
baseHM = in_bst_headmodel(baseHmFile, 1);   % ApplyOrient=1 -> [nCh x nVert]
goodMask = all(isfinite(double(baseHM.Gain)),2);
L = double(baseHM.Gain(goodMask,:));
Surf = in_tess_bst(baseHM.SurfaceFile);
if ~isfield(Surf,'VertConn') || isempty(Surf.VertConn)
    Surf.VertConn = tess_vertconn(Surf.Vertices, Surf.Faces);
end
S = bst_benchmark_sources(struct('Vertices',Surf.Vertices,'VertConn',Surf.VertConn), ...
    'focal', 'nTime', 10, 'Seed', 1);
NC = load(file_fullpath(ncFile)); C = NC.NoiseCov(goodMask,goodMask);
Sim = bst_benchmark_simulate(L, S.Sources, C, 'SNR', 6, 'Seed', 1);

Est = bst_benchmark_inverse(Sim.F, baseHmFile, ncFile, chFile, goodMask, 6);
fn = fieldnames(Est);
assert(~isempty(fn), 'panel must return at least one method estimate.');
for i = 1:numel(fn)
    E = Est.(fn{i});
    assert(size(E,1)==size(L,2), '%s estimate must have nVert rows.', fn{i});
    assert(all(isfinite(E(:))), '%s estimate must be finite.', fn{i});
end
disp('ALL TESTS PASSED');
end
