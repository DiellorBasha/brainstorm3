function test_sa_crossing_profile
% sa_crossing_profile builds crossing paths between facing-wall pairs; along the
% path the normal turns substantially (~pi, descends fundus + back up) while the
% Fiedler frame turns less (smoother). SKIP if no cortex.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'dev', 'benchmarks', 'sign_ambiguity'));
if ~brainstorm('status'), brainstorm nogui; end
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute required: %s', errMsg);
bst_plugin('Load', 'nxr-compute');
SurfaceFile = find_cortex_20484V();
if isempty(SurfaceFile)
    fprintf('SKIP: no 20484-vertex cortex.\n'); return;
end

F = sa_frames(SurfaceFile);
TessMat = in_tess_bst(SurfaceFile);
if isfield(TessMat,'SulciMap') && ~isempty(TessMat.SulciMap)
    SulciMap = TessMat.SulciMap;
else
    SulciMap = tess_sulcimap(TessMat);
end
opts = struct('MaxDist',0.003, 'NormalDot',-0.7, 'Nring',3);
pairs = sa_sulcal_walls(F.Vtx, F.Nv, SulciMap, F.VertConn, opts);

P = sa_crossing_profile(F, pairs, 3);
assert(~isempty(P), 'Expected at least one crossing profile.');

% Number of profiles is bounded by min(nProfiles, size(pairs,1)).
assert(numel(P) <= 3, 'Profile count must not exceed nProfiles=3.');
assert(numel(P) <= size(pairs,1), 'Profile count must not exceed pair count.');
% Each profile s field is empty when no J was passed.
for k = 1:numel(P)
    assert(isempty(P(k).s), 'P(k).s should be [] when J is not passed.');
end

for k = 1:numel(P)
    assert(numel(P(k).path) >= 3, 'Path too short.');
    assert(isequal(size(P(k).arclen), size(P(k).nAngle)), 'arclen/nAngle size mismatch.');
    assert(P(k).nAngle(end) >= P(k).fAngle(end) - 1e-6, ...
        'Normal should turn at least as much as the Fiedler frame across a fold.');
end
fprintf('PASSED: test_sa_crossing_profile (%d profiles, max normal turn %.2f rad).\n', ...
    numel(P), max(arrayfun(@(x) x.nAngle(end), P)));
end

function SurfaceFile = find_cortex_20484V()
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
if isempty(sSubjects), return; end
allSubj = [sSubjects.Subject];
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex'), continue; end
        try, T = load(file_fullpath(surf(iF).FileName), 'Vertices'); catch, continue; end
        if size(T.Vertices,1) == 20484, SurfaceFile = surf(iF).FileName; return; end
    end
end
end
