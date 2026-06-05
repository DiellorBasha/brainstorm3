function test_sa_frames
% Smoke test: sa_frames returns three per-vertex tangent frames + the Fiedler
% decode on the real DB cortex. SKIP if no 20484V cortex.
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
    fprintf('SKIP: no 20484-vertex cortex in the current protocol.\n');
    return;
end

F = sa_frames(SurfaceFile);
nV = size(F.Vtx, 1);
for name = {'fiedler','tang','xyz'}
    fr = F.(name{1});
    assert(isequal(size(fr.e1), [nV 3]) && isequal(size(fr.e2), [nV 3]), ...
        'Frame %s must be [nV x 3].', name{1});
    % Where defined, e1 and e2 are unit and orthogonal.
    nrm = sqrt(sum(fr.e1.^2, 2));
    def = nrm > 1e-6;
    assert(all(abs(nrm(def) - 1) < 1e-3), 'Frame %s e1 not unit.', name{1});
    dt  = abs(sum(fr.e1(def,:) .* fr.e2(def,:), 2));
    assert(all(dt < 1e-3), 'Frame %s e1/e2 not orthogonal.', name{1});
end
assert(isequal(size(F.R.Field), [nV 3]), 'Fiedler field must be [nV x 3].');
assert(isequal(size(F.R.Magnitude), [nV 1]), 'R.Magnitude must be [nV x 1].');
assert(numel(F.R.Phase) == nV, 'R.Phase must have nV elements.');
assert(isstruct(F.R) && isfield(F.R, 'Singularities'), 'R must have Singularities field.');
assert(~isempty(F.ConnEig.ConnLaplacian), 'sa_frames must expose ConnLaplacian K.');
fprintf('PASSED: test_sa_frames (%d vertices).\n', nV);
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
        if size(T.Vertices, 1) == 20484
            SurfaceFile = surf(iF).FileName;
            return;
        end
    end
end
end
