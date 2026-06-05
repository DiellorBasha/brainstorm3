function test_sa_sulcal_walls
% sa_sulcal_walls returns facing sulcal-wall pairs; every pair satisfies the
% sulcal / 3mm / anti-normal / non-adjacent criteria. SKIP if no cortex.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'dev', 'benchmarks', 'sign_ambiguity'));
if ~brainstorm('status'), brainstorm nogui; end
SurfaceFile = find_cortex_20484V();
if isempty(SurfaceFile)
    fprintf('SKIP: no 20484-vertex cortex.\n'); return;
end

TessMat = in_tess_bst(SurfaceFile);
Vtx = TessMat.Vertices; Fcs = double(TessMat.Faces); Nv = TessMat.VertNormals;
if isfield(TessMat,'VertConn') && ~isempty(TessMat.VertConn)
    VertConn = TessMat.VertConn;
else
    VertConn = tess_vertconn(Vtx, Fcs);
end
if isfield(TessMat,'SulciMap') && ~isempty(TessMat.SulciMap)
    SulciMap = TessMat.SulciMap;
else
    SulciMap = tess_sulcimap(TessMat);
end

opts = struct('MaxDist', 0.003, 'NormalDot', -0.7, 'Nring', 3);
pairs = sa_sulcal_walls(Vtx, Nv, SulciMap, VertConn, opts);
assert(~isempty(pairs), 'Expected a non-empty set of sulcal-wall pairs.');
assert(size(pairs,2) == 2, 'pairs must be [nPairs x 2].');

i = pairs(:,1); j = pairs(:,2);
% Both endpoints sulcal.
assert(all(SulciMap(i) > 0 & SulciMap(j) > 0), 'Both endpoints must be sulcal.');
% Within MaxDist.
d = sqrt(sum((Vtx(i,:) - Vtx(j,:)).^2, 2));
assert(all(d <= opts.MaxDist + 1e-9), 'Pairs must be within MaxDist.');
% Anti-aligned normals.
nd = sum(Nv(i,:) .* Nv(j,:), 2);
assert(all(nd < opts.NormalDot + 1e-9), 'Pairs must have anti-aligned normals.');
% Not mesh-adjacent within Nring.
A = double(VertConn > 0);
Aring = A; P = A;
for r = 2:opts.Nring, P = double((P*A) > 0); Aring = double((Aring + P) > 0); end
lin = sub2ind(size(Aring), i, j);
assert(all(Aring(lin) == 0), 'Pairs must not be within the N-ring.');
fprintf('PASSED: test_sa_sulcal_walls (%d pairs).\n', size(pairs,1));
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
