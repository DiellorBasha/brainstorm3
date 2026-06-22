function test_conn_phase
% (A) Gauge cross-check: decoding nxr's smoothVertex raw field via the exported
%     per-vertex frame must align with its world vectors -> the frame is the
%     correct gauge for vertex-domain complex coordinates (validates the decode
%     that bst_conn_phase relies on, incl. the connection eigenmodes which share
%     the same gauge).
% (B) bst_conn_phase outputs on the real cortex: gauge-independent 3D field
%     (tangent), magnitude = |z|, FS-gauge phase that winds, and singularities.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute plugin required: %s', errMsg);
bst_plugin('Load', 'nxr-compute');

srcFile = find_cortex_20484V();
if isempty(srcFile)
    fprintf('SKIP: no 20484-vertex cortex in the current protocol.\n');
    return;
end
% Work on a temp COPY so bst_conn_eigenmodes_ensure does not write ConnEigenmodes
% into the shared DB surface (named tess_cortex_*.mat so file_fullpath accepts it).
tmpFile = fullfile(tempdir, 'tess_cortex_connphase.mat');
copyfile(file_fullpath(srcFile), tmpFile);
cleanup = onCleanup(@() delete(tmpFile));

TessMat = in_tess_bst(tmpFile);
V  = TessMat.Vertices;
F  = double(TessMat.Faces);
N  = TessMat.VertNormals;
nV = size(V, 1);

mctx   = nxr.manifold.context(V, F);
vFrame = nxr.manifold.measure.vertexFrame(mctx);

% ===== (A) Gauge cross-check via smoothVertex (nSym=1) =====
sm    = nxr.manifold.interpolate.smoothVertex(mctx, 1);   % nSym=1: true vector field
rawFlat = sm.vertexFieldRaw;   % [2*nV x 1] real, interleaved: (e1,e2) per vertex
raw   = reshape(rawFlat, 2, nV).';  % [nV x 2]: col1=e1-coord, col2=e2-coord
world = sm.vertexVectors;      % [nV x 3] world-space vectors
decoded = raw(:,1) .* vFrame.e1 + raw(:,2) .* vFrame.e2;   % [nV x 3]
dn = decoded ./ max(sqrt(sum(decoded.^2,2)), eps);
wn = world   ./ max(sqrt(sum(world.^2,2)),   eps);
dp = abs(sum(dn .* wn, 2));    % |cos angle|; allow global nRoSy sign
mask = (sqrt(sum(decoded.^2,2)) > 1e-6) & (sqrt(sum(world.^2,2)) > 1e-6);
assert(median(dp(mask)) > 0.99, ...
    'GAUGE: decoded smoothVertex raw field must align with its world vectors (median |cos|=%.4f).', median(dp(mask)));
fprintf('PASSED (A): vertex-frame gauge validated against smoothVertex (median |cos|=%.4f).\n', median(dp(mask)));

% ===== (B) bst_conn_phase =====
ConnEig = bst_conn_eigenmodes_ensure(tmpFile, 20);   % small, fast; explicit count skips scalar-axis ensure
[Uf, ~] = tess_tangents(tmpFile, 'NoSave', 1);
[Uv, Vv] = bst_tangent_face2vertex(F, Uf, N);
FsFrame = struct('e1', Uv, 'e2', Vv);

R = bst_conn_phase(ConnEig, vFrame, 'Rank', 1, 'FsFrame', FsFrame, 'nSing', 2);

assert(isequal(size(R.Field), [nV 3]), 'Field must be nV x 3.');
supp = find(any(R.Field ~= 0, 2));
assert(~isempty(supp), 'Field must be nonzero on the Fiedler support.');
% Field is tangent to the (nxr) vertex normal by construction.
tang = abs(sum(R.Field(supp,:) .* vFrame.normals(supp,:), 2));
assert(max(tang) < 1e-6, 'Field must be tangent (Field . n ~0).');
% Magnitude equals |z| (nonneg, positive on support).
assert(all(R.Magnitude(supp) > 0), 'Magnitude must be positive on support.');
% FS-gauge phase: finite on support, within [-pi, pi].
assert(all(isfinite(R.Phase(supp))), 'Phase must be finite on support.');
assert(max(abs(R.Phase(supp))) <= pi + 1e-9, 'Phase must lie in [-pi, pi].');
% Phase winds: spans more than half the circle on the support of one component.
comp1 = find(ConnEig.Component(:)==1 & ConnEig.CompRank(:)==1, 1);
idx1 = find(ConnEig.Vectors(:, comp1) ~= 0);
assert((max(R.Phase(idx1)) - min(R.Phase(idx1))) > pi, ...
    'FS-gauge phase should wind over a hemisphere (range > pi).');
% Singularities: nSing per component, with small magnitude (a field dip).
nComp = numel(unique(ConnEig.Component(:)));
assert(numel(R.Singularities) == 2 * nComp, 'Expected 2 singularities per component.');
assert(median(R.Magnitude(R.Singularities)) < median(R.Magnitude(supp)), ...
    'Singularity magnitude should dip below the component median.');

fprintf('PASSED (B): bst_conn_phase field/magnitude/phase/singularities on the cortex (nV=%d).\n', nV);
fprintf('ALL TESTS PASSED: test_conn_phase\n');
end


function SurfaceFile = find_cortex_20484V()
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
if isempty(sSubjects)
    return;
end
allSubj = [sSubjects.Subject];
fallback = '';
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex')
            continue;
        end
        try
            T = load(file_fullpath(surf(iF).FileName), 'Vertices', 'Reg');
        catch
            continue;
        end
        if size(T.Vertices, 1) ~= 20484
            continue;
        end
        hasReg = isfield(T,'Reg') && isstruct(T.Reg) && isfield(T.Reg,'Sphere') ...
                 && isfield(T.Reg.Sphere,'Vertices') && ~isempty(T.Reg.Sphere.Vertices);
        if hasReg
            SurfaceFile = surf(iF).FileName;
            return;
        elseif isempty(fallback)
            fallback = surf(iF).FileName;
        end
    end
end
if isempty(SurfaceFile)
    SurfaceFile = fallback;
end
end
