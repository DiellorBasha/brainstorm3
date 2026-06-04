function test_vertex_frame_plugin
% Verifies nxr.manifold.measure.vertexFrame is callable through the INSTALLED
% Brainstorm plugin (after restaging the rebuilt binary), on the real cortex.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute plugin required: %s', errMsg);
bst_plugin('Load', 'nxr-compute');

SurfaceFile = find_cortex_20484V();
if isempty(SurfaceFile)
    fprintf('SKIP: no 20484-vertex cortex in the current protocol.\n');
    return;
end
TessMat = in_tess_bst(SurfaceFile);
V = TessMat.Vertices;
F = double(TessMat.Faces);
nV = size(V, 1);

mctx = nxr.manifold.context(V, F);
vf = nxr.manifold.measure.vertexFrame(mctx);

assert(isequal(size(vf.e1), [nV 3]) && isequal(size(vf.e2), [nV 3]) && isequal(size(vf.normals), [nV 3]), ...
    'vertexFrame fields must be nV x 3.');
assert(max(abs(sum(vf.e1 .* vf.e2, 2))) < 1e-6, 'e1 . e2 must be ~0.');
assert(max(abs(sqrt(sum(vf.e1.^2,2)) - 1)) < 1e-6, 'e1 must be unit.');
cr = cross(vf.e1, vf.e2, 2);
assert(max(max(abs(cr - vf.normals))) < 1e-6, 'e1 x e2 must equal n.');

fprintf('PASSED: nxr.manifold.measure.vertexFrame works through the installed plugin (nV=%d).\n', nV);
fprintf('ALL TESTS PASSED: test_vertex_frame_plugin\n');
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
