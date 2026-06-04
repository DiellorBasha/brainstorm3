function test_connection_laplacian_nsym
% tess_connection_laplacian assembles for each n-RoSy symmetry (vector/line/cross)
% and the operator stays complex Hermitian. Real 20484-vertex cortex.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute plugin required for this test: %s', errMsg);
bst_plugin('Load', 'nxr-compute');

SurfaceFile = find_cortex_20484V();
if isempty(SurfaceFile)
    fprintf('SKIP: no 20484-vertex cortex in the current protocol (e.g. load TutorialAuditory / Subject01).\n');
    return;
end
TessMat = in_tess_bst(SurfaceFile);
V  = TessMat.Vertices;
F  = double(TessMat.Faces);
nV = size(V, 1);

for nSym = [1 2 4]
    [K, ~, Info] = tess_connection_laplacian(V, F, 'nSym', nSym);
    assert(Info.nSym == nSym, 'Info.nSym should be %d.', nSym);
    assert(isequal(size(K), [nV nV]), 'K must be nV x nV for nSym=%d.', nSym);
    assert(max(max(abs(K - K'))) < 1e-9, 'K must be Hermitian for nSym=%d.', nSym);
    fprintf('PASSED: nSym=%d assembled, Hermitian.\n', nSym);
end

fprintf('ALL TESTS PASSED: test_connection_laplacian_nsym\n');
end


function SurfaceFile = find_cortex_20484V()
% Return a Cortex surface with exactly 20484 vertices in the current protocol,
% preferring one with a FreeSurfer registration sphere; '' if none.
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
