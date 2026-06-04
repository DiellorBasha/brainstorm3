function test_connection_laplacian_smoke
% Smoke test for tess_connection_laplacian on a real 20484-vertex Brainstorm
% cortex (e.g. TutorialAuditory Subject01 tess_cortex_*_low): the complex-
% Hermitian vertex connection Laplacian, its lumped mass, and the Info struct.
% Pure (V,F) operator — no storage, no readout.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% --- Require + load nxr-compute (operator is nxr-only) ---
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute plugin required for this test: %s', errMsg);
bst_plugin('Load', 'nxr-compute');

% --- Resolve a real 20484-vertex cortex from the current protocol ---
SurfaceFile = find_cortex_20484V();
if isempty(SurfaceFile)
    fprintf('SKIP: no 20484-vertex cortex in the current protocol (e.g. load TutorialAuditory / Subject01).\n');
    return;
end
fprintf('Source cortex: %s\n', SurfaceFile);
TessMat = in_tess_bst(SurfaceFile);
V  = TessMat.Vertices;
F  = double(TessMat.Faces);
nV = size(V, 1);

% --- Assemble (defaults: vertex domain, nSym=1, complex) ---
[K, M, Info] = tess_connection_laplacian(V, F);

% --- Operator shape / type / symmetry ---
assert(isequal(size(K), [nV nV]), 'K must be nV x nV.');
assert(issparse(K), 'K must be sparse.');
assert(~isreal(K), 'K must be complex (nonzero imaginary part on a curved cortex).');
assert(max(max(abs(K - K'))) < 1e-9, 'K must be Hermitian (conjugate-symmetric).');

% --- Mass ---
assert(isequal(size(M), [nV nV]), 'M must be nV x nV.');
assert(isdiag(M), 'M must be diagonal (lumped vertex mass).');
assert(all(diag(M) > 0), 'M diagonal entries must be positive.');

% --- Info ---
assert(Info.nSym == 1, 'Info.nSym should default to 1.');
assert(strcmp(Info.Domain, 'vertex'), 'Info.Domain should be ''vertex''.');
assert(strcmp(Info.Format, 'complex'), 'Info.Format should be ''complex''.');
assert(strcmp(Info.Backend, 'nxr'), 'Info.Backend should be ''nxr''.');
assert(Info.baseDim == nV, 'Info.baseDim should equal nV.');

fprintf('PASSED: connection Laplacian assembled (nV=%d): sparse, complex, Hermitian; lumped mass; Info OK.\n', nV);
fprintf('ALL TESTS PASSED: test_connection_laplacian_smoke\n');
end


function SurfaceFile = find_cortex_20484V()
% Return the FileName of a Cortex surface with exactly 20484 vertices ("20484V")
% in the current protocol, preferring one with a FreeSurfer registration sphere;
% '' if none (e.g. no suitable protocol loaded).
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
        if size(T.Vertices, 1) ~= 20484   % "20484V" = vertex count, not a filename
            continue;
        end
        hasReg = isfield(T,'Reg') && isstruct(T.Reg) && isfield(T.Reg,'Sphere') ...
                 && isfield(T.Reg.Sphere,'Vertices') && ~isempty(T.Reg.Sphere.Vertices);
        if hasReg
            SurfaceFile = surf(iF).FileName;   % prefer a registered cortex
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
