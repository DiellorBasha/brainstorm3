function test_conn_eigenmodes_compute
% Compute correctness for tess_conn_eigenmodes on a real 20484-vertex cortex:
% complex eigenvectors, real positive eigenvalues, per-component block structure,
% canonical Order, and M-orthonormality within a component.
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
    fprintf('SKIP: no 20484-vertex cortex in the current protocol (e.g. load TutorialAuditory / Subject01).\n');
    return;
end
TessMat = in_tess_bst(SurfaceFile);
V  = TessMat.Vertices;
F  = double(TessMat.Faces);
nV = size(V, 1);

k = 20;
[ConnEig, K, M] = tess_conn_eigenmodes(V, F, 'nModes', k);

% --- Type / shape ---
assert(~isreal(ConnEig.Vectors), 'Vectors must be complex.');
assert(size(ConnEig.Vectors, 1) == nV, 'Vectors must have nV rows.');
assert(ConnEig.nModes == size(ConnEig.Vectors, 2), 'nModes must equal column count.');
assert(isreal(ConnEig.Values) && all(ConnEig.Values > 0), 'Values must be real and > 0.');

% --- Metadata ---
assert(strcmp(ConnEig.OperatorType, 'Connection-LeviCivita'), 'OperatorType mismatch.');
assert(ConnEig.nSym == 1, 'nSym should be 1.');
assert(ConnEig.nRemoved == 0, 'No DC mode is removed for the connection bundle.');
assert(issparse(ConnEig.ConnLaplacian) && ~isreal(ConnEig.ConnLaplacian), 'ConnLaplacian must be complex sparse.');
assert(isequal(size(K), [nV nV]) && ~isreal(K), 'Returned K must be nV x nV complex.');

% --- Canonical order sorts Values ascending ---
assert(issorted(ConnEig.Values(ConnEig.Order)), 'Order must sort Values ascending.');

% --- Per-component checks (all components) ---
compId = conncomp(graph(tess_vertconn(V, F)));
assert(ConnEig.nComponents == max(compId), 'nComponents mismatch.');
for c = 1:ConnEig.nComponents
    idx  = find(compId == c);
    cols = find(ConnEig.Component == c);
    Uc   = ConnEig.Vectors(idx, cols);
    Mc   = M(idx, idx);
    % Block structure: a component's modes vanish off that component.
    other = setdiff((1:nV)', idx);
    assert(max(max(abs(ConnEig.Vectors(other, cols)))) < 1e-10, 'Component %d modes must vanish off-component.', c);
    % M-orthonormality (Hermitian inner product; '' is conjugate transpose).
    G = Uc' * Mc * Uc;
    assert(max(max(abs(G - eye(size(G))))) < 1e-6, 'Component %d modes must be M-orthonormal.', c);
end

fprintf('PASSED: %d-mode connection eigenbasis (nV=%d, %d components): complex, real>0 eigenvalues, block-structured, M-orthonormal.\n', ...
    ConnEig.nModes, nV, ConnEig.nComponents);
fprintf('ALL TESTS PASSED: test_conn_eigenmodes_compute\n');
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
