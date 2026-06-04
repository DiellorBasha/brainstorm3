function test_connection_laplacian_spectrum
% Spectral sanity for tess_connection_laplacian on a real 20484-vertex cortex.
% nxr's 'solve' is real-only and cannot consume the complex K, so the
% generalized problem K*phi = lambda*M*phi is solved with MATLAB eigs (which
% handles a Hermitian K and a Hermitian-PD M). This is a light operator-level
% sanity check; full spectral validation belongs to the eigenmode milestone.
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
V = TessMat.Vertices;
F = double(TessMat.Faces);

[K, M] = tess_connection_laplacian(V, F);

% Smallest k modes of the complex Hermitian generalized eigenproblem.
% (The whole-mesh operator is block-diagonal across the two hemispheres.)
k = 6;
[Phi, Lam] = eigs(K, M, k, 'smallestabs');
lam = diag(Lam);

assert(isequal(size(Phi), [size(K,1) k]), 'Phi must be N x k.');
assert(max(abs(imag(lam))) < 1e-6 * max(abs(real(lam)) + 1), ...
    'Eigenvalues of a Hermitian operator must be real.');
lam = sort(real(lam), 'ascend');
assert(all(lam > -1e-6), 'Connection-Laplacian eigenvalues must be >= 0.');

fprintf('PASSED: %d smallest eigenvalues real and >= 0; range [%.3g, %.3g].\n', ...
    k, lam(1), lam(end));
fprintf('ALL TESTS PASSED: test_connection_laplacian_spectrum\n');
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
