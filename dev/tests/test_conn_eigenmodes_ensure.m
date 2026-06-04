function test_conn_eigenmodes_ensure
% bst_conn_eigenmodes_ensure: (1) with no count, derives the per-component count
% from the surface's scalar Eigenmodes axis (match-scalar); (2) reuses an existing
% ConnEigenmodes axis idempotently. Works on a temp COPY; uses small mode counts.
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
    fprintf('SKIP: no 20484-vertex cortex in the current protocol (e.g. load TutorialAuditory / Subject01).\n');
    return;
end
tmpFile = fullfile(tempdir, 'tess_cortex_conneig_ensure.mat');
copyfile(file_fullpath(srcFile), tmpFile);
cleanup = onCleanup(@() delete(tmpFile));

TessMat = in_tess_bst(tmpFile);
V = TessMat.Vertices;
F = double(TessMat.Faces);

% --- Seed a SMALL scalar Eigenmodes axis so the match-scalar derivation is fast ---
sEig = tess_eigenmodes(V, F, 'nModes', 15, 'MassType', 'barycentric', 'RemoveDC', 1);
out_tess_eigenmodes(tmpFile, sEig, V, F);
sEigStored = in_tess_eigenmodes(tmpFile);

% --- Ensure with no count: derives + computes + stores ---
ConnEig = bst_conn_eigenmodes_ensure(tmpFile);
assert(ConnEig.nComponents == sEigStored.nComponents, 'Component count should match the mesh.');
% Independent match-scalar check: the connection axis must carry the same
% per-component mode count as the scalar axis it was matched to (compared against
% the scalar axis's ACTUAL per-component counts, not the production formula). The
% production count is the rounded per-mesh average, so allow a 1-mode rounding gap.
for c = 1:ConnEig.nComponents
    nConn   = sum(ConnEig.Component == c);
    nScalar = sum(sEigStored.Component == c);
    assert(abs(nConn - nScalar) <= 1, ...
        'Component %d: connection modes (%d) must match scalar modes (%d) within rounding.', ...
        c, nConn, nScalar);
end
fprintf('PASSED: match-scalar per-component count (component 1: %d conn vs %d scalar).\n', ...
    sum(ConnEig.Component == 1), sum(sEigStored.Component == 1));

% --- The scalar Eigenmodes axis must be untouched by the connection ensure ---
sEigAfter = in_tess_eigenmodes(tmpFile);
assert(sEigAfter.nModes == sEigStored.nModes && sEigAfter.nComponents == sEigStored.nComponents, ...
    'Scalar Eigenmodes axis must remain intact after the connection ensure.');
fprintf('PASSED: scalar Eigenmodes axis intact (%d modes).\n', sEigAfter.nModes);

% --- Second call reuses (idempotent: returns the stored struct, not a recompute) ---
ConnEig2 = bst_conn_eigenmodes_ensure(tmpFile);
assert(ConnEig2.nModes == ConnEig.nModes, 'Reuse must return the same nModes.');
assert(ConnEig2.ComputeTime == ConnEig.ComputeTime, ...
    'Reuse must return the stored struct (identical ComputeTime), not recompute.');
fprintf('PASSED: reuse is idempotent (%d modes, ComputeTime preserved).\n', ConnEig2.nModes);

fprintf('ALL TESTS PASSED: test_conn_eigenmodes_ensure\n');
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
