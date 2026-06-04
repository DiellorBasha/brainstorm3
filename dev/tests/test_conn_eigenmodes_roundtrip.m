function test_conn_eigenmodes_roundtrip
% out_tess_conn_eigenmodes then in_tess_conn_eigenmodes preserves the complex
% eigenvectors (single round-trip) and the complex operator. Works on a temp COPY
% so the DB surface is not mutated.
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
% Temp copy named tess_cortex_*.mat so file_fullpath accepts the absolute path.
tmpFile = fullfile(tempdir, 'tess_cortex_conneig_roundtrip.mat');
copyfile(file_fullpath(srcFile), tmpFile);
cleanup = onCleanup(@() delete(tmpFile));

TessMat = in_tess_bst(tmpFile);
V = TessMat.Vertices;
F = double(TessMat.Faces);

% Not computed yet.
[~, isComputed0] = in_tess_conn_eigenmodes(tmpFile);
assert(~isComputed0, 'ConnEigenmodes should be absent before computing.');

% Compute + store.
ConnEig = tess_conn_eigenmodes(V, F, 'nModes', 20);
out_tess_conn_eigenmodes(tmpFile, ConnEig, V, F);

% Load back.
[ConnEig2, isComputed] = in_tess_conn_eigenmodes(tmpFile);
assert(isComputed, 'ConnEigenmodes should be present after store.');
assert(~isreal(ConnEig2.Vectors), 'Loaded Vectors must be complex.');
assert(isa(ConnEig2.Vectors, 'double'), 'Loaded Vectors must be cast to double.');
assert(isequal(size(ConnEig2.Vectors), size(ConnEig.Vectors)), 'Vector shape must be preserved.');
% Stored as single, loaded as double(single): compare against the single-cast original.
d = max(max(abs(double(single(ConnEig.Vectors)) - ConnEig2.Vectors)));
assert(d < 1e-5, 'Vectors must survive the single round-trip (max diff = %g).', d);
assert(isfield(ConnEig2, 'ConnLaplacian') && issparse(ConnEig2.ConnLaplacian) ...
       && ~isreal(ConnEig2.ConnLaplacian), 'ConnLaplacian must round-trip as complex sparse.');
assert(strcmp(ConnEig2.OperatorType, 'Connection-LeviCivita'), 'OperatorType must round-trip.');
assert(ConnEig2.nModes == ConnEig.nModes, 'nModes must round-trip.');

fprintf('PASSED: ConnEigenmodes round-trip (%d modes; vector max diff %g; complex operator preserved).\n', ...
    ConnEig2.nModes, d);
fprintf('ALL TESTS PASSED: test_conn_eigenmodes_roundtrip\n');
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
