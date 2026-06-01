function test_tess_tangents
% Integration test for tess_tangents on a real FreeSurfer-registered cortex:
% per-face orthonormal frame, right-handedness, determinism, and storage
% round-trip. Runs on a temp COPY so the DB surface is not mutated.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% --- Find a low-res FreeSurfer-registered cortex in the DB ---
srcFile = find_registered_cortex();
assert(~isempty(srcFile), 'No FreeSurfer-registered cortex found in the DB to test against.');
fprintf('Source cortex: %s\n', srcFile);

% --- Work on a temp copy (named tess_cortex_*.mat so file_fullpath accepts it) ---
tmpFile = fullfile(tempdir, 'tess_cortex_tangents_test.mat');
copyfile(file_fullpath(srcFile), tmpFile);
cleanup = onCleanup(@() delete(tmpFile));

% --- Compute + store ---
[U, V] = tess_tangents(tmpFile);
TessMat = in_tess_bst(tmpFile);
Vtx = TessMat.Vertices;  Fcs = double(TessMat.Faces);
nF  = size(Fcs, 1);
assert(isequal(size(U), [nF 3]), 'U must be nF x 3.');
assert(isequal(size(V), [nF 3]), 'V must be nF x 3.');

% --- Orthonormal, right-handed frame ---
assert(all(abs(sqrt(sum(U.^2,2)) - 1) < 1e-4), 'U not unit length.');
assert(all(abs(sqrt(sum(V.^2,2)) - 1) < 1e-4), 'V not unit length.');
assert(max(abs(sum(U.*V, 2))) < 1e-4, 'U and V not orthogonal.');
fn = cross(Vtx(Fcs(:,2),:) - Vtx(Fcs(:,1),:), Vtx(Fcs(:,3),:) - Vtx(Fcs(:,1),:));
fn = fn ./ max(sqrt(sum(fn.^2,2)), eps);
assert(all(sum(cross(U, V, 2) .* fn, 2) > -1e-6), 'Frame not right-handed w.r.t. face normal.');
fprintf('PASSED: per-face orthonormal, right-handed frame (%d faces).\n', nF);

% --- Storage round-trip ---
assert(isfield(TessMat, 'TangentFrame'), 'TangentFrame not stored.');
TF = TessMat.TangentFrame;
assert(strcmp(TF.Domain, 'face'), 'TangentFrame.Domain should be ''face''.');
assert(isequal(single(U), TF.U) && isequal(single(V), TF.V), 'Stored U/V mismatch.');
assert(numel(TF.Singularities.Vertices) == 4, 'Expected 4 singularities (2 per hemisphere).');
assert(all(TF.Singularities.Indices == 1), 'Singularity indices should all be +1.');
fprintf('PASSED: TangentFrame stored + reloaded (4 pole singularities).\n');

% --- Determinism ---
[U2, V2] = tess_tangents(tmpFile, 'NoSave', 1);
assert(isequal(U, U2) && isequal(V, V2), 'tess_tangents is not deterministic.');
fprintf('PASSED: deterministic across runs.\n');

fprintf('ALL TESTS PASSED: test_tess_tangents\n');
end


function SurfaceFile = find_registered_cortex()
% Return a low-res cortex FileName that has Reg.Sphere.Vertices, or '' if none.
SurfaceFile = '';
best = inf;
sSubjects = bst_get('ProtocolSubjects');
allSubj = [sSubjects.Subject];
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex'), continue; end
        try
            T = load(file_fullpath(surf(iF).FileName), 'Reg', 'Vertices');
        catch
            continue;
        end
        if isfield(T,'Reg') && isstruct(T.Reg) && isfield(T.Reg,'Sphere') ...
           && isfield(T.Reg.Sphere,'Vertices') && ~isempty(T.Reg.Sphere.Vertices)
            n = size(T.Vertices, 1);
            if n < best          % prefer the smallest (fastest to solve)
                best = n;
                SurfaceFile = surf(iF).FileName;
            end
        end
    end
end
end
