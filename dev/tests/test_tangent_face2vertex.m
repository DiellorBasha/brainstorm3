function test_tangent_face2vertex
% Transfer the per-face trivial-connection frame (tess_tangents) to a per-vertex
% orthonormal frame and check it is unit, orthonormal, right-handed, and smooth.
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
F  = double(TessMat.Faces);
N  = TessMat.VertNormals;
nV = size(N, 1);

% Per-face trivial-connection frame (e1 = U).
[Uf, ~] = tess_tangents(SurfaceFile, 'NoSave', 1);

[Uv, Vv] = bst_tangent_face2vertex(F, Uf, N);

assert(isequal(size(Uv), [nV 3]) && isequal(size(Vv), [nV 3]), 'Uv/Vv must be nV x 3.');
assert(max(abs(sqrt(sum(Uv.^2,2)) - 1)) < 1e-6, 'Uv must be unit.');
assert(max(abs(sqrt(sum(Vv.^2,2)) - 1)) < 1e-6, 'Vv must be unit.');
assert(max(abs(sum(Uv .* Vv, 2))) < 1e-6, 'Uv . Vv must be ~0.');
Nu = N ./ max(sqrt(sum(N.^2,2)), eps);
assert(max(abs(sum(Uv .* Nu, 2))) < 1e-6, 'Uv must lie in the tangent plane (Uv . n ~0).');
cr = cross(Nu, Uv, 2);
assert(max(max(abs(cr - Vv))) < 1e-6, 'Vv must equal n x Uv (right-handed).');

% Smoothness: the per-vertex e1 should agree with most incident face e1's
% (the trivial-connection field is smooth away from the few singularities).
%
% NOTE: Uf lives in the FACE tangent plane; Uv lives in the VERTEX tangent plane.
% On a curved mesh these planes differ (face vs vertex normals diverge).  Comparing
% across planes introduces a cosine-of-dihedral bias that limits |Uv . Uf| even for
% a perfectly smooth field.  The geometrically correct check projects Uf into the
% vertex tangent plane first, then measures angular agreement in the same space.
Nv1     = Nu(F(:,1), :);                      % outward vertex normal at vertex 1
Uf_proj = Uf - sum(Uf .* Nv1, 2) .* Nv1;    % project face e1 into vertex tangent plane
Uf_plen = sqrt(sum(Uf_proj.^2, 2));
validProj = Uf_plen > 0.1;                    % exclude near-singularity faces
Uf_pnrm = Uf_proj ./ max(Uf_plen, eps);
agree = abs(sum(Uv(F(:,1),:) .* Uf_pnrm, 2));
frac  = mean(agree(validProj) > 0.9);
assert(frac > 0.9, 'Per-vertex frame should align with incident face frames for >90%% of faces (got %.2f).', frac);

fprintf('PASSED: per-vertex FS frame (nV=%d): orthonormal, right-handed, smooth (%.1f%% aligned).\n', nV, 100*frac);
fprintf('ALL TESTS PASSED: test_tangent_face2vertex\n');
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
