function test_laplacian_nxr_parity
% nxr vs MATLAB parity for tess_laplacian.
%
% L (stiffness): nxr must match the MATLAB cotangent assembler to machine
%   precision for every mass type (L is mass-independent).
%
% M (mass): nxr's MATLAB binding exposes one mass — geometry-central's
%   vertexLumpedMassMatrix (diag A/3 per vertex), i.e. the BARYCENTRIC mass.
%   (Older nxr mislabelled this "Voronoi"; it never was a circumcentric/mixed
%   Voronoi mass. Upstream nxr now names it Lumped, matching geometry-central.)
%   tess_laplacian therefore serves nxr's mass for the 'barycentric' request
%   (machine-precision parity) and uses the MATLAB assembler for 'voronoi' and
%   'galerkin'. This test verifies that routing.
%
% Stages + loads the locally-built nxr-compute plugin.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% Stage + load the plugin (same offline install as the lifecycle test)
codeDir  = fileparts(repoRoot);
nxrStage = fullfile(codeDir, 'nxr-compute', 'dist', 'plugin', 'nxr-compute');
assert(isfolder(nxrStage), 'Run nxr-compute/scripts/package-plugin.sh first (missing %s).', nxrStage);
destDir = fullfile(bst_get('UserPluginsDir'), 'nxr-compute');
if isfolder(destDir), file_delete(destDir, 1, 3); end
mkdir(destDir);
copyfile(fullfile(nxrStage, '*'), destDir);
[isOk, errMsg] = bst_plugin('Load', 'nxr-compute');
assert(isOk, 'Load failed: %s', errMsg);
cleanup = onCleanup(@() bst_plugin('Unload', 'nxr-compute'));

[V, F] = tess_sphere(642);

% --- Stiffness L parity: nxr must match MATLAB to machine precision ---
tolL = 1e-6;
for mt = {'barycentric', 'voronoi', 'galerkin'}
    Ln = tess_laplacian(V, F, 'MassType', mt{1}, 'Backend', 'nxr');
    Lm = tess_laplacian(V, F, 'MassType', mt{1}, 'Backend', 'matlab');
    dL = full(max(abs(Ln(:) - Lm(:))));
    fprintf('L parity (%s): max|dL| = %.3e\n', mt{1}, dL);
    assert(dL < tolL, 'Stiffness L parity failed for %s (max|dL|=%.3e).', mt{1}, dL);
end

% --- Mass M routing: nxr serves barycentric; voronoi/galerkin use MATLAB ---
% (1) 'barycentric' via nxr == MATLAB barycentric (machine precision).
[~, Mb_nxr] = tess_laplacian(V, F, 'MassType', 'barycentric', 'Backend', 'nxr');
[~, Mb_mat] = tess_laplacian(V, F, 'MassType', 'barycentric', 'Backend', 'matlab');
dBary = full(max(abs(diag(Mb_nxr) - diag(Mb_mat))));
fprintf('barycentric mass, nxr vs MATLAB: max|d| = %.3e (expect ~0)\n', dBary);
assert(dBary < 1e-9, 'nxr barycentric mass no longer matches MATLAB (max|d|=%.3e).', dBary);

% (2) 'voronoi' via the nxr backend must route to the MATLAB Voronoi mass,
%     NOT nxr's lumped/barycentric mass (the old mislabel bug).
[~, Mv_nxr] = tess_laplacian(V, F, 'MassType', 'voronoi', 'Backend', 'nxr');
[~, Mv_mat] = tess_laplacian(V, F, 'MassType', 'voronoi', 'Backend', 'matlab');
dVor = full(max(abs(diag(Mv_nxr) - diag(Mv_mat))));
fprintf('voronoi mass, nxr-backend vs MATLAB: max|d| = %.3e (expect ~0)\n', dVor);
assert(dVor < 1e-9, ...
    ['voronoi request via the nxr backend does not match the MATLAB Voronoi mass ' ...
     '(max|d|=%.3e); nxr must not serve its lumped mass for ''voronoi''.'], dVor);

% (3) 'galerkin' via the nxr backend must route to the MATLAB Galerkin mass.
[~, Mg_nxr] = tess_laplacian(V, F, 'MassType', 'galerkin', 'Backend', 'nxr');
[~, Mg_mat] = tess_laplacian(V, F, 'MassType', 'galerkin', 'Backend', 'matlab');
dGal = full(max(abs(Mg_nxr(:) - Mg_mat(:))));
fprintf('galerkin mass, nxr-backend vs MATLAB: max|d| = %.3e (expect ~0)\n', dGal);
assert(dGal < 1e-9, 'galerkin request via nxr backend does not match MATLAB Galerkin (max|d|=%.3e).', dGal);

% Sanity: barycentric and Voronoi genuinely differ on this irregular mesh, so
% the routing distinction above is meaningful (not a no-op).
assert(full(max(abs(diag(Mb_mat) - diag(Mv_mat)))) > 1e-6, ...
    'barycentric and Voronoi unexpectedly identical; test mesh too regular.');
% nxr's served mass: total area conserved (vs Voronoi), positive diagonal.
assert(abs(full(sum(diag(Mb_nxr))) - full(sum(diag(Mv_mat)))) < 1e-9, 'nxr mass total area not conserved.');
assert(all(full(diag(Mb_nxr)) > 0), 'nxr mass has non-positive diagonal.');
fprintf('Mass routing verified: nxr->barycentric; MATLAB->voronoi/galerkin.\n');

% --- Sanity on the nxr stiffness itself ---
Ln = tess_laplacian(V, F, 'MassType', 'barycentric', 'Backend', 'nxr');
assert(norm(Ln - Ln', 1) < 1e-9, 'nxr L not symmetric.');
assert(max(abs(sum(Ln, 2))) < 1e-6, 'nxr L row sums not ~0.');
fprintf('ALL TESTS PASSED: test_laplacian_nxr_parity\n');
end
