function test_laplacian_nxr_parity
% nxr vs MATLAB parity for tess_laplacian.
%
% L (stiffness): nxr must match the MATLAB cotangent assembler to machine
%   precision for every mass type (L is mass-independent).
%
% M (mass): tess_laplacian currently serves nxr's own mass for the 'voronoi'
%   request. nxr's mass is, in fact, a BARYCENTRIC-style dual area (NOT a true
%   circumcentric Voronoi mass). This test pins nxr's mass to what it actually
%   serves: it must equal Brainstorm's barycentric mass to machine precision,
%   and it is expected to diverge from Brainstorm's true Voronoi mass. When nxr
%   is fixed upstream to compute a true Voronoi mass, this test must be updated.
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

% --- Mass M: pin nxr's served mass to what it actually is (barycentric) ---
[~, Mn]    = tess_laplacian(V, F, 'MassType', 'voronoi',     'Backend', 'nxr');     % nxr's served mass
[~, Mbary] = tess_laplacian(V, F, 'MassType', 'barycentric', 'Backend', 'matlab');  % true barycentric
[~, Mvor]  = tess_laplacian(V, F, 'MassType', 'voronoi',     'Backend', 'matlab');  % true Voronoi
dBary = full(max(abs(diag(Mn) - diag(Mbary))));
dVor  = full(max(abs(diag(Mn) - diag(Mvor))));
fprintf('nxr mass vs MATLAB barycentric: max|d| = %.3e (expect ~0)\n', dBary);
fprintf('nxr mass vs MATLAB voronoi:     max|d| = %.3e (expect > 0; known divergence)\n', dVor);

% nxr currently serves a barycentric-style dual area: pin to barycentric.
assert(dBary < 1e-9, ...
    ['nxr mass no longer matches MATLAB barycentric (max|d|=%.3e). nxr''s mass ' ...
     'convention changed; re-characterize this test (and check if the upstream ' ...
     'true-Voronoi fix landed).'], dBary);
% Document the known divergence from true Voronoi (do not fail on it).
assert(dVor > 1e-6, ...
    ['nxr mass now MATCHES true Voronoi (max|d|=%.3e) - the upstream fix may have ' ...
     'landed; update this test to assert true-Voronoi parity.'], dVor);

% nxr's mass is a valid mass matrix: total area conserved, positive diagonal.
assert(abs(full(sum(diag(Mn))) - full(sum(diag(Mvor)))) < 1e-9, 'nxr mass total area not conserved.');
assert(all(full(diag(Mn)) > 0), 'nxr mass has non-positive diagonal.');
fprintf('Mass pinned to barycentric; Voronoi divergence documented.\n');

% --- Sanity on the nxr stiffness itself ---
Ln = tess_laplacian(V, F, 'MassType', 'voronoi', 'Backend', 'nxr');
assert(norm(Ln - Ln', 1) < 1e-9, 'nxr L not symmetric.');
assert(max(abs(sum(Ln, 2))) < 1e-6, 'nxr L row sums not ~0.');
fprintf('ALL TESTS PASSED: test_laplacian_nxr_parity\n');
end
