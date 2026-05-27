function test_laplacian_ico
% Verify cotangent Laplacian on an icosphere, and the CheckManifold warning.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

[V, F] = tess_sphere(642);
[L, M] = tess_laplacian(V, F);   % default barycentric mass

% L symmetric
assert(norm(L - L', 1) < 1e-9, 'L is not symmetric.');
% Barycentric M is diagonal and strictly positive
md = full(diag(M));
assert(all(md > 0), 'M has non-positive diagonal entries.');
assert(nnz(M - spdiags(md, 0, numel(md), numel(md))) == 0, 'Barycentric M is not diagonal.');
% L positive semidefinite
lmin = eigs(L, 1, 'smallestreal');
assert(lmin > -1e-6, 'L is not PSD (min eig %.3e).', lmin);
% Constant function is in the null space: row sums ~ 0
assert(max(abs(sum(L, 2))) < 1e-6, 'L row sums are not ~0.');
fprintf('PASSED: L symmetric/PSD, M positive diagonal, constant null space.\n');

% CheckManifold warns on non-manifold input
vExtra = find(~ismember(1:size(V,1), F(1,:)), 1);
Fnm = [F; F(1,1), F(1,2), vExtra];
lastwarn('');
tess_laplacian(V, Fnm, 'CheckManifold', true);
[w, wid] = lastwarn();
assert(~isempty(w), 'CheckManifold=true did not warn on non-manifold input.');
assert(strcmp(wid, 'tess_laplacian:NonManifold'), ...
    'Expected tess_laplacian:NonManifold warning, got: %s', wid);
fprintf('PASSED: CheckManifold warns on non-manifold input.\n');

fprintf('ALL TESTS PASSED: test_laplacian_ico\n');
end
