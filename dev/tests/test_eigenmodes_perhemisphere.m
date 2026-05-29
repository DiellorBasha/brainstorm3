function test_eigenmodes_perhemisphere
% Verify tess_eigenmodes solves per connected component and tags modes with
% Component/CompRank, keeps each mode localized to its component, removes one
% DC mode per component, and stays globally M-orthonormal.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% ----- Two disjoint octahedra (each a closed 2-manifold) -----
V1 = [ 1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
F1 = [5 1 3; 5 3 2; 5 2 4; 5 4 1; 6 3 1; 6 2 3; 6 4 2; 6 1 4];
V  = [V1; bsxfun(@plus, V1, [10 0 0])];   % second component offset in x
F  = [F1; F1 + 6];
nV = size(V,1);

nModes = 3;
[Eig, ~, M] = tess_eigenmodes(V, F, 'nModes', nModes, 'Verbose', 0);

% Two components, nModes each
assert(size(Eig.Vectors,2) == 2*nModes, 'Expected nModes per component (2 components).');
assert(isequal(Eig.Component(:), [1;1;1;2;2;2]), 'Component labels wrong.');
assert(isequal(Eig.CompRank(:),  [1;2;3;1;2;3]), 'CompRank wrong.');
assert(Eig.nRemoved == 2, 'Expected one DC mode removed per component.');

% Localization: component-1 modes zero on component-2 vertices and vice versa
c1 = (Eig.Component==1);  c2 = (Eig.Component==2);
assert(max(max(abs(Eig.Vectors(7:12, c1)))) < 1e-8, 'Comp-1 modes must vanish on comp-2 verts.');
assert(max(max(abs(Eig.Vectors(1:6,  c2)))) < 1e-8, 'Comp-2 modes must vanish on comp-1 verts.');

% Global M-orthonormality (block-diagonal M => per-component orthonormality)
G = Eig.Vectors' * M * Eig.Vectors;
assert(norm(G - eye(size(G)), 'fro') < 1e-6, 'Modes not M-orthonormal.');

% ----- Single-component mesh still works (one octahedron) -----
[EigS, ~, ~] = tess_eigenmodes(V1, F1, 'nModes', 2, 'Verbose', 0);
assert(all(EigS.Component == 1), 'Single component must be labeled 1.');
assert(isequal(EigS.CompRank(:), (1:numel(EigS.Values))'), 'Single-component CompRank must be 1..K.');

disp('ALL TESTS PASSED');
end
