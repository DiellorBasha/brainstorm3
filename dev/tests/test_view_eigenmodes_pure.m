function test_view_eigenmodes_pure
% Verify the paired-grid builder: column k holds each component's rank-k mode,
% summed (disjoint support), so a single display column shows both hemispheres.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% Fabricate a 2-component Eigenmodes struct: comp 1 on verts 1:4, comp 2 on 5:8
nV = 8;
Eig = struct();
Eig.Vectors = zeros(nV, 4);
Eig.Vectors(1:4, 1) = [1;2;3;4];     % comp 1, rank 1
Eig.Vectors(1:4, 2) = [5;6;7;8];     % comp 1, rank 2
Eig.Vectors(5:8, 3) = [9;10;11;12];  % comp 2, rank 1
Eig.Vectors(5:8, 4) = [13;14;15;16]; % comp 2, rank 2
Eig.Values    = [1;2;1;2];
Eig.Component = [1;1;2;2];
Eig.CompRank  = [1;2;1;2];

[Grid, K, Info] = view_eigenmodes('BuildPairedGrid', Eig);
assert(K == 2, 'K must equal the max within-component rank.');
assert(isequal(size(Grid), [nV, 2]), 'Grid must be [nV x K].');
% Column 1 = comp1-rank1 (verts 1:4) + comp2-rank1 (verts 5:8)
assert(isequal(Grid(:,1), [1;2;3;4;9;10;11;12]), 'Paired column 1 wrong.');
% Column 2 = comp1-rank2 + comp2-rank2
assert(isequal(Grid(:,2), [5;6;7;8;13;14;15;16]), 'Paired column 2 wrong.');
assert(Info.K == 2 && numel(Info.Values) == 4, 'Info fields wrong.');

% Single-component fallback (no metadata): each column maps to its own rank
Eig2 = struct('Vectors', magic(4), 'Values', (1:4)');
[Grid2, K2] = view_eigenmodes('BuildPairedGrid', Eig2);
assert(K2 == 4 && isequal(Grid2, magic(4)), 'Single-component grid must pass modes through.');

disp('ALL TESTS PASSED');
end
