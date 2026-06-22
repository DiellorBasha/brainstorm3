function test_benchmark_sources_pure
% Verify the ground-truth generator: focal/patch/distributed regimes on a
% synthetic surface; shapes; seeding reproducibility; spatial profiles.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot, fullfile(repoRoot,'dev','benchmarks'));
if ~brainstorm('status'); brainstorm nogui; end

% Synthetic surface: 6x6 grid of vertices on a plane (deterministic), 4-neighbour conn
n = 6; [X,Y] = ndgrid(1:n, 1:n);
Vertices = [X(:)*0.01, Y(:)*0.01, zeros(n*n,1)];   % 10 mm spacing, in metres
nV = size(Vertices,1);
VertConn = sparse(nV, nV);
idx = reshape(1:nV, n, n);
for i = 1:n
    for j = 1:n
        if i<n; VertConn(idx(i,j), idx(i+1,j)) = 1; VertConn(idx(i+1,j), idx(i,j)) = 1; end
        if j<n; VertConn(idx(i,j), idx(i,j+1)) = 1; VertConn(idx(i,j+1), idx(i,j)) = 1; end
    end
end
Surface = struct('Vertices', Vertices, 'VertConn', VertConn);

% --- focal ---
S = bst_benchmark_sources(Surface, 'focal', 'nTime', 10, 'Seed', 1);
assert(isequal(size(S.GT), [nV 1]), 'GT map must be [nV x 1].');
assert(nnz(S.GT) == 1, 'focal GT must have exactly one active vertex.');
assert(S.GT(S.SeedVertex) == max(S.GT), 'seed vertex must be the active one.');
assert(isequal(size(S.Sources), [nV 10]), 'Sources must be [nV x nTime].');

% --- patch (radius 1 hop -> seed + 4 neighbours = 5 active for an interior seed) ---
P = bst_benchmark_sources(Surface, 'patch', 'Radius', 1, 'Seed', 1, 'SeedVertex', idx(3,3));
assert(P.GT(idx(3,3))>0 && P.GT(idx(2,3))>0 && P.GT(idx(4,3))>0 ...
    && P.GT(idx(3,2))>0 && P.GT(idx(3,4))>0, 'patch must cover seed + 1-hop neighbours.');
assert(P.GT(idx(1,1))==0, 'patch must not reach distant vertices at radius 1.');

% --- distributed (Gaussian falloff, monotone decreasing with distance) ---
D = bst_benchmark_sources(Surface, 'distributed', 'Sigma', 0.02, 'SeedVertex', idx(3,3), 'Seed', 1);
dPeak = D.GT(idx(3,3)); dNear = D.GT(idx(3,4)); dFar = D.GT(idx(1,1));
assert(dPeak >= dNear && dNear > dFar, 'distributed profile must fall off with distance.');
assert(abs(dPeak - max(D.GT)) < 1e-12, 'peak must be at the seed.');

% --- seeding reproducibility (same seed -> same random seed vertex) ---
A = bst_benchmark_sources(Surface, 'focal', 'Seed', 7);
B = bst_benchmark_sources(Surface, 'focal', 'Seed', 7);
assert(A.SeedVertex == B.SeedVertex, 'same Seed must give the same random seed vertex.');

disp('ALL TESTS PASSED');
end
