function S = bst_benchmark_sources(Surface, regime, varargin)
% BST_BENCHMARK_SOURCES: Generate ground-truth cortical sources for benchmarking.
%
% USAGE:  S = bst_benchmark_sources(Surface, regime, 'Param', value, ...)
%
% INPUTS:
%   Surface : struct with .Vertices [nV x 3] (metres) and .VertConn [nV x nV] sparse
%   regime  : 'focal' | 'patch' | 'distributed'
% OPTIONS:
%   'Seed'       : RNG seed for random seed-vertex selection (default 1)
%   'SeedVertex' : explicit seed vertex index (default: random, reproducible via Seed)
%   'Radius'     : patch radius in graph hops (default 2)               [regime 'patch']
%   'Sigma'      : Gaussian spatial scale in metres (default 0.02)      [regime 'distributed']
%   'nTime'      : number of time samples (default 20)
%
% OUTPUT struct S:
%   .GT         [nV x 1]      ground-truth spatial amplitude map (peak normalized to 1)
%   .Sources    [nV x nTime]  GT source matrix = GT * timecourse
%   .SeedVertex scalar        seed vertex index (the GT "true location")
%   .Regime     char          the regime
%
% Authors: Diellor Basha, 2026

Seed = 1; SeedVertex = []; Radius = 2; Sigma = 0.02; nTime = 20;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'seed',       Seed       = varargin{i+1};
        case 'seedvertex', SeedVertex = varargin{i+1};
        case 'radius',     Radius     = varargin{i+1};
        case 'sigma',      Sigma      = varargin{i+1};
        case 'ntime',      nTime      = varargin{i+1};
    end
end
V  = double(Surface.Vertices);
nV = size(V, 1);

% Seed vertex (reproducible)
if isempty(SeedVertex)
    rng(Seed);
    SeedVertex = randi(nV);
end

% Spatial profile
GT = zeros(nV, 1);
switch lower(regime)
    case 'focal'
        GT(SeedVertex) = 1;
    case 'patch'
        keep = graph_ball(Surface.VertConn, SeedVertex, Radius);
        GT(keep) = 1;
    case 'distributed'
        d = sqrt(sum((V - V(SeedVertex,:)).^2, 2));   % Euclidean distance (metres)
        GT = exp(-(d.^2) / (2 * Sigma^2));
    otherwise
        error('bst_benchmark_sources:UnknownRegime', 'Unknown regime: %s', regime);
end
GT = GT / max(GT);   % peak-normalize

% Time course: Gaussian-windowed burst, peak at the centre sample
t  = (1:nTime) - (nTime+1)/2;
tc = exp(-(t.^2) / (2 * (nTime/6)^2));

S = struct();
S.GT         = GT;
S.Sources    = GT * tc;            % [nV x nTime]
S.SeedVertex = SeedVertex;
S.Regime     = lower(regime);
end

function keep = graph_ball(VertConn, seed, radius)
% Vertices within `radius` hops of seed (BFS), inclusive.
keep = false(size(VertConn, 1), 1);
keep(seed) = true;
frontier = seed;
for h = 1:radius
    nbr = any(VertConn(:, frontier) ~= 0, 2);
    new = nbr & ~keep;
    keep = keep | nbr;
    frontier = find(new);
    if isempty(frontier); break; end
end
keep = find(keep);
end
