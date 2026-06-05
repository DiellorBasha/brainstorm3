function pairs = sa_sulcal_walls(Vtx, Nv, SulciMap, VertConn, opts)
% SA_SULCAL_WALLS: Detect facing sulcal-wall vertex pairs.
%
% USAGE:  pairs = sa_sulcal_walls(Vtx, Nv, SulciMap, VertConn, opts)
%
% A pair (i,j) qualifies when:
%   (1) both i and j are sulcal vertices (SulciMap>0),
%   (2) ||Vtx_i - Vtx_j|| <= opts.MaxDist,
%   (3) Nv_i . Nv_j < opts.NormalDot  (anti-aligned -> facing walls),
%   (4) i and j are NOT within each other's opts.Nring mesh neighborhood
%       (so they are across a fold, not neighbors on the same wall).
%
% opts fields: MaxDist (m, e.g. 0.003), NormalDot (e.g. -0.7), Nring (e.g. 3).
% OUTPUT: pairs [nPairs x 2], unique unordered (i<j), sorted.

if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts,'MaxDist'),   opts.MaxDist   = 0.003; end
if ~isfield(opts,'NormalDot'), opts.NormalDot = -0.7;  end
if ~isfield(opts,'Nring'),     opts.Nring     = 3;     end

% (1) Restrict to sulcal vertices.
iSulci = find(SulciMap > 0);
if isempty(iSulci), pairs = zeros(0,2); return; end
Vs = Vtx(iSulci, :);

% (2) Range search among sulcal vertices (3D proximity).
nb = rangesearch(Vs, Vs, opts.MaxDist);

% N-ring adjacency (global indices) for criterion (4).
A = double(VertConn > 0);
Aring = A; P = A;
for r = 2:opts.Nring
    P = double((P * A) > 0);
    Aring = double((Aring + P) > 0);
end

I = []; J = [];
for a = 1:numel(iSulci)
    gi = iSulci(a);
    for b = nb{a}(:)'
        gj = iSulci(b);
        if gj <= gi, continue; end                 % unordered, i<j, drop self
        if sum(Nv(gi,:) .* Nv(gj,:), 2) >= opts.NormalDot, continue; end   % (3)
        if Aring(gi, gj) ~= 0, continue; end        % (4) within N-ring
        I(end+1,1) = gi; J(end+1,1) = gj; %#ok<AGROW>
    end
end
pairs = sortrows([I J]);
end
