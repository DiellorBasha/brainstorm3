function S = sa_smoothness(F, SulciMap)
% SA_SMOOTHNESS: Component 1 (data-free) edge metrics for the sign-ambiguity
% continuity test, partitioned by SulciMap.
%
% USAGE:  S = sa_smoothness(F, SulciMap)
%   F        : struct from sa_frames (uses .VertConn, .Nv, .ConnEig, .R).
%   SulciMap : [nV x 1] binary (1 = sulcal vertex).
%
% Per mesh edge (i,j), upper triangle:
%   dn(i,j) = acos(clamp(n_i . n_j))                      -- normal angular variation
%   df(i,j) = |wrap( arg(zF_i) - arg(zF_j) - arg(R_ij) )| -- Fiedler covariant variation
%             with R_ij = -K(i,j)/|K(i,j)|, zF the per-component Fiedler eigenvector.
% An edge is "sulcal" if either endpoint is sulcal.
%
% OUTPUT S: dn/df medians per partition, the decoupling ratios, and the fraction
% of total df^2 energy on edges touching a singularity neighborhood.

VertConn = F.VertConn;
Nv = F.Nv;
K  = F.ConnEig.ConnLaplacian;
nV = size(Nv, 1);

% Build the combined Fiedler eigenvector zF (per-component Fiedler columns have
% disjoint support, so summing them yields one [nV x 1] complex field).
ConnEig = F.ConnEig;
comps = unique(ConnEig.Component(:))';
zF = zeros(nV, 1);
for c = comps
    col = find(ConnEig.Component(:) == c & ConnEig.CompRank(:) == 1, 1);
    if isempty(col), continue; end
    zc = ConnEig.Vectors(:, col);
    zF(zc ~= 0) = zc(zc ~= 0);
end
assert(~isreal(zF), 'sa_smoothness:realEigenvectors', ...
    'ConnEig.Vectors must be complex (connection-Laplacian eigenvectors are complex-valued).');

% Edge list (upper triangle of the mesh adjacency).
[ii, jj] = find(triu(VertConn > 0, 1));

% dn: normal angular variation.
nd = sum(Nv(ii,:) .* Nv(jj,:), 2);
nd = max(min(nd, 1), -1);
dn = acos(nd);

% df: covariant phase increment using the connection's own transport arg(R_ij).
lin = sub2ind(size(K), ii, jj);
Kij = full(K(lin));
Rij = -Kij ./ max(abs(Kij), eps);            % transport rotation (unit complex)
inc = angle(zF(ii)) - angle(zF(jj)) - angle(Rij);
df  = abs(atan2(sin(inc), cos(inc)));        % wrap to [0, pi], magnitude
% Undefined where either endpoint has no Fiedler support.
defined = (zF(ii) ~= 0) & (zF(jj) ~= 0) & (abs(Kij) > eps);

% Partition edges by SulciMap (sulcal if either endpoint is sulcal).
edgeSulcal = (SulciMap(ii) > 0) | (SulciMap(jj) > 0);

S = struct();
S.dn_sulci_median = median(dn(edgeSulcal));
S.dn_crown_median = median(dn(~edgeSulcal));
S.df_sulci_median = median(df(defined & edgeSulcal));
S.df_crown_median = median(df(defined & ~edgeSulcal));
S.dn_sulci_crown_ratio = S.dn_sulci_median / max(S.dn_crown_median, eps);

% Singularity-energy concentration: fraction of total df^2 on edges within the
% 5-ring of any singularity vertex.
sing = F.R.Singularities(:);
near = false(nV, 1);
if ~isempty(sing)
    A = double(VertConn > 0);
    reach = sparse(sing, 1, 1, nV, 1);
    for r = 1:5, reach = double((A * reach + reach) > 0); end
    near = reach > 0;
end
edgeNearSing = near(ii) | near(jj);
df2 = df(defined).^2;
sel = edgeNearSing(defined);
S.singEnergyFrac = sum(df2(sel)) / max(sum(df2), eps);

% Carry the raw edge data for figures.
S.edges = [ii jj];
S.dn = dn; S.df = df; S.defined = defined; S.edgeSulcal = edgeSulcal;
end
