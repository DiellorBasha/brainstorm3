function M = bst_benchmark_metrics(gt, est, GridLoc, seedVertex)
% BST_BENCHMARK_METRICS: Compare an estimated source map against ground truth.
%
% USAGE:  M = bst_benchmark_metrics(gt, est, GridLoc, seedVertex)
%
% INPUTS:
%   gt         [nV x 1]   ground-truth amplitude map
%   est        [nV x 1]   estimated amplitude map (same vertices)
%   GridLoc    [nV x 3]   vertex positions in metres
%   seedVertex scalar     GT "true location" vertex index
%
% OUTPUT struct M: .LocError (mm) .AUC .NRMSE .Correlation .SpatialDispersion (mm)
%
% Authors: Diellor Basha, 2026
gt  = double(gt(:));  est = double(est(:));
ag  = abs(gt);        ae  = abs(est);
mm  = 1000;

% Localization error: GT seed location vs estimate peak location
[~, iPk] = max(ae);
M.LocError = norm(GridLoc(iPk,:) - GridLoc(seedVertex,:)) * mm;

% Spatial dispersion about the true location, weighted by estimate power
d = sqrt(sum((GridLoc - GridLoc(seedVertex,:)).^2, 2)) * mm;
w = ae.^2;
if sum(w) > 0
    M.SpatialDispersion = sqrt(sum(d.^2 .* w) / sum(w));
else
    M.SpatialDispersion = 0;
end

% NRMSE (norm-matched, GBF definition)
if norm(est) > 0
    est_s = est * (norm(gt) / norm(est));
else
    est_s = est;
end
rng_gt = max(gt) - min(gt);
if rng_gt > 0
    M.NRMSE = sqrt(mean((gt - est_s).^2)) / rng_gt;
else
    M.NRMSE = sqrt(mean((gt - est_s).^2));
end

% Correlation (Pearson, inline — no Statistics Toolbox required)
N = numel(gt);
sg = std(gt); se = std(est);
if sg > 0 && se > 0
    M.Correlation = sum((gt - mean(gt)) .* (est - mean(est))) / ((N-1) * sg * se);
else
    M.Correlation = 0;
end

% AUC: detect GT-active vertices (|gt| > 50% of peak) from estimate magnitude
labels = ag > 0.5 * max(ag);
M.AUC  = local_auc(labels, ae);
end

function auc = local_auc(labels, scores)
% Rank-based ROC AUC (Mann-Whitney). labels logical, scores >= 0.
pos = scores(labels); neg = scores(~labels);
nP = numel(pos); nN = numel(neg);
if nP == 0 || nN == 0; auc = 0.5; return; end
[~, order] = sort(scores); r = zeros(size(scores)); r(order) = 1:numel(scores);
auc = (sum(r(labels)) - nP*(nP+1)/2) / (nP*nN);
end
