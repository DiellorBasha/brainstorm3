function mk = bst_persistence_gate(mk, frac)
% BST_PERSISTENCE_GATE  Per-hemisphere persistence gate for singular-point markers.
% Keeps a core if it is global (Inf persistence) or its persistence >= frac * (max
% FINITE persistence within its OWN hemisphere). frac<=0 keeps everything.
%   mk : struct array with fields .hemi and .persistence
% Author: Diellor Basha, 2026
    if isempty(mk) || frac <= 0, return; end
    hh   = [mk.hemi];
    keep = false(1, numel(mk));
    for h = unique(hh)
        idx = find(hh == h);
        pr  = [mk(idx).persistence];
        mxf = max([pr(isfinite(pr)), 0]);
        if mxf <= 0
            keep(idx) = true;
        else
            keep(idx) = isinf(pr) | (pr >= frac * mxf);
        end
    end
    mk = mk(keep);
end
