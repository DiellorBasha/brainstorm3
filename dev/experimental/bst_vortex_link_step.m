function match = bst_vortex_link_step(prevHeads, curCores, maxJump)
% BST_VORTEX_LINK_STEP  One-frame greedy match of track heads to current cores.
% match(h) = index into curCores matched to prevHeads(h) (0 = none). A pair is
% eligible only if same chirality SIGN and same hemisphere and within maxJump (m).
%   prevHeads, curCores : struct arrays with .pos (1x3) .charge .hemi
% Author: Diellor Basha, 2026
    nH = numel(prevHeads);  nC = numel(curCores);
    match = zeros(1, nH);
    if nH == 0 || nC == 0, return; end
    pairs = zeros(0,3);                       % [h, c, dist]
    for h = 1:nH
        for c = 1:nC
            if sign(prevHeads(h).charge) ~= sign(curCores(c).charge), continue; end
            if prevHeads(h).hemi ~= curCores(c).hemi, continue; end
            d = norm(prevHeads(h).pos - curCores(c).pos);
            if d <= maxJump, pairs(end+1,:) = [h c d]; end %#ok<AGROW>
        end
    end
    if isempty(pairs), return; end
    pairs = sortrows(pairs, 3);
    usedC = false(1, nC);
    for p = 1:size(pairs,1)
        h = pairs(p,1);  c = pairs(p,2);
        if match(h) ~= 0 || usedC(c), continue; end
        match(h) = c;  usedC(c) = true;
    end
end
