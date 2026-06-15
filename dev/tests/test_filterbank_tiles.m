function test_filterbank_tiles()
% Pure-logic test for the spectrum tiling generator.
% Authors: Diellor Basha, 2026
    nFail = 0;
    base = struct('Kernel','mexhat', 'Params',struct(), ...
                  'Direction',[1 0 0], 'Chirality',0, 'Axis',[0 0 1], ...
                  'N',4, 'Spacing','geometric', 'LambdaRange',[1 256], 'Chiralities',0);

    T = bst_filterbank_tiles(base);
    nFail = nFail + chk('4 tiles, no chirality split', numel(T)==4);
    nFail = nFail + chk('each tile carries Kernel', all(arrayfun(@(t) strcmp(t.Kernel,'mexhat'), T)));

    % mexhat peaks at l = 1/t, so the per-tile t must span the range geometrically:
    % t_j = 1/center_j, centers geometric from 1 to 256 -> ratio 4 across 4 tiles.
    centers = arrayfun(@(t) 1./t.Params.t, T);
    ratios  = centers(2:end) ./ centers(1:end-1);
    nFail = nFail + chk('geometric centers (constant ratio)', max(abs(ratios - ratios(1))) < 1e-9);
    nFail = nFail + chk('span covers LambdaRange ends', abs(centers(1)-1)<1e-6 && abs(centers(end)-256)<1e-6);

    % chirality split doubles the bank into +1/-1 with matched scales
    base2 = base; base2.Chiralities = [1 -1];
    T2 = bst_filterbank_tiles(base2);
    nFail = nFail + chk('chirality split doubles count', numel(T2)==8);
    nFail = nFail + chk('both signs present', any([T2.Chirality]==1) && any([T2.Chirality]==-1));

    fprintf('\n==== test_filterbank_tiles: %d failed ====\n', nFail);
    if nFail > 0, error('test_filterbank_tiles FAILED'); end
end

function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
