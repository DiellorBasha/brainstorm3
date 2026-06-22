function test_filterbank_tiles()
% Pure-logic test for the spectrum tiling module: it consumes ONE designed wavelet
% plus tiling options and returns a bank, without ever deciding the wavelet itself.
% Authors: Diellor Basha, 2026
    nFail = 0;
    wavelet = struct('Kernel','mexhat', 'Params',struct('t',0.0123), ...
                     'Direction',[1 0 0], 'Chirality',0, 'Axis',[0 0 1]);
    opts = struct('N',4, 'Spacing','geometric', 'LambdaRange',[1 256], 'Chiralities',[]);

    T = bst_filterbank_tiles(wavelet, opts);
    nFail = nFail + chk('4 tiles, no chirality split', numel(T)==4);
    nFail = nFail + chk('each tile carries Kernel', all(arrayfun(@(t) strcmp(t.Kernel,'mexhat'), T)));

    % mexhat peaks at l = 1/t, so the per-tile t must span the range geometrically:
    % t_j = 1/center_j, centers geometric from 1 to 256 -> ratio 4 across 4 tiles.
    centers = arrayfun(@(t) 1./t.Params.t, T);
    ratios  = centers(2:end) ./ centers(1:end-1);
    nFail = nFail + chk('geometric centers (constant ratio)', max(abs(ratios - ratios(1))) < 1e-9);
    nFail = nFail + chk('span covers LambdaRange ends', abs(centers(1)-1)<1e-6 && abs(centers(end)-256)<1e-6);

    % N==1 keeps the wavelet's designed scale param (does NOT overwrite from LambdaRange)
    opts1 = opts; opts1.N = 1;
    T1 = bst_filterbank_tiles(wavelet, opts1);
    nFail = nFail + chk('N=1 single tile', numel(T1)==1);
    nFail = nFail + chk('N=1 keeps designed scale param', abs(T1(1).Params.t - 0.0123) < 1e-12);

    % chirality split doubles the bank into +1/-1 with matched scales
    opts2 = opts; opts2.Chiralities = [1 -1];
    T2 = bst_filterbank_tiles(wavelet, opts2);
    nFail = nFail + chk('chirality split doubles count', numel(T2)==8);
    nFail = nFail + chk('both signs present', any([T2.Chirality]==1) && any([T2.Chirality]==-1));

    % default opts (no second arg) => single wavelet, kept as-is
    T3 = bst_filterbank_tiles(wavelet);
    nFail = nFail + chk('no opts => 1 tile, scale kept', numel(T3)==1 && abs(T3(1).Params.t-0.0123)<1e-12);

    fprintf('\n==== test_filterbank_tiles: %d failed ====\n', nFail);
    if nFail > 0, error('test_filterbank_tiles FAILED'); end
end

function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
