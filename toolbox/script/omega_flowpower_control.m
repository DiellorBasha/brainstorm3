function out = omega_flowpower_control(subShort, moduleDs, mnFile)
% OMEGA_FLOWPOWER_CONTROL  Apples-to-apples flow!=power control for the Dirac EIGENBASIS vs the
% vanilla wMNE kernel, single-subject and UNSMOOTHED (unlike the group-averaged/ssmoothed maps).
%
% Streams the subject's alpha-band (8-13 Hz) sensor covariance Ca [C x C] once, then for BOTH
% kernels forms the per-vertex flow-POWER maps  p_div(v) = (Kdiv Ca Kdiv')_vv , p_curl likewise,
% and reports their spatial correlation. If flow were just power, div-power and curl-power would
% both track |J|^2 and hence each other; the basis-free wMNE run gave div-curl = 0.33 (decisive).
%
%   moduleDs : nxr module dataset name (for flow.context -> Dirac currentKernel + operators)
%   mnFile   : the subject's wMNE MN kernel result (source of the vanilla kernel + raw link)
%
% Author: Diellor Basha, 2026

    addpath('/Users/diellorbasha/workspace/research/code/nxr-cortical-flow-matlab');

    % --- alpha covariance by streaming the raw recording (verified ReadRawBlock pattern) ---
    R  = in_bst_results(mnFile, 0);
    gc = R.GoodChannel;
    sF = in_bst_data(R.DataFile, 'F');  sFile = sF.F;
    CM = in_bst_channel(bst_get('ChannelFileForStudy', R.DataFile));
    fs = sFile.prop.sfreq;  t0 = sFile.prop.times(1);  t1 = sFile.prop.times(2);
    [b,a] = butter(4, [8 13]/(fs/2), 'bandpass');

    Ca = zeros(numel(gc)); n = 0;  blk = 30;   % 30 s blocks
    for ts = t0:blk:t1
        te = min(ts+blk, t1);  if te - ts < 1, break; end
        F = panel_record('ReadRawBlock', sFile, CM, 1, [ts te], 0, 1, 'all', 1);
        F = F(gc, :);
        Fa = filtfilt(b, a, F')';                 % alpha band
        Ca = Ca + Fa*Fa';  n = n + size(Fa,2);
    end
    Ca = Ca / n;
    fprintf('[%s] alpha covariance streamed: %d chan, %d samples\n', subShort, numel(gc), n);

    % --- kernels: Dirac eigenbasis (module) + vanilla wMNE (the MN result) ---
    ctx    = flow.context(moduleDs, 400);
    Kdiv_d = flow.divergence(ctx).vertexOperator;         % [V x C]
    Kcurl_d= flow.curl(ctx).vertexOperator;
    Kmn    = R.ImagingKernel;                             % [3V x C] wMNE current
    Kdiv_m = differential.divergence(Kmn, ctx.S, ctx.fg);
    Kcurl_m= differential.curl(Kmn, ctx.S, ctx.fg);

    powmap = @(K) sum((K*Ca).*K, 2);                     % diag(K Ca K') per vertex
    pdiv_d = powmap(Kdiv_d);  pcurl_d = powmap(Kcurl_d);
    pdiv_m = powmap(Kdiv_m);  pcurl_m = powmap(Kcurl_m);

    out.dirac_divcurl = corr(pdiv_d, pcurl_d);
    out.mn_divcurl    = corr(pdiv_m, pcurl_m);
    out.cross_div     = corr(pdiv_d, pdiv_m);            % eigenbasis vs wMNE, same quantity
    out.cross_curl    = corr(pcurl_d, pcurl_m);
    fprintf('[%s] div-curl corr (unsmoothed, alpha): dirac=%.3f  wMNE=%.3f  (low = flow!=power)\n', ...
        subShort, out.dirac_divcurl, out.mn_divcurl);
    fprintf('[%s] eigenbasis-vs-wMNE concordance: div=%.3f  curl=%.3f\n', ...
        subShort, out.cross_div, out.cross_curl);
end
