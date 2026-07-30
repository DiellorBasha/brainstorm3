function out = omega_regime_screen(subShort, moduleDs, mnFile)
% OMEGA_REGIME_SCREEN  Standing-vs-traveling regime screen on the DIRAC EIGENBASIS flow, full 600 s.
%
% Cheap linear route: the Laplace-Beltrami coefficients of the divergence flow field come straight
% from sensors via the precomputed eigenbasis operator  c(lambda,t) = coeffOperator_div * b(t)
% (flow.divergence(ctx).coeffOperator, [Ks x C]) -- NO per-frame source reconstruction. Stream the
% full recording, resample to 150 Hz, demean, take the joint (lambda, omega) spectrum
% (filters.jspectrum) and sweep the traveling-wave ridge over phase speeds (filters.travwave).
% A traveling wave concentrates energy on omega = c*sqrt(lambda) (E(c) peaked); a standing
% oscillation is a horizontal band across scales (E(c) flat, peakedness ~2).
%
% Author: Diellor Basha, 2026

    addpath('/Users/diellorbasha/workspace/research/code/nxr-cortical-flow-matlab');
    ctx = flow.context(moduleDs, 400);
    P   = flow.divergence(ctx).coeffOperator;        % [Ks x C]  LBO coeffs of div, from sensors
    lam = ctx.lbo.Lambda(:);

    % --- stream b -> resample 150 Hz -> c = P*b (coefficient stream is small) ---
    R  = in_bst_results(mnFile, 0);
    gc = R.GoodChannel;
    sF = in_bst_data(R.DataFile, 'F');  sFile = sF.F;
    CM = in_bst_channel(bst_get('ChannelFileForStudy', R.DataFile));
    fs = sFile.prop.sfreq;  t0 = sFile.prop.times(1);  t1 = sFile.prop.times(2);
    fs2 = 150;  [pp, qq] = rat(fs2 / fs);
    C = [];  blk = 30;
    for ts = t0:blk:t1
        te = min(ts+blk, t1);  if te - ts < 1, break; end
        F  = panel_record('ReadRawBlock', sFile, CM, 1, [ts te], 0, 1, 'all', 1);
        Fr = resample(F(gc,:)', pp, qq)';            % [C x nT150]
        C  = [C, P*Fr];                              %#ok<AGROW>  [Ks x nT150]
    end
    C = C - mean(C, 2);
    fprintf('[%s] regime: streamed %d coeff-samples @%d Hz, Ks=%d\n', subShort, size(C,2), fs2, size(C,1));

    % --- joint spectrum + traveling-wave speed sweep ---
    [chat, f] = filters.jspectrum(C, fs2);
    jp = abs(chat).^2;
    speeds = linspace(0.05, 1.0, 60);  E = zeros(size(speeds));
    for i = 1:numel(speeds)
        G = filters.travwave(lam, f, speeds(i), 1.5);
        E(i) = sum(jp(:).*G(:)) / sum(jp(:));
    end
    [~, ip] = max(E);
    out.peakedness = max(E) / mean(E);
    out.speedPeak  = speeds(ip);
    out.speedAtEdge = (ip <= 2);                     % peak pinned to c->0 edge = standing
    out.jointPower = jp;  out.f = f;  out.lam = lam;  out.E = E;  out.speeds = speeds;
    fprintf('[%s] speed sweep peaks c=%.3f m/s | peakedness=%.2f (>~3 traveling, ~2 standing) | edge=%d\n', ...
        subShort, out.speedPeak, out.peakedness, out.speedAtEdge);
end
