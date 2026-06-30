function out = bst_eigfilter_controls(action, kernelName, varargin)
% BST_EIGFILTER_CONTROLS: shared per-kernel GUI control spec for the atom designer AND panel_bst_dynamics.
%   S  = bst_eigfilter_controls('Sliders', name, bounds)        -> 1x3 struct(label,unit,lo,hi,def,fmt)
%   kp = bst_eigfilter_controls('ToKernel', name, vals, lmax)   -> kernel param struct (vals=[s1 s2 s3])
% Single source of truth for which params a kernel exposes, their ranges, and how slider values map to
% kernel params. Spectrum-dependent ranges (Scale/Rate) come from bounds (scaleMinMM/scaleMaxMM/rate*).
% Authors: Diellor Basha, 2026
switch lower(action)
    case 'sliders',  out = i_sliders(lower(kernelName), varargin{1});
    case 'tokernel', out = i_tokernel(lower(kernelName), varargin{1}, varargin{2});
    otherwise, error('bst_eigfilter_controls: unknown action ''%s''.', action);
end
end

function S = i_sliders(k, b)
    e   = i_row('',[],[],[],'');
    scD = round((b.scaleMinMM + b.scaleMaxMM)/2);               % spectrum-derived defaults
    rtD = round(((b.scaleMinMM + b.scaleMaxMM)/2)^2);
    sc  = i_row('Scale (mm)', b.scaleMinMM, b.scaleMaxMM, scD, '%.0f');
    switch k
        case 'diffusion',   S = [i_row('Rate (mm^2/s)', b.rateMinMM2, b.rateMaxMM2, rtD, '%.0f'), e, e];
        case {'heat','mexhat','diffgauss','flat','ideal','inverse_heat','log','matern','power','tikhonov'}
            S = [sc, e, e];
        case 'wave',        S = [e, i_row('Speed (m/s)',0.1,10,1,'%.2g'), e];
        case 'kleingordon', S = [e, i_row('Speed (m/s)',0.1,10,1,'%.2g'), e];
        case 'dampedwave',  S = [e, i_row('Speed (m/s)',0.1,10,1,'%.2g'), i_row('Decay (s)',0.05,2,0.5,'%.2g')];
        case 'gabor',       S = [sc, i_row('Freq (Hz)',0,50,10,'%.1f'), i_row('BW (Hz)',0.5,20,2,'%.1f')];
        case 'travwave',    S = [i_row('Speed (m/s)',0.05,3,1,'%.2g'), i_row('RidgeW (Hz)',0.5,20,2,'%.1f'), e];
        case 'resonator',   S = [i_row('Freq (Hz)',0,50,10,'%.1f'), i_row('Q',1,30,6,'%.1f'), e];
        case 'stmatern',    S = [i_row('Corr (mm)', b.scaleMinMM, b.scaleMaxMM, scD, '%.0f'), i_row('nu',0.5,4,1.5,'%.1f'), e];
        otherwise,          S = [sc, e, e];
    end
end

function kp = i_tokernel(k, v, lmax)
    s1 = v(1); s2 = v(2); s3 = v(3);
    kp = struct('lmax', lmax);
    switch k
        case 'wave',        kp.alpha = s2 * sqrt(lmax) / 2;
        case 'kleingordon', kp.alpha = s2 * sqrt(lmax) / 2;  kp.mu = 0.1*lmax;
        case 'dampedwave',  kp.alpha = s2 * sqrt(lmax) / 2;  kp.beta = 1/max(s3,eps);
        case 'diffusion',   kp.tau   = max((s1/1e6) * lmax, eps);
        case 'heat',        lamS = (2*pi/(s1/1000))^2;  kp.t = log(2)/max(lamS,eps);
        case 'mexhat',      lamS = (2*pi/(s1/1000))^2;  kp.t = 1/max(lamS,eps);
        case 'gabor',       kp.k0 = 2*pi/max(s1/1000,eps);  kp.f0 = s2;  kp.sf = max(s3,eps);
        case 'travwave',    kp.c  = s1;  kp.width = max(s2,eps);
        case 'resonator',   kp.f0 = s1;  kp.Q = max(s2,eps);
        case 'stmatern',    kp.kappa = 2*pi/max(s1/1000,eps);  kp.nu = max(s2,eps);
    end
end

function r = i_row(label, lo, hi, def, fmt)
    r = struct('label',label, 'lo',lo, 'hi',hi, 'def',def, 'fmt',fmt);
end
