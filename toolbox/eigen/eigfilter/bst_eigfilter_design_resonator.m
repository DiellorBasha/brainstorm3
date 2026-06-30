function out = bst_eigfilter_design_resonator(params)
% BST_EIGFILTER_DESIGN_RESONATOR: joint-spectral damped harmonic oscillator (Lorentzian), peak at f0 Hz
% with quality Q. Hermitian in signed frequency -> real decaying-oscillation atom. lambda-independent.
% Authors: Diellor Basha, 2026
if nargin >= 1 && ischar(params) && strcmpi(params, 'meta')
    out = struct('name','resonator', 'display','Resonator (f0, Q)', ...
        'params', struct('f0',struct('default',10,'range',[0 Inf]), ...
                         'Q', struct('default',6, 'range',[0 Inf])), ...
        'domain','js', 'separable',true, 'bandpass',true, 'priorAdmissible',false);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'f0') || isempty(params.f0); params.f0 = 10; end
if ~isfield(params,'Q')  || isempty(params.Q);  params.Q  = 6;  end
out = @(l,w) i_reson(l, w, max(params.f0,eps), max(params.Q,eps));   % f0>0: f0=0 makes H=0/0 (NaN) at DC
end
function G = i_reson(l, w, f0, Q)
    ws = i_signed(w);
    H  = (f0.^2) ./ (f0.^2 - ws.^2 + 1i.*ws.*(f0/Q));                         % [1 x N] Hermitian
    G  = ones(numel(l),1) * H;                                               % [K x N], lambda-independent
end
function ws = i_signed(w)
    w = double(w(:)');  Fs = numel(w)*(w(2)-w(1));  ws = w - Fs.*(w >= Fs/2);
end
