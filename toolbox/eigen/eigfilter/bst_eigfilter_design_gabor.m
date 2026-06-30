function out = bst_eigfilter_design_gabor(params)
% BST_EIGFILTER_DESIGN_GABOR: joint-spectral Gabor packet g(l,w) localised at spatial wavenumber k0
% (rad/m) and temporal frequency f0 (Hz, bandwidth sf). Hermitian in temporal frequency -> real atom.
% Authors: Diellor Basha, 2026
if nargin >= 1 && ischar(params) && strcmpi(params, 'meta')
    out = struct('name','gabor', 'display','Gabor packet (scale x freq)', ...
        'params', struct('f0',struct('default',10,'range',[0 Inf]), ...
                         'k0',struct('default',209,'range',[0 Inf]), ...
                         'sf',struct('default',2,'range',[0 Inf])), ...
        'domain','js', 'separable',true, 'bandpass',true, 'priorAdmissible',false);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'f0') || isempty(params.f0); params.f0 = 10;  end
if ~isfield(params,'k0') || isempty(params.k0); params.k0 = 209; end
if ~isfield(params,'sf') || isempty(params.sf); params.sf = 2;   end
out = @(l,w) i_gabor(l, w, params.k0, max(params.k0/3,eps), params.f0, max(params.sf,eps));
end
function G = i_gabor(l, w, k0, sk, f0, sf)
    ws = i_signed(w);
    sl = exp(-((sqrt(double(l(:)))-k0).^2)/(2*sk^2));                          % [K x 1]
    gt = exp(-((ws-f0).^2)/(2*sf^2)) + exp(-((ws+f0).^2)/(2*sf^2));            % [1 x N] Hermitian
    G  = sl * gt;
end
function ws = i_signed(w)
    w = double(w(:)');  Fs = numel(w)*(w(2)-w(1));  ws = w - Fs.*(w >= Fs/2);
end
