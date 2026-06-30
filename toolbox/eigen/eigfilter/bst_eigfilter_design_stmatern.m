function out = bst_eigfilter_design_stmatern(params)
% BST_EIGFILTER_DESIGN_STMATERN: spatiotemporal Whittle-Matern spectral density g = (kappa^2 + lambda +
% (w/v)^2)^(-nu) -- the 1/f aperiodic background / a spatiotemporal prior. v (space-time speed, m/s) fixed.
% Authors: Diellor Basha, 2026
V_SPEED = 1.0;   % m/s, fixed space-time coupling
if nargin >= 1 && ischar(params) && strcmpi(params, 'meta')
    out = struct('name','stmatern', 'display','Spatiotemporal 1/f (Whittle-Matern)', ...
        'params', struct('kappa',struct('default',126,'range',[0 Inf]), ...
                         'nu',   struct('default',1.5,'range',[0 Inf])), ...
        'domain','js', 'separable',false, 'bandpass',false, 'priorAdmissible',true);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'kappa') || isempty(params.kappa); params.kappa = 126; end
if ~isfield(params,'nu')    || isempty(params.nu);    params.nu    = 1.5; end
out = @(l,w) i_stm(l, w, params.kappa, params.nu, V_SPEED);
end
function G = i_stm(l, w, kappa, nu, v)
    ws = i_signed(w);
    base = (kappa.^2 + double(l(:))) * ones(1,numel(ws)) + ones(numel(l),1) * (ws./v).^2;   % [K x N]
    G = base .^ (-nu);
end
function ws = i_signed(w)
    w = double(w(:)');  Fs = numel(w)*(w(2)-w(1));  ws = w - Fs.*(w >= Fs/2);
end
