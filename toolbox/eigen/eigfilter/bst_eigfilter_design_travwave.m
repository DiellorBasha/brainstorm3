function out = bst_eigfilter_design_travwave(params)
% BST_EIGFILTER_DESIGN_TRAVWAVE: joint-spectral traveling wave -- energy on the dispersion ridge
% f = c*sqrt(lambda)/(2*pi) Hz (phase speed c m/s, ridge width Hz). Non-separable; even in freq -> real.
% Authors: Diellor Basha, 2026
if nargin >= 1 && ischar(params) && strcmpi(params, 'meta')
    out = struct('name','travwave', 'display','Traveling wave (speed c)', ...
        'params', struct('c',struct('default',1,'range',[0 Inf]), ...
                         'width',struct('default',2,'range',[0 Inf])), ...
        'domain','js', 'separable',false, 'bandpass',true, 'priorAdmissible',false);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'c')     || isempty(params.c);     params.c = 1;     end
if ~isfield(params,'width') || isempty(params.width); params.width = 2; end
out = @(l,w) i_trav(l, w, params.c, max(params.width,eps));
end
function G = i_trav(l, w, c, sg)
    ws = i_signed(w);
    fr = c .* sqrt(double(l(:))) ./ (2*pi);                                   % [K x 1] ridge freq (Hz)
    G  = exp(-((abs(ws) - fr).^2) / (2*sg^2));                                % [K x N], even in ws
end
function ws = i_signed(w)
    w = double(w(:)');  Fs = numel(w)*(w(2)-w(1));  ws = w - Fs.*(w >= Fs/2);
end
