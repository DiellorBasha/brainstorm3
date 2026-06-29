function out = bst_eigfilter_design_dampedwave(params)
% BST_EIGFILTER_DESIGN_DAMPEDWAVE: Damped-wave kernel g(l,t) = exp(-beta*|t|)*cos(alpha*t*acos(1-2*l/lmax)).
% The wave kernel with temporal damping beta: a seed propagates as a decaying wave (speed alpha, damping
% beta). This is the canonical DYNAMIC eigenwavelet - a traveling, fading wavelet over geometry (LBO) or
% geometry+connectome (LB-Connectome). Domain 'ts' (eigen-time), NON-separable. t = physical time (s).
% (GSPBox gsp_jtv_design_damped_wave analog; stability alpha < 2/Fs.)
%
% USAGE:  g = bst_eigfilter_design_dampedwave(struct('alpha',1,'beta',0.1,'lmax',Lmax)) -> g(l,t) -> [K x N]
%         m = bst_eigfilter_design_dampedwave('meta')                                   -> metadata
%
% Authors: Diellor Basha, 2026

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% Copyright (c) University of Southern California & McGill University
% Distributed under the terms of the GNU General Public License (GPLv3).
% =============================================================================@

if nargin >= 1 && ischar(params) && strcmpi(params, 'meta')
    out = struct('name','dampedwave', 'display','Damped wave (speed alpha, damping beta)', ...
        'params', struct('alpha', struct('default',1,'range',[0 Inf]), ...
                         'beta',  struct('default',0.1,'range',[0 Inf])), ...
        'domain','ts', 'separable',false, 'bandpass',false, 'priorAdmissible',false);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'alpha') || isempty(params.alpha); params.alpha = 1;   end
if ~isfield(params,'beta')  || isempty(params.beta);  params.beta  = 0.1; end
if ~isfield(params,'lmax')  || isempty(params.lmax);  params.lmax  = 1;   end
a = params.alpha;  b = params.beta;  sc = 2/max(params.lmax, eps);
out = @(l,t) exp(-b * abs(double(t(:)'))) .* cos(a * double(t(:))' .* acos(max(min(1 - sc*double(l(:)), 1), -1)));
end
