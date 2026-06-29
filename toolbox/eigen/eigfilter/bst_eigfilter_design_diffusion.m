function out = bst_eigfilter_design_diffusion(params)
% BST_EIGFILTER_DESIGN_DIFFUSION: Physical-time heat kernel g(l,t) = exp(-tau*|t|*l/lmax), with t the
% PHYSICAL time axis (seconds). The time-aware extension of the static heat kernel: tau is a diffusivity
% acting over real time, so the abstract heat scale becomes a physical diffusion (the tau->seconds bridge).
% Domain 'ts' (eigen-time, g(lambda,t)), separable.
%
% USAGE:  g = bst_eigfilter_design_diffusion(struct('tau',1,'lmax',Lmax))   -> handle g(l,t) -> [K x N]
%         m = bst_eigfilter_design_diffusion('meta')                        -> metadata
%
% Authors: Diellor Basha, 2026

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% =============================================================================@

if nargin >= 1 && ischar(params) && strcmpi(params, 'meta')
    out = struct('name','diffusion', 'display','Diffusion (physical-time heat)', ...
        'params', struct('tau', struct('default',1,'range',[0 Inf])), ...
        'domain','ts', 'separable',true, 'bandpass',false, 'priorAdmissible',false);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'tau')  || isempty(params.tau);  params.tau  = 1;  end
if ~isfield(params,'lmax') || isempty(params.lmax); params.lmax = 1;  end
tau = params.tau;  sc = 1/max(params.lmax, eps);
out = @(l,t) exp(-tau * sc * double(l(:)) .* abs(double(t(:))'));   % [K x N]
end
