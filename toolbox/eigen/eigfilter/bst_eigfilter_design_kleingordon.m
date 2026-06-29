function out = bst_eigfilter_design_kleingordon(params)
% BST_EIGFILTER_DESIGN_KLEINGORDON: Klein-Gordon (massive wave) kernel
%     g(l,t) = cos(alpha*t*acos(1 - 2*(l+mu)/lmax)).
% The wave kernel with a mass term mu added to the eigenvalue: a propagation with a low-frequency cutoff
% / minimal oscillation rate set by mu (a dispersive wave that does not propagate the lowest modes
% statically). Domain 'ts' (eigen-time), NON-separable. t = physical time (s).
%
% USAGE:  g = bst_eigfilter_design_kleingordon(struct('alpha',1,'mu',0.1,'lmax',Lmax)) -> g(l,t) -> [K x N]
%         m = bst_eigfilter_design_kleingordon('meta')                                 -> metadata
%
% Authors: Diellor Basha, 2026

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% Copyright (c) University of Southern California & McGill University
% Distributed under the terms of the GNU General Public License (GPLv3).
% =============================================================================@

if nargin >= 1 && ischar(params) && strcmpi(params, 'meta')
    out = struct('name','kleingordon', 'display','Klein-Gordon (massive wave, mass mu)', ...
        'params', struct('alpha', struct('default',1,'range',[0 Inf]), ...
                         'mu',    struct('default',0.1,'range',[0 Inf])), ...
        'domain','ts', 'separable',false, 'bandpass',false, 'priorAdmissible',false);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'alpha') || isempty(params.alpha); params.alpha = 1;   end
if ~isfield(params,'mu')    || isempty(params.mu);    params.mu    = 0.1; end
if ~isfield(params,'lmax')  || isempty(params.lmax);  params.lmax  = 1;   end
a = params.alpha;  mu = params.mu;  sc = 2/max(params.lmax, eps);
out = @(l,t) cos(a * double(t(:))' .* acos(max(min(1 - sc*(double(l(:)) + mu), 1), -1)));   % [K x N]
end
