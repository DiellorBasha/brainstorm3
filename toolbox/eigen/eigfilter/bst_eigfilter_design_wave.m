function out = bst_eigfilter_design_wave(params)
% BST_EIGFILTER_DESIGN_WAVE: Wave-propagation kernel g(l,t) = cos(alpha*t*acos(1 - 2*l/lmax)), the
% impulse response of the wave equation on the eigenbasis. The dispersion acos(1-2l/lmax) makes each
% eigenmode oscillate at a rate set by its eigenvalue l, so a seed propagates as a wave with speed alpha
% (with the mm/Hz calibration, alpha maps to m/s). Domain 'ts' (eigen-time), NON-separable (the cortex
% scale l sets the temporal behaviour). t is the physical time axis (seconds).
%
% USAGE:  g = bst_eigfilter_design_wave(struct('alpha',1,'lmax',Lmax))   -> handle g(l,t) -> [K x N]
%         m = bst_eigfilter_design_wave('meta')                          -> metadata
%
% Authors: Diellor Basha, 2026

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% Copyright (c) University of Southern California & McGill University
% Distributed under the terms of the GNU General Public License (GPLv3).
% =============================================================================@

if nargin >= 1 && ischar(params) && strcmpi(params, 'meta')
    out = struct('name','wave', 'display','Wave propagation (speed alpha)', ...
        'params', struct('alpha', struct('default',1,'range',[0 Inf])), ...
        'domain','ts', 'separable',false, 'bandpass',false, 'priorAdmissible',false);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'alpha') || isempty(params.alpha); params.alpha = 1; end
if ~isfield(params,'lmax')  || isempty(params.lmax);  params.lmax  = 1; end
a = params.alpha;  sc = 2/max(params.lmax, eps);
out = @(l,t) cos(a * double(t(:))' .* acos(max(min(1 - sc*double(l(:)), 1), -1)));   % [K x N]
end
