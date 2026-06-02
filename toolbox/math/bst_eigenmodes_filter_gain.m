function h = bst_eigenmodes_filter_gain(lambdas, FilterType, varargin)
% BST_EIGENMODES_FILTER_GAIN: Per-mode transfer-function gain for eigenmode filtering.
%
% USAGE:  h = bst_eigenmodes_filter_gain(lambdas, 'lowpass',  'CutoffMode', 50)
%         h = bst_eigenmodes_filter_gain(lambdas, 'bandpass', 'ModeRange', [20 80])
%         h = bst_eigenmodes_filter_gain(lambdas, 'heat',     'DiffusionTime', 0.01)
%         h = bst_eigenmodes_filter_gain(lambdas, 'tikhonov', 'RegBeta', 1)
%
% DESCRIPTION:
%     Returns the per-mode gain vector h(lambda_k) for a spatial-spectral filter
%     in the Laplace-Beltrami eigenmode basis. This is the single source of the
%     transfer functions shared by bst_eigenmodes_filter (vertex fields) and
%     process_eigenmodes_coeffsfilter (eigenmode coefficients).
%
%     Types: 'lowpass'/'highpass'/'bandpass' (mode-index cutoffs), 'heat'
%     (exp(-t*lambda)), 'inverse_heat' (exp(+t*lambda) clamped at MaxGain),
%     'tikhonov' (1/(1+beta*lambda)), 'custom' (user TransferFn handle).
%
%     The analytic lambda-kernels ('heat', 'inverse_heat', 'tikhonov') delegate
%     to the eigfilter library via bst_eigfilter_kernel / bst_eigfilter_evaluate;
%     the index masks ('lowpass'/'highpass'/'bandpass') and 'custom' are computed
%     in place here.
%
% INPUTS:
%     lambdas    : [K x 1] eigenvalues (K = number of modes).
%     FilterType : one of the types above.
%
% OPTIONS (name-value):
%     'CutoffMode'    (50)      mode index for low/high-pass
%     'ModeRange'     ([20 80]) [k1 k2] for band-pass
%     'DiffusionTime' (0.01)    t for heat / inverse_heat
%     'MaxGain'       (10)      clamp for inverse_heat
%     'RegBeta'       (1)       beta for tikhonov
%     'TransferFn'    ([])      function handle for custom
%
% OUTPUTS:
%     h : [K x 1] gain vector.
%
% SEE ALSO: bst_eigenmodes_filter, process_eigenmodes_coeffsfilter,
%           bst_eigfilter_kernel, bst_eigfilter_evaluate

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Diellor Basha, 2026

%% ===== PARSE OPTIONS =====
CutoffMode    = 50;
ModeRange     = [20, 80];
DiffusionTime = 0.01;
MaxGain       = 10;
RegBeta       = 1;
TransferFn    = [];
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'cutoffmode',    CutoffMode    = varargin{i+1};
        case 'moderange',     ModeRange     = varargin{i+1};
        case 'diffusiontime', DiffusionTime = varargin{i+1};
        case 'maxgain',       MaxGain       = varargin{i+1};
        case 'regbeta',       RegBeta       = varargin{i+1};
        case 'transferfn',    TransferFn    = varargin{i+1};
    end
end

lambdas = double(lambdas(:));
nModes  = numel(lambdas);
h = zeros(nModes, 1);

%% ===== TRANSFER FUNCTION =====
switch lower(FilterType)
    case 'lowpass'
        c = min(CutoffMode, nModes);
        h(1:c) = 1;
    case 'highpass'
        c = max(1, min(CutoffMode, nModes));
        h(c:end) = 1;
    case 'bandpass'
        k1 = max(1, ModeRange(1));
        k2 = min(nModes, ModeRange(2));
        h(k1:k2) = 1;
    case 'heat'
        if DiffusionTime <= 0
            error('DiffusionTime must be positive (got %g).', DiffusionTime);
        end
        g = bst_eigfilter_kernel('heat', struct('t', DiffusionTime));
        h = bst_eigfilter_evaluate(g, lambdas);
    case 'inverse_heat'
        if DiffusionTime <= 0
            error('DiffusionTime must be positive (got %g).', DiffusionTime);
        end
        g = bst_eigfilter_kernel('inverse_heat', struct('t', DiffusionTime, 'maxgain', MaxGain));
        h = bst_eigfilter_evaluate(g, lambdas);
    case 'tikhonov'
        if RegBeta < 0
            error('RegBeta must be non-negative (got %g).', RegBeta);
        end
        g = bst_eigfilter_kernel('tikhonov', struct('beta', RegBeta));
        h = bst_eigfilter_evaluate(g, lambdas);
    case 'custom'
        if isempty(TransferFn) || ~isa(TransferFn, 'function_handle')
            error('Custom filter requires a TransferFn option (function handle).');
        end
        h = TransferFn(lambdas);
        if numel(h) ~= nModes
            error('TransferFn must return a vector of length %d (got %d).', nModes, numel(h));
        end
        h = h(:);
    otherwise
        error('Unknown filter type: %s. Use lowpass, highpass, bandpass, heat, inverse_heat, tikhonov, or custom.', FilterType);
end
end
