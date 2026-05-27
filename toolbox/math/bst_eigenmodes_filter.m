function Filtered = bst_eigenmodes_filter(Eigenmodes, Data, MassMatrix, FilterType, varargin)
% BST_EIGENMODES_FILTER: Spatial spectral filtering of cortical data via eigenmodes.
%
% USAGE:  Filtered = bst_eigenmodes_filter(Eigenmodes, Data, MassMatrix, 'lowpass',  'CutoffMode', 50)
%         Filtered = bst_eigenmodes_filter(Eigenmodes, Data, MassMatrix, 'bandpass', 'ModeRange', [20 80])
%         Filtered = bst_eigenmodes_filter(Eigenmodes, Data, MassMatrix, 'heat',     'DiffusionTime', 0.01)
%
% INPUTS:
%     Eigenmodes : struct from in_tess_eigenmodes (uses .Vectors, .Values)
%     Data       : [nVertices x nTime] scalar field(s) on the mesh
%     MassMatrix : [nVertices x nVertices] sparse mass matrix (from tess_laplacian)
%     FilterType : 'lowpass','highpass','bandpass','heat','inverse_heat','custom'
%
% OPTIONS:
%     CutoffMode, ModeRange, DiffusionTime, MaxGain, TransferFn (see code)
%
% SEE ALSO: bst_eigenmodes_project, in_tess_eigenmodes, tess_eigenmodes, tess_laplacian
%
% DESCRIPTION:
%     Filters a scalar field on the cortical surface in the spectral domain
%     defined by the Laplace-Beltrami eigenmodes. The data is projected onto
%     the eigenmode basis, multiplied by a filter transfer function, and
%     reconstructed. This is the surface analogue of Fourier-domain filtering.
%
%     Available filter types:
%
%     'lowpass'  — Keep modes 1..CutoffMode, zero out higher modes.
%                  Produces a spatially smooth version of the data.
%
%     'highpass' — Keep modes CutoffMode..end, zero out lower modes.
%                  Isolates fine spatial detail.
%
%     'bandpass' — Keep modes in ModeRange [k1, k2], zero out the rest.
%                  Isolates a spatial frequency band.
%
%     'heat'     — Heat kernel filter: h(lambda_k) = exp(-t * lambda_k).
%                  Smooth, continuous low-pass with diffusion time parameter t.
%                  Small t is approximately identity; large t gives stronger smoothing.
%                  This is equivalent to running heat diffusion on the surface
%                  for time t. Commonly used in neuroimaging for cortical
%                  smoothing that respects surface geometry.
%
%     'inverse_heat' — Inverse heat kernel: h(lambda_k) = exp(+t * lambda_k),
%                  clamped to avoid amplifying high-frequency noise beyond a
%                  MaxGain factor. Used for spatial sharpening / deblurring.
%
%     'custom'   — User-supplied transfer function h(lambda_k) via the
%                  'TransferFn' option. Must be a function handle mapping
%                  eigenvalue vector to gain vector.
%
% MATHEMATICAL NOTES:
%     The filtering operation in matrix form is:
%
%         u_filtered = Phi * diag(h) * Phi' * M * u
%
%     where h is the vector of filter gains h(lambda_k) and M is the mass
%     matrix. For the heat kernel, h_k = exp(-t * lambda_k), which naturally
%     suppresses high spatial frequencies (large lambda_k) while preserving
%     low-frequency structure.

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

%% ===== PARSE INPUTS =====
CutoffMode    = 50;
ModeRange     = [20, 80];
DiffusionTime = 0.01;
MaxGain       = 10;
TransferFn    = [];
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'cutoffmode',    CutoffMode    = varargin{i+1};
        case 'moderange',     ModeRange     = varargin{i+1};
        case 'diffusiontime', DiffusionTime = varargin{i+1};
        case 'maxgain',       MaxGain       = varargin{i+1};
        case 'transferfn',    TransferFn    = varargin{i+1};
    end
end

Phi     = double(Eigenmodes.Vectors);   % [nV x nModes]
lambdas = double(Eigenmodes.Values(:));  % [nModes x 1]
nV      = size(Phi, 1);
nModes  = size(Phi, 2);

%% ===== VALIDATE =====
Data = double(Data);
if size(Data, 1) ~= nV
    error('Data has %d rows but eigenmodes have %d vertices.', size(Data, 1), nV);
end
if (size(MassMatrix, 1) ~= nV) || (size(MassMatrix, 2) ~= nV)
    error('MassMatrix must be %dx%d.', nV, nV);
end

%% ===== BUILD TRANSFER FUNCTION =====
h = zeros(nModes, 1);
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
        h = exp(-DiffusionTime * lambdas);
    case 'inverse_heat'
        if DiffusionTime <= 0
            error('DiffusionTime must be positive (got %g).', DiffusionTime);
        end
        h = min(exp(DiffusionTime * lambdas), MaxGain);
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
        error('Unknown filter type: %s. Use lowpass, highpass, bandpass, heat, inverse_heat, or custom.', FilterType);
end

%% ===== APPLY: u_filtered = Phi * (h .* (Phi' * M * u)) =====
Coeffs   = Phi' * (MassMatrix * Data);
Coeffs   = bsxfun(@times, h, Coeffs);
Filtered = Phi * Coeffs;
end
