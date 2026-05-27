function Filtered = tess_eigenmodes_filter(SurfaceFile, Data, FilterType, varargin)
% TESS_EIGENMODES_FILTER: Spatial frequency filtering of cortical data via eigenmodes.
%
% USAGE:  Filtered = tess_eigenmodes_filter(SurfaceFile, Data, 'lowpass',  'CutoffMode', 50)
%         Filtered = tess_eigenmodes_filter(SurfaceFile, Data, 'bandpass', 'ModeRange', [20 80])
%         Filtered = tess_eigenmodes_filter(SurfaceFile, Data, 'highpass', 'CutoffMode', 30)
%         Filtered = tess_eigenmodes_filter(SurfaceFile, Data, 'heat',     'DiffusionTime', 0.01)
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
%                  Small t ≈ identity; large t → stronger smoothing.
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
% INPUTS:
%     SurfaceFile : Relative or absolute path to the Brainstorm surface file
%                   with precomputed Eigenmodes.
%     Data        : [nVertices x nTime] scalar field(s) on the mesh.
%     FilterType  : One of 'lowpass', 'highpass', 'bandpass', 'heat',
%                   'inverse_heat', 'custom'.
%
% OPTIONS:
%     CutoffMode    : Mode index for lowpass/highpass cutoff (default: 50).
%     ModeRange     : [k1 k2] for bandpass (default: [20 80]).
%     DiffusionTime : Scalar t for heat kernel (default: 0.01).
%     MaxGain       : Maximum gain for inverse_heat (default: 10).
%     TransferFn    : Function handle @(lambdas) -> gains for custom filter.
%     MassMatrix    : Precomputed mass matrix (default: recomputed).
%
% OUTPUTS:
%     Filtered : [nVertices x nTime] filtered data.
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
%
%     The heat kernel at time t is the fundamental solution of the heat
%     equation du/dt = -L*u on the surface. Applying it is equivalent to
%     "blurring" the field by running diffusion for time t.
%
% EXAMPLES:
%     % Smooth source data with heat kernel (t=0.005)
%     Smoothed = tess_eigenmodes_filter(SurfaceFile, Sources, 'heat', 'DiffusionTime', 0.005);
%
%     % Keep only the 30 smoothest spatial patterns
%     LowPass = tess_eigenmodes_filter(SurfaceFile, Sources, 'lowpass', 'CutoffMode', 30);
%
%     % Isolate mid-frequency cortical patterns
%     BandPass = tess_eigenmodes_filter(SurfaceFile, Sources, 'bandpass', 'ModeRange', [10 60]);
%
% SEE ALSO: tess_eigenmodes_project, tess_eigenmodes_load, tess_eigenmodes, tess_laplacian

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
MassMatrix    = [];

for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'cutoffmode',    CutoffMode    = varargin{i+1};
        case 'moderange',     ModeRange     = varargin{i+1};
        case 'diffusiontime', DiffusionTime = varargin{i+1};
        case 'maxgain',       MaxGain       = varargin{i+1};
        case 'transferfn',    TransferFn    = varargin{i+1};
        case 'massmatrix',    MassMatrix    = varargin{i+1};
    end
end

%% ===== LOAD EIGENMODES =====
[Eigenmodes, isComputed] = tess_eigenmodes_load(SurfaceFile);
if ~isComputed
    error('No precomputed eigenmodes found in: %s. Run tess_eigenmodes first.', SurfaceFile);
end

Phi     = Eigenmodes.Vectors;   % [nV x nModes]
lambdas = Eigenmodes.Values;    % [nModes x 1]
nModes  = Eigenmodes.nModes;
nV      = size(Phi, 1);

%% ===== VALIDATE DATA =====
Data = double(Data);
if size(Data, 1) ~= nV
    error('Data has %d rows but surface has %d vertices.', size(Data, 1), nV);
end

%% ===== BUILD TRANSFER FUNCTION =====
h = zeros(nModes, 1);

switch lower(FilterType)
    case 'lowpass'
        % Sharp cutoff at mode CutoffMode
        CutoffMode = min(CutoffMode, nModes);
        h(1:CutoffMode) = 1;

    case 'highpass'
        % Sharp cutoff: keep modes >= CutoffMode
        CutoffMode = max(1, min(CutoffMode, nModes));
        h(CutoffMode:end) = 1;

    case 'bandpass'
        % Keep modes in [k1, k2]
        k1 = max(1, ModeRange(1));
        k2 = min(nModes, ModeRange(2));
        h(k1:k2) = 1;

    case 'heat'
        % Heat kernel: h_k = exp(-t * lambda_k)
        if DiffusionTime <= 0
            error('DiffusionTime must be positive (got %g).', DiffusionTime);
        end
        h = exp(-DiffusionTime * lambdas);

    case 'inverse_heat'
        % Inverse heat kernel with gain clamping
        if DiffusionTime <= 0
            error('DiffusionTime must be positive (got %g).', DiffusionTime);
        end
        h = exp(DiffusionTime * lambdas);
        h = min(h, MaxGain);  % clamp to avoid noise amplification

    case 'custom'
        % User-supplied transfer function
        if isempty(TransferFn) || ~isa(TransferFn, 'function_handle')
            error('Custom filter requires a TransferFn option (function handle).');
        end
        h = TransferFn(lambdas);
        if length(h) ~= nModes
            error('TransferFn must return a vector of length %d (got %d).', nModes, length(h));
        end
        h = h(:);

    otherwise
        error('Unknown filter type: %s. Use lowpass, highpass, bandpass, heat, inverse_heat, or custom.', FilterType);
end

%% ===== GET OR COMPUTE MASS MATRIX =====
if isempty(MassMatrix)
    SurfaceFileFull = file_fullpath(SurfaceFile);
    if isempty(SurfaceFileFull)
        SurfaceFileFull = SurfaceFile;
    end
    TessMat = load(SurfaceFileFull, 'Vertices', 'Faces');
    [~, MassMatrix] = tess_laplacian(double(TessMat.Vertices), double(TessMat.Faces), ...
        'MassType', Eigenmodes.MassType);
end

%% ===== APPLY FILTER =====
% u_filtered = Phi * diag(h) * Phi' * M * u
%
% Step 1: Project to spectral domain  — c = Phi' * M * u       [nModes x nTime]
% Step 2: Apply filter gains          — c_filt = h .* c        [nModes x nTime]
% Step 3: Reconstruct                 — u_filt = Phi * c_filt  [nV x nTime]

Coeffs = Phi' * (MassMatrix * Data);     % [nModes x nTime]
Coeffs = bsxfun(@times, h, Coeffs);      % element-wise multiply by gains
Filtered = Phi * Coeffs;                  % [nV x nTime]

end
