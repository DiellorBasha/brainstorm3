function [Results, errMsg] = tess_eigenmodes_leadfield(HeadModelFile, SurfaceFile, NoiseCovFile, varargin)
% TESS_EIGENMODES_LEADFIELD: Eigenmode-space source mapping via compressed lead field.
%
% USAGE:  [Results, errMsg] = tess_eigenmodes_leadfield(HeadModelFile, SurfaceFile, NoiseCovFile, ...)
%
% DESCRIPTION:
%     Implements eigenmode-space MEG/EEG source mapping by compressing the
%     lead field matrix into the Laplace-Beltrami eigenmode basis. Instead of
%     solving the inverse problem in vertex space (nChannels → nVertices), this
%     method solves it in eigenmode space (nChannels → K modes), yielding a
%     dramatically better-conditioned system.
%
%     The standard pipeline maps:
%         d(t) → u(t) ∈ ℝ^nVertices  (underdetermined: 300 sensors, 15000 sources)
%
%     The eigenmode pipeline maps:
%         d(t) → ĉ(t) ∈ ℝ^K          (nearly determined: 300 sensors, K≈200 modes)
%
%     The compressed forward model is:
%         L̃ = L · Φ  ∈ ℝ^(nChannels × K)
%
%     where L is the constrained lead field [nChannels × nVertices] and
%     Φ is the eigenmode matrix [nVertices × K]. Each column of L̃ is the
%     sensor topography of eigenmode k — the field pattern that spatial mode k
%     produces at the sensor array. This is grounded in Maxwell's equations:
%     the measurement is a spatial low-pass filter that preferentially sees
%     large-scale coherent patterns.
%
%     The eigenmode-space inverse is:
%         M̃ = Σ_d · L̃' · (L̃ · Σ_d · L̃' + λ · Σ_n)^{-1}
%
%     where Σ_d is the source prior covariance in eigenmode space and Σ_n is
%     the noise covariance. The source prior can encode a spectral power law:
%         Σ_d = diag(λ_k^{-α})
%     where α controls the smoothness assumption (α=0: flat, α=1: 1/f, α=2: smooth).
%
%     Available inverse methods:
%
%     'mne'      — Minimum norm estimate in eigenmode space.
%                  Recovers eigenmode coefficients with minimum L2 norm.
%
%     'dspm'     — Dynamic statistical parametric mapping.
%                  Noise-normalized MNE: divides by the expected noise
%                  amplitude per mode. Output is a z-score-like quantity.
%
%     'sloreta'  — Standardized low-resolution electromagnetic tomography.
%                  Normalizes by the resolution kernel diagonal. Produces
%                  unbiased localization (no depth bias in eigenmode space).
%
% INPUTS:
%     HeadModelFile : Path to Brainstorm head model file (relative or absolute)
%     SurfaceFile   : Path to Brainstorm surface file with precomputed eigenmodes
%     NoiseCovFile  : Path to noise covariance file (can be empty for 'mne' without whitening)
%
% OPTIONS (name-value pairs):
%     'Method'        : 'mne', 'dspm', or 'sloreta' (default: 'dspm')
%     'nModes'        : Number of eigenmodes to use (default: all available, clamped to nChannels)
%     'PriorAlpha'    : Spectral power-law exponent for source prior (default: 0)
%                       0 = flat (standard MNE), 1 = 1/f, 2 = diffusion-like
%     'SNR'           : Signal-to-noise ratio for regularization (default: 3)
%     'GoodChannel'   : Logical vector of good channels [nChannels × 1] (default: all true)
%     'ChannelTypes'  : Cell array of channel types to include (default: {'MEG', 'MEG MAG', 'MEG GRAD'})
%
% OUTPUTS:
%     Results : Struct with fields:
%       .ImagingKernel  : [K × nGoodChannels] eigenmode-space imaging kernel
%       .CompressedLF   : [nGoodChannels × K] compressed lead field L̃
%       .EigenGains     : [K × 1] per-mode gain factors (for dSPM/sLORETA normalization)
%       .Whitener       : [nGoodChannels × nGoodChannels] noise whitening matrix
%       .SourcePrior    : [K × 1] diagonal source prior Σ_d
%       .nModes         : Number of modes used
%       .Method         : Inverse method name
%       .SNR            : SNR value used
%       .PriorAlpha     : Prior exponent used
%       .SurfaceFile    : Surface file path
%       .ConditionNumber: Condition number of the whitened compressed system
%       .ConditionNumberFull: Condition number of the full whitened lead field (for comparison)
%     errMsg : Error message string (empty if success)
%
% MATHEMATICAL NOTES:
%     The compressed lead field L̃ = L·Φ has columns that are the sensor
%     topographies of the eigenmodes. Because MEG/EEG sensors act as spatial
%     low-pass filters (Biot-Savart 1/r² decay causes field cancellation for
%     fine spatial patterns), the columns of L̃ have naturally decaying norms:
%     ||L̃(:,k)|| decreases with k. This means the system is naturally
%     well-conditioned for the first K modes.
%
%     The source prior Σ_d = diag(λ_k^{-α}) encodes the assumption that
%     cortical activity has a spatial power spectrum that decays with spatial
%     frequency. This is equivalent to assuming the cortex preferentially
%     expresses large-scale coherent patterns — a physically reasonable prior
%     grounded in the volume conduction physics of cortical columns.
%
%     For α=0 (flat prior), the eigenmode inverse reduces to standard MNE
%     applied to the compressed lead field — equivalent to MNE with a
%     truncated basis constraint.
%
% SEE ALSO: tess_eigenmodes_load, tess_eigenmodes_project, process_eigenmodes_inverse

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

Results = [];
errMsg = '';

%% ===== PARSE OPTIONS =====
Method      = 'dspm';
nModes      = [];
PriorAlpha  = 0;
SNR         = 3;
GoodChannel = [];
ChannelTypes = {};

for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'method',       Method       = lower(varargin{i+1});
        case 'nmodes',       nModes       = varargin{i+1};
        case 'prioralpha',   PriorAlpha   = varargin{i+1};
        case 'snr',          SNR          = varargin{i+1};
        case 'goodchannel',  GoodChannel  = varargin{i+1};
        case 'channeltypes', ChannelTypes = varargin{i+1};
    end
end

% Validate method
if ~ismember(Method, {'mne', 'dspm', 'sloreta'})
    errMsg = ['Unknown method: ' Method '. Use mne, dspm, or sloreta.'];
    return;
end

%% ===== LOAD HEAD MODEL =====
% Load with constrained orientation (ApplyOrient=1)
HeadModel = in_bst_headmodel(HeadModelFile, 1);

% Validate surface head model
if ~strcmpi(HeadModel.HeadModelType, 'surface')
    errMsg = 'Eigenmode inverse requires a surface head model.';
    return;
end

% Get constrained lead field: [nChannels × nVertices]
Gain = double(HeadModel.Gain);
[nAllChannels, nVertices] = size(Gain);

% Apply channel selection
if ~isempty(GoodChannel) && length(GoodChannel) == nAllChannels
    iGoodChan = find(GoodChannel);
else
    iGoodChan = 1:nAllChannels;
end
L = Gain(iGoodChan, :);  % [nGoodChannels × nVertices]
nChannels = size(L, 1);

%% ===== LOAD EIGENMODES =====
[Eigenmodes, isComputed] = tess_eigenmodes_load(SurfaceFile);
if ~isComputed
    errMsg = ['No eigenmodes on surface: ' SurfaceFile '. Run "Compute eigenmodes" first.'];
    return;
end

% Check vertex count
nV_eigen = size(Eigenmodes.Vectors, 1);
if nVertices ~= nV_eigen
    errMsg = sprintf('Vertex mismatch: head model has %d sources, eigenmodes have %d vertices.', ...
        nVertices, nV_eigen);
    return;
end

% Determine number of modes to use
K_available = Eigenmodes.nModes;
if isempty(nModes) || nModes <= 0
    % Default: use min(nChannels, available modes) — beyond nChannels is noise
    K = min(nChannels, K_available);
else
    K = min(nModes, K_available);
end

Phi = double(Eigenmodes.Vectors(:, 1:K));       % [nVertices × K]
lambdas = double(Eigenmodes.Values(1:K));        % [K × 1]
lambdas = lambdas(:);

%% ===== COMPRESS LEAD FIELD =====
% L̃ = L · Φ  : [nChannels × K]
% Each column is the sensor topography of eigenmode k
L_tilde = L * Phi;

%% ===== LOAD AND APPLY NOISE COVARIANCE =====
if ~isempty(NoiseCovFile) && ~isequal(NoiseCovFile, '')
    NoiseCovMat = load(file_fullpath(NoiseCovFile));
    if isfield(NoiseCovMat, 'NoiseCov')
        NoiseCov_full = double(NoiseCovMat.NoiseCov);
    elseif isfield(NoiseCovMat, 'Value')
        NoiseCov_full = double(NoiseCovMat.Value);
    else
        errMsg = 'Could not find noise covariance matrix in file.';
        return;
    end
    % Extract good channels
    NoiseCov = NoiseCov_full(iGoodChan, iGoodChan);

    % Regularize noise covariance (Tikhonov: add small fraction of trace)
    regParam = 0.05;
    NoiseCov = NoiseCov + regParam * (trace(NoiseCov) / nChannels) * eye(nChannels);

    % Compute whitening matrix via SVD
    [Un, Sn2] = svd(NoiseCov);
    sn = sqrt(diag(Sn2));
    % Regularize: floor small singular values
    sn_floor = max(sn) * 1e-6;
    sn(sn < sn_floor) = sn_floor;
    iW = Un * diag(1 ./ sn) * Un';   % [nChannels × nChannels]

    % Whiten the compressed lead field
    L_white = iW * L_tilde;           % [nChannels × K]

    hasWhitener = true;
else
    % No noise covariance — identity whitening
    iW = eye(nChannels);
    L_white = L_tilde;
    hasWhitener = false;
end

%% ===== BUILD SOURCE PRIOR =====
% Spectral power-law prior: Σ_d = diag(λ_k^{-α})
% For α=0: flat prior (standard MNE)
% For α>0: smooth prior (suppresses high spatial frequencies)
if PriorAlpha > 0
    % Avoid division by zero for first eigenvalue (which may be 0 or near-0)
    lambdas_safe = max(lambdas, max(lambdas) * 1e-10);
    SourcePrior = lambdas_safe .^ (-PriorAlpha);
    % Normalize so max prior variance = 1
    SourcePrior = SourcePrior / max(SourcePrior);
else
    SourcePrior = ones(K, 1);
end

%% ===== COMPUTE EIGENMODE-SPACE INVERSE =====
% The system: d = L̃ · c + noise
% Where c ∈ ℝ^K are eigenmode coefficients
%
% Bayesian MAP estimate (weighted MNE):
%   M̃ = Σ_d · L̃' · (L̃ · Σ_d · L̃' + λ · I)^{-1}
%
% In whitened space (L_w = iW · L̃):
%   M̃_w = Σ_d · L_w' · (L_w · Σ_d · L_w' + λ · I)^{-1}
%   M̃ = M̃_w · iW     (de-whiten the data side)

% Apply source prior to whitened lead field
% L_ws = L_w · sqrt(Σ_d) — scale columns by sqrt of prior
sqrtPrior = sqrt(SourcePrior);
L_ws = L_white * diag(sqrtPrior);    % [nChannels × K]

% SVD of scaled whitened lead field
[U, S, V] = svd(L_ws, 'econ');       % U [nCh×r], S [r×r], V [K×r]
s = diag(S);                           % singular values
r = length(s);

% Compute condition numbers (for diagnostics)
ConditionNumber = s(1) / max(s(end), eps);

% Also compute condition of full (uncompressed) whitened lead field for comparison
L_full_white = iW * L;
s_full = svd(L_full_white);
ConditionNumberFull = s_full(1) / max(s_full(end), eps);

% Compute regularization parameter (Brainstorm convention):
%   Lambda = trace(L_ws · L_ws') / (nChannels · SNR²)
% This scales the regularization to the data, so the SNR parameter has
% a consistent meaning regardless of the absolute scale of the lead field.
Lambda = sum(s.^2) / (nChannels * SNR^2);

% Inverse filter in SVD domain:
%   M̃_w = V · diag(s/(s² + λ)) · U'
%
% This is equivalent to the standard Tikhonov/MNE formula but computed
% via SVD which is numerically stable for the K×K system
alpha = s ./ (s.^2 + Lambda);         % [r × 1] — filter factors

% Build imaging kernel
switch Method
    case 'mne'
        % M̃ = diag(sqrtPrior) · V · diag(alpha) · U' · iW
        Kernel = diag(sqrtPrior) * V * diag(alpha) * U' * iW;  % [K × nChannels]
        EigenGains = ones(K, 1);

    case 'dspm'
        % MNE kernel first
        K_mne = diag(sqrtPrior) * V * diag(alpha) * U';  % [K × nChannels] (in whitened space)

        % dSPM normalization: divide each mode by its expected noise amplitude
        % Noise amplitude per mode = sqrt(diag(K_mne · K_mne'))
        % Since we're in whitened space where noise cov = I:
        NoiseNorm = sqrt(sum(K_mne.^2, 2));       % [K × 1]
        NoiseNorm(NoiseNorm < max(NoiseNorm) * 1e-10) = max(NoiseNorm) * 1e-10;

        Kernel = diag(1 ./ NoiseNorm) * K_mne * iW;  % [K × nChannels]
        EigenGains = NoiseNorm;

    case 'sloreta'
        % MNE kernel first
        K_mne = diag(sqrtPrior) * V * diag(alpha) * U';  % [K × nChannels] (in whitened space)

        % sLORETA normalization: divide by sqrt of resolution kernel diagonal
        % Resolution matrix: R = K_mne · L_ws · diag(1/sqrtPrior)
        % For the diagonal: R_kk = K_mne(k,:) · L_white(:,k) · sqrtPrior(k)
        % But in SVD domain this simplifies:
        %   R = diag(sqrtPrior) · V · diag(s · alpha) · V' · diag(1/sqrtPrior)
        % Diagonal of R:
        s_alpha = s .* alpha;                      % [r × 1]
        % R_kk = sum_j V(k,j)^2 · s_alpha(j)     (prior scaling cancels on diagonal)
        ResDiag = sum((V.^2) .* (s_alpha'), 2);   % [K × 1]
        ResDiag(ResDiag < max(ResDiag) * 1e-10) = max(ResDiag) * 1e-10;

        Kernel = diag(1 ./ sqrt(ResDiag)) * K_mne * iW;  % [K × nChannels]
        EigenGains = sqrt(ResDiag);
end

%% ===== PACKAGE OUTPUT =====
Results.ImagingKernel   = Kernel;            % [K × nGoodChannels]
Results.CompressedLF    = L_tilde;           % [nGoodChannels × K]
Results.EigenGains      = EigenGains;        % [K × 1]
Results.Whitener        = iW;                % [nGoodChannels × nGoodChannels]
Results.SourcePrior     = SourcePrior;       % [K × 1]
Results.nModes          = K;
Results.Method          = Method;
Results.SNR             = SNR;
Results.PriorAlpha      = PriorAlpha;
Results.SurfaceFile     = SurfaceFile;
Results.GoodChannel     = iGoodChan;
Results.Eigenvalues     = lambdas;           % [K × 1]
Results.ConditionNumber     = ConditionNumber;
Results.ConditionNumberFull = ConditionNumberFull;

end
