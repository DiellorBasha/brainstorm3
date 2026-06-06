function Results = bst_eigenmode_analytic_inverse(HeadModelFile, DataFile, varargin)
% BST_EIGENMODE_ANALYTIC_INVERSE  Complex analytic inverse in eigenmode space.
%
% USAGE:
%   Results = bst_eigenmode_analytic_inverse(HeadModelFile, DataFile)
%   Results = bst_eigenmode_analytic_inverse(HeadModelFile, DataFile, 'FreqBand', [7 13])
%   Results = bst_eigenmode_analytic_inverse(HeadModelFile, DataFile, 'SNR', 3, 'Prior', 'log')
%
% MOTIVATION:
%   Standard eigenmode MNE solves  θ(t) = K · d(t)  with real d(t) and real θ(t).
%   For traveling waves at a carrier frequency ω₀, each coefficient should be:
%       θ_k(t) = A_k(t) · exp(i(ω₀t + φ_k(t)))
%   Taking the COMPLEX ANALYTIC SIGNAL of the sensor data and passing it through
%   the same linear solver gives complex θ̂(t) directly — no post-processing Hilbert,
%   no frequency doubling, no sign ambiguity.
%
%   The reconstructed source field:
%       s̃(x,t) = Ψ · θ̂(t)               [complex, nVert x nTime]
%   carries:
%       |s̃(x,t)|       amplitude envelope  — always positive, no rectification
%       arg(s̃(x,t))    instantaneous phase — in the eigenmode gauge, comparable
%                       across vertices because eigenmodes span a consistent basis
%
% PIPELINE:
%   1. Load eigenmode headmodel (face-based or vertex-based, isEigenmode=1)
%   2. Load + bandpass sensor data around ω₀
%   3. Compute analytic sensor signal: d̃ = d + i·H[d]  per channel
%   4. Apply eigenmode MNE to d̃:  θ̂(t) = K · d̃(t)   (same K as real case)
%   5. Reconstruct: s̃(x,t) = Ψ[:, ModeIndices] · θ̂(t)
%
%   The kernel K is computed from L̃ (the composed eigenmode leadfield) and the
%   noise covariance — identical to the standard solve, since the Tikhonov
%   problem is linear and extends to complex inputs without modification.
%
% INPUTS:
%   HeadModelFile  Brainstorm relative path to an eigenmode headmodel
%                  (isEigenmode=1; face-based or vertex-based)
%   DataFile       Brainstorm relative path to a data file
%
% OPTIONS (name-value):
%   'FreqBand'   [f_lo f_hi]  bandpass before Hilbert (default: [] = no bandpass)
%   'SNR'        scalar        Tikhonov SNR (default: 3)
%   'Prior'      'log'|'flat'  spectral source prior (default: 'log')
%   'Alpha'      scalar        power-law exponent for 'power' prior (default: 1)
%
% OUTPUT:
%   Results  struct with fields:
%     .ImageGridAmp   [nVert x nTime]  complex source field s̃(x,t)
%                     — Re(.) = real part, angle(.) = instantaneous phase
%     .Amplitude      [nVert x nTime]  abs(s̃) — amplitude envelope
%     .Phase          [nVert x nTime]  angle(s̃) — instantaneous phase [rad]
%     .theta          [K x nTime]      complex eigenmode coefficients θ̂(t)
%     .ImagingKernel  [K x nChan]      mode-space kernel K (for re-application)
%     .Time           [1 x nTime]
%     .Eigenvalues    [K x 1]
%     .nModes         K
%     .SurfaceFile
%
% WAVE DETECTION:
%   The spatial phase map Results.Phase(:, t0) at any time t0 encodes the
%   instantaneous wave structure.  The phase gradient ∇_S Phase(x,t0) points
%   in the propagation direction; its magnitude is the local wavenumber.
%   Because the reconstruction uses a consistent eigenmode basis, phase
%   comparisons between any two vertices are meaningful without further
%   gauge correction.
%
% SEE ALSO: bst_face_eigenmode_leadfield, bst_eigenmode_oscillator_inverse,
%           bst_inverse_eigenmodes
%
% Authors: Diellor Basha, 2026

%% ── Parse options ────────────────────────────────────────────────────────
FreqBand     = [];
SNR          = 3;
Prior        = 'log';
Alpha        = 1;
HelmetFilter = [];   % [nCh x J] helmet eigenmodes for spatial pre-filtering
FaceSpace    = false;  % return s_f = PhiFaces*theta in face space (physically consistent)
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'freqband',     FreqBand     = varargin{k+1};
        case 'snr',          SNR          = varargin{k+1};
        case 'prior',        Prior        = lower(varargin{k+1});
        case 'alpha',        Alpha        = varargin{k+1};
        case 'helmetfilter', HelmetFilter = varargin{k+1};
        case 'facespace',    FaceSpace    = logical(varargin{k+1});
    end
end

%% ── Load eigenmode headmodel ─────────────────────────────────────────────
HM = in_bst_headmodel(HeadModelFile, 0);
if ~isfield(HM,'isEigenmode') || ~HM.isEigenmode
    error('bst_eigenmode_analytic_inverse:NotEigenmode', ...
        'HeadModelFile must be an eigenmode headmodel (isEigenmode=1).');
end
L_tilde  = double(HM.Gain);          % [nCh x K]
lambdas  = double(HM.Eigenvalues(:));
K        = size(L_tilde, 2);
ModeIdx  = HM.ModeIndices;

fprintf('Analytic inverse: K=%d modes, λ=[%.0f…%.0f]\n', K, min(lambdas), max(lambdas));

%% ── Load sensor data ─────────────────────────────────────────────────────
DataMat = in_bst_data(DataFile);
F       = double(DataMat.F);          % [nAllCh x nTime]
Time    = DataMat.Time(:)';
Fs      = 1 / mean(diff(Time));
nTime   = numel(Time);

% Good channels (same subset as the headmodel)
nHMCh = size(L_tilde, 1);
[~, iStudy] = bst_get('DataFile', DataFile);
ChanMat = in_bst_channel(bst_get('ChannelFileForStudy', iStudy));
iMEG = find(strcmpi({ChanMat.Channel.Type}, 'MEG'));
if numel(iMEG) ~= nHMCh
    iMEG = iMEG(1:nHMCh);   % trim to headmodel channel count
end
d = F(iMEG, :);               % [nCh x nTime]

%% ── Bandpass (optional) ─────────────────────────────────────────────────
if ~isempty(FreqBand)
    d = bst_bandpass_hfilter(d, Fs, FreqBand(1), FreqBand(2), 0, 'iir');
    fprintf('  Bandpassed to [%.1f %.1f] Hz\n', FreqBand(1), FreqBand(2));
end

%% ── Helmet spatial pre-filter (optional) ─────────────────────────────────
% Project sensor data onto the J smoothest helmet LBO modes before inversion.
% d_filt = Phi_h * (Phi_h' * d)  removes spatially incoherent sensor noise
% while preserving spatially coherent brain signals.
if ~isempty(HelmetFilter)
    Phi_h  = HelmetFilter;                    % [nCh x J]
    J_helm = size(Phi_h, 2);
    % Proper oblique projector onto col(Phi_h): P = Phi*(Phi'*Phi)^{-1}*Phi'
    % Works correctly even when Phi columns are not Euclidean-orthonormal
    % (as happens after nearest-neighbor interpolation to off-hull sensors).
    AtAinv = inv(Phi_h' * Phi_h);            % [J x J]
    d      = Phi_h * (AtAinv * (Phi_h' * d)); % [nCh x nTime] — smooth modes only
    fprintf('  Helmet filter: J=%d modes  (%.0f%% of channels)\n', ...
        J_helm, 100*J_helm/nHMCh);
end

%% ── Complex analytic signal per channel ──────────────────────────────────
% d̃(t) = d(t) + i · H[d(t)]   where H is the Hilbert transform along time.
% Use explicit regular-transpose (.') throughout: MATLAB's ctranspose (')
% conjugates complex intermediate results and produces the conjugate analytic
% signal (negative instantaneous frequency) when applied to the real-but-
% dimensionally-transposed matrix from hilbert(). See implementation note.
d_analytic = d + 1i * imag(hilbert(d.')).';   % [nCh x nTime]

%% ── Mode-space kernel (standard Tikhonov) ────────────────────────────────
% Same computation as bst_inverse_eigenmodes / solve_modespace, exposed here
% so we can apply it to complex data.
R    = bst_eigenmode_prior(lambdas, K, Prior, Alpha);   % [K x 1]
sP   = sqrt(R(:));                   % column scaling
Lws  = L_tilde .* sP.';             % [nCh x K]
[U, S, V] = svd(Lws, 'econ');
s    = diag(S);
Lam  = sum(s.^2) / (nHMCh * SNR^2);
alph = s ./ (s.^2 + Lam);          % MNE filter factors
% K_mne [K x nCh] maps sensor data to (whitened) eigenmode coefficients
K_mne = diag(sP) * V * diag(alph) * U';

%% ── Apply kernel to complex analytic data ────────────────────────────────
theta_hat = K_mne * d_analytic;     % [K x nTime]  complex eigenmode coefficients

fprintf('  |θ̂| mean=%.3e  max=%.3e\n', mean(abs(theta_hat(:))), max(abs(theta_hat(:))));

%% ── Reconstruct complex source field ─────────────────────────────────────
% Load Phi for vertex-space reconstruction
if isfield(HM,'isFaceBased') && HM.isFaceBased && isfield(HM,'PhiVertices')
    % PhiVertices [nV x K] already has exactly the K selected modes in order
    Phi = double(HM.PhiVertices);
    idx = 1:K;
else
    Eig = in_tess_eigenmodes(HM.SurfaceFile);
    Phi = double(Eig.Vectors);
    idx = ModeIdx(1:K);   % global indices into full Phi
end
s_tilde = Phi(:, idx) * theta_hat;  % [nV x nTime]  complex source field

%% ── Face-space reconstruction (physically consistent) ─────────────────────
% s_f = Ψ·θ̂ [nLHF × nTime] is the primal 2-form: current flux per face.
% This is the quantity the forward model actually encodes — not vertex dipoles.
FaceGridAmp = [];
if FaceSpace && isfield(HM,'PhiFaces') && ~isempty(HM.PhiFaces) && ...
        isfield(HM,'isFaceBased') && HM.isFaceBased
    FaceGridAmp = double(HM.PhiFaces) * theta_hat;   % [nLHF x nTime] complex
    fprintf('  Face-space reconstruction: [%d faces x %d]\n', ...
        size(FaceGridAmp,1), size(FaceGridAmp,2));
end

%% ── Pack results ─────────────────────────────────────────────────────────
Results = struct();
Results.ImageGridAmp  = s_tilde;         % complex [nV x nTime]
Results.Amplitude     = abs(s_tilde);    % [nV x nTime] envelope (no rectification)
Results.Phase         = angle(s_tilde);  % [nV x nTime] instantaneous phase [rad]
Results.theta         = theta_hat;       % [K x nTime]  complex mode coefficients
Results.ImagingKernel = K_mne;           % [K x nCh]    for re-application
Results.PhiVertices   = Phi(:, idx);     % [nV x K]     cached for Poisson sharpen
Results.FaceGridAmp   = FaceGridAmp;     % [nLHF x nTime] complex, or [] if FaceSpace=false
Results.FaceIndices   = [];              % filled below if FaceSpace
if FaceSpace && ~isempty(FaceGridAmp)
    Results.FaceIndices = HM.FaceIndices;
end
Results.Time          = Time;
Results.Eigenvalues   = lambdas;
Results.nModes        = K;
Results.SurfaceFile   = HM.SurfaceFile;
Results.FreqBand      = FreqBand;
Results.SNR           = SNR;
Results.Prior         = Prior;

fprintf('  Done: ImageGridAmp [%d x %d] complex\n', size(s_tilde,1), size(s_tilde,2));
end
