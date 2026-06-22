function Results = bst_eigenmode_oscillator_inverse(HeadModelFile, DataFile, varargin)
% BST_EIGENMODE_OSCILLATOR_INVERSE  Eigenmode inverse with oscillator temporal prior.
%
% USAGE:
%   Results = bst_eigenmode_oscillator_inverse(HeadModelFile, DataFile, 'Freq', 10)
%   Results = bst_eigenmode_oscillator_inverse(..., 'Gamma', 1, 'SNR', 3)
%
% MOTIVATION:
%   Standard eigenmode MNE treats each time point independently — no temporal
%   structure is assumed.  For traveling alpha waves, each eigenmode coefficient
%   should oscillate at the carrier frequency ω₀:
%
%       θ̈_k + 2γ θ̇_k + ω₀² θ_k = u_k(t)
%
%   where γ is the damping rate and u_k(t) is the innovation (unknown input).
%   This is a second-order linear dynamical system — exactly a Kalman filter.
%
%   Encoding this prior means the inverse "knows" the source oscillates at ω₀.
%   It suppresses contributions at other frequencies and produces temporally
%   smooth, oscillatory source trajectories.  For traveling waves, the recovered
%   θ̂_k(t) follow the wave's amplitude and phase envelope across modes.
%
% ALGORITHM:
%   State: x_k(t) = [θ_k(t), θ̇_k(t)]'  ∈ ℝ²  per mode k (total state = 2K)
%   Dynamics (discretised at Ts = 1/Fs):
%       x_k(t+1) = A_osc · x_k(t) + noise_process
%   Observation (all modes together):
%       d(t) = L̃ · [θ₁(t), ..., θ_K(t)]' + noise_sensor
%
%   The Kalman filter propagates the state estimate forward; the RTS smoother
%   then passes backward to incorporate future observations.  The result is the
%   MAP estimate of the oscillatory trajectory given all the data.
%
% INPUTS:
%   HeadModelFile  Brainstorm relative path to an eigenmode headmodel
%   DataFile       Brainstorm relative path to a data file
%
% OPTIONS (name-value):
%   'Freq'      ω₀/(2π) in Hz — carrier frequency of the oscillation (default: 10)
%   'Gamma'     γ — damping rate in Hz; higher = faster amplitude decay (default: 1)
%   'SNR'       sensor SNR for measurement noise variance (default: 3)
%   'FreqBand'  [f_lo f_hi] optional pre-bandpass (default: [])
%   'Prior'     spectral prior type (default: 'log')
%
% OUTPUT:
%   Results  struct:
%     .ImageGridAmp   [nVert x nTime]  real source amplitude (NOT complex here;
%                                      use .AnalyticAmp / .AnalyticPhase below)
%     .AnalyticAmp    [nVert x nTime]  amplitude envelope of θ-reconstructed field
%     .AnalyticPhase  [nVert x nTime]  instantaneous phase [rad]
%     .theta          [K x nTime]      smoothed mode coefficients (real)
%     .theta_dot      [K x nTime]      smoothed mode velocities (real)
%     .Time           [1 x nTime]
%     .Eigenvalues    [K x 1]
%
% INTERPRETATION:
%   The Kalman smoother enforces that each mode oscillates at 'Freq' Hz with
%   damping 'Gamma'.  Modes that are not excited at that frequency are strongly
%   suppressed.  The result is a spatially smooth, temporally constrained source
%   map that is more sensitive to the target oscillation than standard MNE.
%
%   The analytic amplitude and phase are recovered by applying the Hilbert
%   transform to the smoothed real trajectory — valid because the Kalman prior
%   has already restricted the signal to a narrow band.
%
% SEE ALSO: bst_eigenmode_analytic_inverse, bst_face_eigenmode_leadfield
%
% Authors: Diellor Basha, 2026

%% ── Parse options ────────────────────────────────────────────────────────
Freq     = 10;
Gamma    = 1;
SNR      = 3;
FreqBand = [];
Prior    = 'log';
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'freq',     Freq     = varargin{k+1};
        case 'gamma',    Gamma    = varargin{k+1};
        case 'snr',      SNR      = varargin{k+1};
        case 'freqband', FreqBand = varargin{k+1};
        case 'prior',    Prior    = lower(varargin{k+1});
    end
end
omega0 = 2*pi*Freq;

%% ── Load headmodel and data ──────────────────────────────────────────────
HM = in_bst_headmodel(HeadModelFile, 0);
if ~isfield(HM,'isEigenmode') || ~HM.isEigenmode
    error('bst_eigenmode_oscillator_inverse:NotEigenmode', ...
        'HeadModelFile must be an eigenmode headmodel (isEigenmode=1).');
end
L_tilde = double(HM.Gain);
lambdas = double(HM.Eigenvalues(:));
K       = size(L_tilde, 2);
ModeIdx = HM.ModeIndices;

DataMat = in_bst_data(DataFile);
F       = double(DataMat.F);
Time    = DataMat.Time(:)';
Fs      = 1 / mean(diff(Time));
Ts      = 1/Fs;
nTime   = numel(Time);

[~, iStudy] = bst_get('DataFile', DataFile);
ChanMat = in_bst_channel(bst_get('ChannelFileForStudy', iStudy));
iMEG = find(strcmpi({ChanMat.Channel.Type}, 'MEG'));
nHMCh = size(L_tilde,1);
if numel(iMEG) > nHMCh, iMEG = iMEG(1:nHMCh); end
d = F(iMEG, :);

if ~isempty(FreqBand)
    d = bst_bandpass_hfilter(d, Fs, FreqBand(1), FreqBand(2), 0, 'iir');
end

fprintf('Oscillator inverse: K=%d modes, f₀=%.1fHz, γ=%.1fHz\n', K, Freq, Gamma);

%% ── Spectral source prior (same as standard MNE) ────────────────────────
R  = bst_eigenmode_prior(lambdas, K, Prior, 1);   % [K x 1] spectral variance

%% ── Discretised oscillator state matrix ─────────────────────────────────
% Continuous: [θ̈; θ̈] = [0 1; -ω₀² -2γ] [θ; θ̇]  +  input
% Discretise with matrix exponential at Ts
A_c = [0, 1; -omega0^2, -2*pi*Gamma];   % continuous dynamics
A_d = expm(A_c * Ts);                   % discrete state transition [2 x 2]

% State: x = [θ₁, θ̇₁, θ₂, θ̇₂, ..., θ_K, θ̇_K]'  [2K x 1]
% Block-diagonal A matrix
A_full = kron(eye(K), A_d);             % [2K x 2K]

% Observation: d(t) = H · x(t)  where H selects the first component of each block
H_obs = kron(eye(K), [1 0]);            % [K x 2K]  selects θ_k from [θ_k, θ̇_k]
L_obs = L_tilde * H_obs;                % [nCh x 2K]  combined observation matrix

%% ── Noise covariances ────────────────────────────────────────────────────
% Measurement noise: estimated from data variance, scaled by SNR
% Using data variance rather than leadfield power gives physically correct units
sig_y2  = var(d(:)) / SNR^2;
R_noise = sig_y2 * eye(nHMCh);        % [nCh x nCh] sensor noise

% Process noise: drives the oscillation; variance proportional to spectral prior
% Q models "how much can the amplitude change per time step"
q_diag   = reshape([R(:)'; zeros(1,K)], 2*K, 1) * (2*pi*Gamma * Ts);
Q_noise  = diag(max(q_diag, 1e-30));   % [2K x 2K] process noise

%% ── Kalman filter (forward pass) ─────────────────────────────────────────
x_filt  = zeros(2*K, nTime);   % filtered state means
P_filt  = zeros(2*K, 2*K, nTime);
x_pred  = zeros(2*K, nTime);
P_pred  = zeros(2*K, 2*K, nTime);

% Initialise
x_k  = zeros(2*K, 1);
P_k  = Q_noise * 100;          % wide initial uncertainty

for t = 1:nTime
    % Predict
    x_pr     = A_full * x_k;
    P_pr     = A_full * P_k * A_full' + Q_noise;
    x_pred(:,t)   = x_pr;
    P_pred(:,:,t) = P_pr;

    % Update
    S_k = L_obs * P_pr * L_obs' + R_noise;
    K_k = P_pr * L_obs' / S_k;            % Kalman gain [2K x nCh]
    innov = d(:,t) - L_obs * x_pr;
    x_k   = x_pr + K_k * innov;
    P_k   = (eye(2*K) - K_k * L_obs) * P_pr;
    x_filt(:,t)   = x_k;
    P_filt(:,:,t) = P_k;
end

%% ── RTS smoother (backward pass) ─────────────────────────────────────────
x_smooth = zeros(2*K, nTime);
x_smooth(:, nTime) = x_filt(:, nTime);

for t = nTime-1:-1:1
    G_t = P_filt(:,:,t) * A_full' / P_pred(:,:,t+1);
    x_smooth(:,t) = x_filt(:,t) + G_t * (x_smooth(:,t+1) - x_pred(:,t+1));
end

% Extract θ_k(t) — first component of each 2-state block
theta_smooth = x_smooth(1:2:end, :);   % [K x nTime]
tdot_smooth  = x_smooth(2:2:end, :);   % [K x nTime]

%% ── Reconstruct source field ─────────────────────────────────────────────
if isfield(HM,'isFaceBased') && HM.isFaceBased && isfield(HM,'PhiVertices')
    Phi = double(HM.PhiVertices);
    idx = 1:K;
else
    Eig = in_tess_eigenmodes(HM.SurfaceFile);
    Phi = double(Eig.Vectors);
    idx = ModeIdx(1:K);
end
s_real = Phi(:,idx) * theta_smooth;    % [nV x nTime]  real source field

% Analytic amplitude and phase via Hilbert (valid: Kalman already band-limited)
s_analytic  = s_real + 1i*imag(hilbert(s_real')');
amp_map     = abs(s_analytic);
phase_map   = angle(s_analytic);

%% ── Pack results ─────────────────────────────────────────────────────────
Results = struct();
Results.ImageGridAmp   = s_real;         % [nV x nTime] real (standard display)
Results.AnalyticAmp    = amp_map;        % [nV x nTime] amplitude envelope
Results.AnalyticPhase  = phase_map;      % [nV x nTime] instantaneous phase
Results.theta          = theta_smooth;   % [K x nTime]  mode coefficients
Results.theta_dot      = tdot_smooth;    % [K x nTime]  mode velocities
Results.Time           = Time;
Results.Eigenvalues    = lambdas;
Results.nModes         = K;
Results.SurfaceFile    = HM.SurfaceFile;
Results.Freq           = Freq;
Results.Gamma          = Gamma;
Results.SNR            = SNR;

fprintf('  Done: ImageGridAmp [%d x %d]\n', size(s_real,1), size(s_real,2));
end
