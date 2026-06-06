function Results = bst_cwt_fiedler_pipeline(HeadModelFile, DataFile, SurfaceFile, varargin)
% BST_CWT_FIEDLER_PIPELINE  End-to-end CWT/Fiedler traveling-wave pipeline.
%
% USAGE:
%   Results = bst_cwt_fiedler_pipeline(HeadModelFile, DataFile, SurfaceFile)
%   Results = bst_cwt_fiedler_pipeline(..., 'TargetFreq', 10, 'SNR', 3)
%
% DESCRIPTION:
%   Complete parallel pipeline for alpha traveling-wave detection using the
%   CWT/Fiedler approach.  Runs independently of — and comparable to — the
%   existing Hilbert-based pipeline (bst_eigenmode_analytic_inverse →
%   bst_face_sign_correct → bst_face_wavefront_track).
%
%   Stage 1  bst_sensor_cwt           CWT of sensor data
%   Stage 2  bst_eigenmode_cwt_inverse  frequency-resolved eigenmode coefficients
%   Stage 3  bst_lambda_omega_spectrum  (λ,ω) dispersion analysis
%   Stage 4  bst_dispersion_fit         wave speed + scale selection
%   Stage 5  bst_fiedler_joint_field    joint Fiedler torus field Z_joint
%   Stage 6  bst_fiedler_wave_analyze   source, speed, PLV in Fiedler frame
%
%   The Fiedler frame from TessMat.nxr (tess_nxr_populate) is used throughout;
%   the flat-prior MNE inverse keeps mode selection to the (λ,ω) analysis.
%
% INPUTS:
%   HeadModelFile  face-based eigenmode headmodel (isEigenmode=1, isFaceBased=1)
%   DataFile       Brainstorm data file
%   SurfaceFile    Brainstorm cortex surface file
%
% OPTIONS (name-value):
%   'TargetFreq'       Hz for alpha wave (default: 10)
%   'FreqRange'        [f_lo f_hi] for CWT (default: [1 30])
%   'VoicesPerOctave'  CWT frequency resolution (default: 12)
%   'SNR'              Tikhonov SNR for flat-prior inverse (default: 3)
%   'AmpThreshold'     spatial amplitude gate fraction (default: 0.08)
%   'AmpThresholdTime' temporal amplitude gate fraction (default: 0.20)
%   'Verbose'          logical (default: true)
%
% OUTPUT  Results struct:
%   .W_sensor       [nCh × nScales × nTime]  raw CWT coefficients
%   .Theta          [K × nScales × nTime]    eigenmode CWT coefficients
%   .FaceAmp        [nLHF × nScales × nTime] face-space per scale
%   .S_joint        [K × nFreq]              (λ,ω) power spectrum
%   .f_ax, .k_ax                             spectral axes
%   .peak_info                               dispersion ridge summary
%   .DispInfo                                dispersion fit details
%   .Z_joint        [nLHF × nTime]           joint Fiedler field
%   .phi_F          [nLHF × 1]              Fiedler longitude
%   .phi_T          [1 × nTime]             temporal Fiedler phase
%   .confidence     [nLHF × nTime]          joint amplitude confidence
%   .Wave                                    bst_fiedler_wave_analyze output
%   .WTInfo                                  wavelet transform metadata
%   .InvInfo                                 inverse metadata
%
% Authors: Diellor Basha, 2026

t0 = tic;

%% ── Parse options ────────────────────────────────────────────────────────
TargetFreq       = 10;
FreqRange        = [1 30];
VoicesPerOctave  = 12;
SNR              = 3;
AmpThreshold     = 0.08;
AmpThresholdTime = 0.20;
Verbose          = true;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'targetfreq',       TargetFreq       = varargin{k+1};
        case 'freqrange',        FreqRange        = varargin{k+1};
        case 'voicesperoctave',  VoicesPerOctave  = varargin{k+1};
        case 'snr',              SNR              = varargin{k+1};
        case 'ampthreshold',     AmpThreshold     = varargin{k+1};
        case 'ampthresholdtime', AmpThresholdTime = varargin{k+1};
        case 'verbose',          Verbose          = logical(varargin{k+1});
    end
end

vp = @(s) fprintf('[CWT-Fiedler] %s\n', s);

%% ── Stage 1: CWT of sensor data ─────────────────────────────────────────
if Verbose, vp('Stage 1/6 — sensor CWT...'); end
[W_sensor, WTInfo] = bst_sensor_cwt(DataFile, ...
    'FreqRange', FreqRange, 'VoicesPerOctave', VoicesPerOctave, ...
    'Verbose', Verbose);

%% ── Stage 2: Eigenmode CWT inverse (mode coefficients only) ──────────────
% Defer face-space reconstruction: FaceAmp at all scales would be ~60 GB.
% We reconstruct only the selected alpha scale after dispersion fit (Stage 5).
if Verbose, vp('Stage 2/6 — eigenmode CWT inverse...'); end
[Theta, ~, InvInfo] = bst_eigenmode_cwt_inverse(W_sensor, WTInfo, HeadModelFile, ...
    'SNR', SNR, 'FaceSpace', false, 'Verbose', Verbose);

%% ── Stage 3: (λ,ω) joint spectral analysis ───────────────────────────────
if Verbose, vp('Stage 3/6 — (λ,ω) spectrum...'); end
HM = in_bst_headmodel(HeadModelFile, 0);
[S_joint, f_ax, k_ax, lambda_ax, peak_info] = bst_lambda_omega_spectrum( ...
    Theta, double(HM.Eigenvalues(:)), WTInfo, ...
    'FreqRange', FreqRange, 'Normalize', 'mode', 'Verbose', Verbose);

%% ── Stage 4: Dispersion fit → select alpha scale ─────────────────────────
if Verbose, vp('Stage 4/6 — dispersion fit...'); end
[v_est, s_alpha, DispInfo] = bst_dispersion_fit(S_joint, double(HM.Eigenvalues(:)), ...
    f_ax, WTInfo, 'TargetFreq', TargetFreq, 'Verbose', Verbose);

%% ── Stage 5: Joint Fiedler torus field ───────────────────────────────────
% Reconstruct the FULL face field ONLY at the selected alpha scale (memory-safe),
% then demodulate against the Fiedler gauge.
if Verbose, vp('Stage 5/6 — joint Fiedler field...'); end
Phi_f = double(HM.PhiFaces);                       % [nLHF × K]
s_face_alpha = Phi_f * squeeze(Theta(:, s_alpha, :));   % [nLHF × nTime] full reconstruction
[Z_joint, phi_F, phi_T, confidence] = bst_fiedler_joint_field( ...
    s_face_alpha, InvInfo.FaceIndices, SurfaceFile, 'Verbose', Verbose);
FaceAmp = s_face_alpha;   % return the single-scale face field

%% ── Stage 6: Wave analysis in Fiedler coordinates ────────────────────────
if Verbose, vp('Stage 6/6 — Fiedler wave analysis...'); end
Wave = bst_fiedler_wave_analyze(Z_joint, phi_F, InvInfo.FaceIndices, SurfaceFile, WTInfo, ...
    'CenterFreq', TargetFreq, 'AmpThreshold', AmpThreshold, ...
    'AmpThresholdTime', AmpThresholdTime, 'Verbose', Verbose);

%% ── Pack results ─────────────────────────────────────────────────────────
Results.W_sensor   = W_sensor;
Results.Theta      = Theta;
Results.FaceAmp    = FaceAmp;
Results.S_joint    = S_joint;
Results.f_ax       = f_ax;
Results.k_ax       = k_ax;
Results.lambda_ax  = lambda_ax;
Results.peak_info  = peak_info;
Results.DispInfo   = DispInfo;
Results.Z_joint    = Z_joint;
Results.phi_F      = phi_F;
Results.phi_T      = phi_T;
Results.confidence = confidence;
Results.Wave       = Wave;
Results.WTInfo     = WTInfo;
Results.InvInfo    = InvInfo;

if Verbose
    fprintf('\n[CWT-Fiedler] Done in %.1f s\n', toc(t0));
    fprintf('  Dispersion: v=%.2f m/s  R²=%.3f\n', DispInfo.v_linear, DispInfo.R2);
    fprintf('  Source: face %d  [%.1f %.1f %.1f] mm\n', ...
        Wave.SourceFace, Wave.SourcePos*1000);
    fprintf('  PLV: max=%.3f  v_Fiedler=%.2f m/s\n', ...
        max(Wave.PLV,[],'omitnan'), Wave.v_fiedler);
end
end
