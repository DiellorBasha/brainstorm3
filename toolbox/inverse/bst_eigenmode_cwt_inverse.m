function [Theta, FaceAmp, InvInfo] = bst_eigenmode_cwt_inverse(W_sensor, WTInfo, HeadModelFile, varargin)
% BST_EIGENMODE_CWT_INVERSE  Eigenmode MNE inverse applied to CWT sensor data.
%
% USAGE:
%   [Theta, FaceAmp, InvInfo] = bst_eigenmode_cwt_inverse(W_sensor, WTInfo, HeadModelFile)
%   [Theta, FaceAmp, InvInfo] = bst_eigenmode_cwt_inverse(..., 'SNR', 3, 'FaceSpace', true)
%
% Applies the linear MNE kernel K_mne [K x nCh] to CWT coefficients
% W_sensor [nCh x nScales x nTime], giving θ̂(k,s,t) — complex amplitude
% of eigenmode k at wavelet scale s and time t.  Uses a FLAT prior (R=1);
% mode selection from the dispersion spectrum (λ,ω) is done downstream.
%
% INPUTS:
%   W_sensor      [nCh x nScales x nTime] complex — from bst_sensor_cwt
%   WTInfo        struct (.Fs, .f_scales, .nScales, .nTime, .ChanIdx, …)
%   HeadModelFile Brainstorm relative path to a face-based eigenmode headmodel
%
% OPTIONS: 'SNR' (default 3), 'FaceSpace' (default true), 'Verbose' (default true)
%
% OUTPUTS:
%   Theta    [K x nScales x nTime] complex — eigenmode CWT coefficients
%   FaceAmp  [nLHF x nScales x nTime] complex — PhiFaces*Theta per scale; [] if FaceSpace=false
%   InvInfo  struct: .K .lambdas .nScales .f_scales .SNR .FaceIndices .HeadModelFile
%
% SEE ALSO: bst_sensor_cwt, bst_eigenmode_analytic_inverse, bst_face_eigenmode_leadfield
%
% Authors: Diellor Basha, 2026

%% ── Parse options ────────────────────────────────────────────────────────
SNR       = 3;
FaceSpace = true;
Verbose   = true;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'snr',        SNR       = varargin{k+1};
        case 'facespace',  FaceSpace = logical(varargin{k+1});
        case 'verbose',    Verbose   = logical(varargin{k+1});
    end
end

%% ── Load eigenmode headmodel ─────────────────────────────────────────────
HM = in_bst_headmodel(HeadModelFile, 0);
if ~isfield(HM,'isEigenmode') || ~HM.isEigenmode
    error('bst_eigenmode_cwt_inverse:NotEigenmode', ...
        'HeadModelFile must be an eigenmode headmodel (isEigenmode=1).');
end
L_tilde = double(HM.Gain);           % [nCh x K]
lambdas = double(HM.Eigenvalues(:)); % [K x 1]
K       = size(L_tilde, 2);
nHMCh   = size(L_tilde, 1);

%% ── Channel alignment ────────────────────────────────────────────────────
% W_sensor rows correspond to WTInfo.ChanIdx in the original channel list.
% If the row count already matches the headmodel, use directly.
nWTCh = size(W_sensor, 1);
if nWTCh ~= nHMCh
    % Map WTInfo.ChanIdx to the HM's own channel ordering (1-based local index).
    % HM channels are assumed to be the first nHMCh entries of WTInfo.ChanIdx.
    [~, local_ch_idx] = ismember(1:nHMCh, WTInfo.ChanIdx);
    local_ch_idx = local_ch_idx(local_ch_idx > 0);
    L_tilde = L_tilde(local_ch_idx, :);
    W_sensor = W_sensor(local_ch_idx, :, :);
    nHMCh = size(L_tilde, 1);
end

nScales = WTInfo.nScales;
nTime   = WTInfo.nTime;

%% ── Memory estimate ──────────────────────────────────────────────────────
mem_GB = K * nScales * nTime * 16 / 1e9;
if Verbose
    fprintf('CWT inverse: K=%d modes, nScales=%d, nTime=%d, SNR=%.1f\n', ...
        K, nScales, nTime, SNR);
    fprintf('  Flat prior (R=1) — no spectral shaping\n');
    fprintf('  Theta estimated memory: %.2f GB\n', mem_GB);
end
if mem_GB > 2
    warning('bst_eigenmode_cwt_inverse:LargeArray', ...
        'Theta will occupy %.1f GB — consider reducing K or nScales.', mem_GB);
end

%% ── Build MNE kernel (flat prior) ────────────────────────────────────────
% sP = ones(K,1), so Lws = L_tilde and K_mne = V·diag(α)·U' (no diag(sP) wrapper).
[U, S, V] = svd(L_tilde, 'econ');
s    = diag(S);
Lam  = sum(s.^2) / (nHMCh * SNR^2);
alph = s ./ (s.^2 + Lam);
K_mne = V * diag(alph) * U';           % [K x nCh]

%% ── Apply kernel to each CWT scale ───────────────────────────────────────
Theta = zeros(K, nScales, nTime, 'like', 1i);
for si = 1:nScales
    W_s = squeeze(W_sensor(:, si, :));  % [nCh x nTime] complex
    Theta(:, si, :) = K_mne * W_s;     % [K x nTime]  complex
end

if Verbose
    fprintf('  |Θ| mean=%.3e  max=%.3e\n', mean(abs(Theta(:))), max(abs(Theta(:))));
end

%% ── Face-space reconstruction ────────────────────────────────────────────
FaceAmp    = [];
FaceIndices = [];
if FaceSpace && isfield(HM,'PhiFaces') && ~isempty(HM.PhiFaces)
    Phi_f   = double(HM.PhiFaces);     % [nLHF x K]
    nLHF    = size(Phi_f, 1);
    FaceAmp = zeros(nLHF, nScales, nTime, 'like', 1i);
    for si = 1:nScales
        FaceAmp(:, si, :) = Phi_f * squeeze(Theta(:, si, :));
    end
    FaceIndices = HM.FaceIndices;
    if Verbose
        fprintf('  FaceAmp: [%d faces x %d scales x %d]\n', nLHF, nScales, nTime);
    end
end

%% ── Pack InvInfo ─────────────────────────────────────────────────────────
InvInfo.K             = K;
InvInfo.lambdas       = lambdas;
InvInfo.nScales       = nScales;
InvInfo.f_scales      = WTInfo.f_scales(:);
InvInfo.SNR           = SNR;
InvInfo.FaceIndices   = FaceIndices;
InvInfo.HeadModelFile = HeadModelFile;

if Verbose
    fprintf('  Done.\n');
end
end
