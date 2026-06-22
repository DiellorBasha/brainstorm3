% CWT_SINGLE_FREQ  Single-frequency (~9.5 Hz) CWT source mapping via cwtfilterbank.
%
% Pipeline:
%   1. cwtfilterbank over [7 13] Hz (matches the alpha bandpass of the data).
%   2. Pick the filter nearest 9.5 Hz; extract its complex coefficient per sensor.
%   3. Sensor-space amplitude = abs(coeff), phase = angle(coeff)  [nCh x nTime].
%   4. Source-map the COMPLEX coefficient through the face-eigenmode MNE inverse
%      (linear), then take abs (amplitude map) and angle (phase map) on faces.
%      Eigenmode -> face reconstruction is windowed in time to bound memory.
%
% Output: dev/tests/cwt_single_freq_result.mat
%
% Authors: Diellor Basha, 2026

addpath('/Users/diellorbasha/workspace/research/code/brainstorm3');
if ~brainstorm('status'), brainstorm server; end
bst_plugin('Load','nxr-compute');

HM   = 'Subject01/S01_AEF_20131218_01_notch/headmodel_face_eigenmode_260606_1131.mat';
DATA = 'Subject01/S01_AEF_20131218_02_notch/data_block001_band.mat';
SURF = 'Subject01/tess_cortex_pial_low.mat';
TARGET_F = 9.5;            % Hz, pick the filter nearest this
FLIMITS  = [7 13];         % Hz, alpha passband of the data
SNR      = 3;
WIN      = 2000;           % time-window length for face reconstruction

%% ── 1. Load sensor data ──────────────────────────────────────────────────
DataMat = in_bst_data(DATA);
Time = DataMat.Time(:)';  Fs = 1/mean(diff(Time));  nTime = numel(Time);
[~, iStudy] = bst_get('DataFile', DATA);
ChanMat = in_bst_channel(bst_get('ChannelFileForStudy', iStudy));
iCh = find(strcmpi({ChanMat.Channel.Type}, 'MEG'));
d   = double(DataMat.F(iCh, :));   % [nCh x nTime]
nCh = numel(iCh);
fprintf('Data: %d MEG channels, %d samples, Fs=%.1f Hz\n', nCh, nTime, Fs);

%% ── 2. cwtfilterbank over [7 13] Hz ──────────────────────────────────────
fb     = cwtfilterbank('SignalLength', nTime, 'SamplingFrequency', Fs, ...
                       'FrequencyLimits', FLIMITS);
freqs  = centerFrequencies(fb);          % [nScales x 1], descending
[~, idx] = min(abs(freqs - TARGET_F));
f_pick = freqs(idx);
fprintf('cwtfilterbank: %d filters in [%.1f %.1f] Hz; picked filter %d @ %.2f Hz\n', ...
        numel(freqs), FLIMITS(1), FLIMITS(2), idx, f_pick);

%% ── 3. Extract single-frequency complex coefficient per sensor ───────────
C = zeros(nCh, nTime, 'like', 1i);
for ch = 1:nCh
    cf = wt(fb, d(ch,:));     % [nScales x nTime] complex
    C(ch,:) = cf(idx, :);     % single frequency row
end
A_sens = abs(C);              % sensor amplitude  [nCh x nTime]
P_sens = angle(C);            % sensor phase      [nCh x nTime]
fprintf('Sensor coeff @ %.2f Hz: |C| mean=%.3e max=%.3e\n', ...
        f_pick, mean(A_sens(:)), max(A_sens(:)));

%% ── 4a. Eigenmode MNE inverse on the complex coefficient ─────────────────
% Reuse the tested kernel by presenting C as a single-scale CWT tensor.
W3 = reshape(C, [nCh, 1, nTime]);
WTInfo = struct('Fs',Fs, 'f_scales',f_pick, 'nScales',1, 'nTime',nTime, ...
                'nCh',nCh, 'ChanIdx',iCh(:));
[Theta, ~, Inv] = bst_eigenmode_cwt_inverse(W3, WTInfo, HM, ...
                    'SNR', SNR, 'FaceSpace', false, 'Verbose', true);
Theta = squeeze(Theta);       % [K x nTime] complex
K = size(Theta,1);

%% ── 4b. Windowed eigenmode -> face reconstruction; abs/angle on faces ────
Hm     = in_bst_headmodel(HM, 0);
PhiF   = double(Hm.PhiFaces);          % [nLHF x K]
FaceIdx= Hm.FaceIndices(:);            % global face indices of LH faces
nLHF   = size(PhiF, 1);

amp_sum   = zeros(nLHF, 1);            % accumulate |s_face| over time -> mean amp map
amp_peak  = zeros(nLHF, 1);           % |s_face| at global-peak time
gpow      = zeros(1, nTime);          % global source power over time
phase_at_peak = zeros(nLHF, 1);
% first pass: global power to find peak time
for w0 = 1:WIN:nTime
    w1 = min(w0+WIN-1, nTime);
    Sfw = PhiF * Theta(:, w0:w1);      % [nLHF x win] complex
    gpow(w0:w1) = sum(abs(Sfw).^2, 1);
    amp_sum = amp_sum + sum(abs(Sfw), 2);
end
amp_mean = amp_sum / nTime;
[~, t_peak] = max(gpow);
% second pass: grab the single peak-time column for amplitude + phase snapshot
w0 = max(1, t_peak);
Sf_peak = PhiF * Theta(:, w0);         % [nLHF x 1] complex
amp_peak = abs(Sf_peak);
phase_at_peak = angle(Sf_peak);
fprintf('Faces: %d LH faces, K=%d. Peak time idx=%d (%.3f s)\n', ...
        nLHF, K, t_peak, Time(t_peak));

%% ── Source amplitude peak location ───────────────────────────────────────
Tess = in_tess_bst(SURF, 0);
Vtx  = Tess.Vertices;  Fc = double(Tess.Faces);
ctr  = (Vtx(Fc(:,1),:)+Vtx(Fc(:,2),:)+Vtx(Fc(:,3),:))/3;   % [nFaces x 3]
[~, iMaxMean] = max(amp_mean);
[~, iMaxPeak] = max(amp_peak);
src_mean_face = FaceIdx(iMaxMean);  src_mean_pos = ctr(src_mean_face,:);
src_peak_face = FaceIdx(iMaxPeak);  src_peak_pos = ctr(src_peak_face,:);

fprintf('\n=== SINGLE-FREQ (%.2f Hz) SOURCE MAP ===\n', f_pick);
fprintf('Amplitude (time-mean) peak: face %d  [%.1f %.1f %.1f] mm\n', ...
        src_mean_face, src_mean_pos*1000);
fprintf('Amplitude (at t_peak)  peak: face %d  [%.1f %.1f %.1f] mm\n', ...
        src_peak_face, src_peak_pos*1000);
fprintf('Source |s_face|: mean=%.3e  peak-time max=%.3e\n', ...
        mean(amp_mean), max(amp_peak));

%% ── Save compact result ──────────────────────────────────────────────────
R = struct();
R.f_pick      = f_pick;       R.FLIMITS = FLIMITS;   R.freqs = freqs;
R.Fs          = Fs;           R.Time = Time;         R.FaceIndices = FaceIdx;
R.A_sens_mean = mean(A_sens,2);   R.P_sens_peak = P_sens(:, t_peak);
R.amp_mean    = amp_mean;     R.amp_peak = amp_peak; R.phase_at_peak = phase_at_peak;
R.gpow        = gpow;         R.t_peak = t_peak;
R.src_mean_face = src_mean_face;  R.src_mean_pos = src_mean_pos;
R.src_peak_face = src_peak_face;  R.src_peak_pos = src_peak_pos;
R.Theta       = single(Theta);    % [K x nTime] complex, for later reconstruction
outFile = '/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/cwt_single_freq_result.mat';
save(outFile, 'R', '-v7.3');
fprintf('Saved -> %s\n', outFile);
