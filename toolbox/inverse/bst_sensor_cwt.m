function [W, WTInfo] = bst_sensor_cwt(DataFile, varargin)
% BST_SENSOR_CWT  Continuous wavelet transform of MEG sensor time series.
%
% USAGE:
%   [W, WTInfo] = bst_sensor_cwt(DataFile)
%   [W, WTInfo] = bst_sensor_cwt(DataFile, 'FreqRange', [4 30], 'VoicesPerOctave', 16)
%
% INPUTS:
%   DataFile  Brainstorm relative path to a data file
%
% OPTIONS (name-value):
%   'ChanType'        MEG sensor type string (default: 'MEG')
%   'FreqRange'       [f_lo f_hi] Hz (default: [1 30])
%   'VoicesPerOctave' integer, frequency resolution (default: 12)
%   'WaveletName'     'morse' or 'amor' (Morlet) (default: 'morse')
%   'TimeBandwidth'   Morse wavelet time-bandwidth product (default: 5)
%   'Verbose'         logical (default: true)
%
% OUTPUTS:
%   W       [nCh x nScales x nTime] complex double
%             abs(W)   = instantaneous amplitude envelope at each scale
%             angle(W) = instantaneous phase at each scale
%   WTInfo  struct — Fs, f_scales, Time, nScales, nTime, nCh, ChanIdx,
%                    WaveletName, VoicesPerOctave, FreqRange, DataFile
%
% NOTE: Requires Signal Processing Toolbox (cwtfilterbank).
%       Memory: nCh x nScales x nTime x 16 bytes (complex double).
%       For large datasets, reduce FreqRange or VoicesPerOctave.
%
% Authors: Diellor Basha, 2026

%% ── Toolbox check ────────────────────────────────────────────────────────
if ~license('test', 'signal_toolbox') || isempty(which('cwtfilterbank'))
    error('bst_sensor_cwt:NoSignalToolbox', ...
        'Signal Processing Toolbox is required (cwtfilterbank).');
end

%% ── Parse options ────────────────────────────────────────────────────────
ChanType        = 'MEG';
FreqRange       = [1 30];
VoicesPerOctave = 12;
WaveletName     = 'morse';
TimeBandwidth   = 5;
Verbose         = true;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'chantype',         ChanType        = varargin{k+1};
        case 'freqrange',        FreqRange       = varargin{k+1};
        case 'voicesperoctave',  VoicesPerOctave = varargin{k+1};
        case 'waveletname',      WaveletName     = varargin{k+1};
        case 'timebandwidth',    TimeBandwidth   = varargin{k+1};
        case 'verbose',          Verbose         = logical(varargin{k+1});
    end
end

%% ── Load data ────────────────────────────────────────────────────────────
DataMat = in_bst_data(DataFile);
Time    = DataMat.Time(:)';
Fs      = 1 / mean(diff(Time));
nTime   = numel(Time);

[~, iStudy] = bst_get('DataFile', DataFile);
ChanMat     = in_bst_channel(bst_get('ChannelFileForStudy', iStudy));
iCh         = find(strcmpi({ChanMat.Channel.Type}, ChanType));
if isempty(iCh)
    error('bst_sensor_cwt:NoChan', 'No channels of type ''%s'' found.', ChanType);
end
nCh = numel(iCh);
d   = double(DataMat.F(iCh, :));   % [nCh x nTime]

%% ── CWT per channel ──────────────────────────────────────────────────────
% Use cwt() directly with named args (NOT cwtfilterbank → cwt, which ignores
% FrequencyLimits and returns Nyquist-extended scales).
% [wt, f] = cwt(x, wavelet, Fs, 'FrequencyLimits', ...) respects the range.
cwtArgs = {WaveletName, Fs, 'FrequencyLimits', FreqRange, 'VoicesPerOctave', VoicesPerOctave};
if strcmpi(WaveletName, 'morse')
    cwtArgs = [cwtArgs, {'TimeBandwidth', TimeBandwidth}];
end

% Trial run on first channel to discover scale count and frequency axis
[W1, f_scales] = cwt(d(1,:), cwtArgs{:});   % [nScales x nTime], [nScales x 1]
f_scales = f_scales(:);
nScales  = numel(f_scales);

memBytes = nCh * nScales * nTime * 16;   % complex double = 16 bytes
if memBytes > 2e9
    warning('bst_sensor_cwt:LargeMemory', ...
        'W ~%.1f GB. Reduce FreqRange, VoicesPerOctave, or trim nTime.', ...
        memBytes / 1e9);
end

W = zeros(nCh, nScales, nTime, 'like', 1i);
W(1, :, :) = W1;
for ch = 2:nCh
    W(ch, :, :) = cwt(d(ch,:), cwtArgs{:});
end

%% ── Verbose summary ──────────────────────────────────────────────────────
if Verbose
    fprintf('bst_sensor_cwt: %d channels, %d scales [%.1f–%.1f Hz], %d time points\n', ...
        nCh, nScales, f_scales(end), f_scales(1), nTime);
end

%% ── Pack WTInfo ──────────────────────────────────────────────────────────
WTInfo.Fs              = Fs;
WTInfo.f_scales        = f_scales;
WTInfo.Time            = Time;
WTInfo.nScales         = nScales;
WTInfo.nTime           = nTime;
WTInfo.nCh             = nCh;
WTInfo.ChanIdx         = iCh(:);
WTInfo.WaveletName     = WaveletName;
WTInfo.VoicesPerOctave = VoicesPerOctave;
WTInfo.FreqRange       = FreqRange;
WTInfo.DataFile        = DataFile;
end
