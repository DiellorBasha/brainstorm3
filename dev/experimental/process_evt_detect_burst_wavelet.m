function varargout = process_evt_detect_burst_wavelet( varargin )
% PROCESS_EVT_DETECT_BURST_WAVELET: Detect oscillatory bursts using Morlet wavelet SNR.
%
% Detects transient narrowband oscillatory events by computing the ratio of
% Morlet wavelet power in a target frequency band to wavelet power in a
% broadband reference (everything outside the target band).
%
% This is conceptually equivalent to the CFAR (Constant False Alarm Rate)
% approach of Muller et al. (2016, eLife 5:e17267), but replaces Butterworth
% filtering + sliding-window power with Morlet wavelet coefficients. The
% wavelet naturally provides the analytic signal envelope at each frequency,
% and the ratio across bands yields the same amplitude-invariant SNR metric.
%
% Architecture: single broadband CWT decomposition via cwtfilterbank,
% then partition wavelet coefficients into target vs reference frequencies.
% Using one CWT for both ensures consistent power normalization.
%
% Advantages over the Butterworth/CFAR approach:
%   - Time-frequency resolution adapts to frequency (shorter windows for
%     faster oscillations, longer for slower ones).
%   - TimeBandwidth (Morse) controls time-frequency tradeoff directly.
%   - Can detect bursts that drift in peak frequency within the band.
%   - Uses MATLAB's optimized cwtfilterbank for fast computation.
%   - No filter-design stability issues (no SOS/Butterworth needed).
%
% Requires: MATLAB Wavelet Toolbox.
%
% Algorithm:
%   1. Compute CWT (cwtfilterbank) across full broadband range.
%   2. Average wavelet power over target-band frequencies -> P_target(t).
%   3. Average wavelet power over remaining frequencies  -> P_ref(t).
%   4. SNR(t) = 10*log10(P_target / P_ref) per channel.
%      (No PSD normalization needed — CWT power is already per-frequency.)
%   5. Per-channel thresholding at constant SNR level.
%   6. Consensus voting: require N% of channels to exceed threshold.
%   7. Boundary extension to earliest onset / latest offset.
%   8. Merge nearby events and enforce minimum duration.
%
% USAGE:  OutputFiles = process_evt_detect_burst_wavelet('Run', sProcess, sInputs)
%                 evt = process_evt_detect_burst_wavelet('Compute', F, TimeVector, OPTIONS)
%             OPTIONS = process_evt_detect_burst_wavelet('Compute')  % Get default options
%
% Authors: Diellor Basha, 2026

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

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
    % Description
    sProcess.Comment     = 'Detect bursts (Wavelet)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Events';
    sProcess.Index       = 48;
    sProcess.Description = 'https://doi.org/10.7554/eLife.17267';
    % Input/Output
    sProcess.InputTypes  = {'raw', 'data'};
    sProcess.OutputTypes = {'raw', 'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % === EVENT NAME ===
    sProcess.options.eventname.Comment = 'Event name: ';
    sProcess.options.eventname.Type    = 'text';
    sProcess.options.eventname.Value   = 'burst_alpha';
    % Separator
    sProcess.options.sep1.Type    = 'separator';
    sProcess.options.sep1.Comment = ' ';

    % === SENSOR SELECTION ===
    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG';
    sProcess.options.sensortypes.InputTypes = {'data', 'raw'};

    % === TIME WINDOW ===
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];

    % Separator
    sProcess.options.sep2.Type    = 'separator';
    sProcess.options.sep2.Comment = ' ';

    % === TARGET FREQUENCY BAND ===
    sProcess.options.freqband.Comment = {'Delta (2-4 Hz)', 'Theta (4-8 Hz)', ...
        'Alpha (8-13 Hz)', 'Beta (13-30 Hz)', 'Gamma (30-60 Hz)', 'Custom (below)', ...
        '<B>Target frequency band:</B>'; ...
        'delta', 'theta', 'alpha', 'beta', 'gamma', 'custom', ''};
    sProcess.options.freqband.Type    = 'radio_linelabel';
    sProcess.options.freqband.Value   = 'alpha';
    % Custom frequency range
    sProcess.options.freqrange.Comment = 'Custom frequency range:';
    sProcess.options.freqrange.Type    = 'freqrange';
    sProcess.options.freqrange.Value   = {[8, 13], 'Hz', []};

    % === WAVELET PARAMETERS ===
    sProcess.options.sep3.Type    = 'separator';
    sProcess.options.sep3.Comment = ' ';
    sProcess.options.label_wavelet.Comment = '<B>Wavelet parameters:</B>';
    sProcess.options.label_wavelet.Type    = 'label';

    % Wavelet type
    sProcess.options.wavelet.Comment = {'Morse (adjustable time-bandwidth)', ...
        'Morlet (analytic, fixed ~6 cycles)', '<B>Wavelet type:</B>'; ...
        'morse', 'amor', ''};
    sProcess.options.wavelet.Type    = 'radio_linelabel';
    sProcess.options.wavelet.Value   = 'morse';

    % Time-bandwidth product (Morse only: higher = better freq resolution, worse time)
    sProcess.options.timebandwidth.Comment = 'Time-bandwidth product (Morse only): ';
    sProcess.options.timebandwidth.Type    = 'value';
    sProcess.options.timebandwidth.Value   = {60, '', 0};

    % Voices per octave (frequency resolution within CWT)
    sProcess.options.voicesperoctave.Comment = 'Voices per octave: ';
    sProcess.options.voicesperoctave.Type    = 'value';
    sProcess.options.voicesperoctave.Value   = {16, '', 0};

    % === DETECTION PARAMETERS ===
    sProcess.options.sep4.Type    = 'separator';
    sProcess.options.sep4.Comment = ' ';
    sProcess.options.label_detect.Comment = '<B>Detection parameters:</B>';
    sProcess.options.label_detect.Type    = 'label';

    % SNR threshold
    sProcess.options.snrthreshold.Comment = 'SNR threshold: ';
    sProcess.options.snrthreshold.Type    = 'value';
    sProcess.options.snrthreshold.Value   = {14, 'dB', 1};

    % Broadband reference range
    sProcess.options.broadband.Comment = 'Broadband reference range: ';
    sProcess.options.broadband.Type    = 'freqrange';
    sProcess.options.broadband.Value   = {[1, 100], 'Hz', []};

    % NOTE: No PSD normalization option — the CWT produces per-frequency power,
    % so both target and reference are already on the same PSD-like scale.
    % The bandwidth asymmetry that requires PSD normalization in the CFAR
    % process does not apply here.

    % Resample
    sProcess.options.resample.Comment = 'Downsample to (0=off): ';
    sProcess.options.resample.Type    = 'value';
    sProcess.options.resample.Value   = {600, 'Hz', 0};

    % === CONSENSUS / EVENT PARAMETERS ===
    sProcess.options.sep5.Type    = 'separator';
    sProcess.options.sep5.Comment = ' ';
    sProcess.options.label_consensus.Comment = '<B>Consensus and event parameters:</B>';
    sProcess.options.label_consensus.Type    = 'label';

    % Minimum channels percentage
    sProcess.options.minchannels.Comment = 'Min channels for consensus: ';
    sProcess.options.minchannels.Type    = 'value';
    sProcess.options.minchannels.Value   = {5, '%', 0};

    % Minimum burst duration
    sProcess.options.minduration.Comment = 'Min burst duration: ';
    sProcess.options.minduration.Type    = 'value';
    sProcess.options.minduration.Value   = {0.050, 'ms', []};

    % Minimum gap between bursts (merge if closer)
    sProcess.options.mingap.Comment = 'Min gap between bursts (merge if shorter): ';
    sProcess.options.mingap.Type    = 'value';
    sProcess.options.mingap.Value   = {0.050, 'ms', []};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = ['Detect bursts (Wavelet): ', sProcess.options.eventname.Value];
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs)
    % ===== GET OPTIONS =====
    evtName = strtrim(sProcess.options.eventname.Value);
    if isempty(evtName)
        bst_report('Error', sProcess, [], 'Event name must be specified.');
        OutputFiles = {};
        return;
    end

    % Sensor types
    SensorTypes = strtrim(sProcess.options.sensortypes.Value);

    % Frequency band
    freqBands = struct( ...
        'delta', [2 4], 'theta', [4 8], 'alpha', [8 13], ...
        'beta', [13 30], 'gamma', [30 60]);
    bandVal = sProcess.options.freqband.Value;
    if isfield(freqBands, bandVal)
        FreqRange = freqBands.(bandVal);
    else
        FreqRange = sProcess.options.freqrange.Value{1};
    end
    if isempty(FreqRange) || (length(FreqRange) ~= 2)
        bst_report('Error', sProcess, [], 'Invalid frequency range.');
        OutputFiles = {};
        return;
    end

    % Time window
    if isfield(sProcess.options, 'timewindow') && isfield(sProcess.options.timewindow, 'Value') && ...
            iscell(sProcess.options.timewindow.Value) && ~isempty(sProcess.options.timewindow.Value)
        TimeWindow = sProcess.options.timewindow.Value{1};
    else
        TimeWindow = [];
    end

    % Build OPTIONS structure for Compute
    OPTIONS = Compute();  % Get defaults
    OPTIONS.freqRange      = FreqRange;
    OPTIONS.snrThreshold   = sProcess.options.snrthreshold.Value{1};
    OPTIONS.wavelet        = sProcess.options.wavelet.Value;
    OPTIONS.timeBandwidth  = sProcess.options.timebandwidth.Value{1};
    OPTIONS.voicesPerOctave = sProcess.options.voicesperoctave.Value{1};
    OPTIONS.broadband      = sProcess.options.broadband.Value{1};
    OPTIONS.resample       = sProcess.options.resample.Value{1};
    OPTIONS.minChannelsPct = sProcess.options.minchannels.Value{1};
    OPTIONS.minDuration    = sProcess.options.minduration.Value{1};
    OPTIONS.minGap         = sProcess.options.mingap.Value{1};

    % Option structure for in_fread()
    ImportOptions = db_template('ImportOptions');
    ImportOptions.ImportMode      = 'Time';
    ImportOptions.UseCtfComp      = 1;
    ImportOptions.UseSsp          = 1;
    ImportOptions.EventsMode      = 'ignore';
    ImportOptions.DisplayMessages = 0;
    ImportOptions.RemoveBaseline  = 'no';

    % Get current progressbar position
    progressPos = bst_progress('get');

    % ===== PROCESS EACH FILE =====
    iOk = false(1, length(sInputs));
    for iFile = 1:length(sInputs)
        % ===== LOAD DATA =====
        bst_progress('text', sprintf('File %d/%d: Reading data...', iFile, length(sInputs)));
        bst_progress('set', progressPos + round(iFile / length(sInputs) / 3 * 100));

        % Load file descriptor
        isRaw = strcmpi(sInputs(iFile).FileType, 'raw');
        if isRaw
            DataMat = in_bst_data(sInputs(iFile).FileName, 'F', 'Time');
            sFile = DataMat.F;
        else
            DataMat = in_bst_data(sInputs(iFile).FileName, 'Time');
            sFile = in_fopen(sInputs(iFile).FileName, 'BST-DATA');
        end
        % Load channel file
        ChannelMat = in_bst_channel(sInputs(iFile).ChannelFile);

        % Must be continuous
        if ~isempty(sFile.epochs)
            bst_report('Error', sProcess, sInputs(iFile), ...
                'This function can only process continuous recordings (no epochs).');
            continue;
        end

        % Get channel indices for the requested sensor types
        iChannels = channel_find(ChannelMat.Channel, SensorTypes);
        if isempty(iChannels)
            bst_report('Error', sProcess, sInputs(iFile), ...
                ['No channels found matching "' SensorTypes '".']);
            continue;
        end

        % Compute sample bounds
        if ~isempty(TimeWindow)
            SamplesBounds = round(sFile.prop.times(1) .* sFile.prop.sfreq) + ...
                            bst_closest(TimeWindow, DataMat.Time) - 1;
        else
            SamplesBounds = [];
        end

        % Read data
        [F, TimeVector] = in_fread(sFile, ChannelMat, 1, SamplesBounds, iChannels, ImportOptions);
        if isempty(F) || (length(TimeVector) < 2)
            bst_report('Error', sProcess, sInputs(iFile), 'Time window is not valid.');
            continue;
        end

        % ===== DETECT BURSTS =====
        bst_progress('text', sprintf('File %d/%d: Wavelet burst detection...', iFile, length(sInputs)));
        bst_progress('set', progressPos + round(2 * iFile / length(sInputs) / 3 * 100));

        [evt, stats] = Compute(F, TimeVector, OPTIONS);

        % ===== SAVE EVENTS =====
        bst_progress('text', sprintf('File %d/%d: Saving events...', iFile, length(sInputs)));
        bst_progress('set', progressPos + round(3 * iFile / length(sInputs) / 3 * 100));

        if ~isempty(evt) && size(evt, 2) > 0
            % Initialize events structure if needed
            if ~isfield(sFile, 'events') || isempty(sFile.events)
                sFile.events = repmat(db_template('event'), 0);
            end
            % Find or create event category
            iEvt = find(strcmpi({sFile.events.label}, evtName));
            if isempty(iEvt)
                iEvt = length(sFile.events) + 1;
                sEvent = db_template('event');
                sEvent.label = evtName;
                sEvent.color = panel_record('GetNewEventColor', iEvt, sFile.events);
            else
                sEvent = sFile.events(iEvt);
            end
            % Set event times (extended events: 2 x nEvents)
            sEvent.times    = evt;
            sEvent.epochs   = ones(1, size(evt, 2));
            sEvent.channels = [];
            sEvent.notes    = [];
            sFile.events(iEvt) = sEvent;

            % Save back to file
            if isRaw
                DataMat.F = sFile;
            else
                DataMat.Events = sFile.events;
            end
            DataMat = rmfield(DataMat, 'Time');
            bst_save(file_fullpath(sInputs(iFile).FileName), DataMat, 'v6', 1);

            % Report
            nBursts = size(evt, 2);
            meanDur = mean(evt(2,:) - evt(1,:)) * 1000;
            bst_report('Info', sProcess, sInputs(iFile), ...
                sprintf('[Wavelet-%s] %d bursts (mean dur: %.0f ms, SNR thr: %.1f dB, peak: %.1f dB, %d target / %d ref bins, consensus: %d/%d ch)', ...
                OPTIONS.wavelet, nBursts, meanDur, stats.snrThresholdUsed, stats.maxSNR, ...
                stats.nTargetBins, stats.nRefBins, stats.minChannelsUsed, stats.nChannels));
        else
            bst_report('Warning', sProcess, sInputs(iFile), ...
                'No bursts detected. Try lowering the threshold, widening the frequency band, or reducing the consensus percentage.');
        end
        iOk(iFile) = true;
    end

    % Return processed files
    OutputFiles = {sInputs(iOk).FileName};
end


%% ===== COMPUTE =====
% USAGE:  [evt, stats] = Compute(F, TimeVector, OPTIONS)
%                OPTIONS = Compute()    % Return default options
function [evt, stats] = Compute(F, TimeVector, OPTIONS)
    % Default options
    defOptions = struct(...
        'freqRange',      [8, 13], ...   % Target frequency band [low, high] in Hz
        'snrThreshold',   14, ...        % Threshold in dB (higher than CFAR due to log-freq averaging)
        'wavelet',        'morse', ...   % Wavelet type: 'morse' or 'amor'
        'timeBandwidth',  60, ...        % Morse time-bandwidth product (ignored for amor)
        'voicesPerOctave', 16, ...       % Frequency resolution in CWT
        'broadband',      [1, 100], ...  % Broadband reference range [low, high] in Hz
        'resample',       600, ...       % Downsample before wavelet (0=off)
        'minChannelsPct', 5, ...         % Min % of channels for consensus
        'minDuration',    0.050, ...     % Min burst duration in seconds
        'minGap',         0.050);        % Min gap between bursts in seconds

    % Return defaults if no input
    if (nargin == 0)
        evt = defOptions;
        stats = [];
        return;
    end
    OPTIONS = struct_copy_fields(OPTIONS, defOptions, 0);

    % Sampling frequency
    sFreq = 1 / (TimeVector(2) - TimeVector(1));
    nSamplesOrig = length(TimeVector);

    % ===== COMPUTE PER-CHANNEL SNR MASKS =====
    [chanMasks, stats] = ComputeWaveletSNR(F, sFreq, nSamplesOrig, OPTIONS);
    % chanMasks: nChannels x nSamplesOrig (binary, at original sample rate)

    nChannels = size(chanMasks, 1);

    % ===== CONSENSUS VOTING =====
    chanCount = sum(chanMasks, 1);
    minChanAbs = max(1, round(OPTIONS.minChannelsPct / 100 * nChannels));
    consensusMask = chanCount >= minChanAbs;

    stats.minChannelsUsed = minChanAbs;
    stats.nChannels = nChannels;

    % Find consensus intervals
    diffCons = diff([0, consensusMask, 0]);
    consOnsets  = find(diffCons == 1);
    consOffsets = find(diffCons == -1) - 1;

    if isempty(consOnsets)
        stats.nRawBursts = 0;
        evt = [];
        return;
    end

    % ===== EXTEND BOUNDARIES TO EARLIEST/LATEST CHANNEL CROSSING =====
    onsets  = zeros(1, length(consOnsets));
    offsets = zeros(1, length(consOnsets));
    for iEvt = 1:length(consOnsets)
        cOn  = consOnsets(iEvt);
        cOff = consOffsets(iEvt);

        % Channels active during this consensus period
        activeChan = find(any(chanMasks(:, cOn:cOff), 2));

        % Extend to earliest onset / latest offset across active channels
        evtStart = cOn;
        evtEnd   = cOff;
        for iCh = 1:length(activeChan)
            ch = activeChan(iCh);
            s = cOn;
            while s > 1 && chanMasks(ch, s-1)
                s = s - 1;
            end
            e = cOff;
            while e < nSamplesOrig && chanMasks(ch, e+1)
                e = e + 1;
            end
            evtStart = min(evtStart, s);
            evtEnd   = max(evtEnd, e);
        end
        onsets(iEvt)  = evtStart;
        offsets(iEvt) = evtEnd;
    end

    % ===== MERGE CLOSE / OVERLAPPING BURSTS =====
    minGapSamples = round(OPTIONS.minGap * sFreq);
    iBurst = 1;
    while iBurst < length(onsets)
        if (onsets(iBurst + 1) - offsets(iBurst)) < minGapSamples
            offsets(iBurst) = max(offsets(iBurst), offsets(iBurst + 1));
            onsets(iBurst + 1)  = [];
            offsets(iBurst + 1) = [];
        else
            iBurst = iBurst + 1;
        end
    end

    % ===== ENFORCE MINIMUM DURATION =====
    minDurSamples = round(OPTIONS.minDuration * sFreq);
    durations = offsets - onsets + 1;
    keepMask = durations >= minDurSamples;
    onsets  = onsets(keepMask);
    offsets = offsets(keepMask);

    % ===== BUILD EVENT TIMES =====
    if isempty(onsets)
        evt = [];
    else
        evt = [TimeVector(onsets); TimeVector(offsets)];
    end
    stats.nRawBursts = length(onsets);
end


%% ===== COMPUTE WAVELET SNR =====
% Full-CWT power ratio: single broadband cwtfilterbank decomposition,
% then split wavelet coefficients into target vs reference frequencies.
%
% Using one CWT for both target and reference guarantees consistent power
% normalization — no scale mismatch between wavelet and filter-based power.
%
% For each channel and time point:
%   P_target(t) = mean( |W(f,t)|^2 ) for f in target band
%   P_ref(t)    = mean( |W(f,t)|^2 ) for f outside target band
%   SNR(t)      = 10*log10( P_target / P_ref )
%
% Note: CWT bins are log-spaced, so this is a geometric-frequency average.
% Typical SNR values are higher than the CFAR process (default threshold: 14 dB).
%
% Requires: Wavelet Toolbox (cwtfilterbank).
function [chanMasks, stats] = ComputeWaveletSNR(F, sFreq, nSamplesOrig, OPTIONS)
    nChannelsOrig = size(F, 1);

    % --- Check for Wavelet Toolbox ---
    if isempty(which('cwtfilterbank'))
        error(['Wavelet Toolbox is required for the wavelet burst detector. ' ...
               'Use "Detect bursts (CFAR)" as an alternative.']);
    end

    % --- Optional downsampling for speed ---
    targetFs = OPTIONS.resample;
    didResample = false;
    if targetFs > 0 && targetFs < sFreq
        if targetFs <= 2 * OPTIONS.broadband(2)
            error('Resample rate (%.0f Hz) must be > 2x broadband upper limit (%.0f Hz).', ...
                targetFs, OPTIONS.broadband(2));
        end
        [p, q] = rat(targetFs / sFreq, 1e-6);
        nSamplesNew = round(nSamplesOrig * p / q);
        Frs = zeros(nChannelsOrig, nSamplesNew);
        for iChan = 1:nChannelsOrig
            Frs(iChan, :) = resample(F(iChan, :), p, q);
        end
        F = Frs;
        sFreq = targetFs;
        didResample = true;
        clear Frs;
    end

    nChannels = size(F, 1);
    nSamples  = size(F, 2);

    % --- Frequency parameters ---
    fLow  = OPTIONS.freqRange(1);
    fHigh = OPTIONS.freqRange(2);
    bLow  = OPTIONS.broadband(1);
    bHigh = OPTIONS.broadband(2);

    % === BUILD SINGLE BROADBAND CWT FILTER BANK ===
    waveletType = OPTIONS.wavelet;
    vpo = OPTIONS.voicesPerOctave;
    fbArgs = {'SignalLength', nSamples, 'SamplingFrequency', sFreq, ...
              'FrequencyLimits', [bLow, bHigh], 'VoicesPerOctave', vpo};
    if strcmpi(waveletType, 'morse')
        fb = cwtfilterbank(fbArgs{:}, 'Wavelet', 'morse', ...
            'TimeBandwidth', OPTIONS.timeBandwidth);
    else
        fb = cwtfilterbank(fbArgs{:}, 'Wavelet', 'amor');
    end

    % Get center frequencies and partition into target vs reference
    freqsCWT = centerFrequencies(fb);
    isTarget = (freqsCWT >= fLow) & (freqsCWT <= fHigh);
    isRef    = ~isTarget;
    nTargetBins = sum(isTarget);
    nRefBins    = sum(isRef);

    if nTargetBins < 1
        error('No CWT frequency bins in target band [%.1f, %.1f] Hz. Increase VoicesPerOctave.', fLow, fHigh);
    end
    if nRefBins < 1
        error('No CWT frequency bins outside the target band. Widen the broadband range.');
    end

    % === COMPUTE CWT POWER PER CHANNEL ===
    % Average wavelet power across target and reference frequency bins.
    % CWT bins are log-spaced (constant voices per octave), so this is
    % effectively a geometric-frequency average.  The dB threshold should
    % be calibrated accordingly (CWT SNR values are typically higher than
    % the CFAR process at the same bandwidth settings).
    powTarget = zeros(nChannels, nSamples);
    powRef    = zeros(nChannels, nSamples);

    for iChan = 1:nChannels
        coeffs = wt(fb, F(iChan, :));          % [nFreqs x nSamples], complex
        powAll = abs(coeffs).^2;
        powTarget(iChan, :) = mean(powAll(isTarget, :), 1);
        powRef(iChan, :)    = mean(powAll(isRef, :), 1);
    end

    % --- SNR in dB per channel ---
    snrPerChan = -Inf(nChannels, nSamples);
    validMask = (powTarget > 0) & (powRef > 0);
    snrPerChan(validMask) = 10 * log10(powTarget(validMask) ./ powRef(validMask));

    % Per-channel thresholding
    chanMasksDS = snrPerChan >= OPTIONS.snrThreshold;

    % --- Upsample per-channel masks to original sample count ---
    if didResample
        idxOrig = round(linspace(1, nSamples, nSamplesOrig));
        chanMasks = chanMasksDS(:, idxOrig);
    else
        chanMasks = chanMasksDS;
    end

    % Stats (aggregate for reporting)
    maxSNRperTime = max(snrPerChan, [], 1);
    stats = struct(...
        'snrThresholdUsed', OPTIONS.snrThreshold, ...
        'waveletType',      waveletType, ...
        'nTargetBins',      nTargetBins, ...
        'nRefBins',         nRefBins, ...
        'targetFreqs',      freqsCWT(isTarget)', ...
        'meanSNR',          mean(maxSNRperTime(isfinite(maxSNRperTime))), ...
        'maxSNR',           max(maxSNRperTime(isfinite(maxSNRperTime))), ...
        'bwTarget',         fHigh - fLow, ...
        'bwRef',            (bHigh - bLow) - (fHigh - fLow), ...
        'actualSFreq',      sFreq);
end
