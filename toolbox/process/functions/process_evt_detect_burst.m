function varargout = process_evt_detect_burst( varargin )
% PROCESS_EVT_DETECT_BURST: Detect transient oscillatory bursts in continuous recordings.
%
% Detects time intervals where narrowband oscillatory power exceeds a threshold.
% Two detection methods are available:
%
% Method 1: Hilbert envelope + MAD threshold
%   1. Band-pass filter the selected channel(s) in the target frequency band.
%   2. Compute the analytic signal via Hilbert transform; take the envelope.
%   3. Compute threshold = median(envelope) + k * MAD(envelope), where k is
%      user-specified (default 2). MAD-based thresholding is robust to outliers.
%   4. Mark contiguous supra-threshold intervals as burst events.
%
% Method 2: SNR / CFAR (Muller et al. 2016, eLife 5:e17267)
%   1. Bandpass filter in target band (8th-order Butterworth by default).
%   2. Bandstop filter: broadband (1-100 Hz) with target band removed.
%   3. Compute power of bandpass and bandstop in sliding windows (500 ms).
%   4. SNR(t) = 10*log10(P_bandpass / P_bandstop) in dB per channel.
%   5. Threshold at constant SNR level (default 5 dB) — CFAR-like detection.
%   This method is amplitude-invariant: it detects narrowband spectral events
%   regardless of overall signal amplitude.
%
% Both methods then:
%   - Enforce minimum burst duration (reject brief crossings).
%   - Merge nearby events separated by less than the minimum gap.
%
% USAGE:  OutputFiles = process_evt_detect_burst('Run', sProcess, sInputs)
%                 evt = process_evt_detect_burst('Compute', F, TimeVector, OPTIONS)
%             OPTIONS = process_evt_detect_burst('Compute')    % Get default options
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
function sProcess = GetDescription() %#ok<DEFNU>
    % Description
    sProcess.Comment     = 'Detect oscillatory bursts';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Events';
    sProcess.Index       = 46;
    sProcess.Description = '';
    % Input/Output
    sProcess.InputTypes  = {'raw', 'data'};
    sProcess.OutputTypes = {'raw', 'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % === EVENT NAME ===
    sProcess.options.eventname.Comment = 'Event name: ';
    sProcess.options.eventname.Type    = 'text';
    sProcess.options.eventname.Value   = 'burst_beta';
    % Separator
    sProcess.options.sep1.Type    = 'separator';
    sProcess.options.sep1.Comment = ' ';

    % === SENSOR SELECTION ===
    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG';
    sProcess.options.sensortypes.InputTypes = {'raw', 'data'};
    % Sensor selection help
    sProcess.options.sensorhelp.Comment = ['<I><FONT color="#777777">When multiple channels are selected, the RMS envelope<BR>' ...
                                           'across channels is used for detection.</FONT></I>'];
    sProcess.options.sensorhelp.Type    = 'label';

    % === TIME WINDOW ===
    sProcess.options.timewindow.Comment = 'Time window: ';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];

    % Separator
    sProcess.options.sep2.Type    = 'separator';
    sProcess.options.sep2.Comment = ' ';

    % === FREQUENCY BAND ===
    sProcess.options.label_freq.Comment = '<U><B>Frequency band</B></U>:';
    sProcess.options.label_freq.Type    = 'label';
    sProcess.options.freqband.Comment = {'Delta (2-4 Hz)', 'Theta (4-8 Hz)', 'Alpha (8-13 Hz)', ...
                                         'Beta (13-30 Hz)', 'Low gamma (30-60 Hz)', 'Custom', ''};
    sProcess.options.freqband.Type    = 'radio_line';
    sProcess.options.freqband.Value   = 4;  % Default: Beta
    % Custom frequency range
    sProcess.options.freqrange.Comment = 'Custom frequency range: ';
    sProcess.options.freqrange.Type    = 'range';
    sProcess.options.freqrange.Value   = {[13, 30], 'Hz', 2};

    % Separator
    sProcess.options.sep3.Type    = 'separator';
    sProcess.options.sep3.Comment = ' ';

    % === DETECTION METHOD ===
    sProcess.options.label_method.Comment = '<U><B>Detection method</B></U>:';
    sProcess.options.label_method.Type    = 'label';
    sProcess.options.method.Comment = {'Hilbert envelope (MAD threshold)', ...
                                       'SNR / CFAR (Muller et al. 2016)', ''};
    sProcess.options.method.Type    = 'radio_line';
    sProcess.options.method.Value   = 1;

    % Separator
    sProcess.options.sep4.Type    = 'separator';
    sProcess.options.sep4.Comment = ' ';

    % === HILBERT/MAD OPTIONS ===
    sProcess.options.label_hilbert.Comment = '<U><B>Hilbert envelope parameters</B></U>:';
    sProcess.options.label_hilbert.Type    = 'label';
    % Threshold (in units of MAD above median)
    sProcess.options.threshold.Comment = 'Threshold (median + k*MAD): ';
    sProcess.options.threshold.Type    = 'value';
    sProcess.options.threshold.Value   = {2, ' x MAD', 2};

    % Separator
    sProcess.options.sep5.Type    = 'separator';
    sProcess.options.sep5.Comment = ' ';

    % === SNR/CFAR OPTIONS ===
    sProcess.options.label_snr.Comment = '<U><B>SNR/CFAR parameters</B></U>:';
    sProcess.options.label_snr.Type    = 'label';
    sProcess.options.label_snr_help.Comment = ['<I><FONT color="#777777">Bandpass vs bandstop power ratio in sliding windows.<BR>' ...
                                                'Amplitude-invariant detection (Muller et al. 2016, eLife).</FONT></I>'];
    sProcess.options.label_snr_help.Type    = 'label';
    % SNR threshold in dB
    sProcess.options.snrthreshold.Comment = 'SNR threshold: ';
    sProcess.options.snrthreshold.Type    = 'value';
    sProcess.options.snrthreshold.Value   = {5, 'dB', 1};
    % Sliding window length
    sProcess.options.winlength.Comment = 'Sliding window length: ';
    sProcess.options.winlength.Type    = 'value';
    sProcess.options.winlength.Value   = {0.500, 's', 3};
    % Broadband range for bandstop reference
    sProcess.options.broadband.Comment = 'Broadband range (for bandstop reference): ';
    sProcess.options.broadband.Type    = 'range';
    sProcess.options.broadband.Value   = {[1, 100], 'Hz', 1};
    % Filter order
    sProcess.options.filtorder.Comment = 'Butterworth filter order: ';
    sProcess.options.filtorder.Type    = 'value';
    sProcess.options.filtorder.Value   = {8, '', 0};
    % Resample frequency (0 = no resampling)
    sProcess.options.resample.Comment = 'Resample before filtering (0=off): ';
    sProcess.options.resample.Type    = 'value';
    sProcess.options.resample.Value   = {600, 'Hz', 0};
    sProcess.options.resample_help.Comment = ['<I><FONT color="#777777">Downsample to speed up SNR computation. Must be &gt; 2x broadband upper limit.<BR>' ...
                                               'Recommended: 600 Hz for beta, 300 Hz for alpha/theta. 0 = use original rate.</FONT></I>'];
    sProcess.options.resample_help.Type    = 'label';

    % Separator
    sProcess.options.sep6.Type    = 'separator';
    sProcess.options.sep6.Comment = ' ';

    % === CHANNEL CONSENSUS ===
    sProcess.options.label_consensus.Comment = '<U><B>Channel consensus</B></U>:';
    sProcess.options.label_consensus.Type    = 'label';
    sProcess.options.minchannels.Comment = 'Min channels exceeding threshold (%): ';
    sProcess.options.minchannels.Type    = 'value';
    sProcess.options.minchannels.Value   = {5, '% of channels', 1};
    sProcess.options.consensus_help.Comment = ['<I><FONT color="#777777">Per-channel detection with consensus voting. A burst is called<BR>' ...
                                                'when &ge; N channels exceed threshold simultaneously. Event<BR>' ...
                                                'boundaries extend to the earliest onset / latest offset<BR>' ...
                                                'across all participating channels.</FONT></I>'];
    sProcess.options.consensus_help.Type    = 'label';

    % Separator
    sProcess.options.sep7.Type    = 'separator';
    sProcess.options.sep7.Comment = ' ';

    % === COMMON EVENT PARAMETERS ===
    sProcess.options.label_evt.Comment = '<U><B>Event parameters</B></U>:';
    sProcess.options.label_evt.Type    = 'label';
    % Minimum burst duration
    sProcess.options.minduration.Comment = 'Minimum burst duration: ';
    sProcess.options.minduration.Type    = 'value';
    sProcess.options.minduration.Value   = {0.050, 'ms', []};
    % Minimum gap between bursts (merge if closer)
    sProcess.options.mingap.Comment = 'Min gap between bursts (merge if shorter): ';
    sProcess.options.mingap.Type    = 'value';
    sProcess.options.mingap.Value   = {0.050, 'ms', []};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    methodNames = {'Hilbert/MAD', 'SNR/CFAR'};
    iMethod = sProcess.options.method.Value;
    if iMethod < 1 || iMethod > 2
        iMethod = 1;
    end
    Comment = ['Detect bursts (' methodNames{iMethod} '): ', sProcess.options.eventname.Value];
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
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
    freqBands = {[2 4], [4 8], [8 13], [13 30], [30 60]};
    iBand = sProcess.options.freqband.Value;
    if iBand <= 5
        FreqRange = freqBands{iBand};
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

    % Detection method
    iMethod = sProcess.options.method.Value;
    if iMethod == 1
        methodStr = 'hilbert';
    else
        methodStr = 'snr';
    end

    % Build OPTIONS structure for Compute
    OPTIONS = Compute();  % Get defaults
    OPTIONS.method      = methodStr;
    OPTIONS.freqRange   = FreqRange;
    OPTIONS.minDuration = sProcess.options.minduration.Value{1};
    OPTIONS.minGap      = sProcess.options.mingap.Value{1};
    % Hilbert/MAD options
    OPTIONS.threshold   = sProcess.options.threshold.Value{1};
    % SNR/CFAR options
    OPTIONS.snrThreshold = sProcess.options.snrthreshold.Value{1};
    OPTIONS.winLength    = sProcess.options.winlength.Value{1};
    OPTIONS.broadband    = sProcess.options.broadband.Value{1};
    OPTIONS.filtOrder    = sProcess.options.filtorder.Value{1};
    OPTIONS.resample     = sProcess.options.resample.Value{1};
    % Consensus
    OPTIONS.minChannelsPct = sProcess.options.minchannels.Value{1};

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
        bst_progress('text', sprintf('File %d/%d: Detecting bursts...', iFile, length(sInputs)));
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
            consensusStr = sprintf(', consensus: %d/%d channels', stats.minChannelsUsed, stats.nChannels);
            if strcmpi(OPTIONS.method, 'hilbert')
                bst_report('Info', sProcess, sInputs(iFile), ...
                    sprintf('[Hilbert/MAD] %d bursts (mean dur: %.0f ms, thresh: %.2e%s)', ...
                    nBursts, meanDur, stats.thresholdValue, consensusStr));
            else
                bst_report('Info', sProcess, sInputs(iFile), ...
                    sprintf('[SNR/CFAR] %d bursts (mean dur: %.0f ms, SNR: %.1f dB, peak: %.1f dB%s)', ...
                    nBursts, meanDur, stats.snrThresholdUsed, stats.maxSNR, consensusStr));
            end
        else
            bst_report('Warning', sProcess, sInputs(iFile), ...
                'No bursts detected. Try lowering the threshold or widening the frequency band.');
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
        'method',        'hilbert', ... % 'hilbert' or 'snr'
        'freqRange',     [13, 30], ...  % Frequency band [low, high] in Hz
        'threshold',     2, ...         % Hilbert: threshold in units of MAD above median
        'snrThreshold',  5, ...         % SNR: threshold in dB
        'winLength',     0.500, ...     % SNR: sliding window length in seconds
        'broadband',     [1, 100], ...  % SNR: broadband range for bandstop reference
        'filtOrder',     8, ...         % SNR: Butterworth filter order
        'resample',      600, ...       % SNR: resample to this freq before filtering (0=off)
        'minChannelsPct', 5, ...        % Min % of channels for consensus (0=single-channel max)
        'minDuration',   0.050, ...     % Minimum burst duration in seconds
        'minGap',        0.050);        % Minimum gap between bursts in seconds (merge if shorter)

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

    % ===== DISPATCH BY METHOD (returns per-channel masks) =====
    switch lower(OPTIONS.method)
        case 'hilbert'
            [chanMasks, stats] = ComputeHilbert(F, sFreq, OPTIONS);
        case 'snr'
            [chanMasks, stats] = ComputeSNR(F, sFreq, OPTIONS);
        otherwise
            error('Unknown detection method: %s', OPTIONS.method);
    end
    % chanMasks: nChannels x nSamplesOrig (binary, at original sample rate)

    nChannels = size(chanMasks, 1);

    % ===== CONSENSUS VOTING =====
    % Count how many channels exceed threshold at each time point
    chanCount = sum(chanMasks, 1);  % 1 x nSamplesOrig

    % Minimum channels for consensus
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
    % For each consensus interval, find the participating channels and
    % extend the event to the earliest onset and latest offset among them.
    onsets  = zeros(1, length(consOnsets));
    offsets = zeros(1, length(consOnsets));
    for iEvt = 1:length(consOnsets)
        cOn  = consOnsets(iEvt);
        cOff = consOffsets(iEvt);

        % Channels active during this consensus period
        activeChan = find(any(chanMasks(:, cOn:cOff), 2));

        % For each active channel, find the contiguous supra-threshold
        % interval that overlaps with [cOn, cOff] and take the union
        evtStart = cOn;
        evtEnd   = cOff;
        for iCh = 1:length(activeChan)
            ch = activeChan(iCh);
            % Walk backward from cOn to find where this channel's crossing began
            s = cOn;
            while s > 1 && chanMasks(ch, s-1)
                s = s - 1;
            end
            % Walk forward from cOff to find where this channel's crossing ends
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
    % Extension may cause overlap between adjacent consensus events
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


%% ===== COMPUTE: HILBERT ENVELOPE + MAD =====
% Returns per-channel binary masks (nChannels x nSamples)
function [chanMasks, stats] = ComputeHilbert(F, sFreq, OPTIONS)
    % Band-pass filter
    Ffilt = process_bandpass('Compute', F, sFreq, OPTIONS.freqRange(1), OPTIONS.freqRange(2), ...
                             'bst-hfilter-2019', 1);

    % Hilbert envelope per channel
    nChannels = size(Ffilt, 1);
    envelopes = abs(hilbert(Ffilt'))';  % nChannels x nSamples

    % Per-channel robust threshold: median + k * scaled MAD
    chanMasks = false(size(envelopes));
    medEnvs = zeros(1, nChannels);
    madEnvs = zeros(1, nChannels);
    threshValues = zeros(1, nChannels);
    for iCh = 1:nChannels
        env = envelopes(iCh, :);
        medEnvs(iCh) = median(env);
        madEnvs(iCh) = median(abs(env - medEnvs(iCh)));
        threshValues(iCh) = medEnvs(iCh) + OPTIONS.threshold * 1.4826 * madEnvs(iCh);
        chanMasks(iCh, :) = env > threshValues(iCh);
    end

    stats = struct('method', 'hilbert', ...
                   'medianEnv', median(medEnvs), ...
                   'madEnv', median(madEnvs), ...
                   'thresholdValue', median(threshValues));
end


%% ===== COMPUTE: SNR / CFAR (Muller et al. 2016) =====
% Bandpass vs bandstop power ratio in sliding windows.
% Reference: Muller et al. (2016) eLife 5:e17267
%   SNR(t) = 10*log10( P_bandpass(t) / P_bandstop(t) )
%   Detection at constant SNR threshold (CFAR-like).
function [chanMasks, stats] = ComputeSNR(F, sFreq, OPTIONS)
    nChannelsOrig = size(F, 1);
    nSamplesOrig  = size(F, 2);

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

    % --- Design Butterworth filters (SOS for numerical stability) ---
    fLow  = OPTIONS.freqRange(1);
    fHigh = OPTIONS.freqRange(2);
    bLow  = OPTIONS.broadband(1);
    bHigh = OPTIONS.broadband(2);
    filtOrd = OPTIONS.filtOrder;
    nyq = sFreq / 2;

    [zBP, pBP, kBP] = butter(filtOrd, [fLow, fHigh] / nyq, 'bandpass');
    [sosBP, gBP] = zp2sos(zBP, pBP, kBP);

    [zBroad, pBroad, kBroad] = butter(filtOrd, [bLow, bHigh] / nyq, 'bandpass');
    [sosBroad, gBroad] = zp2sos(zBroad, pBroad, kBroad);
    [zBS, pBS, kBS] = butter(filtOrd, [fLow, fHigh] / nyq, 'stop');
    [sosBS, gBS] = zp2sos(zBS, pBS, kBS);

    % --- Filter each channel ---
    F_bp = zeros(nChannels, nSamples);
    F_bs = zeros(nChannels, nSamples);
    for iChan = 1:nChannels
        F_bp(iChan, :) = filtfilt(sosBP, gBP, F(iChan, :));
        tmp = filtfilt(sosBroad, gBroad, F(iChan, :));
        F_bs(iChan, :) = filtfilt(sosBS, gBS, tmp);
    end

    % --- Sliding window power and SNR (vectorized) ---
    winSamples = round(OPTIONS.winLength * sFreq);
    if winSamples < 2
        winSamples = 2;
    end
    powBP = movmean(F_bp.^2, winSamples, 2);
    powBS = movmean(F_bs.^2, winSamples, 2);

    % SNR in dB per channel
    snrPerChan = -Inf(nChannels, nSamples);
    validMask = (powBP > 0) & (powBS > 0);
    snrPerChan(validMask) = 10 * log10(powBP(validMask) ./ powBS(validMask));

    % Per-channel thresholding → binary masks at downsampled rate
    chanMasksDS = snrPerChan >= OPTIONS.snrThreshold;  % nChannels x nSamples

    % --- Upsample per-channel masks to original sample count ---
    if didResample
        idxOrig = round(linspace(1, nSamples, nSamplesOrig));
        chanMasks = chanMasksDS(:, idxOrig);
    else
        chanMasks = chanMasksDS;
    end

    % Stats (aggregate for reporting only)
    maxSNR = max(snrPerChan, [], 1);
    stats = struct('method', 'snr', ...
                   'snrThresholdUsed', OPTIONS.snrThreshold, ...
                   'meanSNR', mean(maxSNR(isfinite(maxSNR))), ...
                   'maxSNR', max(maxSNR(isfinite(maxSNR))), ...
                   'actualSFreq', sFreq);
end
