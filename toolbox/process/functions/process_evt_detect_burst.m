function varargout = process_evt_detect_burst( varargin )
% PROCESS_EVT_DETECT_BURST: Detect transient oscillatory bursts in continuous recordings.
%
% Detects time intervals where the analytic amplitude (Hilbert envelope) in a
% specified frequency band exceeds a threshold relative to the baseline
% distribution. Useful for identifying alpha bursts, beta bursts, and other
% transient spectral events in resting-state or task data.
%
% Algorithm:
%   1. Band-pass filter the selected channel(s) in the target frequency band.
%   2. Compute the analytic signal via Hilbert transform; take the envelope.
%   3. Compute threshold = median(envelope) + k * MAD(envelope), where k is
%      user-specified (default 2). MAD-based thresholding is robust to outliers.
%   4. Mark contiguous supra-threshold intervals as burst events.
%   5. Enforce minimum burst duration (reject brief crossings) and minimum
%      gap between bursts (merge nearby events).
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

    % === THRESHOLD ===
    sProcess.options.label_thresh.Comment = '<U><B>Detection parameters</B></U>:';
    sProcess.options.label_thresh.Type    = 'label';
    % Threshold (in units of MAD above median)
    sProcess.options.threshold.Comment = 'Threshold (median + k*MAD): ';
    sProcess.options.threshold.Type    = 'value';
    sProcess.options.threshold.Value   = {2, ' x MAD', 2};
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
    Comment = ['Detect bursts: ', sProcess.options.eventname.Value];
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

    % Build OPTIONS structure for Compute
    OPTIONS = Compute();  % Get defaults
    OPTIONS.freqRange   = FreqRange;
    OPTIONS.threshold   = sProcess.options.threshold.Value{1};
    OPTIONS.minDuration = sProcess.options.minduration.Value{1};
    OPTIONS.minGap      = sProcess.options.mingap.Value{1};

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
            bst_report('Info', sProcess, sInputs(iFile), ...
                sprintf('%d bursts detected (mean duration: %.0f ms, median envelope: %.2e, threshold: %.2e)', ...
                nBursts, meanDur, stats.medianEnv, stats.thresholdValue));
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
        'freqRange',   [13, 30], ...   % Frequency band [low, high] in Hz
        'threshold',   2, ...          % Threshold in units of MAD above median
        'minDuration', 0.050, ...      % Minimum burst duration in seconds
        'minGap',      0.050);         % Minimum gap between bursts in seconds (merge if shorter)

    % Return defaults if no input
    if (nargin == 0)
        evt = defOptions;
        stats = [];
        return;
    end
    OPTIONS = struct_copy_fields(OPTIONS, defOptions, 0);

    % Sampling frequency
    sFreq = 1 / (TimeVector(2) - TimeVector(1));

    % ===== 1. BAND-PASS FILTER =====
    Ffilt = process_bandpass('Compute', F, sFreq, OPTIONS.freqRange(1), OPTIONS.freqRange(2), ...
                             'bst-hfilter-2019', 1);

    % ===== 2. HILBERT ENVELOPE =====
    % If multiple channels, compute envelope per channel then take RMS across channels
    nChannels = size(Ffilt, 1);
    envelopes = abs(hilbert(Ffilt'))';  % hilbert works on columns, so transpose
    if nChannels > 1
        envelope = sqrt(mean(envelopes.^2, 1));  % RMS across channels
    else
        envelope = envelopes;
    end

    % ===== 3. THRESHOLD =====
    % Robust threshold: median + k * MAD
    medEnv = median(envelope);
    madEnv = median(abs(envelope - medEnv));  % MAD = median absolute deviation
    threshValue = medEnv + OPTIONS.threshold * 1.4826 * madEnv;
    % 1.4826 scales MAD to be consistent with std for normal distributions

    % ===== 4. DETECT SUPRA-THRESHOLD INTERVALS =====
    supraMask = envelope > threshValue;

    % Find onset/offset indices
    diffMask = diff([0, supraMask, 0]);
    onsets  = find(diffMask == 1);
    offsets = find(diffMask == -1) - 1;

    if isempty(onsets)
        evt = [];
        stats = struct('medianEnv', medEnv, 'madEnv', madEnv, ...
                        'thresholdValue', threshValue, 'nRawBursts', 0);
        return;
    end

    % ===== 5. MERGE CLOSE BURSTS =====
    minGapSamples = round(OPTIONS.minGap * sFreq);
    iBurst = 1;
    while iBurst < length(onsets)
        if (onsets(iBurst + 1) - offsets(iBurst)) < minGapSamples
            % Merge: extend current burst to end of next
            offsets(iBurst) = offsets(iBurst + 1);
            onsets(iBurst + 1)  = [];
            offsets(iBurst + 1) = [];
        else
            iBurst = iBurst + 1;
        end
    end

    % ===== 6. ENFORCE MINIMUM DURATION =====
    minDurSamples = round(OPTIONS.minDuration * sFreq);
    durations = offsets - onsets + 1;
    keepMask = durations >= minDurSamples;
    onsets  = onsets(keepMask);
    offsets = offsets(keepMask);

    % ===== 7. BUILD EVENT TIMES =====
    if isempty(onsets)
        evt = [];
    else
        evt = [TimeVector(onsets); TimeVector(offsets)];
    end

    % Statistics for reporting
    stats = struct(...
        'medianEnv',      medEnv, ...
        'madEnv',         madEnv, ...
        'thresholdValue', threshValue, ...
        'nRawBursts',     length(onsets));
end
