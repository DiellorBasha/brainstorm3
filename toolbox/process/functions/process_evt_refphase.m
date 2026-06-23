function varargout = process_evt_refphase( varargin )
% PROCESS_EVT_REFPHASE: Detect periods of high band-limited oscillatory power.
%
% Detects time windows where the amplitude of a target frequency band (e.g.
% alpha) is high, stores them as an extended-event group (onset/offset), and
% optionally marks per-cycle reference-phase point events inside each period.
% This is an AMPLITUDE detector. Unlike spectral-contrast burst detectors that
% ask "is there a narrowband peak vs broadband?" (amplitude-invariant, fires
% continuously on oscillation-rich segments), it asks "is the band STRONG right
% now?" using a relative threshold on the signal's own distribution, then marks
% per-cycle phase-polarity reference events. It is the canonical burst detector;
% the earlier CFAR/Wavelet spectral-contrast variants are retired to
% dev/experimental/.
%
% Detection signal: the band Global Field Power (GFP = std across sensors, the
% "Extra > Show GFP" curve in Butterfly mode; Lehmann & Skrandies 1980). For a
% narrowband signal the GFP is ~82% slow envelope, so SMOOTHING (lowpassing) it
% yields a clean, free amplitude trace for thresholding.
%
% Reference-phase markers: a band-limited GFP is sign-blind, so it oscillates at
% TWICE the target frequency (2f), with no energy at the fundamental. Bandpassing
% the GFP at [2*fLow, 2*fHigh] gives a near-sinusoid whose extrema are robust
% cycle landmarks:
%   - GFP PEAKS  (field-magnitude maxima) = the alpha extrema. Consecutive peaks
%     alternate alpha+/alpha- (every other peak is the same alpha phase); the
%     source field's own sign recovers the polarity downstream.
%   - GFP TROUGHS (field-magnitude minima) = the alpha zero-crossings.
% Each detected period's onset/offset snap to its first/last GFP peak, so a
% period spans an integer number of marked cycles (sign-consistent bounds).
%
% Robustness (ported from the artifact detector process_evt_detect):
%   - Bad segments are excluded from BOTH the threshold estimate and detection.
%   - Filter-transient edges are excluded (avoids spurious boundary events).
%   - Optional per-channel normalization for mixed sensor types (different units).
%
% USAGE:  OutputFiles = process_evt_refphase('Run', sProcess, sInputs)
%   [evt,markers,stats]= process_evt_refphase('Compute', F, TimeVector, OPTIONS, validMask)
%             OPTIONS = process_evt_refphase('Compute')   % Get default options
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
    sProcess.Comment     = 'Detect bursts (phase polarity)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Events';
    sProcess.Index       = 47;
    sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/ArtifactsDetect';
    % Input/Output
    sProcess.InputTypes  = {'raw', 'data'};
    sProcess.OutputTypes = {'raw', 'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % === EVENT NAME ===
    sProcess.options.eventname.Comment = 'Event name: ';
    sProcess.options.eventname.Type    = 'text';
    sProcess.options.eventname.Value   = 'power_alpha';
    % Separator
    sProcess.options.sep1.Type    = 'separator';
    sProcess.options.sep1.Comment = ' ';

    % === SENSOR SELECTION ===
    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG';
    sProcess.options.sensortypes.InputTypes = {'data', 'raw'};
    % Normalize channels (mixed sensor types)
    sProcess.options.normalize.Comment = 'Normalize channels (for mixed sensor types)';
    sProcess.options.normalize.Type    = 'checkbox';
    sProcess.options.normalize.Value   = 0;

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

    % === DETECTION SIGNAL ===
    sProcess.options.sep3.Type    = 'separator';
    sProcess.options.sep3.Comment = ' ';
    sProcess.options.label_signal.Comment = ['<B>Detection signal:</B> smoothed band Global Field Power<BR>' ...
        '<I><FONT color="#777777">GFP = std across sensors ("Show GFP" curve). Smoothing, minimum<BR>' ...
        'duration and merge gap are set automatically from the frequency band.</FONT></I>'];
    sProcess.options.label_signal.Type    = 'label';

    % === RELATIVE THRESHOLD ===
    sProcess.options.sep4.Type    = 'separator';
    sProcess.options.sep4.Comment = ' ';
    sProcess.options.thresholdmode.Comment = {'Percentile of envelope', 'Median + k&middot;MAD', 'Median + k&middot;Std', ...
        '<B>Relative threshold:</B>'; ...
        'percentile', 'mad', 'std', ''};
    sProcess.options.thresholdmode.Type    = 'radio_linelabel';
    sProcess.options.thresholdmode.Value   = 'percentile';
    % Enter / exit thresholds (units depend on mode)
    sProcess.options.enterthresh.Comment = 'Enter threshold: ';
    sProcess.options.enterthresh.Type    = 'value';
    sProcess.options.enterthresh.Value   = {85, '', 1};
    sProcess.options.exitthresh.Comment = 'Exit threshold (hysteresis): ';
    sProcess.options.exitthresh.Type    = 'value';
    sProcess.options.exitthresh.Value   = {75, '', 1};
    sProcess.options.label_thresh.Comment = ['<I><FONT color="#777777">Percentile mode: 0-100 (e.g. 85 / 75 marks the top ~15% of time).<BR>' ...
        'MAD / Std mode: multiplier k above the median (e.g. 2 / 1.5).<BR>' ...
        'Set Enter = Exit to disable hysteresis.</FONT></I>'];
    sProcess.options.label_thresh.Type    = 'label';
    % Phase markers are always emitted: "<name>_peak" (field-magnitude maxima =
    % alpha extrema) and "<name>_trough" (field minima = zero-crossings), from
    % the GFP bandpassed at 2f. Period onset/offset snap to the first/last peak.
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = ['Detect bursts (phase polarity): ', sProcess.options.eventname.Value];
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs)
    OutputFiles = {};
    % ===== GET OPTIONS =====
    evtName = strtrim(sProcess.options.eventname.Value);
    if isempty(evtName)
        bst_report('Error', sProcess, [], 'Event name must be specified.');
        return;
    end
    % Sensor types
    SensorTypes = strtrim(sProcess.options.sensortypes.Value);
    % Frequency band
    freqBands = struct('delta', [2 4], 'theta', [4 8], 'alpha', [8 13], ...
                       'beta', [13 30], 'gamma', [30 60]);
    bandVal = sProcess.options.freqband.Value;
    if isfield(freqBands, bandVal)
        FreqRange = freqBands.(bandVal);
    else
        FreqRange = sProcess.options.freqrange.Value{1};
    end
    if isempty(FreqRange) || (length(FreqRange) ~= 2)
        bst_report('Error', sProcess, [], 'Invalid frequency range.');
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
    OPTIONS.freqRange     = FreqRange;
    OPTIONS.normalize     = sProcess.options.normalize.Value;
    OPTIONS.thresholdMode = sProcess.options.thresholdmode.Value;
    OPTIONS.enterThresh   = sProcess.options.enterthresh.Value{1};
    OPTIONS.exitThresh    = sProcess.options.exitthresh.Value{1};
    % smoothing / minDuration / minGap auto-derived from the band in Compute;
    % phase markers (peak+trough) always emitted.

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
        % Warn on mixed sensor types without normalization (GFP mixes units)
        chTypes = unique({ChannelMat.Channel(iChannels).Type});
        if (numel(chTypes) > 1) && ~OPTIONS.normalize
            bst_report('Warning', sProcess, sInputs(iFile), ...
                sprintf(['Multiple sensor types selected (%s) with different units: GFP will be ' ...
                'dominated by the larger-unit sensors. Enable "Normalize channels" or select one type.'], ...
                strjoin(chTypes, ', ')));
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
        % Build valid mask (exclude bad segments from threshold and detection)
        validMask = i_valid_mask(sFile, TimeWindow, DataMat.Time, length(TimeVector));

        % ===== DETECT PERIODS =====
        bst_progress('text', sprintf('File %d/%d: Detecting band-power periods...', iFile, length(sInputs)));
        bst_progress('set', progressPos + round(2 * iFile / length(sInputs) / 3 * 100));
        [evt, markers, stats] = Compute(F, TimeVector, OPTIONS, validMask);

        % ===== SAVE EVENTS =====
        bst_progress('text', sprintf('File %d/%d: Saving events...', iFile, length(sInputs)));
        bst_progress('set', progressPos + round(3 * iFile / length(sInputs) / 3 * 100));
        if ~isempty(evt) && size(evt, 2) > 0
            if ~isfield(sFile, 'events') || isempty(sFile.events)
                sFile.events = repmat(db_template('event'), 0);
            end
            % Periods (extended) + optional phase-marker trains (simple)
            sFile = i_store_event(sFile, evtName, evt);
            if ~isempty(markers.peak)
                sFile = i_store_event(sFile, [evtName '_peak'], markers.peak);
            end
            if ~isempty(markers.trough)
                sFile = i_store_event(sFile, [evtName '_trough'], markers.trough);
            end
            % Save back to file
            if isRaw
                DataMat.F = sFile;
            else
                DataMat.Events = sFile.events;
            end
            DataMat = rmfield(DataMat, 'Time');
            bst_save(file_fullpath(sInputs(iFile).FileName), DataMat, 'v6', 1);
            % Report
            bst_report('Info', sProcess, sInputs(iFile), ...
                sprintf('[Band-power] %d periods (mean dur: %.0f ms, coverage: %.1f%%), %d peak / %d trough markers (mode: %s, enter/exit: %g/%g)', ...
                size(evt,2), 1000*mean(evt(2,:)-evt(1,:)), stats.coverage, ...
                numel(markers.peak), numel(markers.trough), OPTIONS.thresholdMode, OPTIONS.enterThresh, OPTIONS.exitThresh));
        else
            bst_report('Warning', sProcess, sInputs(iFile), ...
                'No periods detected. Try lowering the enter threshold or the minimum duration.');
        end
        iOk(iFile) = true;
    end
    % Return processed files
    OutputFiles = {sInputs(iOk).FileName};
end


%% ===== COMPUTE =====
% USAGE:  [evt, markers, stats] = Compute(F, TimeVector, OPTIONS, validMask)
%                       OPTIONS = Compute()    % Return default options
%
% F          : [nChannels x nSamples] sensor data
% validMask  : [1 x nSamples] logical, false on samples to exclude (bad segments);
%              [] or omitted => all valid. Filter-transient edges are excluded
%              internally in addition to this mask.
% evt        : [2 x nPeriods] extended event times (onset/offset), cycle-snapped
% markers    : struct('peak',[1 x nPk], 'trough',[1 x nTr]) reference-phase times
% stats      : struct with nWindows, coverage, enterLevel, exitLevel, ...
function [evt, markers, stats] = Compute(F, TimeVector, OPTIONS, validMask)
    % Default options. smoothing/minDuration/minGap default [] => auto-derived
    % from the band center (one knob fewer); set explicitly to override.
    defOptions = struct(...
        'freqRange',     [8, 13], ...     % Target band [low, high] in Hz
        'normalize',     false, ...       % Per-channel scale normalization (mixed sensor types)
        'smoothing',     [], ...          % GFP smoothing window (s); [] = auto = 0.5/fc
        'thresholdMode', 'percentile', ...% 'percentile' | 'mad' | 'std'
        'enterThresh',   85, ...          % Enter level (percentile 0-100, or k for mad/std)
        'exitThresh',    75, ...          % Exit level (hysteresis); clamped <= enter level
        'minDuration',   [], ...          % Min period duration (s); [] = auto = 3/fc
        'minGap',        []);             % Merge gap (s); [] = auto = 2/fc
    % Return defaults if no input
    if (nargin == 0)
        evt = defOptions;
        markers = struct('peak', [], 'trough', []);
        stats = [];
        return;
    end
    if (nargin < 4), validMask = []; end
    OPTIONS = struct_copy_fields(OPTIONS, defOptions, 0);
    markers = struct('peak', [], 'trough', []);

    % Auto-derive timing parameters from the band center (fewer GUI knobs).
    % A band GFP ripples at 2f, so a boxcar of half a period (0.5/fc) nulls that
    % ripple with minimal blur (empirically -23 dB, 99% envelope kept).
    fc = mean(OPTIONS.freqRange);
    if isempty(OPTIONS.smoothing),   OPTIONS.smoothing   = 0.5 / fc; end   % null boxcar at 2f
    if isempty(OPTIONS.minDuration), OPTIONS.minDuration = 3   / fc; end   % >= 3 cycles
    if isempty(OPTIONS.minGap),      OPTIONS.minGap      = 2   / fc; end   % merge if < 2 cycles apart

    % Sampling frequency
    sFreq = 1 / (TimeVector(2) - TimeVector(1));
    nSamples = length(TimeVector);
    % Valid mask
    if isempty(validMask)
        valid = true(1, nSamples);
    else
        valid = logical(validMask(:)');
    end

    % ===== OPTIONAL PER-CHANNEL NORMALIZATION =====
    % Equalize units across sensor types: divide each channel by its robust scale
    % (estimated over valid samples) before the GFP.
    if OPTIONS.normalize && (size(F,1) > 1)
        sc = std(F(:, valid), 0, 2);
        sc(sc <= 0) = 1;
        F = F ./ sc;
    end

    % ===== BUILD DETECTION SIGNAL (smoothed band GFP) =====
    % 1. Bandpass to target band (keep FiltSpec for the transient length)
    [Fbp, FiltSpec] = process_bandpass('Compute', F, sFreq, OPTIONS.freqRange(1), OPTIONS.freqRange(2), 'bst-hfilter-2019', 0);
    % Exclude filter-transient edges from the valid mask
    if isstruct(FiltSpec) && isfield(FiltSpec, 'transient') && ~isempty(FiltSpec.transient)
        smpTrans = round(FiltSpec.transient * sFreq);
        if (smpTrans > 0) && (2*smpTrans < nSamples)
            valid(1:smpTrans) = false;
            valid(end-smpTrans+1:end) = false;
        end
    end
    % 2. Global Field Power across sensors (std across channels; Lehmann & Skrandies)
    if size(Fbp, 1) > 1
        gfp = std(Fbp, 1, 1);             % [1 x nSamples]
    else
        gfp = abs(Fbp);                   % single-channel fallback
    end
    % 3. Lowpass/smooth the GFP -> amplitude detection signal (~82% of GFP energy)
    wSmp = max(1, round(OPTIONS.smoothing * sFreq));
    sig = movmean(gfp, wSmp);

    % ===== RELATIVE THRESHOLD (estimated over valid samples only) =====
    sigValid = sig(valid);
    if isempty(sigValid) || (max(sigValid) <= min(sigValid))
        % No usable samples / flat envelope: nothing to detect
        evt = [];
        stats = i_stats([], [], nSamples, TimeVector, Inf, Inf, OPTIONS);
        return;
    end
    [Thi, Tlo] = i_threshold(sigValid, OPTIONS.thresholdMode, OPTIONS.enterThresh, OPTIONS.exitThresh);

    % ===== HYSTERESIS INTERVAL SCAN (never enter on invalid samples) =====
    mask = false(1, nSamples);
    inWin = false;
    for t = 1:nSamples
        if ~valid(t)
            inWin = false;                % bad/edge sample breaks any open window
            mask(t) = false;
            continue;
        end
        if ~inWin && sig(t) >= Thi
            inWin = true;
        elseif inWin && sig(t) < Tlo
            inWin = false;
        end
        mask(t) = inWin;
    end
    % Candidate periods
    d = diff([0, mask, 0]);
    onsets  = find(d == 1);
    offsets = find(d == -1) - 1;
    % Merge close candidates
    minGapSamples = round(OPTIONS.minGap * sFreq);
    iWin = 1;
    while iWin < length(onsets)
        if (onsets(iWin + 1) - offsets(iWin)) < minGapSamples
            offsets(iWin) = max(offsets(iWin), offsets(iWin + 1));
            onsets(iWin + 1)  = [];
            offsets(iWin + 1) = [];
        else
            iWin = iWin + 1;
        end
    end

    % ===== REFERENCE-PHASE MARKERS FROM GFP AT 2f =====
    % A band-limited GFP oscillates at 2f; bandpass there for a near-sinusoid.
    allPeaks = [];  allTroughs = [];
    fLow2  = 2 * OPTIONS.freqRange(1);
    fHigh2 = min(2 * OPTIONS.freqRange(2), 0.95 * sFreq/2);
    doMarkers = (fLow2 < fHigh2);    % disabled only if 2f reaches Nyquist
    if doMarkers
        gfp2 = process_bandpass('Compute', gfp, sFreq, fLow2, fHigh2, 'bst-hfilter-2019', 0);
        dg = diff(gfp2);
        % Local maxima (peaks, phase 0) and minima (troughs, phase pi)
        pk = find(dg(1:end-1) > 0 & dg(2:end) <= 0) + 1;
        tr = find(dg(1:end-1) < 0 & dg(2:end) >= 0) + 1;
        % Keep only landmarks on valid samples
        pk = pk(valid(pk));
        tr = tr(valid(tr));
    else
        pk = [];  tr = [];
    end

    % ===== SNAP PERIODS TO CYCLES + COLLECT MARKERS =====
    minDurSamples = round(OPTIONS.minDuration * sFreq);
    keepOn = [];  keepOff = [];
    for iWin = 1:length(onsets)
        on = onsets(iWin);  off = offsets(iWin);
        if doMarkers
            pkIn = pk(pk >= on & pk <= off);
            if numel(pkIn) < 2                       % no genuine oscillation: drop
                continue;
            end
            on  = pkIn(1);                           % snap onset/offset to first/last peak
            off = pkIn(end);
        end
        if (off - on + 1) < minDurSamples            % enforce min duration on final bounds
            continue;
        end
        keepOn(end+1)  = on;   %#ok<AGROW>
        keepOff(end+1) = off;  %#ok<AGROW>
        if doMarkers
            allPeaks   = [allPeaks,   pk(pk >= on & pk <= off)];   %#ok<AGROW>
            allTroughs = [allTroughs, tr(tr >= on & tr <= off)];   %#ok<AGROW>
        end
    end

    % ===== BUILD OUTPUTS =====
    if isempty(keepOn)
        evt = [];
    else
        evt = [TimeVector(keepOn); TimeVector(keepOff)];
    end
    if ~isempty(allPeaks),   markers.peak   = TimeVector(sort(allPeaks));   end
    if ~isempty(allTroughs), markers.trough = TimeVector(sort(allTroughs)); end
    stats = i_stats(keepOn, keepOff, nSamples, TimeVector, Thi, Tlo, OPTIONS);
end


%% ===== STATS HELPER =====
function stats = i_stats(onsets, offsets, nSamples, TimeVector, Thi, Tlo, OPTIONS)
    if isempty(onsets)
        cov = 0;  nW = 0;
    else
        cov = 100 * sum(offsets - onsets + 1) / nSamples;
        nW  = numel(onsets);
    end
    stats = struct(...
        'nWindows',    nW, ...
        'coverage',    cov, ...
        'enterLevel',  Thi, ...
        'exitLevel',   Tlo, ...
        'smoothing',   OPTIONS.smoothing, ...    % auto-derived (s)
        'minDuration', OPTIONS.minDuration, ...  % auto-derived (s)
        'minGap',      OPTIONS.minGap, ...       % auto-derived (s)
        'totalDur',    TimeVector(end) - TimeVector(1));
end


%% ===== RELATIVE THRESHOLD =====
% Returns enter level Thi and exit level Tlo (Tlo <= Thi) for a 1xN signal.
function [Thi, Tlo] = i_threshold(sig, mode, enter, exit)
    switch lower(mode)
        case 'percentile'
            % enter/exit are percentiles in 0-100
            Thi = prctile(sig, enter);
            Tlo = prctile(sig, exit);
        case 'mad'
            % enter/exit are k multipliers of the (std-consistent) MAD above the median
            med = median(sig);
            s   = 1.4826 * median(abs(sig - med));
            Thi = med + enter * s;
            Tlo = med + exit  * s;
        case 'std'
            % enter/exit are k multipliers of the std above the median
            med = median(sig);
            s   = std(sig);
            Thi = med + enter * s;
            Tlo = med + exit  * s;
        otherwise
            error('process_evt_refphase:badMode', 'Unknown threshold mode "%s".', mode);
    end
    % Hysteresis sanity: exit level cannot exceed enter level
    if Tlo > Thi
        Tlo = Thi;
    end
end


%% ===== STORE EVENT GROUP =====
% Find-or-create an event group by label and set its times (extended [2xN] or
% simple [1xN]); returns the updated sFile.
function sFile = i_store_event(sFile, label, times)
    if isempty(times)
        return;
    end
    iEvt = find(strcmpi({sFile.events.label}, label));
    if isempty(iEvt)
        iEvt = length(sFile.events) + 1;
        sEvent = db_template('event');
        sEvent.label = label;
        sEvent.color = panel_record('GetNewEventColor', iEvt, sFile.events);
    else
        sEvent = sFile.events(iEvt);
    end
    sEvent.times    = times;
    sEvent.epochs   = ones(1, size(times, 2));
    sEvent.channels = [];
    sEvent.notes    = [];
    sFile.events(iEvt) = sEvent;
end


%% ===== VALID-SAMPLE MASK FROM BAD SEGMENTS =====
% Build a [1 x nRead] logical (false on bad samples) for the read window, mapping
% the file's bad segments into read-window indexing (mirrors process_evt_detect).
function valid = i_valid_mask(sFile, TimeWindow, DataTime, nRead)
    valid = true(1, nRead);
    badSeg = panel_record('GetBadSegments', sFile);
    if isempty(badSeg)
        return;
    end
    % Absolute samples -> read-window indices
    badSeg = badSeg - round(sFile.prop.times(1) .* sFile.prop.sfreq) + 1;
    if ~isempty(TimeWindow)
        badSeg = badSeg - (bst_closest(TimeWindow(1), DataTime) - 1);
    end
    for i = 1:size(badSeg, 2)
        a = max(1, badSeg(1,i));
        b = min(nRead, badSeg(2,i));
        if a <= b
            valid(a:b) = false;
        end
    end
end
