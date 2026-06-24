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
% Reference-phase markers (4 per cycle): a band-limited GFP is sign-blind, so it
% oscillates at TWICE the target frequency (2f), and its ANALYTIC PHASE is monotonic
% -- so its landmarks never skip. phi=0 (upward) = a field-magnitude maximum (alpha
% extremum), phi wrapping +pi->-pi = a field minimum (alpha zero-crossing). The GFP
% cannot tell peak from trough or rising from falling, so the polarity is resolved
% ONCE per period: the sign of the highest-(band-)power sensor at its strongest field
% maximum (where that sensor sits at its own extremum, so the sign is unambiguous)
% says whether that maximum is an alpha peak or trough; the rest follow by enforced
% alternation -- maxima alternate peak/trough, the single crossing after a peak is
% "_falling", after a trough is "_rising". Reading the polarity at one reliable
% extremum (not the noisy zero-crossing slope, per marker) is what makes the four
% phases un-skippable. Each period's onset/offset snap to its first/last alpha
% extremum (integer cycles). PhaseValue() maps each type to its numeric phase [rad].
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

    % === EVENT NAMING (automatic) ===
    sProcess.options.label_naming.Comment = ['<I><FONT color="#777777">Events are named from the band: "&lt;band&gt; (lo-hi Hz)" for the periods, plus 4<BR>' ...
        'phase markers "&lt;band&gt;_peak / _trough / _rising / _falling". Band color runs warm<BR>' ...
        '(low f) to cool (high f); markers: peak=red, trough=blue, rising=green, falling=orange.</FONT></I>'];
    sProcess.options.label_naming.Type    = 'label';
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
    Comment = ['Detect bursts (phase polarity): ', sProcess.options.freqband.Value];
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs)
    OutputFiles = {};
    % ===== GET OPTIONS =====
    % Sensor types
    SensorTypes = strtrim(sProcess.options.sensortypes.Value);
    % Frequency band
    freqBands = struct('delta', [2 4], 'theta', [4 8], 'alpha', [8 13], ...
                       'beta', [13 30], 'gamma', [30 60]);
    bandVal = sProcess.options.freqband.Value;
    if isfield(freqBands, bandVal)
        FreqRange = freqBands.(bandVal);
        bandTok   = bandVal;                               % 'alpha', 'beta', ...
    else
        FreqRange = sProcess.options.freqrange.Value{1};
        bandTok   = 'custom';
    end
    if isempty(FreqRange) || (length(FreqRange) ~= 2)
        bst_report('Error', sProcess, [], 'Invalid frequency range.');
        return;
    end
    % Standard, band-derived event labels and colors:
    %   main   "<band> (lo-hi Hz)"  -> warm(low f) .. cool(high f) fill color
    %   peaks  "<band>_peaks"       -> red    (same across all bands)
    %   markers "<band>_peak/_trough/_rising/_falling" -> fixed phase colors
    evtName     = sprintf('%s (%g-%g Hz)', bandTok, FreqRange(1), FreqRange(2));
    peakName    = sprintf('%s_peak',    bandTok);
    troughName  = sprintf('%s_trough',  bandTok);
    risingName  = sprintf('%s_rising',  bandTok);
    fallingName = sprintf('%s_falling', bandTok);
    bandColor    = i_band_color(mean(FreqRange));
    peakColor    = [0.90 0.10 0.10];   % red    (alpha maximum)
    troughColor  = [0.10 0.25 0.90];   % blue   (alpha minimum)
    risingColor  = [0.10 0.70 0.20];   % green  (rising zero-crossing)
    fallingColor = [0.95 0.55 0.10];   % orange (falling zero-crossing)
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
            % Periods (extended) + 4 phase-marker trains (simple), with standard colors
            sFile = i_store_event(sFile, evtName,     evt,             bandColor);
            sFile = i_store_event(sFile, peakName,    markers.peak,    peakColor);
            sFile = i_store_event(sFile, troughName,  markers.trough,  troughColor);
            sFile = i_store_event(sFile, risingName,  markers.rising,  risingColor);
            sFile = i_store_event(sFile, fallingName, markers.falling, fallingColor);
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
                sprintf('[%s] %d periods (%.1f%% coverage); markers peak/trough/rising/falling = %d/%d/%d/%d (%s %g/%g)', ...
                evtName, size(evt,2), stats.coverage, ...
                numel(markers.peak), numel(markers.trough), numel(markers.rising), numel(markers.falling), ...
                OPTIONS.thresholdMode, OPTIONS.enterThresh, OPTIONS.exitThresh));
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
        markers = struct('peak', [], 'trough', [], 'rising', [], 'falling', []);
        stats = [];
        return;
    end
    if (nargin < 4), validMask = []; end
    OPTIONS = struct_copy_fields(OPTIONS, defOptions, 0);
    markers = struct('peak', [], 'trough', [], 'rising', [], 'falling', []);

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

    % ===== PHASE-MARKER TIMING FROM THE MONOTONIC GFP PHASE =====
    % The band GFP at 2f has a monotonic analytic phase, so its landmarks never
    % skip: phi=0 (upward) = a field-magnitude maximum (an alpha extremum, peak OR
    % trough); phi wrapping +pi->-pi = a field minimum (an alpha zero-crossing,
    % rising OR falling). The GFP is sign-blind, so within a burst its maxima
    % alternate peak/trough and its minima alternate falling/rising -- we resolve
    % the polarity ONCE per period (at the highest-|ref| maximum, where the
    % reference channel is at its own alpha extremum and its sign is unambiguous)
    % and propagate every other label by alternation -- no per-marker sign test.
    fLow2  = 2 * OPTIONS.freqRange(1);
    fHigh2 = min(2 * OPTIONS.freqRange(2), 0.95 * sFreq/2);
    doMarkers = (fLow2 < fHigh2);    % disabled only if 2f reaches Nyquist
    if doMarkers
        gfp2 = process_bandpass('Compute', gfp, sFreq, fLow2, fHigh2, 'bst-hfilter-2019', 0);
        phi  = angle(hilbert(gfp2));                                  % monotonic GFP phase, (-pi,pi]
        gpk  = find(phi(1:end-1) <  0    & phi(2:end) >= 0)    + 1;   % phi=0 upward  -> field extremum
        gtr  = find(phi(1:end-1) >  pi/2 & phi(2:end) < -pi/2) + 1;   % +pi->-pi wrap -> zero-crossing
        gpk  = gpk(valid(gpk));
        gtr  = gtr(valid(gtr));
    else
        gpk = [];  gtr = [];
    end

    % Reference channel for the one-per-period polarity seed = highest band-power
    % sensor (at a field maximum it sits at its own alpha extremum, so its sign is
    % unambiguous -- unlike the noisy zero-crossing slope the old code relied on).
    if size(Fbp, 1) > 1
        [~, iRef] = max(var(Fbp(:, valid), 0, 2));
    else
        iRef = 1;
    end
    ref = Fbp(iRef, :);

    % ===== SNAP PERIODS TO CYCLES + LABEL 4 PHASES BY ALTERNATION =====
    minDurSamples = round(OPTIONS.minDuration * sFreq);
    keepOn = [];  keepOff = [];
    aPk = []; aTr = []; aRis = []; aFal = [];
    for iWin = 1:length(onsets)
        on = onsets(iWin);  off = offsets(iWin);
        if doMarkers
            pkIn = gpk(gpk >= on & gpk <= off);
            if numel(pkIn) < 2                       % no genuine oscillation: drop
                continue;
            end
            on  = pkIn(1);                           % snap onset/offset to first/last field extremum
            off = pkIn(end);
        end
        if (off - on + 1) < minDurSamples            % enforce min duration on final bounds
            continue;
        end
        keepOn(end+1)  = on;   %#ok<AGROW>
        keepOff(end+1) = off;  %#ok<AGROW>
        if ~doMarkers, continue; end
        % --- polarity seed: the field maximum with the largest |ref| (most reliable sign) ---
        [~, iSeed] = max(abs(ref(pkIn)));
        seedIsPeak = (ref(pkIn(iSeed)) >= 0);        % that maximum is an alpha peak (else trough)
        sameAsSeed = (mod((1:numel(pkIn)) - iSeed, 2) == 0);   % maxima alternate peak/trough
        if seedIsPeak
            peakMax   = pkIn(sameAsSeed);   troughMax = pkIn(~sameAsSeed);
        else
            troughMax = pkIn(sameAsSeed);   peakMax   = pkIn(~sameAsSeed);
        end
        isPeakMax = sameAsSeed == seedIsPeak;        % true where maximum k is an alpha peak
        % --- exactly ONE zero-crossing between each consecutive pair of maxima ---
        % Structural alternation: max,min,max,min cannot break even if phase noise
        % adds/drops a wrap. A crossing right after a peak is falling, after a trough rising.
        risC = [];  falC = [];
        for k = 1:numel(pkIn)-1
            a = pkIn(k);  b = pkIn(k+1);
            seg = gtr(gtr > a & gtr < b);
            if isempty(seg)
                [~, rel] = min(gfp2(a:b));  m = a + rel - 1;   % phase-noise gap: use the gfp2 minimum
            else
                m = seg(ceil(numel(seg)/2));                   % the wrap (ignore any noise extras)
            end
            if isPeakMax(k), falC(end+1) = m; else, risC(end+1) = m; end %#ok<AGROW>
        end
        aPk  = [aPk,  peakMax];    aTr  = [aTr,  troughMax]; %#ok<AGROW>
        aRis = [aRis, risC];       aFal = [aFal, falC];      %#ok<AGROW>
    end

    % ===== BUILD OUTPUTS =====
    if isempty(keepOn)
        evt = [];
    else
        evt = [TimeVector(keepOn); TimeVector(keepOff)];
    end
    if ~isempty(aPk),  markers.peak    = TimeVector(sort(aPk));  end
    if ~isempty(aTr),  markers.trough  = TimeVector(sort(aTr));  end
    if ~isempty(aRis), markers.rising  = TimeVector(sort(aRis)); end
    if ~isempty(aFal), markers.falling = TimeVector(sort(aFal)); end
    stats = i_stats(keepOn, keepOff, nSamples, TimeVector, Thi, Tlo, OPTIONS);
end


%% ===== PHASE VALUE (radians) FOR A PHASE-MARKER TYPE =====
% The numeric alpha phase each marker type lands on (analytic-phase convention,
% phi = angle(hilbert(reference alpha))): peak=0, falling=+pi/2, trough=+/-pi,
% rising=-pi/2. Stored on the atom so the four phases stay distinguishable as a
% value, not only as a label. Returns NaN for an unknown name.
function v = PhaseValue(name)
    switch lower(name)
        case 'peak',    v =  0;
        case 'falling', v =  pi/2;
        case 'trough',  v =  pi;
        case 'rising',  v = -pi/2;
        otherwise,      v =  NaN;
    end
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
% simple [1xN]) and color (RGB triplet; [] => auto). Returns the updated sFile.
function sFile = i_store_event(sFile, label, times, color)
    if isempty(times)
        return;
    end
    if (nargin < 4), color = []; end
    iEvt = find(strcmpi({sFile.events.label}, label));
    if isempty(iEvt)
        iEvt = length(sFile.events) + 1;
        sEvent = db_template('event');
        sEvent.label = label;
        if ~isempty(color)
            sEvent.color = color;
        else
            sEvent.color = panel_record('GetNewEventColor', iEvt, sFile.events);
        end
    else
        sEvent = sFile.events(iEvt);
        if ~isempty(color)
            sEvent.color = color;    % keep the standard color on re-run
        end
    end
    sEvent.times    = times;
    sEvent.epochs   = ones(1, size(times, 2));
    sEvent.channels = [];
    sEvent.notes    = [];
    sFile.events(iEvt) = sEvent;
end


%% ===== BAND COLOR (warm low-f -> cool high-f) =====
% Maps a band center frequency to a fill color along an HSV ramp from orange
% (low frequency) to violet (high frequency), deliberately avoiding pure red
% (reserved for peak markers) and pure blue (reserved for trough markers).
% The center range [3,45] Hz spans the delta..gamma band centers (log scale).
function rgb = i_band_color(fc)
    fLo = 3;  fHi = 45;
    t   = (log(min(max(fc, fLo), fHi)) - log(fLo)) / (log(fHi) - log(fLo));
    hue = (30 + t*240) / 360;        % 30deg (orange) .. 270deg (violet)
    rgb = hsv2rgb([hue, 0.75, 0.85]);
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
