function varargout = process_helmholtz_events( varargin )
% PROCESS_HELMHOLTZ_EVENTS: Helmholtz-Hodge maps of the source field at events.
%
% Batch version of the view_helmholtz GUI. At each timepoint of a phase-marker
% event group (e.g. "alpha_peak" from process_evt_refphase), it reconstructs the
% UNCONSTRAINED source vector field from a kernel-link results file and the
% (band-passed) sensor data, runs the Helmholtz-Hodge decomposition, and stores
% the resulting static source maps as "sources" (results) files:
%
%   <event> | Source vector J   : the source vector field      (nComponents=3)
%   <event> | Total field |J|   : field magnitude  (Ht.Fmag)   (nComponents=1)
%   <event> | Potential Phi     : irrotational scalar (Ht.Phi) (nComponents=1, signed)
%   <event> | Stream Psi        : solenoidal scalar  (Ht.Psi)  (nComponents=1, signed)
%
% Pipeline (mirrors view_helmholtz):
%   J(t) = ImagingKernel * Fbp(GoodChannel, t),  Fbp = bandpass(sensors, band)
%   Op   = bst_helmholtz('Prepare', {Dirac,LBO}, Manifold, Surf, 'Domain','vertex')
%   Ht   = bst_helmholtz('Frame', Op, J(t), false)   % cores off (detection is a later step)
%
% Input  : a kernel-link results file (UNCONSTRAINED, nComponents==3).
% Output : 4 results files, one frame per event (Time = the event times).
%
% USAGE:  OutputFiles = process_helmholtz_events('Run', sProcess, sInputs)
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
    sProcess.Comment     = 'Helmholtz maps at events';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 337;
    sProcess.Description = 'https://neuroimage.usc.edu/brainstorm';
    % Input/Output: an unconstrained source kernel link
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % === EVENT NAME ===
    sProcess.options.eventname.Comment = 'Phase-marker event (timepoints): ';
    sProcess.options.eventname.Type    = 'text';
    sProcess.options.eventname.Value   = 'alpha_peak';
    sProcess.options.label_evt.Comment = ['<I><FONT color="#777777">A simple-event group from "Detect bursts (phase polarity)", e.g.<BR>' ...
        'alpha_peak / alpha_trough / alpha_rising / alpha_falling.</FONT></I>'];
    sProcess.options.label_evt.Type    = 'label';
    % Separator
    sProcess.options.sep1.Type    = 'separator';
    sProcess.options.sep1.Comment = ' ';

    % === FREQUENCY BAND (must match the event's band) ===
    sProcess.options.freqband.Comment = {'Delta (2-4 Hz)', 'Theta (4-8 Hz)', ...
        'Alpha (8-13 Hz)', 'Beta (13-30 Hz)', 'Gamma (30-60 Hz)', 'Custom (below)', ...
        '<B>Bandpass applied to the sensors (match the event band):</B>'; ...
        'delta', 'theta', 'alpha', 'beta', 'gamma', 'custom', ''};
    sProcess.options.freqband.Type    = 'radio_linelabel';
    sProcess.options.freqband.Value   = 'alpha';
    sProcess.options.freqrange.Comment = 'Custom frequency range:';
    sProcess.options.freqrange.Type    = 'freqrange';
    sProcess.options.freqrange.Value   = {[8, 13], 'Hz', []};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = ['Helmholtz maps at events: ', sProcess.options.eventname.Value];
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs)
    OutputFiles = {};
    % ===== OPTIONS =====
    evtName = strtrim(sProcess.options.eventname.Value);
    if isempty(evtName)
        bst_report('Error', sProcess, [], 'Event name must be specified.');
        return;
    end
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

    progressPos = bst_progress('get');
    for iFile = 1:length(sInputs)
        ResultsFile = sInputs(iFile).FileName;
        bst_progress('text', sprintf('File %d/%d: Loading kernel and data...', iFile, length(sInputs)));

        % ===== LOAD KERNEL LINK =====
        R = in_bst_results(ResultsFile, 0, 'ImagingKernel', 'GoodChannel', ...
            'SurfaceFile', 'DataFile', 'nComponents');
        if isempty(R.ImagingKernel)
            bst_report('Error', sProcess, sInputs(iFile), ...
                'This process needs a kernel-link results file (ImagingKernel), not a full source map.');
            continue;
        end
        if (R.nComponents ~= 3)
            bst_report('Error', sProcess, sInputs(iFile), ...
                sprintf('Source model must be UNCONSTRAINED (nComponents=3); got %d.', R.nComponents));
            continue;
        end
        if isempty(R.DataFile)
            bst_report('Error', sProcess, sInputs(iFile), 'Results file has no linked recording (DataFile).');
            continue;
        end

        % ===== LOAD RECORDING + EVENTS =====
        ChannelFile = bst_get('ChannelFileForStudy', R.DataFile);
        ChannelMat  = in_bst_channel(ChannelFile);
        [F, TimeVector, events] = i_load_recording(R.DataFile, ChannelMat);
        if isempty(F)
            bst_report('Error', sProcess, sInputs(iFile), 'Could not read the linked recording.');
            continue;
        end
        % Find the requested event group
        if isempty(events)
            bst_report('Error', sProcess, sInputs(iFile), 'The recording has no events.');
            continue;
        end
        iEvt = find(strcmpi({events.label}, evtName), 1);
        if isempty(iEvt)
            bst_report('Error', sProcess, sInputs(iFile), ...
                sprintf('Event group "%s" not found in the recording.', evtName));
            continue;
        end
        evtTimes = events(iEvt).times(1, :);     % markers are simple events ([1 x nE])
        iT = bst_closest(evtTimes, TimeVector);  % nearest sample index per event
        nE = numel(iT);
        if (nE == 0)
            bst_report('Error', sProcess, sInputs(iFile), sprintf('Event group "%s" is empty.', evtName));
            continue;
        end

        % ===== BANDPASS SENSORS (good channels only) =====
        bst_progress('text', sprintf('File %d/%d: Band-passing sensors...', iFile, length(sInputs)));
        sFreq = 1 / (TimeVector(2) - TimeVector(1));
        Fgood = F(R.GoodChannel, :);
        Fbp   = process_bandpass('Compute', Fgood, sFreq, FreqRange(1), FreqRange(2), 'bst-hfilter-2019', 0);

        % ===== PREPARE HELMHOLTZ OPERATOR (once) =====
        bst_progress('text', sprintf('File %d/%d: Preparing Helmholtz operator...', iFile, length(sInputs)));
        Dirac = bst_get_operator_node(R.SurfaceFile, 'Dirac');
        LBO   = bst_get_operator_node(R.SurfaceFile, 'Laplace-Beltrami');
        Surf  = in_tess_bst(R.SurfaceFile, 0);
        Mani  = tess_manifold(R.SurfaceFile);
        Op    = bst_helmholtz('Prepare', {Dirac, LBO}, Mani, Surf, 'Domain', 'vertex');
        nV    = Op.nVtot;

        % ===== DECOMPOSE AT EACH EVENT =====
        Jvec = zeros(3*nV, nE);
        Fmag = zeros(nV, nE);
        Phi  = zeros(nV, nE);
        Psi  = zeros(nV, nE);
        for k = 1:nE
            bst_progress('text', sprintf('File %d/%d: Helmholtz frame %d/%d...', iFile, length(sInputs), k, nE));
            bst_progress('set', progressPos + round(100 * k / nE));
            J  = R.ImagingKernel * Fbp(:, iT(k));     % [3nV x 1] unconstrained source vector
            Ht = bst_helmholtz('Frame', Op, J, false);% cores off (detection is a later step)
            Fmag(:, k)  = Ht.Fmag;
            Phi(:, k)   = Ht.Phi;
            Psi(:, k)   = Ht.Psi;
            Jvec(:, k)  = reshape(Ht.Vtot', [], 1);   % [nV x 3] -> [3nV x 1] interleaved
        end

        % ===== SAVE THE 4 SOURCE MAPS =====
        bst_progress('text', sprintf('File %d/%d: Saving source maps...', iFile, length(sInputs)));
        T = evtTimes;
        OutputFiles{end+1} = i_save_results(Jvec, 3, [evtName ' | Source vector J'], '',      T, R, sInputs(iFile)); %#ok<AGROW>
        OutputFiles{end+1} = i_save_results(Fmag, 1, [evtName ' | Total field |J|'], 'source', T, R, sInputs(iFile)); %#ok<AGROW>
        OutputFiles{end+1} = i_save_results(Phi,  1, [evtName ' | Potential Phi'],   'stat2',  T, R, sInputs(iFile)); %#ok<AGROW>
        OutputFiles{end+1} = i_save_results(Psi,  1, [evtName ' | Stream Psi'],      'stat2',  T, R, sInputs(iFile)); %#ok<AGROW>

        bst_report('Info', sProcess, sInputs(iFile), ...
            sprintf('Helmholtz at %d "%s" events (band %g-%g Hz): saved Source vector J, |J|, Phi, Psi.', ...
            nE, evtName, FreqRange(1), FreqRange(2)));
    end
end


%% ===== LOAD RECORDING (data or raw link) =====
function [F, TimeVector, events] = i_load_recording(DataFile, ChannelMat)
    F = [];  TimeVector = [];  events = [];
    DataMat = in_bst_data(DataFile, 'F', 'Time', 'Events');
    if isstruct(DataMat.F)
        % Raw link: F holds the sFile descriptor -> read full continuous block
        sFile = DataMat.F;
        ImportOptions = db_template('ImportOptions');
        ImportOptions.ImportMode='Time'; ImportOptions.UseCtfComp=1; ImportOptions.UseSsp=1;
        ImportOptions.EventsMode='ignore'; ImportOptions.DisplayMessages=0; ImportOptions.RemoveBaseline='no';
        [F, TimeVector] = in_fread(sFile, ChannelMat, 1, [], [], ImportOptions);
        events = sFile.events;
    else
        F = DataMat.F;
        TimeVector = DataMat.Time;
        if isfield(DataMat, 'Events'), events = DataMat.Events; end
    end
end


%% ===== SAVE A SOURCE (RESULTS) MAP =====
function OutFile = i_save_results(Field, nComp, Comment, ColormapType, TimeVec, R, sInput)
    % Duplicate a single frame so the viewer always has a time axis
    if size(Field, 2) < 2
        Field   = [Field, Field];
        TimeVec = [TimeVec, TimeVec + 1e-3];
    end
    ResultsMat = db_template('resultsmat');
    ResultsMat.ImageGridAmp  = Field;
    ResultsMat.Time          = TimeVec;
    ResultsMat.nComponents   = nComp;
    ResultsMat.Comment       = Comment;
    ResultsMat.Function      = 'process_helmholtz_events';
    ResultsMat.SurfaceFile   = R.SurfaceFile;
    ResultsMat.DataFile      = R.DataFile;
    ResultsMat.HeadModelType = 'surface';
    if ~isempty(ColormapType)
        ResultsMat.ColormapType = ColormapType;
    end
    ResultsMat = bst_history('add', ResultsMat, 'compute', ...
        sprintf('Helmholtz-Hodge map at events (process_helmholtz_events): %s', Comment));
    % Save + register (output folder comes from the real recording, not the
    % "link|kernel|data" input filename whose bst_fileparts would be bogus)
    OutFile = bst_process('GetNewFilename', bst_fileparts(R.DataFile), 'results_helmholtz');
    bst_save(OutFile, ResultsMat, 'v6');
    db_add_data(sInput.iStudy, OutFile, ResultsMat);
end
