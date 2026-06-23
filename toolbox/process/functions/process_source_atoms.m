function varargout = process_source_atoms( varargin )
% PROCESS_SOURCE_ATOMS: Populate a spatiotemporal "atom" table from source peaks.
%
% MILESTONE-1 TEST POPULATE (deliberately trivial): at each timepoint of a
% phase-marker event group, reconstruct the unconstrained source magnitude |J|
% and take its strongest local maxima as "atoms" (db_template('atom')). This is
% only to generate data for the atom system (data model + GUI); the principled
% source/sink and vortex detection (bst_helmholtz cores) comes later.
%
% Input  : an unconstrained source kernel link (results, nComponents==3).
% Output : a dynamics_*.mat table (bst_dynamics), loadable with view_dynamics.
%
% USAGE:  OutputFiles = process_source_atoms('Run', sProcess, sInputs)
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
    sProcess.Comment     = 'Detect source atoms (peaks)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 338;
    sProcess.Description = 'https://neuroimage.usc.edu/brainstorm';
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % Event group
    sProcess.options.eventname.Comment = 'Phase-marker event (timepoints): ';
    sProcess.options.eventname.Type    = 'text';
    sProcess.options.eventname.Value   = 'alpha_peak';
    % Separator
    sProcess.options.sep1.Type    = 'separator';
    sProcess.options.sep1.Comment = ' ';
    % Frequency band
    sProcess.options.freqband.Comment = {'Delta (2-4 Hz)', 'Theta (4-8 Hz)', ...
        'Alpha (8-13 Hz)', 'Beta (13-30 Hz)', 'Gamma (30-60 Hz)', 'Custom (below)', ...
        '<B>Bandpass applied to the sensors (match the event band):</B>'; ...
        'delta', 'theta', 'alpha', 'beta', 'gamma', 'custom', ''};
    sProcess.options.freqband.Type    = 'radio_linelabel';
    sProcess.options.freqband.Value   = 'alpha';
    sProcess.options.freqrange.Comment = 'Custom frequency range:';
    sProcess.options.freqrange.Type    = 'freqrange';
    sProcess.options.freqrange.Value   = {[8, 13], 'Hz', []};
    % Number of peaks per event
    sProcess.options.npeaks.Comment = 'Peaks (local |J| maxima) per event: ';
    sProcess.options.npeaks.Type    = 'value';
    sProcess.options.npeaks.Value   = {3, '', 0};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = ['Detect source atoms: ', sProcess.options.eventname.Value];
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs)
    OutputFiles = {};
    evtName = strtrim(sProcess.options.eventname.Value);
    if isempty(evtName)
        bst_report('Error', sProcess, [], 'Event name must be specified.');
        return;
    end
    freqBands = struct('delta', [2 4], 'theta', [4 8], 'alpha', [8 13], 'beta', [13 30], 'gamma', [30 60]);
    bandVal = sProcess.options.freqband.Value;
    if isfield(freqBands, bandVal)
        FreqRange = freqBands.(bandVal);
    else
        FreqRange = sProcess.options.freqrange.Value{1};
    end
    nPeaks = sProcess.options.npeaks.Value{1};

    for iFile = 1:length(sInputs)
        ResultsFile = sInputs(iFile).FileName;
        % --- Load kernel link ---
        R = in_bst_results(ResultsFile, 0, 'ImagingKernel', 'GoodChannel', 'SurfaceFile', 'DataFile', 'nComponents');
        if isempty(R.ImagingKernel) || (R.nComponents ~= 3)
            bst_report('Error', sProcess, sInputs(iFile), 'Need an UNCONSTRAINED source kernel link (nComponents=3, ImagingKernel).');
            continue;
        end
        % --- Recording + events ---
        ChannelMat = in_bst_channel(bst_get('ChannelFileForStudy', R.DataFile));
        [F, TimeVector, events] = i_load_recording(R.DataFile, ChannelMat);
        if isempty(F), bst_report('Error', sProcess, sInputs(iFile), 'Could not read the recording.'); continue; end
        iEvt = find(strcmpi({events.label}, evtName), 1);
        if isempty(iEvt), bst_report('Error', sProcess, sInputs(iFile), sprintf('Event "%s" not found.', evtName)); continue; end
        evtTimes = events(iEvt).times(1, :);
        iT = bst_closest(evtTimes, TimeVector);
        nE = numel(iT);
        % --- Bandpass + surface ---
        sFreq = 1 / (TimeVector(2) - TimeVector(1));
        Fbp   = process_bandpass('Compute', F(R.GoodChannel, :), sFreq, FreqRange(1), FreqRange(2), 'bst-hfilter-2019', 0);
        Surf  = in_tess_bst(R.SurfaceFile, 0);
        nV    = size(Surf.Vertices, 1);

        % --- Build the table ---
        T = bst_dynamics('New', sprintf('%s atoms (%g-%g Hz)', evtName, FreqRange(1), FreqRange(2)));
        T.DataFile = R.DataFile;  T.SurfaceFile = R.SurfaceFile;
        T.Options = struct('method','source-peak', 'eventname',evtName, 'band',FreqRange, 'nPeaks',nPeaks);
        for k = 1:nE
            J    = R.ImagingKernel * Fbp(:, iT(k));               % [3nV x 1]
            Jmag = sqrt(sum(reshape(J, 3, []).^2, 1))';           % [nV x 1] source magnitude
            vPk  = i_local_maxima(Jmag, Surf.VertConn, nPeaks);   % strongest local maxima
            for v = vPk(:)'
                A = bst_dynamics('NewAtom');
                A.label       = 'peak';
                A.category    = 'peak';
                A.color       = [0.85 0.20 0.75];
                A.time        = evtTimes(k);
                A.sample      = iT(k);
                A.sourceEvent = evtName;
                A.vertex      = v;
                A.pos         = Surf.Vertices(v, :);
                A.hemi        = uint8(1 + (Surf.Vertices(v,2) < 0));   % SCS: Y>0 = left
                A.band        = FreqRange;
                A.bandName    = bandVal;
                A.strength    = Jmag(v);
                A.DataFile    = R.DataFile;
                A.ResultsFile = ResultsFile;
                A.SurfaceFile = R.SurfaceFile;
                T = bst_dynamics('Add', T, A);
            end
        end

        % --- Save (by path; Phase-1 has no tree node yet, so NOT returned as a
        %     tracked OutputFile -- bst_process can't resolve the dynamics type) ---
        OutFile = bst_process('GetNewFilename', bst_fileparts(R.DataFile), 'dynamics');
        bst_dynamics('Save', OutFile, T);
        bst_report('Info', sProcess, sInputs(iFile), ...
            sprintf('%d source atoms from %d "%s" events. Saved: %s  (open with view_dynamics).', ...
            T.nAtoms, nE, evtName, OutFile));
    end
end


%% ===== LOAD RECORDING (data or raw) =====
function [F, TimeVector, events] = i_load_recording(DataFile, ChannelMat)
    F = [];  TimeVector = [];  events = [];
    DataMat = in_bst_data(DataFile, 'F', 'Time', 'Events');
    if isstruct(DataMat.F)
        sFile = DataMat.F;
        ImportOptions = db_template('ImportOptions');
        ImportOptions.ImportMode='Time'; ImportOptions.UseCtfComp=1; ImportOptions.UseSsp=1;
        ImportOptions.EventsMode='ignore'; ImportOptions.DisplayMessages=0; ImportOptions.RemoveBaseline='no';
        [F, TimeVector] = in_fread(sFile, ChannelMat, 1, [], [], ImportOptions);
        events = sFile.events;
    else
        F = DataMat.F;  TimeVector = DataMat.Time;
        if isfield(DataMat, 'Events'), events = DataMat.Events; end
    end
end


%% ===== TOP-N LOCAL MAXIMA ON THE SURFACE GRAPH =====
function vPk = i_local_maxima(field, VertConn, nPeaks)
    vPk = [];
    [~, order] = sort(field, 'descend');
    for v = order'
        nb = find(VertConn(v, :));
        if isempty(nb) || all(field(v) >= field(nb))
            vPk(end+1) = v; %#ok<AGROW>
            if numel(vPk) >= nPeaks, break; end
        end
    end
end
