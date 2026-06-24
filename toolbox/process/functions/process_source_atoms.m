function varargout = process_source_atoms( varargin )
% PROCESS_SOURCE_ATOMS: Populate a nested spatiotemporal atom table from source peaks.
%
% MILESTONE-1 TEST POPULATE (deliberately trivial detection): builds the atom
% hierarchy that mirrors the refphase output --
%
%   <band> (lo-hi Hz)          extended window  [parent]
%   |- <band>_peak             simple points    [child]
%   |- <band>_trough           simple points
%   |- <band>_rising           simple points
%   '- <band>_falling          simple points
%
% The parent window group carries the period time-windows (no spatial markers).
% Each child phase group carries, per phase-marker event time, the strongest
% local maxima of the unconstrained source magnitude |J| as occurrences (vertex,
% pos, hemi, strength). The principled source/sink and vortex detection
% (bst_helmholtz cores) replaces this peak populate in a later phase.
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

    % Frequency band (selects the refphase groups: "<band> (..)" window + "<band>_<phase>")
    sProcess.options.freqband.Comment = {'Delta (2-4 Hz)', 'Theta (4-8 Hz)', ...
        'Alpha (8-13 Hz)', 'Beta (13-30 Hz)', 'Gamma (30-60 Hz)', 'Custom (below)', ...
        '<B>Band (matches the refphase event groups):</B>'; ...
        'delta', 'theta', 'alpha', 'beta', 'gamma', 'custom', ''};
    sProcess.options.freqband.Type    = 'radio_linelabel';
    sProcess.options.freqband.Value   = 'alpha';
    sProcess.options.freqrange.Comment = 'Custom frequency range:';
    sProcess.options.freqrange.Type    = 'freqrange';
    sProcess.options.freqrange.Value   = {[8, 13], 'Hz', []};
    % Peaks per phase-marker event
    sProcess.options.npeaks.Comment = 'Peaks (local |J| maxima) per event: ';
    sProcess.options.npeaks.Type    = 'value';
    sProcess.options.npeaks.Value   = {3, '', 0};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = ['Detect source atoms: ', sProcess.options.freqband.Value];
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs)
    OutputFiles = {};
    freqBands = struct('delta', [2 4], 'theta', [4 8], 'alpha', [8 13], 'beta', [13 30], 'gamma', [30 60]);
    bandVal = sProcess.options.freqband.Value;
    if isfield(freqBands, bandVal)
        FreqRange = freqBands.(bandVal);  bandTok = bandVal;
    else
        FreqRange = sProcess.options.freqrange.Value{1};  bandTok = 'custom';
    end
    nPeaks = sProcess.options.npeaks.Value{1};
    phases = {'peak','trough','rising','falling'};
    phaseColors = struct('peak',[0.90 0.10 0.10], 'trough',[0.10 0.25 0.90], ...
                         'rising',[0.10 0.70 0.20], 'falling',[0.95 0.55 0.10]);

    for iFile = 1:length(sInputs)
        ResultsFile = sInputs(iFile).FileName;
        R = in_bst_results(ResultsFile, 0, 'ImagingKernel', 'GoodChannel', 'SurfaceFile', 'DataFile', 'nComponents');
        if isempty(R.ImagingKernel) || (R.nComponents ~= 3)
            bst_report('Error', sProcess, sInputs(iFile), 'Need an UNCONSTRAINED source kernel link (nComponents=3, ImagingKernel).');
            continue;
        end
        ChannelMat = in_bst_channel(bst_get('ChannelFileForStudy', R.DataFile));
        [F, TimeVector, events] = i_load_recording(R.DataFile, ChannelMat);
        if isempty(F) || isempty(events)
            bst_report('Error', sProcess, sInputs(iFile), 'Could not read the recording / no events.');
            continue;
        end
        sFreq = 1 / (TimeVector(2) - TimeVector(1));
        Fbp   = process_bandpass('Compute', F(R.GoodChannel, :), sFreq, FreqRange(1), FreqRange(2), 'bst-hfilter-2019', 0);
        Surf  = in_tess_bst(R.SurfaceFile, 0);

        T = bst_dynamics('New', sprintf('%s atoms (%g-%g Hz)', bandTok, FreqRange(1), FreqRange(2)));
        T.DataFile = R.DataFile;  T.SurfaceFile = R.SurfaceFile;
        T.Options  = struct('method','source-peak', 'band',FreqRange, 'bandName',bandTok, 'nPeaks',nPeaks);

        % ===== PARENT: extended window group "<band> (lo-hi Hz)" =====
        winLabel = '';
        iWin = find(~cellfun(@isempty, regexp({events.label}, ['^' regexptranslate('escape', bandTok) ' \('], 'once')), 1);
        if ~isempty(iWin)
            w = events(iWin);
            winLabel = w.label;
            G0 = bst_dynamics('NewGroup', winLabel);
            G0.type = 'extended';  G0.times = w.times;  G0.epochs = ones(1, size(w.times,2));
            G0.color = w.color;  if isempty(G0.color), G0.color = [0.5 0.5 0.5]; end
            G0.band = FreqRange;  G0.bandName = bandTok;
            G0.DataFile = R.DataFile;  G0.ResultsFile = ResultsFile;  G0.SurfaceFile = R.SurfaceFile;
            T = bst_dynamics('AddGroup', T, G0);
        end

        % ===== CHILDREN: simple phase groups "<band>_<phase>" =====
        for ip = 1:numel(phases)
            ph = phases{ip};
            iE = find(strcmpi({events.label}, [bandTok '_' ph]), 1);
            if isempty(iE), continue; end
            et = events(iE).times(1, :);
            iT = bst_closest(et, TimeVector);
            accT = [];  accV = [];  accP = zeros(0,3);  accH = [];  accS = [];
            for k = 1:numel(iT)
                J    = R.ImagingKernel * Fbp(:, iT(k));
                Jmag = sqrt(sum(reshape(J, 3, []).^2, 1))';
                vPk  = i_local_maxima(Jmag, Surf.VertConn, nPeaks);
                for v = vPk(:)'
                    accT(end+1) = et(k);                              %#ok<AGROW>
                    accV(end+1) = v;                                  %#ok<AGROW>
                    accP(end+1,:) = Surf.Vertices(v, :);              %#ok<AGROW>
                    accH(end+1) = 1 + (Surf.Vertices(v,2) < 0);       %#ok<AGROW>  SCS Y>0=left
                    accS(end+1) = Jmag(v);                            %#ok<AGROW>
                end
            end
            if isempty(accT), continue; end
            Gp = bst_dynamics('NewGroup', [bandTok '_' ph]);
            % phase = NUMERIC alpha phase (rad); readable type is the label suffix.
            Gp.type = 'simple';  Gp.parent = winLabel;  Gp.phase = process_evt_refphase('PhaseValue', ph);
            Gp.times = accT;  Gp.epochs = ones(1, numel(accT));
            Gp.band = FreqRange;  Gp.bandName = bandTok;  Gp.Function = 'magnitude';
            Gp.vertices = accV;  Gp.pos = accP;  Gp.hemi = accH;  Gp.strength = accS;
            Gp.color = phaseColors.(ph);
            Gp.DataFile = R.DataFile;  Gp.ResultsFile = ResultsFile;  Gp.SurfaceFile = R.SurfaceFile;
            T = bst_dynamics('AddGroup', T, Gp);
        end

        OutFile = bst_process('GetNewFilename', bst_fileparts(R.DataFile), 'dynamics');
        bst_dynamics('Save', OutFile, T);
        nOcc = sum(arrayfun(@(g) numel(g.vertices), T.Groups));
        bst_report('Info', sProcess, sInputs(iFile), ...
            sprintf('%d atom groups (%d spatial occurrences), parent "%s". Saved: %s  (open with view_dynamics).', ...
            T.nGroups, nOcc, winLabel, OutFile));
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
