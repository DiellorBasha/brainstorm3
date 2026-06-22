function varargout = process_eigenmodes_wiener( varargin )
% PROCESS_EIGENMODES_WIENER: Wiener (noise-floor) spectral filter of eigenmode coefficients.
%
% USAGE:  sProcess = process_eigenmodes_wiener('GetDescription')
%       OutputFiles = process_eigenmodes_wiener('Run', sProcess, sInputsA, sInputsB)
%
% DESCRIPTION:
%     Files A = data recording(s) (surface head model + eigenmodes required).
%     Files B = empty-room recording(s). Builds the data's eigenmode coefficient
%     time series and the data/noise Welch PSDs (shared with the denoise
%     process), forms a per-(k,f) Wiener gain G = max(max(P-Alpha*N,0)/P, Gmin)
%     via bst_eigenmodes_noisefloor, and applies it to the coefficient time
%     series in the frequency domain (bst_eigenmodes_wiener). Both inputs must
%     be imported (not raw).
%
% SEE ALSO: bst_eigenmodes_wiener, bst_eigenmodes_noisefloor, process_eigenmodes_denoise

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
%
% Authors: Diellor Basha, 2026

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'Eigenmode Wiener filter';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 336.95;   % after the coefficient filter (336.9)
    sProcess.Description = '';
    sProcess.InputTypes  = {'data', 'raw'};
    sProcess.OutputTypes = {'matrix'};
    sProcess.nInputs     = 2;        % Files A = data, Files B = empty-room
    sProcess.nMinFiles   = 1;

    sProcess.options.nmodes.Comment   = 'Number of eigenmodes (0 = auto, min of channels and available): ';
    sProcess.options.nmodes.Type      = 'value';
    sProcess.options.nmodes.Value     = {0, '', 0};

    sProcess.options.noisewin.Comment = 'Welch window length: ';
    sProcess.options.noisewin.Type    = 'value';
    sProcess.options.noisewin.Value   = {2, 's', 2};

    sProcess.options.alpha.Comment    = 'Over-subtraction factor alpha (>=0, >1 over-subtracts): ';
    sProcess.options.alpha.Type       = 'value';
    sProcess.options.alpha.Value      = {1, '', 2};

    sProcess.options.gainfloor.Comment = 'Gain floor Gmin (0..1): ';
    sProcess.options.gainfloor.Type    = 'value';
    sProcess.options.gainfloor.Value   = {0, '', 2};

    sProcess.options.domirror.Comment = 'Mirror signal edges (reduce FFT edge artifacts)';
    sProcess.options.domirror.Type    = 'checkbox';
    sProcess.options.domirror.Value   = 1;

    sProcess.options.dorecon.Comment  = 'Also reconstruct filtered vertex sources (Q = Phi * theta_filt)';
    sProcess.options.dorecon.Type     = 'checkbox';
    sProcess.options.dorecon.Value    = 0;

    sProcess.options.savegain.Comment = 'Also save the Wiener gain spectrum G(k,f)';
    sProcess.options.savegain.Type    = 'checkbox';
    sProcess.options.savegain.Value   = 0;

    sProcess.options.label_info.Comment = ['<FONT color="#777777">Files A = data, Files B = empty-room. ' ...
        'Estimates a per-mode Wiener gain G(&lambda;,&omega;) from the noise floor and applies it to the ' ...
        'eigenmode coefficient time series in the frequency domain.<BR>Both inputs must be imported; requires ' ...
        'a surface head model + eigenmodes.</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    a = sProcess.options.alpha.Value{1};
    g = sProcess.options.gainfloor.Value{1};
    Comment = sprintf('Eigenmode Wiener filter (a=%.1f, Gmin=%.2f)', a, g);
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputsA, sInputsB) %#ok<DEFNU>
    OutputFiles = {};
    if isempty(sInputsB)
        bst_report('Error', sProcess, sInputsA, 'Select the empty-room recording(s) as Files B.');
        return;
    end
    nModesOpt = sProcess.options.nmodes.Value{1};
    WinLen    = sProcess.options.noisewin.Value{1};
    Alpha     = sProcess.options.alpha.Value{1};
    GainFloor = sProcess.options.gainfloor.Value{1};
    DoMirror  = sProcess.options.domirror.Value;
    DoRecon   = sProcess.options.dorecon.Value;
    SaveGain  = sProcess.options.savegain.Value;

    % Shared setup: coefficients + data/noise PSDs (errors already reported inside)
    S = process_eigenmodes_denoise('GetCoeffsAndPSD', sProcess, sInputsA, sInputsB, nModesOpt, WinLen);
    if isempty(S), return; end

    % Wiener gain from the noise floor, applied in the frequency domain
    NF = bst_eigenmodes_noisefloor(S.Pdata, S.Nnoise, 'Alpha', Alpha, 'GainFloor', GainFloor);
    G  = NF.Gain;   % [K x nFreq] in [0,1]
    Filtered = bst_eigenmodes_wiener(S.Coeffs, S.sfreq, G, S.Freqs, 'Mirror', logical(DoMirror));

    RowNames = cell(S.K,1);
    for k = 1:S.K, RowNames{k} = sprintf('Mode %d (lam=%.1f)', k, S.lambdas(k)); end
    [sStudyOut, iStudyOut] = bst_get('Study', sInputsA(1).iStudy);
    StudyDir = bst_fileparts(file_fullpath(sStudyOut.FileName));

    % --- Filtered coefficients (matrix) ---
    Mout = db_template('matrixmat');
    Mout.Value       = Filtered;
    Mout.Time        = S.Time;
    Mout.Description  = RowNames;
    Mout.SurfaceFile = S.SurfaceFile;
    Mout.Comment     = sprintf('EigenWiener (%d modes, a=%.1f, Gmin=%.2f) | %s', S.K, Alpha, GainFloor, sInputsA(1).Comment);
    Mout.nAvg        = 1;
    Mout = bst_history('add', Mout, 'eigenmodes_wiener', Mout.Comment);
    OutFile = bst_process('GetNewFilename', StudyDir, 'matrix_eigenwiener');
    bst_save(OutFile, Mout, 'v6');
    db_add_data(iStudyOut, OutFile, Mout);
    OutputFiles{end+1} = file_short(OutFile);

    % --- Optional Wiener gain spectrum (timefreq sidecar; not in OutputFiles) ---
    if SaveGain
        TFmat = db_template('timefreqmat');
        TFmat.TF          = reshape(G, [S.K, 1, numel(S.Freqs)]);
        TFmat.Freqs       = S.Freqs;
        TFmat.Time        = [S.Time(1), S.Time(end)];
        TFmat.RowNames    = RowNames;
        TFmat.Measure     = 'power';
        TFmat.Method      = 'psd';
        TFmat.DataType    = 'matrix';
        TFmat.SurfaceFile = S.SurfaceFile;
        TFmat.Comment     = sprintf('EigenWienerGain (%d modes, a=%.1f, Gmin=%.2f) | %s', S.K, Alpha, GainFloor, sInputsA(1).Comment);
        TFmat.nAvg        = 1;
        TFmat = bst_history('add', TFmat, 'eigenmodes_wiener', TFmat.Comment);
        GainFile = bst_process('GetNewFilename', StudyDir, 'timefreq_eigenwienergain');
        bst_save(GainFile, TFmat, 'v6');
        db_add_data(iStudyOut, GainFile, TFmat);
    end

    % --- Optional filtered vertex reconstruction (results sidecar; not in OutputFiles) ---
    if DoRecon
        Q = S.Phi * Filtered;   % [nVert x nTime]
        ResMat = db_template('resultsmat');
        ResMat.ImageGridAmp  = Q;
        ResMat.ImagingKernel = [];
        ResMat.nComponents   = 1;
        ResMat.Time          = S.Time;
        ResMat.SurfaceFile   = S.SurfaceFile;
        ResMat.HeadModelType = 'surface';
        ResMat.ColormapType  = 'source';
        ResMat.Comment       = sprintf('EigenWiener recon (%d modes) | %s', S.K, sInputsA(1).Comment);
        ResMat.nAvg          = 1;
        ResMat = bst_history('add', ResMat, 'eigenmodes_wiener', ResMat.Comment);
        ReconFile = bst_process('GetNewFilename', StudyDir, 'results_eigenwiener');
        bst_save(ReconFile, ResMat, 'v6');
        db_add_data(iStudyOut, ReconFile, ResMat);
    end

    bst_report('Info', sProcess, sInputsA, sprintf('Wiener-filtered %d eigenmode coefficients (alpha=%.1f, Gmin=%.2f).', S.K, Alpha, GainFloor));
end
