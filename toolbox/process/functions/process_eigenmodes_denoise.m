function varargout = process_eigenmodes_denoise( varargin )
% PROCESS_EIGENMODES_DENOISE: Joint (lambda,omega) noise floor from empty-room recordings.
%
% USAGE:  sProcess = process_eigenmodes_denoise('GetDescription')
%       OutputFiles = process_eigenmodes_denoise('Run', sProcess, sInputsA, sInputsB)
%                 S = process_eigenmodes_denoise('GetCoeffsAndPSD', sProcess, sInputsA, sInputsB, nModesOpt, WinLen)
%
% DESCRIPTION:
%     Files A = data recording(s) (must have a surface head model + eigenmodes).
%     Files B = empty-room recording(s). Builds the data's eigenmode transform
%     kernel A = pinv(L*Phi), applies it to BOTH recordings on their common good
%     channels, Welch-PSDs both, then computes SNR(k,f), a power-spectral-
%     subtraction cleaned spectrum, and a reliable-mode cutoff K*(f).
%     Subtraction is on POWER (averaged PSD), never complex coefficients.
%     Both recordings must be imported (not raw). Requires precomputed eigenmodes.
%
% SEE ALSO: bst_eigenmodes_noisefloor, bst_eigenmodes_transform, in_tess_eigenmodes

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
    sProcess.Comment     = 'Eigenmode noise-floor denoising';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 336.6;   % just after the eigenmode transform (336.5)
    sProcess.Description = '';
    sProcess.InputTypes  = {'data', 'raw'};
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 2;       % Files A = data, Files B = empty-room
    sProcess.nMinFiles   = 1;

    sProcess.options.nmodes.Comment    = 'Number of eigenmodes (0 = auto, min of channels and available): ';
    sProcess.options.nmodes.Type       = 'value';
    sProcess.options.nmodes.Value      = {0, '', 0};

    sProcess.options.noisewin.Comment  = 'Welch window length: ';
    sProcess.options.noisewin.Type     = 'value';
    sProcess.options.noisewin.Value    = {2, 's', 2};

    sProcess.options.alpha.Comment     = 'Over-subtraction factor alpha (>=0, >1 over-subtracts): ';
    sProcess.options.alpha.Type        = 'value';
    sProcess.options.alpha.Value       = {1, '', 2};

    sProcess.options.snrthresh.Comment = 'Reliable-mode SNR threshold (linear): ';
    sProcess.options.snrthresh.Type    = 'value';
    sProcess.options.snrthresh.Value   = {1, '', 2};

    sProcess.options.floorfrac.Comment = 'Spectral floor (fraction of noise): ';
    sProcess.options.floorfrac.Type    = 'value';
    sProcess.options.floorfrac.Value   = {0, '', 2};

    sProcess.options.label_info.Comment = ['<FONT color="#777777">Files A = data, Files B = empty-room. ' ...
        'Subtraction is on Welch-averaged power (not complex coefficients).<BR>' ...
        'Outputs an SNR(&lambda;,&omega;) spectrum + cleaned power spectrum. Both inputs must be imported.</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    a = sProcess.options.alpha.Value{1};
    Comment = sprintf('Eigenmode noise-floor denoising (alpha=%.1f)', a);
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
    SnrThresh = sProcess.options.snrthresh.Value{1};
    FloorFrac = sProcess.options.floorfrac.Value{1};

    S = process_eigenmodes_denoise('GetCoeffsAndPSD', sProcess, sInputsA, sInputsB, nModesOpt, WinLen);
    if isempty(S), return; end

    % ===== COMBINE =====
    Out = bst_eigenmodes_noisefloor(S.Pdata, S.Nnoise, 'Alpha', Alpha, 'Floor', FloorFrac, 'SnrThresh', SnrThresh);

    % ===== SAVE =====
    RowNames = cell(S.K,1);
    for k = 1:S.K, RowNames{k} = sprintf('Mode %d (lam=%.1f)', k, S.lambdas(k)); end
    [sStudyOut, iStudyOut] = bst_get('Study', sInputsA(1).iStudy);
    StudyDir = bst_fileparts(file_fullpath(sStudyOut.FileName));

    OutputFiles{end+1} = SaveTF(Out.SNR, S.Freqs, RowNames, S.Time, StudyDir, iStudyOut, ...
        sprintf('EigenSNR (%d modes) | %s', S.K, sInputsA(1).Comment), 'timefreq_eigensnr', S.SurfaceFile); %#ok<AGROW>
    OutputFiles{end+1} = SaveTF(Out.CleanPSD, S.Freqs, RowNames, S.Time, StudyDir, iStudyOut, ...
        sprintf('EigenCleanPSD (%d modes, a=%.1f) | %s', S.K, Alpha, sInputsA(1).Comment), 'timefreq_eigencleanpsd', S.SurfaceFile); %#ok<AGROW>

    bst_report('Info', sProcess, sInputsA, sprintf('Denoise: %d modes (condition %.1f), median reliable-mode cutoff K*=%d at SNR>=%.1f.', ...
        S.K, S.Info.ConditionNumber, round(median(Out.Kstar)), SnrThresh));
end


%% ===== SHARED: COMMON-CHANNEL KERNEL, COEFFICIENTS, DATA/NOISE PSDs =====
function S = GetCoeffsAndPSD(sProcess, sInputsA, sInputsB, nModesOpt, WinLen) %#ok<DEFNU>
    S = [];
    % ===== DATA STUDY: head model, surface, eigenmodes, gain =====
    [sStudyA,~,~,~] = bst_get('Study', sInputsA(1).iStudy);
    if isempty(sStudyA.iHeadModel) || sStudyA.iHeadModel < 1
        bst_report('Error', sProcess, sInputsA, 'No head model for the data study (Files A).');
        return;
    end
    HeadModelFile = sStudyA.HeadModel(sStudyA.iHeadModel).FileName;
    HMmeta = in_bst_headmodel(HeadModelFile, 0, 'HeadModelType', 'SurfaceFile');
    if ~strcmpi(HMmeta.HeadModelType, 'surface')
        bst_report('Error', sProcess, sInputsA, 'Eigenmode denoise requires a surface head model.');
        return;
    end
    SurfaceFile = HMmeta.SurfaceFile;
    [Em, isC] = in_tess_eigenmodes(SurfaceFile);
    if ~isC
        bst_report('Error', sProcess, sInputsA, ['No eigenmodes on surface: ' SurfaceFile '. Run "Compute eigenmodes" first.']);
        return;
    end
    HM = in_bst_headmodel(HeadModelFile, 1);   % ApplyOrient=1 -> [nch x nVert]
    Gain = double(HM.Gain);
    if size(Gain,2) ~= size(Em.Vectors,1)
        bst_report('Error', sProcess, sInputsA, sprintf('Head model (%d) / eigenmode (%d) vertex mismatch.', size(Gain,2), size(Em.Vectors,1)));
        return;
    end

    % ===== CHANNELS: common good channels (by name) between data and noise =====
    ChanFileA = bst_get('ChannelFileForStudy', sStudyA.FileName);
    if isempty(ChanFileA)
        bst_report('Error', sProcess, sInputsA, 'No channel file for the data study (Files A).');
        return;
    end
    ChA = in_bst_channel(ChanFileA);
    DA  = in_bst_data(sInputsA(1).FileName);
    if isstruct(DA.F)
        bst_report('Error', sProcess, sInputsA, 'Files A must be imported data (not raw). Import a block first.');
        return;
    end
    iA = good_channel(ChA.Channel, DA.ChannelFlag, 'MEG');
    if isempty(iA), iA = good_channel(ChA.Channel, DA.ChannelFlag, 'EEG'); end
    if isempty(iA)
        bst_report('Error', sProcess, sInputsA, 'No good MEG or EEG channels found in data (Files A).');
        return;
    end

    [sStudyB,~,~,~] = bst_get('Study', sInputsB(1).iStudy);
    ChanFileB = bst_get('ChannelFileForStudy', sStudyB.FileName);
    if isempty(ChanFileB)
        bst_report('Error', sProcess, sInputsB, 'No channel file for the empty-room study (Files B).');
        return;
    end
    ChB = in_bst_channel(ChanFileB);
    DB  = in_bst_data(sInputsB(1).FileName);
    if isstruct(DB.F)
        bst_report('Error', sProcess, sInputsB, 'Files B (empty-room) must be imported data (not raw). Import a block first.');
        return;
    end
    iB = good_channel(ChB.Channel, DB.ChannelFlag, 'MEG');
    if isempty(iB), iB = good_channel(ChB.Channel, DB.ChannelFlag, 'EEG'); end
    if isempty(iB)
        bst_report('Error', sProcess, sInputsB, 'No good MEG or EEG channels found in empty-room (Files B).');
        return;
    end

    [~, ia, ib] = intersect({ChA.Channel(iA).Name}, {ChB.Channel(iB).Name}, 'stable');
    if isempty(ia)
        bst_report('Error', sProcess, sInputsA, 'No common good channels between data and empty-room.');
        return;
    end
    iCommonA = iA(ia);
    iCommonB = iB(ib);

    % ===== KERNEL ON COMMON CHANNELS =====
    nCh = numel(iCommonA);
    if isempty(nModesOpt) || nModesOpt <= 0
        K = min(nCh, Em.nModes);
    else
        K = min(nModesOpt, Em.nModes);
    end
    Phi     = double(Em.Vectors(:, 1:K));
    lambdas = double(Em.Values(1:K));
    lambdas = lambdas(:);
    [A, Info] = bst_eigenmodes_transform(Gain(iCommonA, :), Phi);   % [K x nCh]

    % ===== COEFFICIENTS + WELCH PSDs (density, same window) =====
    thD = A * double(DA.F(iCommonA, :));
    thN = A * double(DB.F(iCommonB, :));
    sfA = 1 / (DA.Time(2) - DA.Time(1));
    sfB = 1 / (DB.Time(2) - DB.Time(1));
    [TFd, Fv]  = bst_psd(thD, sfA, WinLen, 50, [], [], [], 'physical');
    [TFn, Fvn] = bst_psd(thN, sfB, WinLen, 50, [], [], [], 'physical');
    if numel(Fv) ~= numel(Fvn) || max(abs(Fv(:) - Fvn(:))) > 1e-6
        bst_report('Error', sProcess, sInputsA, 'Data and empty-room PSD frequency grids differ (different sampling rate or window).');
        return;
    end

    % ===== PACK =====
    S = struct();
    S.Coeffs      = thD;                  % [K x nTime] data coefficient time series
    S.Pdata       = reshape(TFd, K, []);  % [K x nFreq]
    S.Nnoise      = reshape(TFn, K, []);  % [K x nFreq]
    S.Freqs       = Fv(:)';               % [1 x nFreq]
    S.lambdas     = lambdas;              % [K x 1]
    S.Phi         = Phi;                  % [nVert x K]
    S.SurfaceFile = SurfaceFile;
    S.K           = K;
    S.Time        = DA.Time;              % [1 x nTime]
    S.sfreq       = sfA;
    S.Info        = Info;
end


%% ===== SAVE ONE (mode x freq) MATRIX AS A TIMEFREQ FILE =====
function OutFile = SaveTF(M2d, Freqs, RowNames, Time, StudyDir, iStudyOut, Comment, prefix, SurfaceFile)
    [K, nFreq] = size(M2d);
    TFmat = db_template('timefreqmat');
    TFmat.TF          = reshape(M2d, [K, 1, nFreq]);
    TFmat.Freqs       = Freqs(:)';
    TFmat.Time        = [Time(1), Time(end)];
    TFmat.RowNames    = RowNames;
    TFmat.Measure     = 'power';
    TFmat.Method      = 'psd';
    TFmat.DataType    = 'matrix';
    TFmat.SurfaceFile = SurfaceFile;
    TFmat.Comment     = Comment;
    TFmat.nAvg        = 1;
    TFmat = bst_history('add', TFmat, 'eigenmodes_denoise', Comment);
    FullFile = bst_process('GetNewFilename', StudyDir, prefix);
    bst_save(FullFile, TFmat, 'v6');
    db_add_data(iStudyOut, FullFile, TFmat);
    OutFile = file_short(FullFile);
end
