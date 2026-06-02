function varargout = process_eigenmodes_inverse( varargin )
% PROCESS_EIGENMODES_INVERSE: Mode-space source mapping on a composed eigenmode leadfield.
%
% Consumes a composed eigenmode head model (headmodel_eigenmode_*.mat, from
% process_eigenmode_leadfield) and produces eigenmode coefficients (matrix) and/or
% reconstructed cortex sources (kernel-only results: ImagingKernel = Phi*M~).
%
% Authors: Diellor Basha, 2026
eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'Eigenmode source mapping';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 339;
    sProcess.InputTypes  = {'data', 'raw'};
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.options.method.Comment = {'MNE (minimum norm)', 'dSPM (noise-normalized)', ...
        'sLORETA (standardized)', 'Inverse method:'; 'mne', 'dspm', 'sloreta', ''};
    sProcess.options.method.Type  = 'radio_linelabel';
    sProcess.options.method.Value = 'dspm';
    sProcess.options.prior.Comment = {'Log (2026)', 'Flat (none)', 'Power (1/f)', 'Spectral prior:'; ...
        'log', 'flat', 'power', ''};
    sProcess.options.prior.Type  = 'radio_linelabel';
    sProcess.options.prior.Value = 'log';
    sProcess.options.snr.Comment = 'Signal-to-noise ratio: ';
    sProcess.options.snr.Type    = 'value';
    sProcess.options.snr.Value   = {3, '', 1};
    sProcess.options.nmodes.Comment = 'Number of eigenmodes (0 = all in head model): ';
    sProcess.options.nmodes.Type    = 'value';
    sProcess.options.nmodes.Value   = {0, '', 0};
    sProcess.options.outputtype.Comment = {'Coefficients (matrix)', 'Sources (results)', 'Both', 'Output:'; ...
        'coefficients', 'sources', 'both', ''};
    sProcess.options.outputtype.Type  = 'radio_linelabel';
    sProcess.options.outputtype.Value = 'both';
end

function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sprintf('Eigenmode %s (%s prior)', upper(sProcess.options.method.Value), ...
        sProcess.options.prior.Value);
end

function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    Method     = lower(sProcess.options.method.Value);
    Prior      = lower(sProcess.options.prior.Value);
    SNR        = sProcess.options.snr.Value{1};
    nModes     = sProcess.options.nmodes.Value{1};
    OutputType = lower(sProcess.options.outputtype.Value);

    [sStudy, ~] = bst_get('Study', sInputs(1).iStudy);
    if isempty(sStudy.iHeadModel) || sStudy.iHeadModel < 1
        bst_report('Error', sProcess, sInputs, 'No head model for this study.'); return;
    end
    HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;
    HM = in_bst_headmodel(HeadModelFile, 0);
    if ~isfield(HM,'isEigenmode') || ~HM.isEigenmode
        bst_report('Error', sProcess, sInputs, ...
            'Active head model is not an eigenmode leadfield. Run "Compute eigenmode leadfield" first.'); return;
    end

    NoiseCovFile = '';
    if ~isempty(sStudy.NoiseCov) && ~isempty(sStudy.NoiseCov(1).FileName)
        NoiseCovFile = sStudy.NoiseCov(1).FileName;
    else
        bst_report('Warning', sProcess, sInputs, 'No noise covariance: using identity whitening.');
    end
    ChannelFile = bst_get('ChannelFileForStudy', sStudy.FileName);
    if isempty(ChannelFile)
        bst_report('Error', sProcess, sInputs, 'No channel file found for this study.'); return;
    end
    ChannelMat  = in_bst_channel(ChannelFile);

    % Good channels from the first input's flags (MEG, else EEG)
    DataFlag = in_bst_data(sInputs(1).FileName, 'ChannelFlag');
    ChannelFlag = ones(numel(ChannelMat.Channel), 1);
    if isfield(DataFlag,'ChannelFlag') && ~isempty(DataFlag.ChannelFlag)
        ChannelFlag = DataFlag.ChannelFlag;
    end
    iSel = good_channel(ChannelMat.Channel, ChannelFlag, 'MEG');
    if isempty(iSel); iSel = good_channel(ChannelMat.Channel, ChannelFlag, 'EEG'); end
    if isempty(iSel); bst_report('Error', sProcess, sInputs, 'No good MEG/EEG channels.'); return; end
    GoodChannel = false(numel(ChannelMat.Channel), 1); GoodChannel(iSel) = true;

    % Solve
    [Inv, errMsg] = bst_inverse_eigenmodes(HeadModelFile, NoiseCovFile, ChannelFile, GoodChannel, ...
        'Method', Method, 'Prior', Prior, 'SNR', SNR, 'nModes', nModes);
    if ~isempty(errMsg); bst_report('Error', sProcess, sInputs, errMsg); return; end
    K = Inv.nModes;

    % Eigenmodes (Phi) for reconstruction
    [Eig, isComputed] = in_tess_eigenmodes(HM.SurfaceFile);
    if ~isComputed
        bst_report('Error', sProcess, sInputs, ...
            ['No eigenmodes on surface: ' HM.SurfaceFile '. Run "Compute eigenmodes" first.']); return;
    end
    Phi = double(Eig.Vectors(:, 1:K));               % [nVert x K]
    lambdas = Inv.Eigenvalues;

    for iInput = 1:numel(sInputs)
        sInput = sInputs(iInput);
        if sInput.iStudy ~= sInputs(1).iStudy
            bst_report('Warning', sProcess, sInput, ...
                'Skipped: input is from a different study than the head model. Run the process per study.');
            continue;
        end
        DataMat = in_bst_data(sInput.FileName);
        isRaw = isstruct(DataMat.F);
        [sStudyOut, iStudyOut] = bst_get('Study', sInput.iStudy);
        StudyDir = bst_fileparts(file_fullpath(sStudyOut.FileName));

        % Cortex results node (kernel-only): ImagingKernel = Phi * M~
        if ismember(OutputType, {'sources','both'})
            ResMat = db_template('resultsmat');
            ResMat.ImagingKernel = Phi * Inv.ImagingKernel;     % [nVert x nGoodCh]
            ResMat.ImageGridAmp  = [];
            ResMat.nComponents   = 1;
            ResMat.Comment       = sprintf('Eigenmode %s (%d modes, %s) | %s', ...
                upper(Method), K, Prior, sInput.Comment);
            ResMat.Function      = ['eigenmode_' Method];
            ResMat.Time          = DataMat.Time;
            if isRaw; ResMat.Time = []; end   % raw kernel: viewer fetches time from the data file
            ResMat.DataFile      = sInput.FileName;
            ResMat.HeadModelFile = HeadModelFile;
            ResMat.HeadModelType = 'surface';
            ResMat.SurfaceFile   = HM.SurfaceFile;
            ResMat.GoodChannel   = iSel;
            ResMat.ChannelFlag   = ChannelFlag;
            ResMat.nAvg          = DataMat.nAvg; ResMat.Leff = DataMat.Leff;
            ResMat = bst_history('add', ResMat, 'eigenmodes_inverse', ...
                sprintf('Eigenmode %s, %d modes, prior=%s, SNR=%.1f', Method, K, Prior, SNR));
            OutFile = bst_process('GetNewFilename', StudyDir, 'results_eigeninverse');
            bst_save(OutFile, ResMat, 'v6');
            db_add_data(iStudyOut, OutFile, ResMat);
            OutputFiles{end+1} = file_short(OutFile); %#ok<AGROW>
        end

        % Coefficients matrix node: theta = M~ * d (imported data only)
        if ismember(OutputType, {'coefficients','both'}) && ~isRaw
            theta = Inv.ImagingKernel * double(DataMat.F(iSel, :));   % [K x nTime]
            MatMat = db_template('matrixmat');
            MatMat.Value       = theta;
            MatMat.Time        = DataMat.Time;
            MatMat.nAvg        = DataMat.nAvg; MatMat.Leff = DataMat.Leff;
            MatMat.SurfaceFile = HM.SurfaceFile;
            MatMat.Comment     = sprintf('EigenCoeffs %s (%d modes, %s) | %s', ...
                upper(Method), K, Prior, sInput.Comment);
            RowNames = cell(K,1);
            for k = 1:K; RowNames{k} = sprintf('Mode %d (lam=%.3g)', k, lambdas(k)); end
            MatMat.Description = RowNames;
            MatMat = bst_history('add', MatMat, 'eigenmodes_inverse', ...
                sprintf('Eigenmode coefficients %s, %d modes, prior=%s', Method, K, Prior));
            OutFile = bst_process('GetNewFilename', StudyDir, 'matrix_eigencoeffs');
            bst_save(OutFile, MatMat, 'v6');
            db_add_data(iStudyOut, OutFile, MatMat);
            OutputFiles{end+1} = file_short(OutFile); %#ok<AGROW>
        elseif ismember(OutputType, {'coefficients','both'}) && isRaw
            bst_report('Warning', sProcess, sInput, 'Coefficients require imported data; skipped for raw.');
        end
    end
end
