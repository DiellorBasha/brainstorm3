function varargout = process_eigenmodes_inverse( varargin )
% PROCESS_EIGENMODES_INVERSE: Source mapping via compressed eigenmode lead field.
%
% USAGE:  sProcess = process_eigenmodes_inverse('GetDescription')
%       OutputFiles = process_eigenmodes_inverse('Run', sProcess, sInputs)
%
% DESCRIPTION:
%     Computes MEG/EEG source estimates directly in eigenmode space by
%     compressing the lead field into the Laplace-Beltrami eigenmode basis.
%     Instead of mapping 300 sensors to 15000 vertices (underdetermined),
%     maps sensors to K eigenmodes (nearly determined: 300 sensors, K~200 modes).
%
%     The compressed forward model is L_tilde = L * Phi, where each column
%     is the sensor topography of eigenmode k. The eigenmode-space inverse
%     yields coefficients c_k(t) that describe the spatial frequency content
%     of cortical activity at each time point.
%
%     Outputs can be:
%     - Eigenmode coefficients (matrix file: rows = modes, cols = time)
%     - Reconstructed vertex sources (results file: u = Phi * c for visualization)
%     - Both (a results file + associated matrix file)
%
%     Requires precomputed eigenmodes on the cortex surface.
%
% SEE ALSO: bst_inverse_eigenmodes, in_tess_eigenmodes, process_eigenmodes

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
    sProcess.Comment     = 'Eigenmode source mapping';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 339;
    sProcess.Description = '';
    % Input / output
    sProcess.InputTypes  = {'data', 'raw'};
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % === INVERSE METHOD ===
    sProcess.options.method.Comment = {'MNE (minimum norm)', 'dSPM (noise-normalized)', ...
                                       'sLORETA (standardized)', ...
                                       'Inverse method:'; ...
                                       'mne', 'dspm', 'sloreta', ''};
    sProcess.options.method.Type    = 'radio_linelabel';
    sProcess.options.method.Value   = 'dspm';

    % === SEPARATOR ===
    sProcess.options.sep1.Type = 'separator';

    % === NUMBER OF MODES ===
    sProcess.options.nmodes.Comment = 'Number of eigenmodes (0 = auto, min of channels and available): ';
    sProcess.options.nmodes.Type    = 'value';
    sProcess.options.nmodes.Value   = {0, '', 0};

    % === SOURCE PRIOR ===
    sProcess.options.prioralpha.Comment = 'Source prior exponent alpha (0=flat, 1=1/f, 2=smooth): ';
    sProcess.options.prioralpha.Type    = 'value';
    sProcess.options.prioralpha.Value   = {0, '', 2};

    % === SNR ===
    sProcess.options.snr.Comment = 'Signal-to-noise ratio: ';
    sProcess.options.snr.Type    = 'value';
    sProcess.options.snr.Value   = {3, '', 1};

    % === SEPARATOR ===
    sProcess.options.sep2.Type = 'separator';

    % === OUTPUT TYPE ===
    sProcess.options.outputtype.Comment = {'Eigenmode coefficients only (matrix file)', ...
                                           'Reconstructed sources only (results file with kernel)', ...
                                           'Both (results file + matrix file)', ...
                                           'Output:'; ...
                                           'coefficients', 'sources', 'both', ''};
    sProcess.options.outputtype.Type    = 'radio_linelabel';
    sProcess.options.outputtype.Value   = 'both';

    % === INFO LABEL ===
    sProcess.options.label_info.Comment = ['<FONT color="#777777">Maps sensors directly to eigenmode space ' ...
        'via compressed lead field.<BR>Requires precomputed eigenmodes on the cortex surface.</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Method = upper(sProcess.options.method.Value);
    nModes = sProcess.options.nmodes.Value{1};
    alpha  = sProcess.options.prioralpha.Value{1};
    if nModes > 0
        modeStr = sprintf('%d modes', nModes);
    else
        modeStr = 'auto modes';
    end
    if alpha > 0
        Comment = sprintf('Eigenmode %s (%s, alpha=%.1f)', Method, modeStr, alpha);
    else
        Comment = sprintf('Eigenmode %s (%s)', Method, modeStr);
    end
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};

    % Get options
    Method     = lower(sProcess.options.method.Value);
    nModes     = sProcess.options.nmodes.Value{1};
    PriorAlpha = sProcess.options.prioralpha.Value{1};
    SNR        = sProcess.options.snr.Value{1};
    OutputType = lower(sProcess.options.outputtype.Value);

    % ===== GET STUDY AND HEAD MODEL =====
    % Get the study associated with the first input
    [sStudy, iStudy, ~, ~] = bst_get('Study', sInputs(1).iStudy);

    % Get head model
    if isempty(sStudy.iHeadModel) || sStudy.iHeadModel < 1
        bst_report('Error', sProcess, sInputs, 'No head model available for this study.');
        return;
    end
    HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;

    % Load head model to get surface file
    HeadModelMat = in_bst_headmodel(HeadModelFile, 0, 'HeadModelType', 'SurfaceFile');
    if ~strcmpi(HeadModelMat.HeadModelType, 'surface')
        bst_report('Error', sProcess, sInputs, 'Eigenmode inverse requires a surface head model.');
        return;
    end
    SurfaceFile = HeadModelMat.SurfaceFile;

    % ===== CHECK EIGENMODES =====
    [Eigenmodes, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if ~isComputed
        bst_report('Error', sProcess, sInputs, ...
            ['No eigenmodes on surface: ' SurfaceFile '. Run "Compute eigenmodes" first.']);
        return;
    end

    % Check vertex count consistency between head model and eigenmodes
    nVertEigen = size(Eigenmodes.Vectors, 1);
    HeadModelFull = in_bst_headmodel(HeadModelFile, 0, 'Gain');
    nVertHM = size(HeadModelFull.Gain, 2) / 3;  % unconstrained: 3 orientations per vertex
    if nVertEigen ~= nVertHM
        bst_report('Error', sProcess, sInputs, ...
            sprintf(['Head model has %d vertices but eigenmodes have %d vertices.\n' ...
            'This typically happens because computing eigenmodes repaired the surface mesh.\n' ...
            'Please recompute the head model (right-click study > Compute head model).'], ...
            nVertHM, nVertEigen));
        return;
    end

    % ===== GET NOISE COVARIANCE =====
    NoiseCovFile = '';
    if ~isempty(sStudy.NoiseCov) && ~isempty(sStudy.NoiseCov(1).FileName)
        NoiseCovFile = sStudy.NoiseCov(1).FileName;
    else
        bst_report('Warning', sProcess, sInputs, ...
            'No noise covariance found. Using identity noise model (not recommended for dSPM/sLORETA).');
    end

    % ===== GET CHANNEL FILE =====
    ChannelFile = bst_get('ChannelFileForStudy', sStudy.FileName);
    if isempty(ChannelFile)
        bst_report('Error', sProcess, sInputs, 'No channel file found.');
        return;
    end
    ChannelMat = in_bst_channel(ChannelFile);

    % Determine good channels (MEG only for now)
    % Load first data file to get channel flags
    DataMat = in_bst_data(sInputs(1).FileName, 'ChannelFlag');
    nAllChannels = length(ChannelMat.Channel);
    if isfield(DataMat, 'ChannelFlag') && ~isempty(DataMat.ChannelFlag)
        ChannelFlag = DataMat.ChannelFlag;
    else
        ChannelFlag = ones(nAllChannels, 1);
    end

    % Select MEG channels that are good
    iMEG = good_channel(ChannelMat.Channel, ChannelFlag, 'MEG');
    if isempty(iMEG)
        % Try EEG
        iMEG = good_channel(ChannelMat.Channel, ChannelFlag, 'EEG');
    end
    if isempty(iMEG)
        bst_report('Error', sProcess, sInputs, 'No good MEG or EEG channels found.');
        return;
    end

    GoodChannel = false(nAllChannels, 1);
    GoodChannel(iMEG) = true;

    % ===== COMPUTE EIGENMODE INVERSE =====
    bst_progress('start', 'Eigenmode inverse', 'Computing compressed lead field...', 0, 100);
    bst_progress('set', 5);

    [InvResults, errMsg] = bst_inverse_eigenmodes(HeadModelFile, SurfaceFile, NoiseCovFile, ...
        'Method', Method, ...
        'nModes', nModes, ...
        'PriorAlpha', PriorAlpha, ...
        'SNR', SNR, ...
        'GoodChannel', GoodChannel);

    if ~isempty(errMsg)
        bst_report('Error', sProcess, sInputs, errMsg);
        return;
    end

    K = InvResults.nModes;
    bst_progress('set', 30);

    % Report condition number improvement
    bst_report('Info', sProcess, sInputs, ...
        sprintf('Compressed lead field: %d channels -> %d modes (condition: %.1f vs full: %.1e)', ...
        sum(GoodChannel), K, InvResults.ConditionNumber, InvResults.ConditionNumberFull));

    % ===== LOAD EIGENMODES FOR RECONSTRUCTION =====
    Phi = double(Eigenmodes.Vectors(:, 1:K));   % [nVertices × K]
    lambdas = double(Eigenmodes.Values(1:K));

    % ===== PROCESS EACH INPUT FILE =====
    nInputs = length(sInputs);
    for iInput = 1:nInputs
        sInput = sInputs(iInput);
        bst_progress('set', round(30 + 70 * (iInput - 1) / nInputs));

        % Load data file
        DataMat = in_bst_data(sInput.FileName);
        isRaw = isstruct(DataMat.F);

        % Comment suffix for alpha
        alphaStr = '';
        if PriorAlpha > 0
            alphaStr = sprintf(', a=%.1f', PriorAlpha);
        end

        % ===== RAW DATA: kernel-only mode =====
        if isRaw
            % For raw files, we can only save kernels (like standard Brainstorm inverse).
            % Eigenmode coefficients cannot be precomputed for raw — user must import first.
            if strcmpi(OutputType, 'coefficients')
                bst_report('Warning', sProcess, sInput, ...
                    'Cannot compute eigenmode coefficients for raw files. Import the data first, or use "sources" or "both" output mode.');
                continue;
            end

            % Get study output directory
            [sStudyOut, iStudyOut] = bst_get('Study', sInput.iStudy);
            StudyDir = bst_fileparts(file_fullpath(sStudyOut.FileName));

            % --- Reconstructed sources (kernel-only for raw) ---
            VertexKernel = Phi * InvResults.ImagingKernel;  % [nVertices × nGoodChannels]

            ResMat = db_template('resultsmat');
            ResMat.ImagingKernel = VertexKernel;
            ResMat.ImageGridAmp  = [];
            ResMat.nComponents   = 1;
            ResMat.Comment       = sprintf('EigenInv %s (%d modes%s) | %s', ...
                upper(Method), K, alphaStr, sInput.Comment);
            ResMat.Function      = Method;
            ResMat.Time          = [];   % No time for raw kernel
            ResMat.DataFile      = sInput.FileName;
            ResMat.HeadModelFile = HeadModelFile;
            ResMat.HeadModelType = 'surface';
            ResMat.SurfaceFile   = SurfaceFile;
            ResMat.GoodChannel   = iMEG;
            ResMat.ChannelFlag   = ChannelFlag;
            ResMat.nAvg          = 1;
            ResMat.Leff          = 1;

            % History
            ResMat = bst_history('add', ResMat, 'eigenmodes_inverse', ...
                sprintf('Eigenmode %s inverse: %d modes, SNR=%.1f, alpha=%.1f', ...
                Method, K, SNR, PriorAlpha));
            ResMat = bst_history('add', ResMat, 'eigenmodes_inverse', ...
                sprintf('Condition: %.1f (full: %.1e). Prior: lambda^(-%.1f)', ...
                InvResults.ConditionNumber, InvResults.ConditionNumberFull, PriorAlpha));

            % Save
            OutputFile = bst_process('GetNewFilename', StudyDir, 'results_eigeninverse');
            bst_save(OutputFile, ResMat, 'v6');
            db_add_data(iStudyOut, OutputFile, ResMat);
            OutputFiles{end+1} = file_short(OutputFile);

            % Also save eigenmode kernel as matrix if "both" was requested
            if strcmpi(OutputType, 'both')
                MatrixMat = db_template('matrixmat');
                % Store the eigenmode-space kernel [K × nGoodChannels] so it can be
                % multiplied with imported data segments later
                MatrixMat.Value       = InvResults.ImagingKernel;
                MatrixMat.Time        = 1:size(InvResults.ImagingKernel, 2);
                MatrixMat.nAvg        = 1;
                MatrixMat.Leff        = 1;
                MatrixMat.SurfaceFile = SurfaceFile;
                MatrixMat.DisplayUnits = '';
                MatrixMat.Comment     = sprintf('EigenKernel %s (%d modes%s) | %s', ...
                    upper(Method), K, alphaStr, sInput.Comment);

                % Row descriptions
                RowNames = cell(K, 1);
                for k = 1:K
                    RowNames{k} = sprintf('Mode %d (lam=%.1f)', k, lambdas(k));
                end
                MatrixMat.Description = RowNames;

                % History
                MatrixMat = bst_history('add', MatrixMat, 'eigenmodes_inverse', ...
                    sprintf('Eigenmode kernel: %s, %d modes, SNR=%.1f, alpha=%.1f', ...
                    Method, K, SNR, PriorAlpha));

                OutputFile = bst_process('GetNewFilename', StudyDir, 'matrix_eigenkernel');
                bst_save(OutputFile, MatrixMat, 'v6');
                db_add_data(iStudyOut, OutputFile, MatrixMat);
                OutputFiles{end+1} = file_short(OutputFile);
            end

        % ===== IMPORTED DATA: full computation =====
        else
            F = double(DataMat.F(iMEG, :));    % [nGoodChannels × nTime]
            TimeVector = DataMat.Time;
            nTime = size(F, 2);

            % Apply eigenmode inverse: ĉ(t) = M̃ · d(t)
            Coeffs = InvResults.ImagingKernel * F;    % [K × nTime]

            % Get study output directory
            [sStudyOut, iStudyOut] = bst_get('Study', sInput.iStudy);
            StudyDir = bst_fileparts(file_fullpath(sStudyOut.FileName));

            % --- Eigenmode coefficients (matrix file) ---
            if ismember(OutputType, {'coefficients', 'both'})
                MatrixMat = db_template('matrixmat');
                MatrixMat.Value       = Coeffs;
                MatrixMat.Time        = TimeVector;
                MatrixMat.nAvg        = DataMat.nAvg;
                MatrixMat.Leff        = DataMat.Leff;
                MatrixMat.SurfaceFile = SurfaceFile;
                MatrixMat.DisplayUnits = '';
                MatrixMat.Comment     = sprintf('EigenInv %s (%d modes%s) | %s', ...
                    upper(Method), K, alphaStr, sInput.Comment);

                % Row descriptions
                RowNames = cell(K, 1);
                for k = 1:K
                    RowNames{k} = sprintf('Mode %d (lam=%.1f)', k, lambdas(k));
                end
                MatrixMat.Description = RowNames;

                % History
                MatrixMat = bst_history('add', MatrixMat, 'eigenmodes_inverse', ...
                    sprintf('Eigenmode %s inverse: %d modes, SNR=%.1f, alpha=%.1f', ...
                    Method, K, SNR, PriorAlpha));
                MatrixMat = bst_history('add', MatrixMat, 'eigenmodes_inverse', ...
                    sprintf('Condition number: %.1f (full: %.1e)', ...
                    InvResults.ConditionNumber, InvResults.ConditionNumberFull));
                MatrixMat = bst_history('add', MatrixMat, 'eigenmodes_inverse', ...
                    sprintf('Input: %s', sInput.FileName));

                % Save
                OutputFile = bst_process('GetNewFilename', StudyDir, 'matrix_eigeninverse');
                bst_save(OutputFile, MatrixMat, 'v6');
                db_add_data(iStudyOut, OutputFile, MatrixMat);
                OutputFiles{end+1} = file_short(OutputFile);
            end

            % --- Reconstructed sources (results file with kernel) ---
            if ismember(OutputType, {'sources', 'both'})
                VertexKernel = Phi * InvResults.ImagingKernel;  % [nVertices × nGoodChannels]

                ResMat = db_template('resultsmat');
                ResMat.ImagingKernel = VertexKernel;
                ResMat.ImageGridAmp  = [];
                ResMat.nComponents   = 1;
                ResMat.Comment       = sprintf('EigenInv %s (%d modes%s) | %s', ...
                    upper(Method), K, alphaStr, sInput.Comment);
                ResMat.Function      = Method;
                ResMat.Time          = TimeVector;
                ResMat.DataFile      = sInput.FileName;
                ResMat.HeadModelFile = HeadModelFile;
                ResMat.HeadModelType = 'surface';
                ResMat.SurfaceFile   = SurfaceFile;
                ResMat.GoodChannel   = iMEG;
                ResMat.ChannelFlag   = ChannelFlag;
                ResMat.nAvg          = DataMat.nAvg;
                ResMat.Leff          = DataMat.Leff;

                % History
                ResMat = bst_history('add', ResMat, 'eigenmodes_inverse', ...
                    sprintf('Eigenmode %s inverse: %d modes, SNR=%.1f, alpha=%.1f', ...
                    Method, K, SNR, PriorAlpha));
                ResMat = bst_history('add', ResMat, 'eigenmodes_inverse', ...
                    sprintf('Condition: %.1f (full: %.1e). Prior: lambda^(-%.1f)', ...
                    InvResults.ConditionNumber, InvResults.ConditionNumberFull, PriorAlpha));

                % Save
                OutputFile = bst_process('GetNewFilename', StudyDir, 'results_eigeninverse');
                bst_save(OutputFile, ResMat, 'v6');
                db_add_data(iStudyOut, OutputFile, ResMat);
                OutputFiles{end+1} = file_short(OutputFile);
            end
        end
    end

    bst_progress('set', 100);
    bst_progress('stop');
end
