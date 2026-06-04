function varargout = process_eigenmodes_transform( varargin )
% PROCESS_EIGENMODES_TRANSFORM: Unregularized sensor->eigenmode spatial transform.
%
% USAGE:  sProcess = process_eigenmodes_transform('GetDescription')
%       OutputFiles = process_eigenmodes_transform('Run', sProcess, sInputs)
%
% DESCRIPTION:
%     Maps sensor recordings directly to eigenmode coefficients via the
%     unregularized transform A = pinv(L*Phi) (see bst_eigenmodes_transform).
%     The output is the raw eigenmode coefficient time series Theta [K x nTime]
%     as a Brainstorm matrix file -- the spatial analogue of a time series ready
%     for FFT. Feed that matrix file into "Frequency > FFT" (process_fft) to see
%     the joint (lambda, omega) spectrum. No regularization is applied; high-mode
%     coefficients are noisy by design.
%
%     Optionally also reconstructs raw vertex sources Q = Phi * Theta.
%     Requires precomputed eigenmodes on the cortex surface.
%
% SEE ALSO: bst_eigenmodes_transform, in_tess_eigenmodes, process_eigenmodes_inverse

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
    sProcess.Comment     = 'Eigenmode transform (spatial FFT)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 336.5;   % leads the eigenmode Sources cluster (filter 337, spectrum 338, inverse 339)
    sProcess.Description = '';
    sProcess.InputTypes  = {'data', 'raw'};
    % Primary output is a matrix file (matrix_eigentransform); OutputTypes mirrors the
    % input convention used by the sibling process_eigenmodes_inverse.
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % === NUMBER OF MODES ===
    sProcess.options.nmodes.Comment = 'Number of eigenmodes (0 = auto, min of channels and available): ';
    sProcess.options.nmodes.Type    = 'value';
    sProcess.options.nmodes.Value   = {0, '', 0};

    % === OPTIONAL VERTEX RECONSTRUCTION ===
    sProcess.options.dorecon.Comment = 'Also reconstruct raw vertex sources (Q = Phi * Theta)';
    sProcess.options.dorecon.Type    = 'checkbox';
    sProcess.options.dorecon.Value   = 0;

    % === INFO LABEL ===
    sProcess.options.label_info.Comment = ['<FONT color="#777777">Unregularized sensor&rarr;eigenmode transform ' ...
        '(A = pinv(L&middot;&Phi;)).<BR>Coefficients are raw and noisy at high modes by design &mdash; ' ...
        'run FFT on the matrix output to see the (&lambda;,&omega;) spectrum.<BR>' ...
        'Requires precomputed eigenmodes on the cortex surface.</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    nModes = sProcess.options.nmodes.Value{1};
    if nModes > 0
        Comment = sprintf('Eigenmode transform (%d modes)', nModes);
    else
        Comment = 'Eigenmode transform (auto modes)';
    end
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};

    nModesOpt = sProcess.options.nmodes.Value{1};
    DoRecon   = sProcess.options.dorecon.Value;

    % ===== STUDY + HEAD MODEL + SURFACE =====
    [sStudy, ~, ~, ~] = bst_get('Study', sInputs(1).iStudy);
    if isempty(sStudy.iHeadModel) || sStudy.iHeadModel < 1
        bst_report('Error', sProcess, sInputs, 'No head model available for this study.');
        return;
    end
    HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;
    HeadModelMat  = in_bst_headmodel(HeadModelFile, 0, 'HeadModelType', 'SurfaceFile');
    if ~strcmpi(HeadModelMat.HeadModelType, 'surface')
        bst_report('Error', sProcess, sInputs, 'Eigenmode transform requires a surface head model.');
        return;
    end
    SurfaceFile = HeadModelMat.SurfaceFile;

    % ===== EIGENMODES =====
    [Eigenmodes, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if ~isComputed
        bst_report('Error', sProcess, sInputs, ...
            ['No eigenmodes on surface: ' SurfaceFile '. Run "Compute eigenmodes" first.']);
        return;
    end
    nVertEigen = size(Eigenmodes.Vectors, 1);

    % ===== CONSTRAINED GAIN (fixed orientation) =====
    % HeadModelType is guaranteed 'surface' above, so ApplyOrient=1 collapses the
    % stored [nch x 3*nVert] gain to the constrained [nch x nVert]; size(Gain,2)==nVert.
    HM   = in_bst_headmodel(HeadModelFile, 1);   % ApplyOrient=1 -> [nch x nVert]
    Gain = double(HM.Gain);
    nVertHM = size(Gain, 2);
    if nVertEigen ~= nVertHM
        bst_report('Error', sProcess, sInputs, ...
            sprintf(['Head model has %d vertices but eigenmodes have %d vertices.\n' ...
            'Recompute the head model (right-click study > Compute head model).'], ...
            nVertHM, nVertEigen));
        return;
    end

    % ===== CHANNELS =====
    ChannelFile = bst_get('ChannelFileForStudy', sStudy.FileName);
    if isempty(ChannelFile)
        bst_report('Error', sProcess, sInputs, 'No channel file found.');
        return;
    end
    ChannelMat   = in_bst_channel(ChannelFile);
    % ChannelFlag is read from the first input; the resulting kernel is reused for
    % every input, so a consistent good-channel set across inputs is assumed.
    DataMat0     = in_bst_data(sInputs(1).FileName, 'ChannelFlag');
    nAllChannels = length(ChannelMat.Channel);
    if isfield(DataMat0, 'ChannelFlag') && ~isempty(DataMat0.ChannelFlag)
        ChannelFlag = DataMat0.ChannelFlag;
    else
        ChannelFlag = ones(nAllChannels, 1);
    end
    iMEG = good_channel(ChannelMat.Channel, ChannelFlag, 'MEG');
    if isempty(iMEG)
        iMEG = good_channel(ChannelMat.Channel, ChannelFlag, 'EEG');
    end
    if isempty(iMEG)
        bst_report('Error', sProcess, sInputs, 'No good MEG or EEG channels found.');
        return;
    end

    % ===== BUILD TRANSFORM KERNEL =====
    nCh         = numel(iMEG);
    K_available = Eigenmodes.nModes;
    if isempty(nModesOpt) || nModesOpt <= 0
        K = min(nCh, K_available);
    else
        K = min(nModesOpt, K_available);
    end
    % Canonical selection (whole-brain lowest spatial frequencies, never one hemisphere).
    if isfield(Eigenmodes,'Order') && ~isempty(Eigenmodes.Order)
        order = double(Eigenmodes.Order(:));
    else
        [~, order] = sort(double(Eigenmodes.Values(:)), 'ascend');
    end
    sel     = order(1:K);
    Phi     = double(Eigenmodes.Vectors(:, sel));
    lambdas = double(Eigenmodes.Values(sel));

    [Kernel, Info] = bst_eigenmodes_transform(Gain(iMEG, :), Phi);   % [K x nCh]

    bst_report('Info', sProcess, sInputs, ...
        sprintf('Eigenmode transform: %d channels -> %d modes (rank %d, condition %.1f)', ...
        nCh, K, Info.Rank, Info.ConditionNumber));

    % Row descriptions (reused across matrix outputs)
    RowNames = cell(K, 1);
    for k = 1:K
        RowNames{k} = sprintf('Mode %d (lam=%.1f)', k, lambdas(k));
    end

    % ===== PER-INPUT =====
    nInputs = numel(sInputs);
    for iInput = 1:nInputs
        sInput = sInputs(iInput);
        [sStudyOut, iStudyOut] = bst_get('Study', sInput.iStudy);
        StudyDir = bst_fileparts(file_fullpath(sStudyOut.FileName));

        DataMat = in_bst_data(sInput.FileName);
        isRaw   = isstruct(DataMat.F);

        if isRaw
            % Raw: cannot precompute coefficients; save kernel-only results.
            ResMat = db_template('resultsmat');
            ResMat.ImagingKernel = Phi * Kernel;          % [nVert x nCh]
            ResMat.ImageGridAmp  = [];
            ResMat.nComponents   = 1;
            ResMat.Comment       = sprintf('EigenTransform (%d modes) | %s', K, sInput.Comment);
            ResMat.Function      = 'eigentransform';
            ResMat.Time          = [];
            ResMat.DataFile      = sInput.FileName;
            ResMat.HeadModelFile = HeadModelFile;
            ResMat.HeadModelType = 'surface';
            ResMat.SurfaceFile   = SurfaceFile;
            ResMat.GoodChannel   = iMEG;
            ResMat.ChannelFlag   = ChannelFlag;
            ResMat.nAvg          = 1;
            ResMat.Leff          = 1;
            ResMat = bst_history('add', ResMat, 'eigenmodes_transform', ...
                sprintf('Unregularized eigenmode transform: %d modes, rank %d', K, Info.Rank));

            OutputFile = bst_process('GetNewFilename', StudyDir, 'results_eigentransform');
            bst_save(OutputFile, ResMat, 'v6');
            db_add_data(iStudyOut, OutputFile, ResMat);
            OutputFiles{end+1} = file_short(OutputFile); %#ok<AGROW>

            bst_report('Warning', sProcess, sInput, ...
                'Raw file: saved kernel only. Import the data to compute eigenmode coefficients and their FFT.');
        else
            F     = double(DataMat.F(iMEG, :));   % [nCh x nTime]
            Theta = Kernel * F;                   % [K x nTime]

            % --- Coefficients matrix file (this is the FFT input) ---
            MatrixMat = db_template('matrixmat');
            MatrixMat.Value        = Theta;
            MatrixMat.Time         = DataMat.Time;
            MatrixMat.nAvg         = DataMat.nAvg;
            MatrixMat.Leff         = DataMat.Leff;
            MatrixMat.SurfaceFile  = SurfaceFile;
            MatrixMat.DisplayUnits = '';
            MatrixMat.Description  = RowNames;
            MatrixMat.Comment      = sprintf('EigenTransform (%d modes) | %s', K, sInput.Comment);
            MatrixMat = bst_history('add', MatrixMat, 'eigenmodes_transform', ...
                sprintf('Unregularized eigenmode transform: %d modes, rank %d, condition %.1f', ...
                K, Info.Rank, Info.ConditionNumber));
            MatrixMat = bst_history('add', MatrixMat, 'eigenmodes_transform', ...
                sprintf('Input: %s', sInput.FileName));

            OutputFile = bst_process('GetNewFilename', StudyDir, 'matrix_eigentransform');
            bst_save(OutputFile, MatrixMat, 'v6');
            db_add_data(iStudyOut, OutputFile, MatrixMat);
            OutputFiles{end+1} = file_short(OutputFile); %#ok<AGROW>

            % --- Optional raw vertex reconstruction ---
            if DoRecon
                ResMat = db_template('resultsmat');
                ResMat.ImagingKernel = Phi * Kernel;       % [nVert x nCh]
                ResMat.ImageGridAmp  = [];
                ResMat.nComponents   = 1;
                ResMat.Comment       = sprintf('EigenTransform recon (%d modes) | %s', K, sInput.Comment);
                ResMat.Function      = 'eigentransform';
                ResMat.Time          = DataMat.Time;
                ResMat.DataFile      = sInput.FileName;
                ResMat.HeadModelFile = HeadModelFile;
                ResMat.HeadModelType = 'surface';
                ResMat.SurfaceFile   = SurfaceFile;
                ResMat.GoodChannel   = iMEG;
                ResMat.ChannelFlag   = ChannelFlag;
                ResMat.nAvg          = DataMat.nAvg;
                ResMat.Leff          = DataMat.Leff;
                ResMat = bst_history('add', ResMat, 'eigenmodes_transform', ...
                    sprintf('Raw vertex reconstruction from %d-mode transform', K));

                OutputFile = bst_process('GetNewFilename', StudyDir, 'results_eigentransform');
                bst_save(OutputFile, ResMat, 'v6');
                db_add_data(iStudyOut, OutputFile, ResMat);
                OutputFiles{end+1} = file_short(OutputFile); %#ok<AGROW>
            end
        end
    end
end
