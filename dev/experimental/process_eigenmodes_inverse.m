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
    sProcess.OutputTypes = {'results', 'results'};
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
    sProcess.options.whiten.Comment = 'Apply noise whitening (recommended)';
    sProcess.options.whiten.Type    = 'checkbox';
    sProcess.options.whiten.Value   = 1;
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
    % Map process options to the shared Compute() OPTIONS
    OPTIONS = struct();
    OPTIONS.InverseMethod    = 'eigenmode';
    OPTIONS.InverseMeasure   = lower(sProcess.options.method.Value);     % mne|dspm|sloreta
    OPTIONS.EigenmodePrior   = lower(sProcess.options.prior.Value);      % log|flat|power
    OPTIONS.SnrFixed         = sProcess.options.snr.Value{1};
    OPTIONS.SnrMethod        = 'fixed';
    OPTIONS.nModes           = sProcess.options.nmodes.Value{1};
    OPTIONS.NoiseMethod      = 'reg';
    OPTIONS.NoiseReg         = 0.1;
    OPTIONS.SourceOrient     = {'fixed'};
    OPTIONS.ComputeKernel    = 1;
    OPTIONS.DisplayMessages  = 0;
    OPTIONS.DataTypes        = [];     % auto-detect MEG else EEG in Compute
    OPTIONS.EigenmodeWhiten  = sProcess.options.whiten.Value;            % default 1
    OPTIONS.SaveCoefficients = ismember(lower(sProcess.options.outputtype.Value), {'coefficients','both'});
    % Resolve study/data indices for the shared core, then delegate
    % (kernel-only, one per file: replicate process_inverse_2018.Run case 2)
    iStudies = [sInputs.iStudy];
    iDatas   = [sInputs.iItem];
    [OutputFiles, errMessage] = process_inverse_2018('Compute', iStudies, iDatas, OPTIONS);
    if ~isempty(errMessage)
        bst_report('Error', sProcess, sInputs, errMessage);
    end
end
