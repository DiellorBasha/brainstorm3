function varargout = process_eigenmodes_coeffsfilter( varargin )
% PROCESS_EIGENMODES_COEFFSFILTER: Spatial-spectral filter of eigenmode coefficients.
%
% USAGE:  sProcess = process_eigenmodes_coeffsfilter('GetDescription')
%       OutputFiles = process_eigenmodes_coeffsfilter('Run', sProcess, sInputs)
%
% DESCRIPTION:
%     Applies a per-mode spatial-spectral gain h(lambda_k) (see
%     bst_eigenmodes_filter_gain) directly to an eigenmode-coefficient matrix
%     (matrix_eigentransform, [K x nTime]): theta_filt = h .* theta. This is the
%     "filter" step of the transform-first workflow. Optionally reconstructs the
%     filtered vertex sources Q = Phi * theta_filt.
%
% SEE ALSO: bst_eigenmodes_filter_gain, process_eigenmodes_transform

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
    sProcess.Comment     = 'Eigenmode coefficient filter';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 336.9;   % after the dispersion analysis (336.8)
    sProcess.Description = '';
    sProcess.InputTypes  = {'matrix'};
    sProcess.OutputTypes = {'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    sProcess.options.filtertype.Comment = {'Low-pass', 'High-pass', 'Band-pass', 'Heat (smooth)', 'Inverse heat (sharpen)', 'Tikhonov', 'Filter:'; ...
                                           'lowpass', 'highpass', 'bandpass', 'heat', 'inverse_heat', 'tikhonov', ''};
    sProcess.options.filtertype.Type    = 'radio_linelabel';
    sProcess.options.filtertype.Value   = 'heat';

    sProcess.options.cutoffmode.Comment = 'Cutoff mode index (low/high-pass): ';
    sProcess.options.cutoffmode.Type    = 'value';
    sProcess.options.cutoffmode.Value   = {50, '', 0};

    sProcess.options.moderange_low.Comment  = 'Band-pass: lower mode index: ';
    sProcess.options.moderange_low.Type     = 'value';
    sProcess.options.moderange_low.Value    = {20, '', 0};

    sProcess.options.moderange_high.Comment = 'Band-pass: upper mode index: ';
    sProcess.options.moderange_high.Type    = 'value';
    sProcess.options.moderange_high.Value   = {80, '', 0};

    sProcess.options.diffusiontime.Comment  = 'Diffusion time (heat / inverse heat): ';
    sProcess.options.diffusiontime.Type     = 'value';
    sProcess.options.diffusiontime.Value    = {0.005, '', 4};

    sProcess.options.regbeta.Comment        = 'Tikhonov beta (h = 1/(1+beta*lambda)): ';
    sProcess.options.regbeta.Type           = 'value';
    sProcess.options.regbeta.Value          = {1, '', 4};

    sProcess.options.dorecon.Comment        = 'Also reconstruct filtered vertex sources (Q = Phi * theta_filt)';
    sProcess.options.dorecon.Type           = 'checkbox';
    sProcess.options.dorecon.Value          = 0;

    sProcess.options.label_info.Comment = ['<FONT color="#777777">Applies a per-mode gain h(&lambda;) to the ' ...
        'eigenmode coefficients (the "filter" step of transform-first).<BR>Input: an eigenmode-coefficient matrix ' ...
        '(matrix_eigentransform); requires eigenmodes on the surface.</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = ['Eigenmode coefficient filter (' sProcess.options.filtertype.Value ')'];
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    FilterType = lower(sProcess.options.filtertype.Value);
    CutoffMode = sProcess.options.cutoffmode.Value{1};
    ModeLow    = sProcess.options.moderange_low.Value{1};
    ModeHigh   = sProcess.options.moderange_high.Value{1};
    DiffTime   = sProcess.options.diffusiontime.Value{1};
    RegBeta    = sProcess.options.regbeta.Value{1};
    DoRecon    = sProcess.options.dorecon.Value;

    for iInput = 1:numel(sInputs)
        sInput = sInputs(iInput);
        M = in_bst_matrix(sInput.FileName);
        Coeffs = double(M.Value);                       % [K x nTime]
        K = size(Coeffs, 1);

        if ~isfield(M, 'SurfaceFile') || isempty(M.SurfaceFile)
            bst_report('Error', sProcess, sInput, 'Matrix has no SurfaceFile; cannot obtain eigenmodes.');
            continue;
        end
        [Em, isComputed] = in_tess_eigenmodes(M.SurfaceFile);
        if ~isComputed
            bst_report('Error', sProcess, sInput, ['No eigenmodes on surface: ' M.SurfaceFile]);
            continue;
        end
        if size(Em.Values, 1) < K
            bst_report('Error', sProcess, sInput, sprintf('Coefficients have %d modes but surface has only %d eigenvalues.', K, size(Em.Values,1)));
            continue;
        end
        lambdas = double(Em.Values(1:K));

        h = bst_eigenmodes_filter_gain(lambdas, FilterType, ...
            'CutoffMode', CutoffMode, 'ModeRange', [ModeLow, ModeHigh], ...
            'DiffusionTime', DiffTime, 'RegBeta', RegBeta);
        CoeffsFilt = bsxfun(@times, h, Coeffs);         % [K x nTime]

        [sStudyOut, iStudyOut] = bst_get('Study', sInput.iStudy);
        StudyDir = bst_fileparts(file_fullpath(sStudyOut.FileName));

        % --- Filtered coefficients (matrix) ---
        Mout = db_template('matrixmat');
        Mout.Value        = CoeffsFilt;
        Mout.Time         = M.Time;
        Mout.Description  = M.Description;
        Mout.SurfaceFile  = M.SurfaceFile;
        Mout.Comment      = sprintf('EigenFilt [%s] | %s', FilterType, sInput.Comment);
        if isfield(M, 'nAvg') && ~isempty(M.nAvg), Mout.nAvg = M.nAvg; else, Mout.nAvg = 1; end
        if isfield(M, 'Leff') && ~isempty(M.Leff), Mout.Leff = M.Leff; end
        Mout = bst_history('add', Mout, 'eigenmodes_coeffsfilter', Mout.Comment);
        OutFile = bst_process('GetNewFilename', StudyDir, 'matrix_eigenfilt');
        bst_save(OutFile, Mout, 'v6');
        db_add_data(iStudyOut, OutFile, Mout);
        OutputFiles{end+1} = file_short(OutFile); %#ok<AGROW>

        % --- Optional filtered vertex reconstruction (results) ---
        if DoRecon
            Phi = double(Em.Vectors(:, 1:K));
            Q   = Phi * CoeffsFilt;                     % [nVert x nTime]
            ResMat = db_template('resultsmat');
            ResMat.ImageGridAmp  = Q;
            ResMat.ImagingKernel = [];
            ResMat.nComponents   = 1;
            ResMat.Time          = M.Time;
            ResMat.SurfaceFile   = M.SurfaceFile;
            ResMat.HeadModelType = 'surface';
            ResMat.ColormapType  = 'source';
            ResMat.Comment       = sprintf('EigenFilt recon [%s] | %s', FilterType, sInput.Comment);
            if isfield(M, 'nAvg') && ~isempty(M.nAvg), ResMat.nAvg = M.nAvg; else, ResMat.nAvg = 1; end
            if isfield(M, 'Leff') && ~isempty(M.Leff), ResMat.Leff = M.Leff; end
            ResMat = bst_history('add', ResMat, 'eigenmodes_coeffsfilter', ResMat.Comment);
            OutFile2 = bst_process('GetNewFilename', StudyDir, 'results_eigenfilt');
            bst_save(OutFile2, ResMat, 'v6');
            db_add_data(iStudyOut, OutFile2, ResMat);
            % Recon is a DB sidecar (visible in the tree); not added to OutputFiles
            % because OutputTypes is {'matrix'} (mixing types breaks pipeline chaining).
        end

        bst_report('Info', sProcess, sInput, sprintf('Filtered %d eigenmode coefficients with %s.', K, FilterType));
    end
end
