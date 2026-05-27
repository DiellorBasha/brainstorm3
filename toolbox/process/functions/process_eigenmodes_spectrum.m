function varargout = process_eigenmodes_spectrum( varargin )
% PROCESS_EIGENMODES_SPECTRUM: Spectral decomposition of source data onto eigenmodes.
%
% USAGE:  sProcess = process_eigenmodes_spectrum('GetDescription')
%       OutputFiles = process_eigenmodes_spectrum('Run', sProcess, sInputs)
%
% DESCRIPTION:
%     Projects source-mapped brain activity onto the Laplace-Beltrami
%     eigenmode basis and saves the spectral coefficients as a Brainstorm
%     matrix file. This reveals the spatial frequency content of cortical
%     activity over time:
%
%         c_k(t) = phi_k' * M * u(t)
%
%     The output matrix has rows = modes and columns = time, viewable as
%     a mode-vs-time image (analogous to a spectrogram). Additional outputs:
%
%     - Power spectrum:  |c_k(t)|^2  — spatial power per mode per time
%     - Mean power:      mean(|c_k|^2, t) — time-averaged spatial spectrum
%     - Cumulative energy: sum(|c_k|^2, k=1..K) / sum(|c_k|^2) — fraction
%       of total energy captured by the first K modes
%
%     Requires precomputed eigenmodes (from process_eigenmodes).
%
% SEE ALSO: bst_eigenmodes_project, in_tess_eigenmodes, process_eigenmodes

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
    sProcess.Comment     = 'Eigenmode spectral decomposition';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 338;
    sProcess.Description = '';
    % Input / output
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % === OUTPUT TYPE ===
    sProcess.options.outputtype.Comment = {'Spectral coefficients c_k(t)', ...
                                           'Power spectrum |c_k(t)|^2', ...
                                           'Time-averaged power spectrum', ...
                                           'Cumulative energy ratio', ...
                                           'Output:'; ...
                                           'coefficients', 'power', 'meanpower', 'cumulative', ''};
    sProcess.options.outputtype.Type    = 'radio_linelabel';
    sProcess.options.outputtype.Value   = 'coefficients';

    % === INFO ===
    sProcess.options.label_info.Comment = ['<FONT color="#777777">Output is a matrix file: rows = eigenmode index, columns = time.<BR>' ...
                                           'Requires precomputed eigenmodes on the source surface.</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    OutputType = sProcess.options.outputtype.Value;
    switch lower(OutputType)
        case 'coefficients'
            Comment = 'Eigenmode coefficients c_k(t)';
        case 'power'
            Comment = 'Eigenmode power |c_k(t)|^2';
        case 'meanpower'
            Comment = 'Eigenmode mean power spectrum';
        case 'cumulative'
            Comment = 'Eigenmode cumulative energy';
        otherwise
            Comment = 'Eigenmode spectrum';
    end
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};

    % Get option
    OutputType = lower(sProcess.options.outputtype.Value);

    % Process each input file
    for iInput = 1:length(sInputs)
        sInput = sInputs(iInput);

        % ===== LOAD SOURCE DATA =====
        FileMat = in_bst_results(sInput.FileName, 1, ...
            'ImageGridAmp', 'ImagingKernel', 'SurfaceFile', ...
            'nComponents', 'HeadModelType', 'Atlas', 'Time', 'DataFile');

        % Validate
        if isfield(FileMat, 'HeadModelType') && ~isempty(FileMat.HeadModelType) && ~strcmpi(FileMat.HeadModelType, 'surface')
            bst_report('Error', sProcess, sInput, 'Only surface source models are supported.');
            continue;
        end
        if isfield(FileMat, 'Atlas') && ~isempty(FileMat.Atlas)
            bst_report('Error', sProcess, sInput, 'Atlas-based source models are not supported.');
            continue;
        end
        if FileMat.nComponents ~= 1
            bst_report('Error', sProcess, sInput, 'Constrained source orientation required.');
            continue;
        end

        % Get full source matrix
        if ~isempty(FileMat.ImagingKernel)
            % Kernel-based: expand
            if isempty(FileMat.DataFile)
                bst_report('Error', sProcess, sInput, 'Kernel result requires an associated data file.');
                continue;
            end
            DataMat = in_bst_data(FileMat.DataFile);
            Sources = FileMat.ImagingKernel * DataMat.F(FileMat.GoodChannel, :);
            TimeVector = DataMat.Time;
        else
            Sources = FileMat.ImageGridAmp;
            TimeVector = FileMat.Time;
        end

        SurfaceFile = FileMat.SurfaceFile;
        if isempty(SurfaceFile)
            bst_report('Error', sProcess, sInput, 'No surface file associated with source map.');
            continue;
        end

        % ===== CHECK EIGENMODES =====
        [Eigenmodes, isComputed] = in_tess_eigenmodes(SurfaceFile);
        if ~isComputed
            bst_report('Error', sProcess, sInput, ...
                ['No eigenmodes on: ' SurfaceFile '. Run "Compute eigenmodes" first.']);
            continue;
        end

        nV_eigen = size(Eigenmodes.Vectors, 1);
        nV_data  = size(Sources, 1);
        if nV_data ~= nV_eigen
            bst_report('Error', sProcess, sInput, ...
                sprintf('Vertex mismatch: sources %d, eigenmodes %d. Recompute sources on updated surface.', ...
                nV_data, nV_eigen));
            continue;
        end

        % Mass matrix for the M-weighted projection (computed once per surface).
        sSurf = in_tess_bst(SurfaceFile, 0);
        [~, MassMatrix] = tess_laplacian(sSurf.Vertices, sSurf.Faces, 'MassType', Eigenmodes.MassType);

        % ===== PROJECT ONTO EIGENMODES =====
        bst_progress('start', 'Eigenmode spectrum', 'Computing spectral projection...', 0, 100);
        bst_progress('set', 10);

        Coeffs = bst_eigenmodes_project(Eigenmodes, double(Sources), MassMatrix);
        nModes = size(Coeffs, 1);
        nTime  = size(Coeffs, 2);
        bst_progress('set', 70);

        % ===== COMPUTE OUTPUT =====
        switch OutputType
            case 'coefficients'
                OutputData = Coeffs;
                OutputTime = TimeVector;
                Units = '';
                CommentStr = sprintf('Eigenmode c_k(t) [%d modes]', nModes);

            case 'power'
                OutputData = Coeffs.^2;
                OutputTime = TimeVector;
                Units = '';
                CommentStr = sprintf('Eigenmode power [%d modes]', nModes);

            case 'meanpower'
                OutputData = mean(Coeffs.^2, 2);   % [nModes x 1]
                OutputTime = 0;                     % single time point
                Units = '';
                CommentStr = sprintf('Eigenmode mean power [%d modes]', nModes);

            case 'cumulative'
                MeanPower = mean(Coeffs.^2, 2);
                TotalPower = sum(MeanPower);
                if TotalPower > 0
                    OutputData = cumsum(MeanPower) / TotalPower;
                else
                    OutputData = zeros(nModes, 1);
                end
                OutputTime = 0;
                Units = '';
                CommentStr = sprintf('Eigenmode cumulative energy [%d modes]', nModes);
        end

        % ===== BUILD MATRIX FILE =====
        MatrixMat = db_template('matrixmat');
        MatrixMat.Value       = OutputData;
        MatrixMat.Time        = OutputTime;
        MatrixMat.Comment     = CommentStr;
        MatrixMat.SurfaceFile = SurfaceFile;
        MatrixMat.nAvg        = 1;
        MatrixMat.Leff        = 1;
        MatrixMat.DisplayUnits = Units;

        % Row descriptions (mode labels with eigenvalue)
        lambdas = Eigenmodes.Values;
        RowNames = cell(nModes, 1);
        for k = 1:nModes
            RowNames{k} = sprintf('Mode %d (lam=%.1f)', k, lambdas(k));
        end
        MatrixMat.Description = RowNames;

        % History
        MatrixMat = bst_history('add', MatrixMat, 'eigenmodes_spectrum', ...
            sprintf('Spectral decomposition (%s) of: %s', OutputType, sInput.FileName));
        MatrixMat = bst_history('add', MatrixMat, 'eigenmodes_spectrum', ...
            sprintf('%d modes, %s mass, eigenvalue range [%.2f, %.2f]', ...
            nModes, Eigenmodes.MassType, lambdas(1), lambdas(end)));

        % ===== SAVE =====
        [sStudy, iStudy] = bst_get('Study', sInput.iStudy);
        StudyDir = bst_fileparts(file_fullpath(sStudy.FileName));
        OutputFile = bst_process('GetNewFilename', StudyDir, 'matrix_eigenspectrum');
        bst_save(OutputFile, MatrixMat, 'v6');
        db_add_data(iStudy, OutputFile, MatrixMat);

        bst_progress('set', 100);
        bst_progress('stop');

        OutputFiles{end+1} = file_short(OutputFile);
    end
end
