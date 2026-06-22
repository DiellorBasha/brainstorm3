function [OutputFiles, Messages, isError] = bst_eigen(Data, OPTIONS)
% BST_EIGEN: Spatial-spectral / differential-geometry analysis of surface-mapped data
%            in an operator eigenbasis. The spatial analogue of BST_TIMEFREQ.
%
% USAGE:  [OutputFiles, Messages, isError] = bst_eigen(Data, OPTIONS)
%                                   OPTIONS = bst_eigen();
%
% DESCRIPTION:
%     bst_eigen is the single ORCHESTRATOR for eigen-domain ANALYSIS, built to mirror
%     bst_timefreq end-to-end (read -> dispatch -> post-process -> write database):
%
%       bst_timefreq : recordings + TIME axis (sfreq, TimeVector)        -> temporal transform
%       bst_eigen    : source map  + SPATIAL-SPECTRAL axis (Lambda, Phi)  -> spatial transform
%
%     Where bst_timefreq reads recordings and obtains the time axis (from which the DFT is
%     computed), bst_eigen reads a surface-mapped data series and obtains the SPATIAL axis
%     from an eigen_ file: the operator eigenvalues Lambda (the spatial-frequency axis) and
%     eigenvectors Phi (the spatial modes). It is general over the operator family carried
%     in the eigen_ metadata -- Laplace-Beltrami (scalar), Connection Laplacian (tangent
%     vector), Dirac (3D embedded vector) -- dispatching on EigenMat.Variant.
%
%     Module boundary (see also the eigen-module-reorg memory):
%       - tess_eigen / tess_operators : structural PRE-COMPUTATION of the eigenbasis
%             (the "spatial axis"), analogous to importing/preprocessing recordings
%             (the "time axis"). Produces the eigen_ / operator_ / manifold_ nodes.
%       - bst_dirac / bst_inverse_dirac : SOURCE MAPPING that happens to consume eigen_/
%             operator_ files (a specific use of the eigen system). NOT analysis.
%       - bst_eigen (this file) : the ANALYSIS orchestrator over (eigen_ file + source-
%             mapped data [+ operator_/manifold_ if needed]). Owns read/dispatch/post/write.
%             Output is either a result_ file (analyzed source map) or an eigen-spectrum
%             (eigen coefficients vs Lambda, the spatial analogue of a PSD) stored as a
%             TimefreqMat with Freqs = sqrt(Lambda).
%
% INPUTS:
%     - Data    : one of
%          - String filename (a results_ / matrix_ source-mapped file)
%          - Cell-array of filenames
%          - Matrix of source time series [nSource x nTime]
%          - Cell-array of such matrices
%     - OPTIONS : struct (call bst_eigen() with no args to get the defaults), fields below.
%
% OUTPUTS:
%     - OutputFiles : cell-array of created files (or file contents when not saved)
%     - Messages    : string, errors/warnings
%     - isError     : 1 if an error happened
%
% NOTE: The 'spectrum' method is wired end-to-end for the Laplace-Beltrami (scalar) variant;
%       the 'project'/'filter' methods and the Dirac/Connection data layouts remain stubs.
%       LH and RH are ALWAYS computed SEPARATELY (the basis is B-orthonormal per hemisphere).

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

% ===== DEFAULT OPTIONS =====
Def_OPTIONS.Comment       = '';
Def_OPTIONS.Method        = 'spectrum';  % {'spectrum'(wired),'project','filter'(stubs)}
Def_OPTIONS.EigenFile     = [];          % eigen_ node (the spatial axis); [] => resolve from SurfaceFile
Def_OPTIONS.Variant       = [];          % operator family hint ('Laplace-Beltrami'|'Connection Laplacian'|'Dirac'|...)
Def_OPTIONS.nModes        = [];          % use a subset of modes (the spatial "band"); [] => all
Def_OPTIONS.Measure       = 'power';     % spectrum measure {'power','magnitude'}
Def_OPTIONS.WinLength     = [];          % windowing: window length in SECONDS ([] => one window)
Def_OPTIONS.WinOverlap    = 50;          % windowing: overlap in percent
Def_OPTIONS.WinFunc       = 'mean';      % windowing: aggregate across windows {'mean','std','mean+std'}
Def_OPTIONS.TimeWindow    = [];          % restrict the input time series
Def_OPTIONS.iTargetStudy  = [];          % output study ('NoSave' => return contents, do not save)

% Return the default options
if (nargin == 0)
    OutputFiles = Def_OPTIONS;
    return;
end
% Copy default options to OPTIONS structure (do not replace defined values)
OPTIONS = struct_copy_fields(OPTIONS, Def_OPTIONS, 0);


% ===== PARSE INPUTS =====
OutputFiles = {};
Messages    = [];
isError     = 0;
% Data: list of source-mapped blocks/files to analyze
if ~iscell(Data)
    Data = {Data};
end


% ===== LOOP ON FILES =====
bst_progress('start', 'Eigen analysis', 'Analyzing in the eigen basis...', 0, length(Data));
for iData = 1:length(Data)
    % ===== GET INITIAL DATA FILE =====
    isFile = ischar(Data{iData});
    if isFile
        InitFile = Data{iData};
        [~, iStudy, ~, DataType] = bst_get('AnyFile', InitFile);
        if strcmpi(DataType, 'link')
            DataType = 'results';
        end
    else
        InitFile = [];
        DataType = 'matrix';
        iStudy   = [];
    end
    % Output study
    if isequal(OPTIONS.iTargetStudy, 'NoSave')
        iTargetStudy = [];
    elseif ~isempty(OPTIONS.iTargetStudy)
        iTargetStudy = OPTIONS.iTargetStudy;
    else
        iTargetStudy = iStudy;
    end

    % ===== READ DATA (the surface-mapped source series) =====
    % Analogue of bst_timefreq's READ DATA: obtain F [nSource x nTime], the time vector,
    % and -- crucially -- the SurfaceFile that ties this data to its eigenbasis.
    SurfaceFile = [];
    TimeVector  = [];
    if isFile
        switch (DataType)
            case 'results'
                ResultsMat  = in_bst_results(InitFile, 1, 'ImageGridAmp', 'Time', 'SurfaceFile', 'nComponents');
                F           = ResultsMat.ImageGridAmp;     % [nSource x nTime] (full source map)
                TimeVector  = ResultsMat.Time;
                SurfaceFile = ResultsMat.SurfaceFile;
                clear ResultsMat;
            case 'matrix'
                sMat        = in_bst_matrix(InitFile);
                F           = sMat.Value;
                TimeVector  = sMat.Time;
                clear sMat;
            otherwise
                Messages = ['Unsupported data type for eigen analysis: ' DataType];
                isError  = 1;
                return;
        end
        % Restrict to the requested time window
        if ~isempty(OPTIONS.TimeWindow) && ~isempty(TimeVector)
            iTime      = bst_closest(OPTIONS.TimeWindow, TimeVector);
            iTime      = iTime(1):iTime(2);
            TimeVector = TimeVector(iTime);
            F          = F(:, iTime);
        end
    else
        % ===== PROCESS DATA BLOCKS (raw matrix in) =====
        F = Data{iData};
    end
    % Sampling frequency for seconds->samples windowing (single sample => 1 Hz placeholder)
    if ~isempty(TimeVector) && (numel(TimeVector) >= 2)
        sfreq = 1 ./ (TimeVector(2) - TimeVector(1));
    else
        sfreq = 1;
    end

    % ===== READ EIGEN BASIS (the spatial-spectral axis) =====
    % Analogue of bst_timefreq deriving sfreq/TimeVector: obtain Lambda (the spatial-
    % frequency axis), Phi (the spatial modes) and the per-hemisphere mass B. EigenMat.Variant
    % tells the dispatch which operator family (and therefore which data layout) we analyze.
    [EigenMat, OperatorMat] = GetEigenBasis(OPTIONS, SurfaceFile);

    % Guard: the dispatch engines consume the source-mapped data F in the EigenMat basis.
    if isempty(F)
        Messages = 'bst_eigen: empty source-mapped data; nothing to analyze.';
        isError  = 1;
        return;
    end

    % ===== COMPUTE TRANSFORM (dispatch) =====
    % Mirrors bst_timefreq's method switch. Engines are the bst_psd/morlet tier of this module.
    Result = [];
    switch lower(OPTIONS.Method)
        case 'spectrum'
            % Windowed eigenspectrum (the bst_psd analogue), computed SEPARATELY per hemisphere.
            [Result, Messages, isError] = ComputeEigenspectrum(F, EigenMat, OperatorMat, sfreq, OPTIONS);
            if isError, return; end
        case 'project'
            % TODO: Coef = Phi' * B * F  (manifold Fourier transform onto the eigenbasis).
        case 'filter'
            % TODO: eigen-domain spectral filter (project -> apply h(Lambda) -> reconstruct).
        otherwise
            Messages = ['Unknown eigen method: ' OPTIONS.Method];
            isError  = 1;
            return;
    end

    % ===== POST-PROCESS -- PLACEHOLDER =====
    Result = PostprocessEigen(Result, OPTIONS, EigenMat);

    % ===== SAVE FILE =====
    SaveFile(iTargetStudy, InitFile, DataType, Result, OPTIONS, EigenMat, TimeVector, SurfaceFile);
    bst_progress('inc', 1);
end


%% ===== POST-PROCESS (placeholder) =====
    function Result = PostprocessEigen(Result, OPTIONS, EigenMat) %#ok<INUSD>
        % TODO: mode-band selection (OPTIONS.nModes), normalization, relative spectra.
        % Pass-through for now.
    end


%% ===== SAVE FILE (assemble + write the output) =====
    function SaveFile(iTargetStudy, DataFile, DataType, Result, OPTIONS, EigenMat, TimeVector, SurfaceFile)
        if isempty(Result)
            return;
        end
        % Build the output struct. Currently only the eigen-spectrum output is assembled,
        % reusing TimefreqMat (Freqs = sqrt(Lambda), the PSD analogue). A 'result_' output
        % (reconstructed/filtered source map) is a TODO for the project/filter methods.
        switch lower(Result.Type)
            case 'spectrum'
                FileMat    = BuildSpectrumTimefreq(Result, OPTIONS, EigenMat, DataType, DataFile, TimeVector, SurfaceFile);
                filePrefix = 'timefreq_eigenspectrum';
            otherwise
                error('bst_eigen:OutputType', 'Unknown result type: %s', Result.Type);
        end
        % Write to the database, or return contents when no target study (mirrors bst_timefreq).
        if ~isempty(iTargetStudy)
            sTargetStudy = bst_get('Study', iTargetStudy);
            FileName     = bst_process('GetNewFilename', bst_fileparts(sTargetStudy.FileName), filePrefix);
            bst_save(FileName, FileMat, 'v6');
            db_add_data(iTargetStudy, FileName, FileMat);
            OutputFiles{end+1} = FileName;
        else
            OutputFiles{end+1} = FileMat;
        end
    end

end


%% ===== GET EIGEN BASIS (resolve + load the spatial axis and the mass) =====
function [EigenMat, OperatorMat] = GetEigenBasis(OPTIONS, SurfaceFile)
    % Resolve the eigen_ node: explicit OPTIONS.EigenFile, else find the one on this surface.
    % Finding is bst_get's job; loading is in_bst_eigen's / in_bst_operator's.
    EigenFile = OPTIONS.EigenFile;
    if isempty(EigenFile) && ~isempty(SurfaceFile)
        % TODO: EigenFile = bst_get('EigenFileForSurface', SurfaceFile, OPTIONS.Variant);
        error('bst_eigen:AutoResolveTODO', ...
            'Auto-resolving the eigen_ basis from the surface is not implemented yet; pass OPTIONS.EigenFile.');
    end
    if isempty(EigenFile)
        error('bst_eigen:NoEigenFile', 'No eigen_ basis specified (OPTIONS.EigenFile).');
    end
    % Load the eigenbasis (the spatial axis: Lambda + Phi + Variant) and the operator (mass B).
    EigenMat = in_bst_eigen(EigenFile);
    if ~isfield(EigenMat, 'OperatorFile') || isempty(EigenMat.OperatorFile)
        error('bst_eigen:NoOperator', 'Eigen node has no linked OperatorFile (mass B unavailable).');
    end
    OperatorMat = in_bst_operator(EigenMat.OperatorFile);   % carries Mass(1x2), the per-hemisphere B
end


%% ===== COMPUTE EIGENSPECTRUM (per-hemisphere, separate) =====
function [Result, Messages, isError] = ComputeEigenspectrum(F, EigenMat, OperatorMat, sfreq, OPTIONS)
    Result   = [];
    Messages = '';
    isError  = 0;
    % First slice: Laplace-Beltrami (scalar) layout only -- F is [nVertices x nTime] and a
    % hemisphere's field is F(GlobalVertices{h}, :). Dirac/Connection vector layouts (3nV /
    % complex-tangent) need a layout map and are a TODO.
    if ~strcmpi(EigenMat.Variant, 'Laplace-Beltrami')
        Messages = sprintf(['bst_eigen: ''spectrum'' is wired for the Laplace-Beltrami (scalar) variant only; ' ...
                            'got ''%s'' (Dirac/Connection vector layouts are a TODO).'], EigenMat.Variant);
        isError  = 1;
        return;
    end
    % Window length: seconds -> samples (the engine windows in samples).
    if ~isempty(OPTIONS.WinLength)
        WinSamples = round(OPTIONS.WinLength * sfreq);
    else
        WinSamples = [];
    end
    nHemi = numel(EigenMat.Phi);
    Hemi  = struct('Tag', {}, 'Spectrum', {}, 'Std', {}, 'kAxis', {}, 'Lambda', {});
    HemiTags = {'LH', 'RH'};
    Nwin = [];
    for h = 1:nHemi
        Phi_h = EigenMat.Phi{h};
        if isempty(Phi_h)
            continue;
        end
        idxH = EigenMat.GlobalVertices{h}(:);
        % Layout check: scalar field rows must index into F.
        if max(idxH) > size(F, 1)
            Messages = sprintf(['bst_eigen: hemisphere %d indexes vertex %d but the source map has %d rows. ' ...
                                'The data layout does not match a scalar Laplace-Beltrami basis ' ...
                                '(is this an unconstrained/vector source map?).'], h, max(idxH), size(F,1));
            isError  = 1;
            return;
        end
        U_h   = F(idxH, :);
        B_h   = OperatorMat.Mass{h};            % per-hemisphere mass (B-orthonormal with Phi_h)
        Lam_h = EigenMat.Lambda{h}(:);
        [S_h, k_h, Nwin, msg, Sstd_h] = bst_eigenspectrum(U_h, Phi_h, B_h, Lam_h, ...
            WinSamples, OPTIONS.WinOverlap, OPTIONS.WinFunc, OPTIONS.Measure);
        if isempty(S_h)
            Messages = ['bst_eigen: ' msg];
            isError  = 1;
            return;
        end
        Hemi(end+1) = struct('Tag', HemiTags{min(h,2)}, 'Spectrum', S_h, 'Std', Sstd_h, ...
                             'kAxis', k_h, 'Lambda', Lam_h);   %#ok<AGROW>
    end
    if isempty(Hemi)
        Messages = 'bst_eigen: no hemispheres to analyze.';
        isError  = 1;
        return;
    end
    Result.Type    = 'spectrum';
    Result.Measure = OPTIONS.Measure;
    Result.Nwin    = Nwin;
    Result.Hemi    = Hemi;
end


%% ===== BUILD THE EIGEN-SPECTRUM TIMEFREQMAT =====
function FileMat = BuildSpectrumTimefreq(Result, OPTIONS, EigenMat, DataType, DataFile, TimeVector, SurfaceFile)
    % LH and RH are stored as SEPARATE ROWS (always computed separately). Because the two
    % hemispheres have distinct -- but rank-aligned and near-identical -- eigenvalues, the
    % shared TimefreqMat.Freqs axis is the per-rank MEAN of sqrt(Lambda); the EXACT per-
    % hemisphere axes/eigenvalues are preserved in Options for an interleaved/exact view.
    nH = numel(Result.Hemi);
    K  = min(arrayfun(@(s) numel(s.Spectrum), Result.Hemi));   % common mode count (truncate if unequal)
    TF       = zeros(nH, 1, K);
    Std      = [];
    RowNames = cell(nH, 1);
    kAll     = zeros(nH, K);
    LambdaAll = cell(1, nH);
    kAxisAll  = cell(1, nH);
    hasStd = ~isempty(Result.Hemi(1).Std);
    if hasStd, Std = zeros(nH, 1, K); end
    for h = 1:nH
        s = Result.Hemi(h);
        TF(h, 1, :) = reshape(s.Spectrum(1:K), 1, 1, K);
        if hasStd && ~isempty(s.Std)
            Std(h, 1, :) = reshape(s.Std(1:K), 1, 1, K); %#ok<AGROW>
        end
        RowNames{h}  = ['eigen ' s.Tag];
        kAll(h, :)   = s.kAxis(1:K)';
        kAxisAll{h}  = s.kAxis;
        LambdaAll{h} = s.Lambda;
    end
    Freqs = mean(kAll, 1)';   % nominal shared spatial-frequency axis = sqrt(Lambda), per-rank mean

    FileMat = db_template('timefreqmat');
    FileMat.TF        = TF;
    FileMat.Std       = Std;
    % Static measure: store the window span [first last] (placeholder when no time axis)
    if isempty(TimeVector)
        FileMat.Time = [0, 1];
    elseif numel(TimeVector) >= 2
        FileMat.Time = TimeVector([1, end]);
    else
        FileMat.Time = [TimeVector(1), TimeVector(1)];
    end
    FileMat.Freqs     = Freqs;
    FileMat.RowNames  = RowNames;
    FileMat.Measure   = Result.Measure;
    FileMat.Method    = 'eigenspectrum';
    FileMat.DataType  = DataType;
    FileMat.SurfaceFile = SurfaceFile;
    FileMat.nAvg      = 1;
    FileMat.Leff      = 1;
    if ~isempty(DataFile)
        FileMat.DataFile = file_short(DataFile);
    end
    if isempty(OPTIONS.Comment)
        FileMat.Comment = sprintf('Eigenspectrum (%s, %d modes)', Result.Measure, K);
    else
        FileMat.Comment = OPTIONS.Comment;
    end
    % Provenance + exact per-hemisphere axes (the shared Freqs is only a nominal display axis)
    FileMat.Options.Method      = 'eigenspectrum';
    FileMat.Options.Measure     = Result.Measure;
    FileMat.Options.Variant     = EigenMat.Variant;
    FileMat.Options.EigenFile   = OPTIONS.EigenFile;
    FileMat.Options.Nwin        = Result.Nwin;
    FileMat.Options.kAxis       = kAxisAll;     % exact sqrt(Lambda) per hemisphere
    FileMat.Options.Lambda      = LambdaAll;    % exact eigenvalues per hemisphere
    FileMat.Options.RowHemi     = {Result.Hemi.Tag};
    FileMat = bst_history('add', FileMat, 'compute', 'Eigen-spectrum decomposition (bst_eigen)');
end
