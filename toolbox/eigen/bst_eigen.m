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
% NOTE: The 'spectrum' method is wired end-to-end for Laplace-Beltrami (scalar [nV x nT]),
%       Dirac (3D-vector map [3nV x nT] embedded as a pure-imaginary quaternion [4nV]),
%       face Dirac/Hodge (3-vector face map [3nF x nT]), and Connection Laplacian (a COMPLEX
%       tangent field [nV x nT] in the operator frame -- real-3D embedding needs the not-yet-
%       persisted per-vertex frame). The 'project'/'filter' methods remain stubs. LH and RH
%       are ALWAYS computed SEPARATELY (the basis is B-orthonormal per hemisphere).

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
try
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
                break;
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
        break;
    end

    % ===== COMPUTE TRANSFORM (dispatch) =====
    % Mirrors bst_timefreq's method switch. Engines are the bst_psd/morlet tier of this module.
    Result = [];
    switch lower(OPTIONS.Method)
        case 'spectrum'
            % Windowed eigenspectrum (the bst_psd analogue), computed SEPARATELY per hemisphere.
            [Result, Messages, isError] = ComputeEigenspectrum(F, EigenMat, OperatorMat, sfreq, OPTIONS);
            if isError, break; end
        case 'project'
            % TODO: Coef = Phi' * B * F  (manifold Fourier transform onto the eigenbasis).
        case 'filter'
            % TODO: eigen-domain spectral filter (project -> apply h(Lambda) -> reconstruct).
        otherwise
            Messages = ['Unknown eigen method: ' OPTIONS.Method];
            isError  = 1;
            break;
    end

    % ===== POST-PROCESS -- PLACEHOLDER =====
    Result = PostprocessEigen(Result, OPTIONS, EigenMat);

    % ===== SAVE FILE =====
    SaveFile(iTargetStudy, InitFile, DataType, Result, OPTIONS, EigenMat, TimeVector, SurfaceFile);
    bst_progress('inc', 1);
end
catch ME
    bst_progress('stop');
    rethrow(ME);
end
bst_progress('stop');


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
        if isempty(EigenMat.Phi{h})
            continue;
        end
        % Map the source map F into THIS hemisphere's native row layout (variant-specific).
        [U_h, Phi_h, B_h, Lam_h, msg] = ExtractHemiField(F, EigenMat, OperatorMat, h);
        if ~isempty(msg)
            Messages = ['bst_eigen: ' msg];
            isError  = 1;
            return;
        end
        % Same engine for all variants (C = Phi'*(B*U) is layout-agnostic; ' is conjugate
        % transpose, so the complex connection Hermitian product is handled automatically).
        [S_h, k_h, Nwin, msg2, Sstd_h] = bst_eigenspectrum(U_h, Phi_h, B_h, Lam_h, ...
            WinSamples, OPTIONS.WinOverlap, OPTIONS.WinFunc, OPTIONS.Measure);
        if isempty(S_h)
            Messages = ['bst_eigen: ' msg2];
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


%% ===== EXTRACT THE HEMISPHERE FIELD IN THE BASIS'S NATIVE LAYOUT =====
function [U_h, Phi_h, B_h, Lam_h, msg] = ExtractHemiField(F, EigenMat, OperatorMat, h)
    % Returns this hemisphere's field U_h in the row layout that matches Phi_h / B_h, plus
    % the basis and eigenvalues. msg is '' on success, else an error string. The transform
    % C = Phi_h'*(B_h*U_h) is identical across variants; only U_h's construction differs.
    msg   = '';
    U_h   = [];
    Phi_h = EigenMat.Phi{h};
    Lam_h = EigenMat.Lambda{h}(:);
    B_h   = OperatorMat.Mass{h};
    nT    = size(F, 2);
    switch EigenMat.Variant
        case 'Laplace-Beltrami'
            % Scalar field: one row per vertex. F is [nVertices x nTime].
            idx = EigenMat.GlobalVertices{h}(:);
            if max(idx) > size(F, 1)
                msg = sprintf(['hemisphere %d indexes vertex %d but the source map has %d rows ' ...
                    '(expected a scalar [nVertices x nTime] map for Laplace-Beltrami).'], h, max(idx), size(F,1));
                return;
            end
            U_h = F(idx, :);

        case 'Dirac'
            % 3D vector source map [3*nVertices x nTime], embedded as a pure-imaginary
            % quaternion field [4*nVh x nTime]: per vertex block rows = [w=0; x; y; z]
            % (matching bst_dirac.m). Mass B_h = kron(Mass_vertex, I4).
            vH = EigenMat.GlobalVertices{h}(:);
            if 3*max(vH) > size(F, 1)
                msg = sprintf(['Dirac needs an unconstrained 3-vector source map [3*nVertices x nTime]; ' ...
                    'hemisphere %d needs row %d but the map has %d rows.'], h, 3*max(vH), size(F,1));
                return;
            end
            nVh = numel(vH);
            U_h = zeros(4*nVh, nT);
            U_h(2:4:end, :) = F((vH-1)*3 + 1, :);   % x  (quaternion i)
            U_h(3:4:end, :) = F((vH-1)*3 + 2, :);   % y  (quaternion j)
            U_h(4:4:end, :) = F((vH-1)*3 + 3, :);   % z  (quaternion k); w rows stay 0

        case 'Connection Laplacian'
            % Complex tangent field, one (complex) row per vertex, expressed in the operator's
            % intrinsic frame. The real-3D-tangent -> complex embedding needs that per-vertex
            % frame {e1,e2}, which is NOT persisted in the operator/eigen node yet, so we
            % require the caller to pass the field already complex (in the operator frame).
            idx = EigenMat.GlobalVertices{h}(:);
            if max(idx) > size(F, 1)
                msg = sprintf('hemisphere %d indexes vertex %d but the field has %d rows.', h, max(idx), size(F,1));
                return;
            end
            if isreal(F)
                msg = ['Connection Laplacian needs a COMPLEX tangent field [nVertices x nTime] in the ' ...
                       'operator frame; embedding a real 3-vector field requires the per-vertex tangent ' ...
                       'frame, which is not persisted in the operator node yet (TODO: store the frame).'];
                return;
            end
            U_h = F(idx, :);

        case {'Dirac-Face', 'Hodge-Face'}
            % Face-domain quaternion field: 3-vector face map [3*nFaces x nTime] embedded as
            % pure-imaginary quaternion [4*nFh x nTime], indexed by GlobalFaces.
            fH = EigenMat.GlobalFaces{h}(:);
            if isempty(fH)
                msg = sprintf('variant ''%s'' hemisphere %d has no GlobalFaces.', EigenMat.Variant, h);
                return;
            end
            if 3*max(fH) > size(F, 1)
                msg = sprintf(['%s needs a 3-vector face map [3*nFaces x nTime]; hemisphere %d needs row %d ' ...
                    'but the map has %d rows.'], EigenMat.Variant, h, 3*max(fH), size(F,1));
                return;
            end
            nFh = numel(fH);
            U_h = zeros(4*nFh, nT);
            U_h(2:4:end, :) = F((fH-1)*3 + 1, :);
            U_h(3:4:end, :) = F((fH-1)*3 + 2, :);
            U_h(4:4:end, :) = F((fH-1)*3 + 3, :);

        otherwise
            msg = sprintf('unsupported eigen variant ''%s''.', EigenMat.Variant);
    end
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
