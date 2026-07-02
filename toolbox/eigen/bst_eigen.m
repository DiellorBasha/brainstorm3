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

% ===== VERB DISPATCH (string first arg) =====
if ischar(Data) && strcmpi(Data, 'Axes')          % ax = bst_eigen('Axes', OPTIONS)
    OutputFiles = BuildAxes(OPTIONS);  Messages = '';  isError = 0;
    return;
end
if (nargin >= 1) && ischar(Data) && strcmpi(Data, 'FieldSpec')
    OutputFiles = i_field_spec(OPTIONS);   % OPTIONS = ax
    return;
end

% ===== DEFAULT OPTIONS =====
Def_OPTIONS.Comment       = '';
Def_OPTIONS.Method        = 'spectrum';  % {'spectrum','filter','wavelet'} wired; 'project' stub
Def_OPTIONS.EigenFile     = [];          % eigen_ node (the spatial axis); [] => resolve from SurfaceFile
Def_OPTIONS.Variant       = [];          % operator family hint ('Laplace-Beltrami'|'Connection Laplacian'|'Dirac'|...)
Def_OPTIONS.nModes        = [];          % use a subset of modes (the spatial "band"); [] => all
Def_OPTIONS.Tau           = [];          % Dirac-type smoothing for basis resolution; [] => any
Def_OPTIONS.Measure       = 'power';     % spectrum measure {'power','magnitude'}
Def_OPTIONS.KernelName    = 'flat';      % 'filter' method: eigfilter kernel name (or handle / [K x 1] gain)
Def_OPTIONS.KernelParams  = struct();    % 'filter' method: kernel parameters
Def_OPTIONS.Nf            = 6;           % 'wavelet' method: number of frame members (per family)
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
        case 'filter'
            % Spatial eigen filter: project -> scale by h(Lambda) -> reconstruct, per hemisphere.
            [Ffilt, Messages, isError] = bst_eigenfilter('Analysis', F, EigenMat, OperatorMat, OPTIONS.KernelName, OPTIONS.KernelParams);
            if isError, break; end
            Result = struct('Type', 'filter', 'Field', Ffilt);
        case 'wavelet'
            % Spatial graph-wavelet transform: design a frame, analyze -> scalogram.
            % The frame must span the GLOBAL spectrum (max over hemispheres), else a hemi
            % whose lambda exceeds the design lmax would have those modes zeroed.
            lams = EigenMat.Lambda(~cellfun(@isempty, EigenMat.Lambda));
            if isempty(lams)
                Messages = 'bst_eigen: eigen node has no eigenvalues for the wavelet frame.';
                isError = 1; break;
            end
            lrange = [min(cellfun(@min, lams)), max(cellfun(@max, lams))];
            frame = bst_eigenwavelet('Design', OPTIONS.KernelName, OPTIONS.Nf, lrange);
            if strcmp(EigenMat.Variant, 'Connection Laplacian')
                % The connection operator lives on the TANGENT bundle: a complex number per
                % vertex in the canonical frame. Convert the physical 3-vector source to that
                % complex field (z = <v,e1> + i<v,e2>) via the operator's own frame, analyze,
                % then decode the complex scalogram back to ambient 3-D tangent vectors for
                % the saved scalogram. The normal component of the input is dropped (it is not
                % representable on the tangent bundle).
                [E1, E2] = GetConnectionFrame(EigenMat, OperatorMat);
                Fc = TangentToComplex(F, E1, E2);
                [Wc, Messages, isError] = bst_eigenwavelet('Analysis', Fc, EigenMat, OperatorMat, frame);
                if isError, break; end
                W = ComplexToTangent(Wc, E1, E2);          % [3nV x nT x M] ambient 3-vector
            else
                [W, Messages, isError] = bst_eigenwavelet('Analysis', F, EigenMat, OperatorMat, frame);
                if isError, break; end
                % Dirac family: Analysis carries the FULL quaternion [4nV x nT x M]; project to
                % the physical 3-vector for the saved source scalogram (the cortical viewer is
                % 1/3-component). The full-quaternion transform stays available via bst_eigenwavelet.
                if any(strcmp(EigenMat.Variant, {'Dirac', 'Dirac-Face', 'Hodge-Face'}))
                    W = bst_eigenwavelet('ToVec', W, EigenMat);
                end
            end
            Result = struct('Type', 'wavelet', 'W', W, 'Frame', frame);
        case 'project'
            % TODO: Coef = Phi' * B * F  (manifold Fourier transform onto the eigenbasis).
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
        % Build the output struct: an eigen-spectrum reuses TimefreqMat (Freqs = sqrt(Lambda),
        % the PSD analogue); a filtered source map is written as a results_ file (ResultsMat).
        switch lower(Result.Type)
            case 'spectrum'
                FileMat    = BuildSpectrumTimefreq(Result, OPTIONS, EigenMat, DataType, DataFile, TimeVector, SurfaceFile);
                filePrefix = 'timefreq_eigenspectrum';
            case 'filter'
                FileMat    = BuildFilterResult(Result, OPTIONS, EigenMat, DataFile, TimeVector, SurfaceFile);
                filePrefix = 'results_eigenfilter';
            case 'wavelet'
                FileMat    = BuildWaveletTimefreq(Result, OPTIONS, EigenMat, DataType, DataFile, TimeVector, SurfaceFile);
                filePrefix = 'timefreq_eigenwavelet';
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
        % Resolve the eigen_ node implicitly from the surface + operator family (the spatial
        % analogue of bst_timefreq deriving its axis from the recordings). Finding is bst_get's
        % job; a single surface can host several eigenbases, so Variant selects among them.
        Variant = OPTIONS.Variant;
        if isempty(Variant); Variant = 'Laplace-Beltrami'; end
        [sSubject, ~, iSurface, iEigen] = bst_get('EigenFileForSurface', SurfaceFile, Variant, OPTIONS.nModes, OPTIONS.Tau);
        if isempty(iEigen)
            error('bst_eigen:NoEigenForSurface', ...
                'No ''%s'' eigenbasis found for this surface — compute it first.', Variant);
        end
        EigenFile = sSubject.Surface(iSurface).Eigen(iEigen).FileName;
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


%% ===== BUILD THE EIGEN-FILTER RESULTS FILE =====
function FileMat = BuildFilterResult(Result, OPTIONS, EigenMat, DataFile, TimeVector, SurfaceFile)
    % Filtered source map -> results_ file. Unconstrained (3-vector, nComponents=3) for the
    % Dirac/face variants; scalar (nComponents=1) otherwise (Connection stays complex).
    if any(strcmp(EigenMat.Variant, {'Dirac', 'Dirac-Face', 'Hodge-Face'}))
        nComp = 3;
    else
        nComp = 1;
    end
    nT = size(Result.Field, 2);
    if ischar(OPTIONS.KernelName)
        kstr = OPTIONS.KernelName;
    else
        kstr = 'custom';
    end

    FileMat = db_template('resultsmat');
    FileMat.ImagingKernel = [];
    FileMat.ImageGridAmp  = Result.Field;          % [nSource x nTime] filtered source map
    FileMat.nComponents   = nComp;
    FileMat.HeadModelType = 'surface';
    FileMat.SurfaceFile   = SurfaceFile;
    FileMat.Function      = 'eigenfilter';
    if ~isempty(TimeVector) && (numel(TimeVector) == nT)
        FileMat.Time = TimeVector;
    else
        FileMat.Time = 0:(nT-1);
    end
    if ~isempty(DataFile)
        FileMat.DataFile = file_short(DataFile);
    end
    if isempty(OPTIONS.Comment)
        FileMat.Comment = sprintf('Eigenfilter (%s, %s)', EigenMat.Variant, kstr);
    else
        FileMat.Comment = OPTIONS.Comment;
    end
    FileMat.Options.Method       = 'eigenfilter';
    FileMat.Options.Variant      = EigenMat.Variant;
    FileMat.Options.EigenFile    = OPTIONS.EigenFile;
    FileMat.Options.KernelName   = kstr;
    FileMat.Options.KernelParams = OPTIONS.KernelParams;
    FileMat = bst_history('add', FileMat, 'compute', 'Eigen-domain spatial filter (bst_eigen)');
end


%% ===== BUILD THE EIGEN-WAVELET (SCALOGRAM) TIMEFREQMAT =====
function FileMat = BuildWaveletTimefreq(Result, OPTIONS, EigenMat, DataType, DataFile, TimeVector, SurfaceFile)
    % Spatial graph-wavelet scalogram W [nSrc x nTime x M] -> source TimefreqMat, with the
    % per-member centroid sqrt(lambda) as the (spatial-)frequency axis (reuses the cortical
    % time-freq viewer). 3-vector for Dirac/face + Connection Laplacian (decoded to ambient
    % tangent vectors); scalar otherwise.
    W = Result.W;                          % [nSrc x nTime x M]
    [~, nT, M] = size(W);
    if any(strcmp(EigenMat.Variant, {'Dirac', 'Dirac-Face', 'Hodge-Face', 'Connection Laplacian'}))
        nComp = 3;
    else
        nComp = 1;
    end
    FileMat = db_template('timefreqmat');
    FileMat.TF          = W;
    FileMat.Freqs       = Result.Frame.Centers(:)';   % per-member centroid sqrt(lambda)
    if ~isempty(TimeVector) && (numel(TimeVector) == nT)
        FileMat.Time = TimeVector;
    else
        FileMat.Time = 0:(nT-1);
    end
    FileMat.Measure     = 'other';
    FileMat.Method      = 'eigenwavelet';
    FileMat.DataType    = DataType;
    FileMat.SurfaceFile = SurfaceFile;
    FileMat.nComponents = nComp;
    FileMat.RowNames    = [];                          % source rows
    FileMat.nAvg = 1; FileMat.Leff = 1;
    if ~isempty(DataFile)
        FileMat.DataFile = file_short(DataFile);
    end
    if isempty(OPTIONS.Comment)
        FileMat.Comment = sprintf('Eigenwavelet (%s, %s, %d scales)', EigenMat.Variant, Result.Frame.Family, M);
    else
        FileMat.Comment = OPTIONS.Comment;
    end
    FileMat.Options.Method    = 'eigenwavelet';
    FileMat.Options.Variant   = EigenMat.Variant;
    FileMat.Options.EigenFile = OPTIONS.EigenFile;
    FileMat.Options.Family    = Result.Frame.Family;
    FileMat.Options.Nf        = Result.Frame.Nf;
    FileMat.Options.Centers   = Result.Frame.Centers;
    FileMat = bst_history('add', FileMat, 'compute', 'Spatial graph-wavelet transform (bst_eigen)');
end


%% ===== CONNECTION-LAPLACIAN TANGENT <-> COMPLEX BRIDGE =====
% The connection operator acts on a complex tangent field (one z per vertex in the
% canonical frame). These convert a physical 3-vector source map to/from that field using
% the operator's OWN frame (bst_operator_frame -> provably the eigenmode decoding basis).

function [E1, E2] = GetConnectionFrame(EigenMat, OperatorMat)
    % Global per-vertex tangent axes [nVg x 3], scattered from each hemisphere's frame.
    nVg = 0;
    for h = 1:numel(EigenMat.GlobalVertices)
        if ~isempty(EigenMat.GlobalVertices{h})
            nVg = max(nVg, max(EigenMat.GlobalVertices{h}(:)));
        end
    end
    E1 = zeros(nVg, 3);  E2 = zeros(nVg, 3);
    for h = 1:numel(EigenMat.GlobalVertices)
        gv = EigenMat.GlobalVertices{h}(:);
        if isempty(gv); continue; end
        Fr = bst_operator_frame(OperatorMat, h);
        E1(gv, :) = Fr.e1;  E2(gv, :) = Fr.e2;
    end
end

function Z = TangentToComplex(F3, E1, E2)
    % Physical 3-vector field [3nVg x nT] -> complex tangent field [nVg x nT]:
    % z(v) = <u(v), e1(v)> + i <u(v), e2(v)>  (the normal component of u is dropped).
    nVg = size(E1, 1);  nT = size(F3, 2);
    F3r = reshape(F3, 3, nVg, nT);
    re = sum(F3r .* reshape(E1.', 3, nVg), 1);   % [1 x nVg x nT]
    im = sum(F3r .* reshape(E2.', 3, nVg), 1);
    Z  = reshape(re, nVg, nT) + 1i * reshape(im, nVg, nT);
end

function V = ComplexToTangent(Wc, E1, E2)
    % Complex tangent scalogram [nVg x nT x M] -> ambient 3-vector [3nVg x nT x M]:
    % u(v) = real(z(v)) * e1(v) + imag(z(v)) * e2(v).
    [nVg, nT, M] = size(Wc);
    V = zeros(3 * nVg, nT, M);
    for m = 1:M
        Z = Wc(:, :, m);                          % [nVg x nT]
        Vk = zeros(3, nVg, nT);
        for k = 1:3
            Vk(k, :, :) = real(Z) .* E1(:, k) + imag(Z) .* E2(:, k);   % [nVg x nT]
        end
        V(:, :, m) = reshape(Vk, 3 * nVg, nT);
    end
end


%% ===== AXES: canonical eigenbasis x time x temporal-frequency =====
function [Time, tlag, omega, nT, NFFT, Fs] = BuildTimeAxis(OPTIONS)
    % Time = actual time vector (s); tlag = lags from 0 (for ts kernels); omega = temporal-freq grid (Hz).
    if isfield(OPTIONS,'Time') && ~isempty(OPTIONS.Time)
        Time = OPTIONS.Time(:)';  Fs = 1 / (Time(2) - Time(1));
    else
        Fs = OPTIONS.SampleRate;  tw = OPTIONS.TimeWindow;
        Time = tw(1) : 1/Fs : tw(2);
    end
    nT = numel(Time);  NFFT = nT;
    omega = (0:NFFT-1) * (Fs / NFFT);            % temporal-frequency grid (Hz)
    tlag  = (0:nT-1) / Fs;                        % time-lag axis (s)
end

function ax = BuildAxes(OPTIONS)
    % ax = bst_eigen('Axes', OPTIONS): assemble the JOINT (cortex, time, temporal-frequency) axes a JTV
    % atom is realised on. The cortex axis is an eigenbasis of OPTIONS.SurfaceFile (find-or-create via
    % tess_eigen); the time axis is OPTIONS.Time, or OPTIONS.TimeWindow + OPTIONS.SampleRate; omega is the
    % temporal FFT grid. Single source-of-truth shared by view_atom_designer and bst_dynamics('Axes').
    if ~isfield(OPTIONS,'SurfaceFile') || isempty(OPTIONS.SurfaceFile)
        error('bst_eigen(''Axes''): OPTIONS.SurfaceFile is required.');
    end
    Variant = 'Laplace-Beltrami'; if isfield(OPTIONS,'Variant') && ~isempty(OPTIONS.Variant), Variant = OPTIONS.Variant; end
    nModes  = 200;                if isfield(OPTIONS,'nModes')  && ~isempty(OPTIONS.nModes),  nModes  = OPTIONS.nModes;  end
    % --- spatial axis: eigenbasis (find-or-create) ---
    EigenMat = tess_eigen(OPTIONS.SurfaceFile, Variant, 'nModes', nModes);
    Op       = in_bst_operator(EigenMat.OperatorFile);
    % --- temporal + spectral axis ---
    [Time, tlag, omega, nT, NFFT, Fs] = BuildTimeAxis(OPTIONS);
    % --- assemble (cell fields after struct() to avoid struct-array widening) ---
    ax = struct('Variant',EigenMat.Variant, 'Time',Time, 'Fs',Fs, 'nT',nT, 'NFFT',NFFT, ...
                'omega',omega, 'tlag',tlag, 'SurfaceFile',OPTIONS.SurfaceFile);
    ax.EigenMat = EigenMat;  ax.Operator = Op;
    ax.Phi = EigenMat.Phi;   ax.Lambda = EigenMat.Lambda;  ax.Mass = Op.Mass;
    ax.GlobalVertices = EigenMat.GlobalVertices;
end

%% ===== FIELD SPEC: (field_type, domain) -> layout, from operator metadata =====
function spec = i_field_spec(ax)
    ft = '';  dom = '';
    if isfield(ax,'Operator') && isstruct(ax.Operator) && isfield(ax.Operator,'Registry') ...
            && ~isempty(ax.Operator.Registry) && isfield(ax.Operator.Registry,'Primary') ...
            && ~isempty(ax.Operator.Registry.Primary)
        P = ax.Operator.Registry.Primary;
        if isfield(P,'field_type'), ft  = P.field_type; end
        if isfield(P,'domain'),     dom = P.domain;     end
    end
    switch lower(ft)
        case 'real',       C = 1;
        case 'complex',    C = 2;
        case 'quaternion', C = 4;
        otherwise,         C = [];      % pre-registry -> infer below
    end
    if isempty(C)                        % guarded fallback: derive from the Phi row layout
        nV = numel(ax.GlobalVertices{1});
        C  = round(size(ax.Phi{1},1) / max(nV,1));
        switch C, case 1, ft='real'; case 2, ft='complex'; case 4, ft='quaternion'; otherwise, ft='real'; C=1; end
    end
    if isempty(dom), dom = 'vertex'; end               % completeness default (vertex is universal)
    switch C
        case 1, kind='scalar';     width=1; nComp=1;
        case 2, kind='tangent';    width=1; nComp=1;   % complex row per element
        case 4, kind='quaternion'; width=4; nComp=3;
        otherwise, kind='scalar';  width=1; nComp=1;
    end
    spec = struct('field_type',ft, 'domain',dom, 'width',width, 'nComponents',nComp, 'C',C, 'kind',kind);
end
