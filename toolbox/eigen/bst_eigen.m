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
%             (eigen coefficients vs Lambda, the spatial analogue of a PSD).
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
% NOTE: SKELETON. The READ and WRITE tiers are wired to the file-based system; the
%       DISPATCH (compute engines) and POST-PROCESS tiers are PLACEHOLDERS, to be filled
%       in incrementally (delegating to the eigfilter library and the eigen analysis
%       engines, the bst_psd/morlet tier of this module).

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
Def_OPTIONS.Method        = 'project';   % {'project','filter','spectrum'} (dispatch placeholders)
Def_OPTIONS.EigenFile     = [];          % eigen_ node (the spatial axis); [] => resolve from SurfaceFile
Def_OPTIONS.Variant       = [];          % operator family hint ('Laplace-Beltrami'|'Connection Laplacian'|'Dirac'|...)
Def_OPTIONS.nModes        = [];          % use a subset of modes (the spatial "band"); [] => all
Def_OPTIONS.Measure       = 'power';     % spectrum measure {'none','power','magnitude'}
Def_OPTIONS.OutputType    = 'result';    % {'result','spectrum'} : analyzed source map vs eigen-spectrum
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

    % ===== READ EIGEN BASIS (the spatial-spectral axis) =====
    % Analogue of bst_timefreq deriving sfreq/TimeVector: here we obtain Lambda (the
    % spatial-frequency axis) and Phi (the spatial modes) from the eigen_ node. The
    % EigenMat.Variant field tells the dispatch which operator family (and therefore which
    % data layout) we are analyzing.
    EigenMat = GetEigenBasis(OPTIONS, SurfaceFile);

    % ===== COMPUTE TRANSFORM (dispatch) -- PLACEHOLDER =====
    % Mirrors bst_timefreq's method switch. Each case will delegate to a domain compute
    % engine (the bst_psd/morlet tier of this module) and MUST honour EigenMat.Variant
    % (scalar vs tangent-vector vs 3D-vector data layout). Stubs for now; F holds the
    % source-mapped data and EigenMat the basis.
    % Guard: the dispatch engines consume the source-mapped data F in the EigenMat basis.
    if isempty(F)
        Messages = 'bst_eigen: empty source-mapped data; nothing to analyze.';
        isError  = 1;
        return;
    end
    Result = [];   % engine output (coefficients / filtered map / spectrum)
    switch lower(OPTIONS.Method)
        case 'project'
            % TODO: Coef = Phi' * B * F  (manifold Fourier transform onto the eigenbasis).
            %       Delegate to the project engine; respect Variant data layout.
        case 'filter'
            % TODO: eigen-domain spectral filter (project -> apply h(Lambda) -> reconstruct).
            %       Delegate to the eigfilter library (bst_eigfilter_*).
        case 'spectrum'
            % TODO: eigen power spectrum vs Lambda (the spatial analogue of a PSD).
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
        % TODO: measure (power/magnitude), mode-band selection (OPTIONS.nModes), unconstrained
        % orientation handling for Dirac results, normalization. Pass-through for now.
    end


%% ===== SAVE FILE (placeholder skeleton) =====
    function SaveFile(iTargetStudy, DataFile, DataType, Result, OPTIONS, EigenMat, TimeVector, SurfaceFile)
        % Build the output struct. Two output kinds, mirroring how bst_timefreq always
        % writes a timefreq_*.mat: here we write either a result_ (analyzed source map) or
        % an eigen-spectrum file (eigen coefficients vs Lambda, the PSD analogue).
        switch lower(OPTIONS.OutputType)
            case 'result'
                FileMat      = struct();
                FileMat.Type = 'result';        % placeholder
                filePrefix   = 'results_eigen';
            case 'spectrum'
                FileMat      = struct();
                FileMat.Type = 'eigenspectrum'; % placeholder
                filePrefix   = 'timefreq_eigen';
            otherwise
                error('bst_eigen:OutputType', 'Unknown OutputType: %s', OPTIONS.OutputType);
        end
        % Common metadata (placeholder)
        FileMat.Comment      = OPTIONS.Comment;
        FileMat.DataType     = DataType;
        FileMat.Method       = OPTIONS.Method;
        FileMat.Time         = TimeVector;
        FileMat.SurfaceFile  = SurfaceFile;
        FileMat.EigenFile    = OPTIONS.EigenFile;
        FileMat.Variant      = EigenMat.Variant;
        FileMat.Result       = Result;            % placeholder payload
        FileMat.FilePrefix   = filePrefix;        % intended output filename prefix (real save uses it)
        FileMat.iTargetStudy = iTargetStudy;      % intended output study (real save uses it)
        if ~isempty(DataFile)
            FileMat.DataFile = file_short(DataFile);
        end
        % TODO: when iTargetStudy is set and the templates/payloads are defined, write to disk:
        %   sTargetStudy = bst_get('Study', iTargetStudy);
        %   FileName = bst_process('GetNewFilename', bst_fileparts(sTargetStudy.FileName), filePrefix);
        %   bst_save(FileName, FileMat, 'v6');  db_add_data(iTargetStudy, FileName, FileMat);
        %   OutputFiles{end+1} = FileName;  return;
        OutputFiles{end+1} = FileMat;   % skeleton: return contents (no disk write yet)
    end

end


%% ===== GET EIGEN BASIS (resolve + load the spatial axis) =====
function EigenMat = GetEigenBasis(OPTIONS, SurfaceFile)
    % Resolve the eigen_ node: explicit OPTIONS.EigenFile, else find the one attached to
    % this data's surface. Finding is bst_get's job; loading is in_bst_eigen's.
    EigenFile = OPTIONS.EigenFile;
    if isempty(EigenFile) && ~isempty(SurfaceFile)
        % TODO (auto-resolve): EigenFile = bst_get('EigenFileForSurface', SurfaceFile, OPTIONS.Variant);
        error('bst_eigen:AutoResolveTODO', ...
            'Auto-resolving the eigen_ basis from the surface is not implemented yet; pass OPTIONS.EigenFile.');
    end
    if isempty(EigenFile)
        error('bst_eigen:NoEigenFile', 'No eigen_ basis specified (OPTIONS.EigenFile).');
    end
    % Load the eigenbasis (the spatial axis: Lambda + Phi + Variant).
    EigenMat = in_bst_eigen(EigenFile);
    % NOTE: the operator node (EigenMat.OperatorFile) carries the mass B; the compute engines
    % load it on demand via in_bst_operator(EigenMat.OperatorFile) when M-weighting is needed.
end
