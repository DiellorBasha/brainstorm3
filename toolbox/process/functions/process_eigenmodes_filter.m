function varargout = process_eigenmodes_filter( varargin )
% PROCESS_EIGENMODES_FILTER: Spatial frequency filtering of sources via Laplace-Beltrami eigenmodes.
%
% USAGE:  sProcess = process_eigenmodes_filter('GetDescription')
%           sInput = process_eigenmodes_filter('Run', sProcess, sInput)
%
% DESCRIPTION:
%     Filters source-mapped brain activity in the spatial frequency domain
%     defined by the Laplace-Beltrami eigenmodes of the cortical surface.
%     This is the surface analogue of Fourier-domain filtering: the data is
%     projected onto the eigenmode basis, multiplied by a transfer function,
%     and reconstructed.
%
%     Unlike Gaussian smoothing (process_ssmooth), eigenmode filtering:
%       - Respects the intrinsic geometry of the cortical surface
%       - Provides precise frequency cutoffs in the spectral domain
%       - Allows band-pass isolation of specific spatial scales
%       - Implements heat kernel smoothing (equivalent to surface diffusion)
%
%     Requires precomputed eigenmodes on the surface (from process_eigenmodes).
%
% SEE ALSO: bst_eigenmodes_filter, bst_eigenmodes_project, tess_eigenmodes, process_eigenmodes

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
    % Description of the process
    sProcess.Comment     = 'Spatial spectral filter (eigenmodes)';
    sProcess.FileTag     = 'eigenfilter';
    sProcess.Category    = 'Filter';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 337;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'results', 'timefreq'};
    sProcess.OutputTypes = {'results', 'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Process all vertices at once (eigenmode projection needs full vertex dim)
    sProcess.processDim  = [];

    % === FILTER TYPE ===
    sProcess.options.filtertype.Comment = {'Low-pass (smooth)', 'High-pass (detail)', ...
                                           'Band-pass', 'Heat kernel (diffusion)', ...
                                           'Filter type:'; ...
                                           'lowpass', 'highpass', 'bandpass', 'heat', ''};
    sProcess.options.filtertype.Type    = 'radio_linelabel';
    sProcess.options.filtertype.Value   = 'lowpass';

    % === SEPARATOR ===
    sProcess.options.sep1.Type = 'separator';

    % === CUTOFF MODE (for lowpass/highpass) ===
    sProcess.options.cutoffmode.Comment = 'Cutoff mode index (for low/high-pass): ';
    sProcess.options.cutoffmode.Type    = 'value';
    sProcess.options.cutoffmode.Value   = {50, '', 0};

    % === MODE RANGE (for bandpass) ===
    sProcess.options.moderange_low.Comment = 'Band-pass: lower mode index: ';
    sProcess.options.moderange_low.Type    = 'value';
    sProcess.options.moderange_low.Value   = {20, '', 0};

    sProcess.options.moderange_high.Comment = 'Band-pass: upper mode index: ';
    sProcess.options.moderange_high.Type    = 'value';
    sProcess.options.moderange_high.Value   = {80, '', 0};

    % === DIFFUSION TIME (for heat kernel) ===
    sProcess.options.diffusiontime.Comment = 'Diffusion time (for heat kernel): ';
    sProcess.options.diffusiontime.Type    = 'value';
    sProcess.options.diffusiontime.Value   = {0.005, '', 4};

    % === INFO LABEL ===
    sProcess.options.label_info.Comment = ['<FONT color="#777777">Requires precomputed eigenmodes on the surface<BR>' ...
                                           '(use process "Compute eigenmodes" first).</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    FilterType = sProcess.options.filtertype.Value;
    switch lower(FilterType)
        case 'lowpass'
            Comment = sprintf('Eigenmode low-pass (k<=%d)', ...
                sProcess.options.cutoffmode.Value{1});
        case 'highpass'
            Comment = sprintf('Eigenmode high-pass (k>=%d)', ...
                sProcess.options.cutoffmode.Value{1});
        case 'bandpass'
            Comment = sprintf('Eigenmode band-pass [%d, %d]', ...
                sProcess.options.moderange_low.Value{1}, ...
                sProcess.options.moderange_high.Value{1});
        case 'heat'
            Comment = sprintf('Eigenmode heat kernel (t=%.4f)', ...
                sProcess.options.diffusiontime.Value{1});
        otherwise
            Comment = 'Eigenmode spatial filter';
    end
end


%% ===== RUN =====
function sInput = Run(sProcess, sInput) %#ok<DEFNU>
    % ===== GET OPTIONS =====
    FilterType    = lower(sProcess.options.filtertype.Value);
    CutoffMode    = sProcess.options.cutoffmode.Value{1};
    ModeRangeLow  = sProcess.options.moderange_low.Value{1};
    ModeRangeHigh = sProcess.options.moderange_high.Value{1};
    DiffusionTime = sProcess.options.diffusiontime.Value{1};

    % ===== LOAD SOURCE FILE METADATA =====
    switch (file_gettype(sInput.FileName))
        case {'results', 'link'}
            FileMat = in_bst_results(sInput.FileName, 0, 'SurfaceFile', 'GridLoc', 'Atlas', 'nComponents', 'HeadModelType');
            nComponents = FileMat.nComponents;
        case 'timefreq'
            FileMat = in_bst_timefreq(sInput.FileName, 0, 'SurfaceFile', 'GridLoc', 'Atlas', 'HeadModelType');
            nComponents = 1;
        otherwise
            bst_report('Error', sProcess, sInput, 'Unsupported file format.');
            sInput = [];
            return;
    end

    % Check: surface head model only
    if isfield(FileMat, 'HeadModelType') && ~isempty(FileMat.HeadModelType) && ~strcmpi(FileMat.HeadModelType, 'surface')
        bst_report('Error', sProcess, sInput, 'Eigenmode filtering is only supported for surface source models.');
        sInput = [];
        return;
    end
    % Check: not atlas-based
    if isfield(FileMat, 'Atlas') && ~isempty(FileMat.Atlas)
        bst_report('Error', sProcess, sInput, 'Eigenmode filtering is not supported for atlas-based source models.');
        sInput = [];
        return;
    end
    % Check: constrained orientation only
    if nComponents ~= 1
        bst_report('Error', sProcess, sInput, ...
            ['Eigenmode filtering requires constrained source orientation (1 component per vertex).' 10 ...
             'With unconstrained orientations, project to constrained first.']);
        sInput = [];
        return;
    end

    % ===== GET SURFACE FILE =====
    SurfaceFile = FileMat.SurfaceFile;
    if isempty(SurfaceFile)
        bst_report('Error', sProcess, sInput, 'No surface file associated with this source map.');
        sInput = [];
        return;
    end

    % ===== CHECK EIGENMODES =====
    [Eigenmodes, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if ~isComputed
        bst_report('Error', sProcess, sInput, ...
            ['No precomputed eigenmodes found on surface: ' SurfaceFile 10 ...
             'Run process "Compute eigenmodes" on this subject first.']);
        sInput = [];
        return;
    end

    % Check vertex count
    nV_eigen = size(Eigenmodes.Vectors, 1);
    nV_data  = size(sInput.A, 1);
    if nV_data ~= nV_eigen
        bst_report('Error', sProcess, sInput, ...
            sprintf(['Vertex count mismatch: source data has %d vertices but eigenmodes have %d.\n' ...
                     'This can happen if the mesh was repaired during eigenmode computation.\n' ...
                     'Recompute the head model and source estimates on the updated surface.'], ...
                     nV_data, nV_eigen));
        sInput = [];
        return;
    end

    % Mass matrix for the M-weighted filter (computed once, reused across freqs).
    sSurf = in_tess_bst(SurfaceFile, 0);
    [~, MassMatrix] = tess_laplacian(sSurf.Vertices, sSurf.Faces, 'MassType', Eigenmodes.MassType);

    % ===== APPLY FILTER =====
    % Build option arguments for bst_eigenmodes_filter
    switch FilterType
        case 'lowpass'
            filterArgs = {'CutoffMode', CutoffMode};
            infoStr = sprintf('low-pass (k<=%d)', CutoffMode);
        case 'highpass'
            filterArgs = {'CutoffMode', CutoffMode};
            infoStr = sprintf('high-pass (k>=%d)', CutoffMode);
        case 'bandpass'
            filterArgs = {'ModeRange', [ModeRangeLow, ModeRangeHigh]};
            infoStr = sprintf('band-pass [%d, %d]', ModeRangeLow, ModeRangeHigh);
        case 'heat'
            if DiffusionTime <= 0
                bst_report('Error', sProcess, sInput, 'Diffusion time must be positive.');
                sInput = [];
                return;
            end
            filterArgs = {'DiffusionTime', DiffusionTime};
            infoStr = sprintf('heat kernel (t=%.4f)', DiffusionTime);
        otherwise
            bst_report('Error', sProcess, sInput, ['Unknown filter type: ' FilterType]);
            sInput = [];
            return;
    end

    % Apply the filter — handles [nV x nTime] and [nV x nTime x nFreq]
    nFreqs = size(sInput.A, 3);
    if nFreqs > 1
        % Time-frequency data: filter each frequency band separately
        for iFreq = 1:nFreqs
            sInput.A(:, :, iFreq) = bst_eigenmodes_filter(Eigenmodes, ...
                sInput.A(:, :, iFreq), MassMatrix, FilterType, filterArgs{:});
        end
    else
        sInput.A = bst_eigenmodes_filter(Eigenmodes, sInput.A, MassMatrix, FilterType, filterArgs{:});
    end

    % ===== UPDATE HISTORY =====
    sInput.CommentTag = FormatComment(sProcess);
    sInput.HistoryComment = sprintf('Eigenmode spatial filter: %s (%d modes, %s mass)', ...
        infoStr, Eigenmodes.nModes, Eigenmodes.MassType);
end
