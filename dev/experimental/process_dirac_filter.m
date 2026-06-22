function varargout = process_dirac_filter( varargin )
% PROCESS_DIRAC_FILTER: Spatial spectral filter of an UNCONSTRAINED source vector
% field in the Dirac eigenmode domain (the vector analogue of
% process_eigenmodes_filter).
%
% Projects the 3-vector source field onto the curvature-aware Dirac eigenbasis,
% applies a scale kernel g(lambda) on the Dirac eigenvalues (shared eigfilter
% library), and reconstructs the filtered vector field. Operates on unconstrained
% (3 components/vertex) surface source maps.
%
% SEE ALSO: bst_dirac_eigenmodes_filter, bst_eigenmodes_filter, tess_eigen, bst_dirac

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
    sProcess.Comment     = 'Spatial spectral filter (Dirac eigenmodes)';
    sProcess.FileTag     = 'diracfilter';
    sProcess.Category    = 'Filter';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 338;
    sProcess.Description = '';
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Process all vertices at once (eigenmode projection needs the full vertex dim)
    sProcess.processDim  = [];

    % === FILTER TYPE ===
    sProcess.options.filtertype.Comment = {'Low-pass (smooth)', 'High-pass (detail)', ...
                                           'Band-pass', 'Heat kernel (diffusion)', ...
                                           'Filter type:'; ...
                                           'lowpass', 'highpass', 'bandpass', 'heat', ''};
    sProcess.options.filtertype.Type    = 'radio_linelabel';
    sProcess.options.filtertype.Value   = 'lowpass';
    sProcess.options.sep1.Type = 'separator';
    % === CUTOFF MODE (lowpass/highpass) ===
    sProcess.options.cutoffmode.Comment = 'Cutoff mode index (for low/high-pass): ';
    sProcess.options.cutoffmode.Type    = 'value';
    sProcess.options.cutoffmode.Value   = {100, '', 0};
    % === MODE RANGE (bandpass) ===
    sProcess.options.moderange_low.Comment = 'Band-pass: lower mode index: ';
    sProcess.options.moderange_low.Type    = 'value';
    sProcess.options.moderange_low.Value   = {40, '', 0};
    sProcess.options.moderange_high.Comment = 'Band-pass: upper mode index: ';
    sProcess.options.moderange_high.Type    = 'value';
    sProcess.options.moderange_high.Value   = {400, '', 0};
    % === DIFFUSION TIME (heat; in Dirac-eigenvalue units) ===
    sProcess.options.diffusiontime.Comment = 'Diffusion time (heat; Dirac-&lambda; units): ';
    sProcess.options.diffusiontime.Type    = 'value';
    sProcess.options.diffusiontime.Value   = {2e5, '', 0};
    % === INFO ===
    sProcess.options.label_info.Comment = ['<FONT color="#777777">Requires an UNCONSTRAINED (3 comp/vertex) surface<BR>' ...
                                           'source map; the Dirac eigenbasis is found-or-created on its surface.</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    switch lower(sProcess.options.filtertype.Value)
        case 'lowpass',  Comment = sprintf('| dirac lowpass(k<=%d)',  sProcess.options.cutoffmode.Value{1});
        case 'highpass', Comment = sprintf('| dirac highpass(k>=%d)', sProcess.options.cutoffmode.Value{1});
        case 'bandpass', Comment = sprintf('| dirac bandpass[%d,%d]', sProcess.options.moderange_low.Value{1}, sProcess.options.moderange_high.Value{1});
        case 'heat',     Comment = sprintf('| dirac heat(t=%g)', sProcess.options.diffusiontime.Value{1});
        otherwise,       Comment = '| dirac filter';
    end
end


%% ===== RUN =====
function sInput = Run(sProcess, sInput) %#ok<DEFNU>
    FilterType    = lower(sProcess.options.filtertype.Value);
    CutoffMode    = sProcess.options.cutoffmode.Value{1};
    ModeRangeLow  = sProcess.options.moderange_low.Value{1};
    ModeRangeHigh = sProcess.options.moderange_high.Value{1};
    DiffusionTime = sProcess.options.diffusiontime.Value{1};

    % ----- source metadata -----
    FileMat = in_bst_results(sInput.FileName, 0, 'SurfaceFile', 'Atlas', 'nComponents', 'HeadModelType');
    if isfield(FileMat,'HeadModelType') && ~isempty(FileMat.HeadModelType) && ~strcmpi(FileMat.HeadModelType,'surface')
        bst_report('Error', sProcess, sInput, 'Dirac filtering is only supported for surface source models.'); sInput=[]; return;
    end
    if isfield(FileMat,'Atlas') && ~isempty(FileMat.Atlas)
        bst_report('Error', sProcess, sInput, 'Dirac filtering is not supported for atlas-based source models.'); sInput=[]; return;
    end
    if isempty(FileMat.nComponents) || (FileMat.nComponents ~= 3)
        bst_report('Error', sProcess, sInput, ...
            'Dirac filtering requires UNCONSTRAINED source orientation (3 components/vertex).'); sInput=[]; return;
    end
    SurfaceFile = FileMat.SurfaceFile;
    if isempty(SurfaceFile)
        bst_report('Error', sProcess, sInput, 'No surface file associated with this source map.'); sInput=[]; return;
    end

    % ----- Dirac eigenbasis + operator mass (found-or-created) -----
    bst_progress('text', 'Loading Dirac eigenbasis...');
    EigenMat = tess_eigen(SurfaceFile, 'Dirac');
    OpMat    = in_bst_operator(EigenMat.OperatorFile);
    MassCell = OpMat.Mass;
    nVert = double(max(cellfun(@(x) max(x(:)), EigenMat.GlobalVertices)));
    if size(sInput.A,1) ~= 3*nVert
        bst_report('Error', sProcess, sInput, sprintf( ...
            'Vertex mismatch: source data is %d rows but the Dirac basis spans 3*%d=%d.', ...
            size(sInput.A,1), nVert, 3*nVert)); sInput=[]; return;
    end

    % ----- filter arguments -----
    switch FilterType
        case 'lowpass',  filterArgs = {'CutoffMode', CutoffMode}; infoStr = sprintf('low-pass (k<=%d)', CutoffMode);
        case 'highpass', filterArgs = {'CutoffMode', CutoffMode}; infoStr = sprintf('high-pass (k>=%d)', CutoffMode);
        case 'bandpass', filterArgs = {'ModeRange', [ModeRangeLow, ModeRangeHigh]}; infoStr = sprintf('band-pass [%d, %d]', ModeRangeLow, ModeRangeHigh);
        case 'heat'
            if DiffusionTime <= 0
                bst_report('Error', sProcess, sInput, 'Diffusion time must be positive.'); sInput=[]; return;
            end
            filterArgs = {'DiffusionTime', DiffusionTime}; infoStr = sprintf('heat (t=%g)', DiffusionTime);
        otherwise
            bst_report('Error', sProcess, sInput, ['Unknown filter type: ' FilterType]); sInput=[]; return;
    end

    % ----- apply (scale filter; real vector in -> real vector out) -----
    bst_progress('text', 'Dirac eigenmode filtering...');
    sInput.A = real(bst_dirac_eigenmodes_filter(EigenMat, MassCell, sInput.A, FilterType, filterArgs{:}));

    sInput.CommentTag = FormatComment(sProcess);
    sInput.HistoryComment = sprintf('Dirac eigenmode spatial filter: %s', infoStr);
end
