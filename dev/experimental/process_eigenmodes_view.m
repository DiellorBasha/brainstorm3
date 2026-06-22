function varargout = process_eigenmodes_view( varargin )
% PROCESS_EIGENMODES_VIEW: Export eigenmodes as a source map for browsing in the 3D viewer.
%
% USAGE:  sProcess = process_eigenmodes_view('GetDescription')
%       OutputFiles = process_eigenmodes_view('Run', sProcess, sInputs)
%
% DESCRIPTION:
%     Creates a Brainstorm source results file where each "time point"
%     corresponds to one Laplace-Beltrami eigenmode. The standard 3D viewer
%     and time slider let users browse through eigenmodes and see their
%     spatial patterns on the cortical surface.
%
%     The Time axis is set to the eigenmode indices (1, 2, ..., nModes),
%     so the time slider position shows which mode is currently displayed.
%     The eigenvalues (spatial frequencies squared) are stored in the file
%     history and comment for reference.
%
%     Requires precomputed eigenmodes (from process_eigenmodes).
%
% SEE ALSO: in_tess_eigenmodes, process_eigenmodes

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
    sProcess.Comment     = 'View eigenmodes on surface';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'Import', 'Import anatomy'};
    sProcess.Index       = 26;
    sProcess.Description = '';
    % Input / output
    sProcess.InputTypes  = {'import'};
    sProcess.OutputTypes = {'import'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 0;

    % === SUBJECT NAME ===
    sProcess.options.subjectname.Comment = 'Subject name:';
    sProcess.options.subjectname.Type    = 'subjectname';
    sProcess.options.subjectname.Value   = '';

    % === SURFACE TYPE ===
    sProcess.options.surfacetype.Comment = {'Cortex', 'InnerSkull', 'OuterSkull', 'Scalp', 'Other', 'Surface type:'; ...
                                            'cortex', 'innerskull', 'outerskull', 'scalp', 'other', ''};
    sProcess.options.surfacetype.Type    = 'radio_linelabel';
    sProcess.options.surfacetype.Value   = 'cortex';

    % === SURFACE INDEX ===
    sProcess.options.isurfidx.Comment = 'Surface index (for Other type): ';
    sProcess.options.isurfidx.Type    = 'value';
    sProcess.options.isurfidx.Value   = {1, '', 0};

    % === SEPARATOR ===
    sProcess.options.sep1.Type = 'separator';

    % === MODE RANGE ===
    sProcess.options.moderange_low.Comment = 'First mode to export: ';
    sProcess.options.moderange_low.Type    = 'value';
    sProcess.options.moderange_low.Value   = {1, '', 0};

    sProcess.options.moderange_high.Comment = 'Last mode to export: ';
    sProcess.options.moderange_high.Type    = 'value';
    sProcess.options.moderange_high.Value   = {200, '', 0};

    % === NORMALIZE ===
    sProcess.options.normalize.Comment = 'Normalize each mode to unit max (better colormap contrast)';
    sProcess.options.normalize.Type    = 'checkbox';
    sProcess.options.normalize.Value   = 0;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sprintf('View eigenmodes [%d-%d]', ...
        sProcess.options.moderange_low.Value{1}, ...
        sProcess.options.moderange_high.Value{1});
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};

    % ===== GET OPTIONS =====
    SubjectName   = file_standardize(sProcess.options.subjectname.Value);
    SurfType      = sProcess.options.surfacetype.Value;
    ModeRangeLow  = sProcess.options.moderange_low.Value{1};
    ModeRangeHigh = sProcess.options.moderange_high.Value{1};
    DoNormalize   = sProcess.options.normalize.Value;

    if isempty(SubjectName)
        bst_report('Error', sProcess, [], 'Subject name is empty.');
        return;
    end

    % ===== GET SUBJECT AND SURFACE =====
    [sSubject, iSubject] = bst_get('Subject', SubjectName);
    if isempty(iSubject)
        bst_report('Error', sProcess, [], ['Subject "' SubjectName '" does not exist.']);
        return;
    end
    if isempty(sSubject.Surface)
        bst_report('Error', sProcess, [], ['No surfaces for subject "' SubjectName '".']);
        return;
    end

    [SurfaceFile, errMsg] = GetSurfaceFile(sSubject, SurfType, sProcess);
    if ~isempty(errMsg)
        bst_report('Error', sProcess, [], errMsg);
        return;
    end

    % ===== LOAD EIGENMODES =====
    [Eigenmodes, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if ~isComputed
        bst_report('Error', sProcess, [], ...
            ['No eigenmodes on surface: ' SurfaceFile '. Run "Compute eigenmodes" first.']);
        return;
    end

    % Clamp range
    ModeRangeLow  = max(1, ModeRangeLow);
    ModeRangeHigh = min(Eigenmodes.nModes, ModeRangeHigh);
    iModes = ModeRangeLow:ModeRangeHigh;
    nExport = length(iModes);
    if nExport < 1
        bst_report('Error', sProcess, [], 'Empty mode range.');
        return;
    end

    % ===== BUILD RESULTS FILE =====
    Phi = Eigenmodes.Vectors(:, iModes);   % [nV x nExport]
    lambdas = Eigenmodes.Values(iModes);

    % Optionally normalize each mode for better colormap contrast
    if DoNormalize
        for k = 1:nExport
            maxVal = max(abs(Phi(:, k)));
            if maxVal > 0
                Phi(:, k) = Phi(:, k) / maxVal;
            end
        end
    end

    % Create results structure
    ResMat = db_template('resultsmat');
    ResMat.ImageGridAmp  = Phi;
    ResMat.ImagingKernel = [];
    ResMat.nComponents   = 1;
    ResMat.Time          = iModes;   % mode indices as "time" axis
    ResMat.SurfaceFile   = SurfaceFile;
    ResMat.HeadModelType = 'surface';
    ResMat.GoodChannel   = [];
    ResMat.ChannelFlag   = [];
    ResMat.DataFile      = [];
    ResMat.nAvg          = 1;
    ResMat.Leff          = 1;
    ResMat.ColormapType  = 'source';

    % Build descriptive comment
    if DoNormalize
        normStr = ', normalized';
    else
        normStr = '';
    end
    ResMat.Comment = sprintf('Eigenmodes %d-%d (%s mass, %d total%s)', ...
        ModeRangeLow, ModeRangeHigh, Eigenmodes.MassType, Eigenmodes.nModes, normStr);

    % Add history
    ResMat = bst_history('add', ResMat, 'eigenmodes_view', ...
        sprintf('Exported %d eigenmodes [%d-%d] from %s', ...
        nExport, ModeRangeLow, ModeRangeHigh, SurfaceFile));
    ResMat = bst_history('add', ResMat, 'eigenmodes_view', ...
        sprintf('Eigenvalue range: [%.2f, %.2f]', lambdas(1), lambdas(end)));

    % ===== SAVE TO INTRA-SUBJECT STUDY =====
    % Get the intra-subject study (where anatomy results go)
    [sStudy, iStudy] = bst_get('AnalysisIntraStudy', iSubject);
    if isempty(iStudy)
        bst_report('Error', sProcess, [], 'Could not find intra-subject study.');
        return;
    end

    % Generate output filename
    StudyDir = bst_fileparts(file_fullpath(sStudy.FileName));
    OutputFile = bst_process('GetNewFilename', StudyDir, 'results_eigenmodes');
    bst_save(OutputFile, ResMat, 'v6');
    db_add_data(iStudy, OutputFile, ResMat);

    bst_report('Info', sProcess, [], ...
        sprintf('Exported %d eigenmodes to: %s', nExport, file_short(OutputFile)));

    OutputFiles = {file_short(OutputFile)};
end


%% ===== HELPER: GET SURFACE FILE =====
function [SurfaceFile, errMsg] = GetSurfaceFile(sSubject, SurfType, sProcess)
    SurfaceFile = '';
    errMsg = '';
    switch lower(SurfType)
        case 'cortex'
            if isempty(sSubject.iCortex) || sSubject.iCortex < 1 || sSubject.iCortex > length(sSubject.Surface)
                errMsg = 'No default cortex surface.';
                return;
            end
            SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
        case 'innerskull'
            if isempty(sSubject.iInnerSkull) || sSubject.iInnerSkull < 1 || sSubject.iInnerSkull > length(sSubject.Surface)
                errMsg = 'No inner skull surface.';
                return;
            end
            SurfaceFile = sSubject.Surface(sSubject.iInnerSkull).FileName;
        case 'outerskull'
            if isempty(sSubject.iOuterSkull) || sSubject.iOuterSkull < 1 || sSubject.iOuterSkull > length(sSubject.Surface)
                errMsg = 'No outer skull surface.';
                return;
            end
            SurfaceFile = sSubject.Surface(sSubject.iOuterSkull).FileName;
        case 'scalp'
            if isempty(sSubject.iScalp) || sSubject.iScalp < 1 || sSubject.iScalp > length(sSubject.Surface)
                errMsg = 'No scalp surface.';
                return;
            end
            SurfaceFile = sSubject.Surface(sSubject.iScalp).FileName;
        case 'other'
            iSurf = sProcess.options.isurfidx.Value{1};
            if isempty(iSurf) || iSurf < 1 || iSurf > length(sSubject.Surface)
                errMsg = sprintf('Surface index %d out of range.', iSurf);
                return;
            end
            SurfaceFile = sSubject.Surface(iSurf).FileName;
        otherwise
            errMsg = ['Unknown surface type: ' SurfType];
    end
end
