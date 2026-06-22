function varargout = process_eigenmodes( varargin )
% PROCESS_EIGENMODES: Compute Laplace-Beltrami eigenmodes on a cortical surface.
%
% USAGE:  sProcess = process_eigenmodes('GetDescription')
%       OutputFiles = process_eigenmodes('Run', sProcess, sInputs)
%
% DESCRIPTION:
%     Computes the Laplace-Beltrami eigenmodes (eigenvectors and eigenvalues)
%     of a triangulated cortical surface and stores them in the surface file.
%     These eigenmodes form an orthonormal basis for scalar functions on the
%     surface, analogous to Fourier modes on flat domains. They enable:
%       - Spectral analysis of source-mapped brain activity
%       - Spatial frequency filtering (low-pass, band-pass)
%       - Compact representation of cortical patterns
%       - Heat kernel smoothing and diffusion
%
%     The pipeline validates the mesh (repair is opt-in via the Repair option),
%     assembles the cotangent Laplacian and mass matrix, solves the generalized
%     eigenvalue problem L*phi = lambda*M*phi, and stores the result.
%
% SEE ALSO: tess_eigenmodes, out_tess_eigenmodes, in_tess_eigenmodes, tess_laplacian

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
    sProcess.Comment     = 'Compute eigenmodes';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'Import', 'Import anatomy'};
    sProcess.Index       = 25;
    sProcess.Description = '';
    % Definition of the input accepted by this process
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

    % === SURFACE INDEX (for 'other' type) ===
    sProcess.options.isurfidx.Comment = 'Surface index (for Other type): ';
    sProcess.options.isurfidx.Type    = 'value';
    sProcess.options.isurfidx.Value   = {1, '', 0};

    % === SEPARATOR ===
    sProcess.options.sep1.Type = 'separator';

    % === NUMBER OF MODES ===
    sProcess.options.nmodes.Comment = 'Number of eigenmodes: ';
    sProcess.options.nmodes.Type    = 'value';
    sProcess.options.nmodes.Value   = {300, '', 0};

    % === MASS MATRIX TYPE ===
    sProcess.options.masstype.Comment = {'Barycentric (fast, robust)', 'Voronoi (accurate, Meyer 2003)', 'Galerkin (consistent FEM)', 'Mass matrix:'; ...
                                         'barycentric', 'voronoi', 'galerkin', ''};
    sProcess.options.masstype.Type    = 'radio_linelabel';
    sProcess.options.masstype.Value   = 'barycentric';

    % === REMOVE DC MODES ===
    sProcess.options.removedc.Comment = 'Remove DC (constant) modes';
    sProcess.options.removedc.Type    = 'checkbox';
    sProcess.options.removedc.Value   = 1;

    % === REPAIR (RISKY, OPT-IN) ===
    sProcess.options.repair.Comment = 'Attempt repair if surface is non-manifold (risky: changes vertex count)';
    sProcess.options.repair.Type    = 'checkbox';
    sProcess.options.repair.Value   = 0;

    % === OVERWRITE ===
    sProcess.options.overwrite.Comment = 'Overwrite existing eigenmodes';
    sProcess.options.overwrite.Type    = 'checkbox';
    sProcess.options.overwrite.Value   = 0;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sprintf('Compute eigenmodes (%d modes, %s mass)', ...
        sProcess.options.nmodes.Value{1}, ...
        sProcess.options.masstype.Value);
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};

    % ===== GET OPTIONS =====
    SubjectName = file_standardize(sProcess.options.subjectname.Value);
    if isempty(SubjectName)
        bst_report('Error', sProcess, [], 'Subject name is empty.');
        return;
    end
    nModes    = sProcess.options.nmodes.Value{1};
    MassType  = sProcess.options.masstype.Value;
    RemoveDC  = sProcess.options.removedc.Value;
    Repair    = sProcess.options.repair.Value;
    Overwrite = sProcess.options.overwrite.Value;
    SurfType  = sProcess.options.surfacetype.Value;

    % Validate nModes
    if isempty(nModes) || nModes < 1
        bst_report('Error', sProcess, [], 'Number of eigenmodes must be >= 1.');
        return;
    end

    % ===== GET SUBJECT =====
    [sSubject, iSubject] = bst_get('Subject', SubjectName);
    if isempty(iSubject)
        bst_report('Error', sProcess, [], ['Subject "' SubjectName '" does not exist.']);
        return;
    end
    if isempty(sSubject.Surface)
        bst_report('Error', sProcess, [], ['No surfaces available for subject "' SubjectName '".']);
        return;
    end

    % ===== SELECT SURFACE =====
    [SurfaceFile, errMsg] = GetSurfaceFile(sProcess, sSubject, SurfType);
    if ~isempty(errMsg)
        bst_report('Error', sProcess, [], errMsg);
        return;
    end

    % ===== CHECK EXISTING =====
    if ~Overwrite
        [~, isComputed] = in_tess_eigenmodes(SurfaceFile);
        if isComputed
            bst_report('Info', sProcess, [], ...
                ['Eigenmodes already computed for: ' SurfaceFile '. Skipping (set Overwrite to recompute).']);
            OutputFiles = {'import'};
            return;
        end
    end

    % ===== COMPUTE =====
    bst_report('Info', sProcess, [], ...
        sprintf('Computing %d eigenmodes on %s (%s mass)...', nModes, SurfaceFile, MassType));
    errMsg = Compute(SurfaceFile, nModes, MassType, RemoveDC, Repair, false);
    if ~isempty(errMsg)
        bst_report('Error', sProcess, [], errMsg);
        return;
    end

    % ===== UPDATE TREE =====
    % Reload the surface in the database to reflect changes
    db_reload_subjects(iSubject);

    OutputFiles = {'import'};
end


%% ===== COMPUTE =====
function errMsg = Compute(SurfaceFile, nModes, MassType, RemoveDC, Repair, isInteractive)
    errMsg = '';

    % Load surface
    SurfaceFileFull = file_fullpath(SurfaceFile);
    if ~exist(SurfaceFileFull, 'file')
        errMsg = ['Surface file not found: ' SurfaceFile];
        return;
    end
    TessMat = load(SurfaceFileFull, 'Vertices', 'Faces');
    if ~isfield(TessMat, 'Vertices') || ~isfield(TessMat, 'Faces')
        errMsg = 'Surface file does not contain Vertices and Faces.';
        return;
    end

    Vertices = double(TessMat.Vertices);
    Faces    = double(TessMat.Faces);

    % ===== MANIFOLD CHECK (never auto-repair) =====
    [~, ~, isManifold, report] = tess_manifold(Vertices, Faces, 'Repair', 0, 'Verbose', 0);
    if ~isManifold && ~Repair
        if isfield(report, 'summary') && ~isempty(report.summary)
            defects = strjoin(report.summary, '; ');
        else
            defects = 'non-manifold mesh';
        end
        if isInteractive
            resp = java_dialog('question', ...
                ['Surface is not a clean 2-manifold (' defects ').' 10 10 ...
                 'Attempt repair? This changes the vertex count and may desync the ' ...
                 'head model, lead field, and atlas built on this surface.'], ...
                'Compute eigenmodes', [], {'Repair', 'Cancel'}, 'Cancel');
            if ~strcmpi(resp, 'Repair')
                errMsg = 'Cancelled: surface is non-manifold and repair was declined.';
                return;
            end
            Repair = true;
        else
            errMsg = ['Surface is not a clean 2-manifold (' defects '). ' ...
                      'Re-mesh with icosphere downsampling, or enable the repair option.'];
            return;
        end
    end

    % Progress bar: the eigensolver gives no progress feedback, so use an
    % indeterminate (marquee) bar -- the 3-arg 'start' form omits the bounds.
    bst_progress('start', 'Eigenmodes', 'Computing Laplace-Beltrami eigenmodes...');

    % Compute eigenmodes
    try
        [Eigenmodes, ~, ~, Vertices, Faces] = tess_eigenmodes(Vertices, Faces, ...
            'nModes',    nModes, ...
            'MassType',  MassType, ...
            'RemoveDC',  RemoveDC, ...
            'Repair',    Repair, ...
            'Verbose',   1);
    catch ME
        bst_progress('stop');
        errMsg = ['Eigenmode computation failed: ' ME.message];
        return;
    end

    % Save to surface file
    bst_progress('text', 'Saving eigenmodes to surface file...');
    try
        out_tess_eigenmodes(SurfaceFile, Eigenmodes, Vertices, Faces, true);
    catch ME
        bst_progress('stop');
        errMsg = ['Failed to save eigenmodes: ' ME.message];
        return;
    end

    bst_progress('stop');
end


%% ===== COMPUTE INTERACTIVE (called from the cortex right-click menu) =====
function ComputeInteractive(iSubject, SurfaceFile) %#ok<DEFNU>
    % One consolidated dialog: number of modes + mass type
    opts = gui_show_dialog('Compute eigenmodes', @panel_eigenmodes_compute);
    if isempty(opts)
        return;   % cancelled
    end
    nModes   = opts.nModes;
    MassType = opts.MassType;
    % Confirm overwrite if eigenmodes already exist
    [~, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if isComputed
        isOverwrite = java_dialog('confirm', ...
            'Eigenmodes already exist for this surface. Overwrite?', 'Compute eigenmodes');
        if ~isOverwrite
            return;
        end
    end
    % Compute (RemoveDC=true, Repair=false, isInteractive=true)
    errMsg = Compute(SurfaceFile, nModes, MassType, true, false, true);
    if ~isempty(errMsg)
        bst_error(errMsg, 'Compute eigenmodes', 0);
        return;
    end
    % Reflect the new eigenmodes in the database
    db_reload_subjects(iSubject);
end


%% ===== HELPER: GET SURFACE FILE =====
function [SurfaceFile, errMsg] = GetSurfaceFile(sProcess, sSubject, SurfType)
    SurfaceFile = '';
    errMsg = '';

    switch lower(SurfType)
        case 'cortex'
            if isempty(sSubject.iCortex) || sSubject.iCortex < 1 || sSubject.iCortex > length(sSubject.Surface)
                errMsg = 'No default cortex surface set for this subject.';
                return;
            end
            SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;

        case 'innerskull'
            if isempty(sSubject.iInnerSkull) || sSubject.iInnerSkull < 1 || sSubject.iInnerSkull > length(sSubject.Surface)
                errMsg = 'No inner skull surface set for this subject.';
                return;
            end
            SurfaceFile = sSubject.Surface(sSubject.iInnerSkull).FileName;

        case 'outerskull'
            if isempty(sSubject.iOuterSkull) || sSubject.iOuterSkull < 1 || sSubject.iOuterSkull > length(sSubject.Surface)
                errMsg = 'No outer skull surface set for this subject.';
                return;
            end
            SurfaceFile = sSubject.Surface(sSubject.iOuterSkull).FileName;

        case 'scalp'
            if isempty(sSubject.iScalp) || sSubject.iScalp < 1 || sSubject.iScalp > length(sSubject.Surface)
                errMsg = 'No scalp surface set for this subject.';
                return;
            end
            SurfaceFile = sSubject.Surface(sSubject.iScalp).FileName;

        case 'other'
            iSurf = sProcess.options.isurfidx.Value{1};
            if isempty(iSurf) || iSurf < 1 || iSurf > length(sSubject.Surface)
                errMsg = sprintf('Surface index %d out of range (subject has %d surfaces).', ...
                    iSurf, length(sSubject.Surface));
                return;
            end
            SurfaceFile = sSubject.Surface(iSurf).FileName;

        otherwise
            errMsg = ['Unknown surface type: ' SurfType];
    end
end
