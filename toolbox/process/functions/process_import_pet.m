function varargout = process_import_pet( varargin )
% PROCESS_IMPORT_PET: Import a PET volume and run the recommended PET pipeline.
%
% The recommended pipeline (all in volume space) mirrors the validated flow and
% produces two anatomy nodes:
%   1. Imported PET   : realign -> mean-aggregate -> coregister + reslice to the MRI.
%   2. Fully processed: PVC -> SUVR (tracer-aware reference) from node 1.
% Surface projection is a separate step, not part of the recommended import.
%
% Interactive entry:  process_import_pet('ComputeInteractive', iSubject)
% Batch entry:        run from the Process tab (GetDescription options).
% Both paths gather options and call the single engine Compute().

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
    sProcess.Comment     = 'Import & process PET';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'Import', 'Import anatomy'};
    sProcess.Index       = 3;
    sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/PetImport';
    sProcess.InputTypes  = {'import'};
    sProcess.OutputTypes = {'import'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 0;
    % Option: Subject name
    sProcess.options.subjectname.Comment = 'Subject name:';
    sProcess.options.subjectname.Type    = 'subjectname';
    sProcess.options.subjectname.Value   = 'NewSubject';
    % Option: PET file
    SelectOptions = {...
        '', ...                            % Filename
        '', ...                            % FileFormat
        'open', ...                        % Dialog type
        'Import PET volume...', ...        % Window title
        'ImportAnat', ...                  % LastUsedDir
        'single', ...                      % Selection mode
        'files', ...                       % Selection type
        bst_get('FileFilters', 'mri'), ... % File formats
        'MriIn'};                          % DefaultFormats
    sProcess.options.petfile.Comment = 'PET volume:';
    sProcess.options.petfile.Type    = 'filename';
    sProcess.options.petfile.Value   = SelectOptions;
    % Option: Volume comment
    sProcess.options.comment.Comment = 'Volume name (empty = default name):';
    sProcess.options.comment.Type    = 'text';
    sProcess.options.comment.Value   = '';
    % Option: Processing mode (Recommended reveals nothing; Advanced reveals the knobs)
    sProcess.options.mode.Comment = {'Recommended', 'Advanced', 'Processing:'; 'recommended', 'advanced', ''};
    sProcess.options.mode.Type    = 'radio_linelabel';
    sProcess.options.mode.Value   = 'recommended';
    sProcess.options.mode.Controller.advanced = 'advanced';
    % ===== ADVANCED OPTIONS (Class 'advanced') =====
    % Run pipeline (off = raw import only)
    sProcess.options.doprocess.Comment = 'Run recommended processing (off = raw import only)';
    sProcess.options.doprocess.Type    = 'checkbox';
    sProcess.options.doprocess.Value   = 1;
    sProcess.options.doprocess.Class   = 'advanced';
    % Aggregate frames
    sProcess.options.aggregate.Comment = 'Aggregate frames: ';
    sProcess.options.aggregate.Type    = 'combobox_label';
    sProcess.options.aggregate.Value   = {'mean', {'Mean','Sum','Median','Max','Min','First','Last'; ...
                                                   'mean','sum','median','max','min','first','last'}};
    sProcess.options.aggregate.Class   = 'advanced';
    % Registration method
    sProcess.options.register.Comment = 'Register to MRI: ';
    sProcess.options.register.Type    = 'combobox_label';
    sProcess.options.register.Value   = {'spm', {'SPM','MNI','Ignore'; 'spm','mni','ignore'}};
    sProcess.options.register.Class   = 'advanced';
    % Reslice
    sProcess.options.reslice.Comment = 'Reslice on the MRI grid';
    sProcess.options.reslice.Type    = 'checkbox';
    sProcess.options.reslice.Value   = 1;
    sProcess.options.reslice.Class   = 'advanced';
    % PVC method
    sProcess.options.pvcmethod.Comment = 'Partial volume correction: ';
    sProcess.options.pvcmethod.Type    = 'combobox_label';
    sProcess.options.pvcmethod.Value   = {'gtm', {'GTM (Rousset)','Muller-Gartner','None'; 'gtm','mg','none'}};
    sProcess.options.pvcmethod.Class   = 'advanced';
    % PVC FWHM (0 = auto from scanner metadata)
    sProcess.options.pvcfwhm.Comment = 'PVC PSF FWHM (0 = auto from scanner): ';
    sProcess.options.pvcfwhm.Type    = 'value';
    sProcess.options.pvcfwhm.Value   = {0, 'mm', 1};
    sProcess.options.pvcfwhm.Class   = 'advanced';
    % SUVR reference (empty = tracer-aware default)
    sProcess.options.suvrref.Comment = 'SUVR reference ROI (empty = tracer default): ';
    sProcess.options.suvrref.Type    = 'text';
    sProcess.options.suvrref.Value   = '';
    sProcess.options.suvrref.Class   = 'advanced';
    % Project to surface
    sProcess.options.doproject.Comment = 'Project SUVR to cortical surface';
    sProcess.options.doproject.Type    = 'checkbox';
    sProcess.options.doproject.Value   = 0;
    sProcess.options.doproject.Class   = 'advanced';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN (batch) =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    % Subject
    SubjectName = file_standardize(sProcess.options.subjectname.Value);
    if isempty(SubjectName)
        bst_report('Error', sProcess, [], 'Subject name is empty.');
        return;
    end
    % PET file
    PetFile = sProcess.options.petfile.Value{1};
    if (length(sProcess.options.petfile.Value) < 2) || isempty(sProcess.options.petfile.Value{2})
        FileFormat = 'ALL';
    else
        FileFormat = sProcess.options.petfile.Value{2};
    end
    if isempty(PetFile)
        bst_report('Error', sProcess, [], 'PET file not selected.');
        return;
    end
    % Comment
    Comment = '';
    if isfield(sProcess.options, 'comment') && ~isempty(sProcess.options.comment.Value)
        Comment = sProcess.options.comment.Value;
    end
    % Get/create subject
    [sSubject, iSubject] = bst_get('Subject', SubjectName);
    if isempty(sSubject)
        [~, iSubject] = db_add_subject(SubjectName);
    end
    if isempty(iSubject)
        bst_report('Error', sProcess, [], ['Cannot create subject "' SubjectName '".']);
        return;
    end
    % Assemble options struct from the process options
    Options = GetOptionsFromProcess(sProcess);
    % Run the engine
    [OutputFiles, errMsg] = Compute(iSubject, PetFile, FileFormat, Comment, Options);
    if ~isempty(errMsg)
        bst_report('Error', sProcess, [], errMsg);
    end
end


%% ===== COMPUTE INTERACTIVE =====
function ComputeInteractive(iSubject) %#ok<DEFNU>
    % Ask for the PET file
    LastUsedDirs = bst_get('LastUsedDirs');
    DefaultFormats = bst_get('DefaultFormats');
    [PetFile, FileFormat] = java_getfile('open', 'Import PET volume...', LastUsedDirs.ImportAnat, ...
        'single', 'files', bst_get('FileFilters', 'mri'), DefaultFormats.MriIn);
    if isempty(PetFile)
        return;
    end
    LastUsedDirs.ImportAnat = bst_fileparts(PetFile);
    bst_set('LastUsedDirs', LastUsedDirs);
    % Show the options dialog (repurposed panel_import_pet); returns an Options struct
    Options = gui_show_dialog('Import & process PET', @panel_import_pet, 1, []);
    if isempty(Options)   % user cancelled
        return;
    end
    % Run the engine
    [~, errMsg] = Compute(iSubject, PetFile, FileFormat, '', Options);
    if ~isempty(errMsg)
        bst_error(errMsg, 'Import & process PET', 0);
    end
    % Refresh the tree
    panel_protocols('UpdateNode', 'Subject', iSubject);
    panel_protocols('SelectNode', [], 'subject', iSubject, -1);
end


%% ===== ENGINE: import + two-node pipeline (matches the validated pipeline) =====
% Produces exactly two anatomy nodes:
%   1. Imported PET   - realigned + mean-aggregated + coregistered to the subject MRI
%                       (== the validated PetAggCoreg volume), done in-memory.
%   2. Fully processed - PVC -> SUVR from node 1. The transient PVC volume that
%                       pet_process creates is removed so only these two nodes remain.
% Import-only mode (RunPipeline=false) saves just the raw volume.
function [OutputFiles, errMsg] = Compute(iSubject, PetFile, FileFormat, Comment, Options)
    OutputFiles = {};
    errMsg = '';
    if isempty(FileFormat)
        FileFormat = 'ALL';
    end
    bst_progress('start', 'Import & process PET', 'Loading PET volume...');

    % ----- Load raw PET into memory (no DB node yet) + capture metadata -----
    sMri = in_mri(PetFile, FileFormat, 0, 0, 1);   % isPet = 1
    if isempty(sMri)
        errMsg = ['Cannot read PET file: "' PetFile '".'];
        bst_progress('stop');
        return;
    end
    petMeta = import_pet('ReadMetadata', PetFile, sMri);   % kept aside: realign/coregister drop sMri.PET
    % Base comment (strip .nii / .gz)
    if isempty(Comment)
        [~, fBase, fExt] = bst_fileparts(PetFile);
        if strcmpi(fExt, '.gz')
            [~, fBase] = bst_fileparts(fBase);
        end
        Comment = fBase;
    end

    % ----- Import-only: save the raw volume as a single node -----
    if ~Options.RunPipeline
        bst_progress('text', 'Saving PET volume...');
        sMri.PET = petMeta;
        OutputFiles{end+1} = SaveVolume(iSubject, sMri, ['PET ' Comment]);
        bst_progress('stop');
        return;
    end

    % ----- Prerequisite: reference MRI -----
    sSubject = bst_get('Subject', iSubject);
    if isempty(sSubject.Anatomy) || isempty(sSubject.iAnatomy)
        errMsg = 'No reference MRI in this subject: cannot process the PET volume.';
        bst_progress('stop');
        return;
    end
    sMriRef = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);

    % ===== NODE 1: realign + mean-aggregate + coregister (validated, in-memory) =====
    if (size(sMri.Cube, 4) > 1)
        bst_progress('text', 'Realigning and aggregating frames...');
        sMri = mri_realign(sMri, 'spm_realign', 0, Options.Aggregate);   % realign + aggregate to 3D
    end
    if ~strcmpi(Options.Register, 'ignore')
        bst_progress('text', 'Coregistering PET to MRI...');
        [sMri, errCoreg] = mri_coregister(sMri, sMriRef, Options.Register, Options.Reslice, 0);
        if ~isempty(errCoreg)
            errMsg = ['Coregistration failed: ' errCoreg];
            bst_progress('stop');
            return;
        end
    end
    sMri.PET = petMeta;
    importedFile = SaveVolume(iSubject, sMri, ['PET ' Comment]);
    OutputFiles{end+1} = importedFile;

    % ===== NODE 2: fully processed = PVC -> SUVR from node 1 =====
    % PVC PSF FWHM from the scanner metadata (derived volumes drop sMri.PET)
    pvcFwhm = Options.PvcFwhm;
    if isempty(pvcFwhm) || (pvcFwhm <= 0)
        pvcFwhm = import_pet('ScannerFwhm', petMeta);
    end
    % Tracer-aware SUVR reference (unless overridden)
    SuvrRef = Options.SuvrRef;
    if isempty(SuvrRef)
        tracerName = '';
        if isstruct(petMeta) && isfield(petMeta, 'Tracer')
            tracerName = petMeta.Tracer.Name;
        end
        SuvrRef = GetTracerReference(tracerName);
    end
    pvcOpts = [];
    if ~strcmpi(Options.PvcMethod, 'none')
        pvcOpts = struct('method', Options.PvcMethod, 'fwhm', pvcFwhm);
    end
    % Snapshot so the transient PVC volume can be removed afterwards
    sSubject = bst_get('Subject', iSubject);
    anatBefore = {sSubject.Anatomy.FileName};
    bst_progress('text', 'Partial volume correction + SUVR...');
    [suvrFile, errSuvr, surfFile] = pet_process(importedFile, 'ASEG', SuvrRef, 'Brainmask', 1, Options.DoProject, pvcOpts);

    % ----- Keep only node 1 (imported) + SUVR; drop the transient PVC volume -----
    sSubject = bst_get('Subject', iSubject);
    newFiles = setdiff({sSubject.Anatomy.FileName}, anatBefore);
    for iNew = 1:numel(newFiles)
        if isempty(suvrFile) || ~strcmp(newFiles{iNew}, file_short(suvrFile))
            DeleteAnatomy(iSubject, newFiles{iNew});   % transient PVC volume
        end
    end
    panel_protocols('UpdateNode', 'Subject', iSubject);

    if isempty(suvrFile)
        errMsg = ['SUVR step failed: ' errSuvr];
        bst_progress('stop');
        return;
    end
    OutputFiles{end+1} = suvrFile;
    if ~isempty(surfFile)
        OutputFiles{end+1} = surfFile;   % optional (Advanced): SUVR projected on the cortex
    end
    if ~isempty(errSuvr)
        disp(['BST> ' errSuvr]);   % non-fatal
    end
    bst_progress('stop');
end


%% ===== save an in-memory volume struct as a new anatomy node =====
function BstFile = SaveVolume(iSubject, sMri, Comment)
    sSubject = bst_get('Subject', iSubject);
    sMri.Comment = Comment;
    % Anatomy folder for this subject
    if ~isempty(sSubject.Anatomy)
        anatDir = bst_fileparts(file_fullpath(sSubject.Anatomy(1).FileName));
    else
        ProtocolInfo = bst_get('ProtocolInfo');
        anatDir = bst_fullfile(ProtocolInfo.SUBJECTS, bst_fileparts(sSubject.FileName));
    end
    MriFileFull = file_unique(bst_fullfile(anatDir, ['subjectimage_' file_standardize(Comment) '.mat']));
    out_mri_bst(sMri, MriFileFull);
    % Register the new volume with the subject
    iAnatomy = length(sSubject.Anatomy) + 1;
    sSubject.Anatomy(iAnatomy) = db_template('Anatomy');
    sSubject.Anatomy(iAnatomy).FileName = file_short(MriFileFull);
    sSubject.Anatomy(iAnatomy).Comment  = Comment;
    bst_set('Subject', iSubject, sSubject);
    db_save();
    BstFile = file_short(MriFileFull);
end


%% ===== delete an anatomy volume node (file + DB entry) =====
function DeleteAnatomy(iSubject, BstFile)
    sSubject = bst_get('Subject', iSubject);
    if isempty(sSubject.Anatomy)
        return;
    end
    iAnat = find(strcmp({sSubject.Anatomy.FileName}, file_short(BstFile)), 1);
    if isempty(iAnat)
        return;
    end
    file_delete(file_fullpath(BstFile), 1);
    sSubject.Anatomy(iAnat) = [];
    % Keep iAnatomy valid after the removal
    if isequal(sSubject.iAnatomy, iAnat)
        sSubject.iAnatomy = [];
    elseif ~isempty(sSubject.iAnatomy) && (sSubject.iAnatomy > iAnat)
        sSubject.iAnatomy = sSubject.iAnatomy - 1;
    end
    bst_set('Subject', iSubject, sSubject);
    db_save();
end


%% ===== assemble options from a process sProcess =====
function Options = GetOptionsFromProcess(sProcess)
    Options = DefaultOptions();
    o = sProcess.options;
    % Recommended mode = defaults; Advanced mode = read the fields
    if isfield(o, 'mode') && strcmpi(o.mode.Value, 'advanced')
        Options.RunPipeline = logical(o.doprocess.Value);
        Options.Aggregate   = o.aggregate.Value{1};
        Options.Register    = o.register.Value{1};
        Options.Reslice     = logical(o.reslice.Value);
        Options.PvcMethod   = o.pvcmethod.Value{1};
        Options.PvcFwhm     = o.pvcfwhm.Value{1};
        Options.SuvrRef     = o.suvrref.Value;
        Options.DoProject   = logical(o.doproject.Value);
    end
end


%% ===== default (recommended) options =====
function Options = DefaultOptions()
    Options.RunPipeline = 1;
    Options.Aggregate   = 'mean';
    Options.Register    = 'spm';
    Options.Reslice     = 1;
    Options.PvcMethod   = 'gtm';   % native (no PETPVE12); MG available in Advanced
    Options.PvcFwhm     = [];      % auto from scanner metadata
    Options.SuvrRef     = '';      % tracer-aware default
    Options.DoProject   = 0;       % projection is a separate step (Advanced opt-in)
end


%% ===== tracer-aware default SUVR reference ROI =====
function refName = GetTracerReference(tracerName)
    % Map a PET tracer to its recommended default SUVR reference ROI (an ASEG
    % atlas region name understood by pet_process). Empty override -> this default.
    %
    % >>> EXPAND THIS REGISTRY <<<  add tracer patterns -> reference ROI as needed,
    % and align the ROI names with the actual ASEG atlas labels in the database.
    refName = 'Cerebellum';   % fallback: cerebellar reference
    if isempty(tracerName) || ~ischar(tracerName)
        return;
    end
    t = lower(tracerName);
    amyloid = {'pib','pittsburgh','av45','av-45','florbetapir','fbp','fbb','florbetaben','fmm','flutemetamol','nav4694','azd4694'};
    tau     = {'av1451','av-1451','flortaucipir','ftp','mk6240','mk-6240','pi2620','pi-2620','ro948','gtp1'};
    if any(cellfun(@(k) ~isempty(strfind(t, k)), amyloid)) %#ok<STREMP>
        refName = 'Cerebellum';         % amyloid: whole cerebellum
    elseif any(cellfun(@(k) ~isempty(strfind(t, k)), tau)) %#ok<STREMP>
        refName = 'Cerebellum';         % tau: cerebellar GM (refine label when available)
    end
end
