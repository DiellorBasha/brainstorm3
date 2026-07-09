function varargout = process_import_pet( varargin )
% PROCESS_IMPORT_PET: Import a PET volume and run the recommended PET pipeline.
%
% The recommended pipeline (all in volume space, benchmarked r=0.99 vs VLPP):
%   raw import -> realign+aggregate -> coregister+reslice to T1 -> PVC ->
%   SUVR (tracer-aware reference) -> project to cortical surface.
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
    sProcess.options.pvcmethod.Value   = {'mg', {'Muller-Gartner','GTM (Rousset)','None'; 'mg','gtm','none'}};
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
    sProcess.options.doproject.Value   = 1;
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


%% ===== ENGINE: import + pipeline =====
function [OutputFiles, errMsg] = Compute(iSubject, PetFile, FileFormat, Comment, Options)
    OutputFiles = {};
    errMsg = '';
    % ----- 1. Raw import (import_pet owns load + metadata) -----
    bst_progress('start', 'Import & process PET', 'Importing PET volume...');
    [DbPetFile, sMri] = import_pet(iSubject, PetFile, FileFormat, 0, 0, Comment);
    if isempty(DbPetFile) || iscell(DbPetFile)
        errMsg = ['Cannot import PET file: "' PetFile '".'];
        bst_progress('stop');
        return;
    end
    OutputFiles{end+1} = DbPetFile;
    % Raw import only?
    if ~Options.RunPipeline
        bst_progress('stop');
        return;
    end
    % ----- Prerequisite: reference MRI (T1) -----
    sSubject = bst_get('Subject', iSubject);
    if isempty(sSubject.Anatomy) || isempty(sSubject.iAnatomy)
        errMsg = 'No reference MRI in this subject: cannot run the PET pipeline. Raw volume imported.';
        bst_progress('stop');
        return;
    end
    MriFile = sSubject.Anatomy(sSubject.iAnatomy).FileName;
    % ----- 2. Realign + aggregate (skip if already static 3D) -----
    if isfield(sMri, 'Cube') && (size(sMri.Cube, 4) > 1)
        bst_progress('text', 'Realigning and aggregating frames...');
        PetAggFile = mri_realign(DbPetFile, 'spm_realign', 0, Options.Aggregate);
    else
        PetAggFile = DbPetFile;   % single-volume PET: nothing to realign/aggregate
    end
    % ----- 3. Coregister + reslice to the T1 -----
    if ~strcmpi(Options.Register, 'ignore')
        bst_progress('text', 'Coregistering PET to MRI...');
        PetCoregFile = mri_coregister(PetAggFile, MriFile, Options.Register, Options.Reslice);
    else
        PetCoregFile = PetAggFile;
    end
    OutputFiles{end+1} = PetCoregFile;
    % ----- 4. Partial volume correction -----
    PetPvcFile = PetCoregFile;
    if ~strcmpi(Options.PvcMethod, 'none')
        bst_progress('text', 'Partial volume correction...');
        pvcFwhm = Options.PvcFwhm;
        if isempty(pvcFwhm) || (pvcFwhm <= 0)
            pvcFwhm = [];   % auto: derive from scanner metadata inside pet_pvc
        end
        [PetPvcFile, errPvc] = pet_pvc(PetCoregFile, MriFile, pvcFwhm, struct('method', Options.PvcMethod));
        if ~isempty(errPvc) || isempty(PetPvcFile)
            % Non-fatal: continue SUVR on the uncorrected (coregistered) volume
            disp(['BST> PET PVC failed: ' errPvc '. Continuing without PVC.']);
            PetPvcFile = PetCoregFile;
        else
            OutputFiles{end+1} = PetPvcFile;
        end
    end
    % ----- 5. SUVR (tracer-aware reference) + optional surface projection -----
    SuvrRef = Options.SuvrRef;
    if isempty(SuvrRef)
        tracerName = '';
        if isfield(sMri, 'PET') && isstruct(sMri.PET) && isfield(sMri.PET, 'Tracer')
            tracerName = sMri.PET.Tracer.Name;
        end
        SuvrRef = GetTracerReference(tracerName);
    end
    bst_progress('text', 'Computing SUVR...');
    [PetSuvrFile, errSuvr, SurfFile] = pet_process(PetPvcFile, 'ASEG', SuvrRef, 'Brainmask', 1, Options.DoProject);
    if ~isempty(errSuvr)
        errMsg = ['SUVR step failed: ' errSuvr];
        bst_progress('stop');
        return;
    end
    if ~isempty(PetSuvrFile)
        OutputFiles{end+1} = PetSuvrFile;
    end
    if ~isempty(SurfFile)
        OutputFiles{end+1} = SurfFile;
    end
    bst_progress('stop');
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
    Options.PvcMethod   = 'mg';
    Options.PvcFwhm     = [];      % auto from scanner metadata
    Options.SuvrRef     = '';      % tracer-aware default
    Options.DoProject   = 1;
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
