function Out = preventad_pet_import(BidsPetDir, SubjectName, Opts)
% PREVENTAD_PET_IMPORT: Import + process PREVENT-AD PET volumes for ONE subject into
% the EXISTING Brainstorm subject of the same name (registering PET to the SAME
% FreeSurfer/ico5 anatomy already imported by the MEG pipeline).
%
% For the matched subject this:
%   1. resolves the existing subject (must already exist from the MEG import) and
%      its structural T1 + ASEG volume atlas,
%   2. for each PET tracer volume (BIDS trc-* under ses-*/pet/):
%        import_mri (as PET) -> mri_realign (frames -> mean) ->
%        mri_coregister (to the subject T1) -> pet_process (SUVR + surface).
%
% PET volumes attach to the subject's anatomy, so the surface-projected SUVR lands
% on the SAME ico5 cortex as the MEG/Dirac source maps (per-subject MEG<->amyloid<->tau
% fusion on a shared mesh). No anatomy is (re)imported.
%
% USAGE:  Out = preventad_pet_import(BidsPetDir, 'sub-MTL0002')
%         Out = preventad_pet_import(BidsPetDir, 'sub-MTL0002', Opts)
%
% INPUTS:
%   BidsPetDir : PREVENT-AD PET BIDS root (default: the SpikeData-2 pet dataset).
%   SubjectName: BIDS subject id, e.g. 'sub-MTL0002' (must already exist in the
%                'preventad' protocol with FreeSurfer anatomy).
%   Opts (all optional; defaults shown):
%       .AtlasName   'ASEG'        anatomy Comment of the volume atlas for SUVR/mask
%       .RefROI      'Cerebellum'  SUVR reference region (amyloid & tau: cerebellar)
%       .MaskROI     'Brainmask'   masking ROI (used when ApplyMask=1)
%       .ApplyMask   1             apply the brain mask
%       .DoProject   1             project SUVR to the cortex surface
%       .DoSUVR      1             run pet_process (SUVR); 0 = stop after coregister
%       .Aggregation 'mean'        dynamic-frame aggregation (mri_realign)
%       .CoregMethod 'spm'         mri_coregister method
%
% OUTPUT:
%   Out : struct array (one per processed tracer) with fields
%         .Subject .Tracer .Imported .Aggregated .Coregistered .Suvr .Surface .Skipped
%
% NOTE (verify once, after the MEG batch): SUVR needs an ASEG *_volatlas on the
% subject. The MEG import used process_import_bids (icosphere); confirm that path
% set isVolumeAtlas so the ASEG atlas exists, else pet_process errors here.
%
% SEE ALSO: tutorial_pet_processing, import_mri, mri_realign, mri_coregister, pet_process
%
% Author: Diellor Basha, 2026

    Out = struct('Subject',{},'Tracer',{},'Imported',{},'Aggregated',{}, ...
                 'Coregistered',{},'Suvr',{},'Surface',{},'Skipped',{});

    % ---- defaults ----
    if (nargin < 1) || isempty(BidsPetDir)
        BidsPetDir = '/Volumes/SpikeData-2/workspace/library/datasets/preventad/pet';
    end
    if (nargin < 2) || isempty(SubjectName)
        error('SubjectName is required (e.g. ''sub-MTL0002'').');
    end
    if (nargin < 3) || isempty(Opts), Opts = struct(); end
    Def = struct('AtlasName','ASEG', 'RefROI','Cerebellum', 'MaskROI','Brainmask', ...
                 'ApplyMask',1, 'DoProject',1, 'DoSUVR',1, 'Aggregation','mean', 'CoregMethod','spm');
    fn = fieldnames(Def);
    for i = 1:numel(fn)
        if ~isfield(Opts, fn{i}) || isempty(Opts.(fn{i})), Opts.(fn{i}) = Def.(fn{i}); end
    end
    if ~file_exist(BidsPetDir)
        error('PET BIDS directory not found: %s', BidsPetDir);
    end

    % ---- select existing protocol ----
    ProtocolName = 'preventad';
    if ~brainstorm('status'), brainstorm nogui; end
    iProtocol = bst_get('Protocol', ProtocolName);
    if isempty(iProtocol)
        error('Unknown protocol: %s (run the MEG import first).', ProtocolName);
    end
    if (iProtocol ~= bst_get('iProtocol'))
        gui_brainstorm('SetCurrentProtocol', iProtocol);
    end

    % ---- resolve the EXISTING subject (matched by name) ----
    [sSubject, iSubject] = bst_get('Subject', SubjectName);
    if isempty(iSubject) || isempty(sSubject) || isempty(sSubject.Anatomy)
        error(['Subject "%s" not found (or has no anatomy) in protocol "%s". ' ...
               'PET registers to the EXISTING MEG anatomy; import MEG first.'], SubjectName, ProtocolName);
    end

    % ---- structural T1 (reference for coregistration) ----
    MriFile = local_find_t1(sSubject);
    if isempty(MriFile)
        error('Could not identify the structural T1 anatomy for "%s".', SubjectName);
    end

    % ---- ASEG volume atlas presence (needed by pet_process) ----
    if Opts.DoSUVR
        iAtlas = find(strcmpi({sSubject.Anatomy.Comment}, Opts.AtlasName), 1);
        if isempty(iAtlas)
            error(['Volume atlas "%s" not found for "%s" (needed for SUVR). ' ...
                   'Confirm the FreeSurfer ASEG atlas was imported (isVolumeAtlas), ' ...
                   'or set Opts.DoSUVR=0.'], Opts.AtlasName, SubjectName);
        end
    end

    % ---- discover PET tracer volumes (BIDS ses-*/pet/*_pet.nii.gz) ----
    petFiles = local_list_pet(BidsPetDir, SubjectName);
    if isempty(petFiles)
        error('No PET volumes found under %s/%s/ses-*/pet/.', BidsPetDir, SubjectName);
    end

    % ---- process each tracer ----
    for k = 1:numel(petFiles)
        petFile = petFiles{k};
        tracer  = local_tracer_tag(petFile);

        % idempotency: skip if a SUVR output for this tracer already exists
        [sSubject] = bst_get('Subject', iSubject);   % refresh (anatomy grows per tracer)
        if local_already_done(sSubject, tracer)
            o = local_rec(SubjectName, tracer); o.Skipped = true;
            Out(end+1) = o; %#ok<AGROW>
            bst_progress('text', sprintf('PET %s / %s: already processed, skipping', SubjectName, tracer));
            continue;
        end

        bst_progress('text', sprintf('PET %s / %s: import', SubjectName, tracer));
        % import_mri needs 'PET' in the Comment to set volType=PET (display); also carries the tracer
        Comment = ['PET ' tracer];
        impFile = import_mri(iSubject, petFile, [], 0, 0, Comment);

        bst_progress('text', sprintf('PET %s / %s: realign+aggregate', SubjectName, tracer));
        aggFile = mri_realign(impFile, 'spm_realign', 0, Opts.Aggregation);

        bst_progress('text', sprintf('PET %s / %s: coregister to T1', SubjectName, tracer));
        coregFile = mri_coregister(aggFile, MriFile, Opts.CoregMethod, 1);

        suvrFile = ''; surfFile = '';
        if Opts.DoSUVR
            bst_progress('text', sprintf('PET %s / %s: SUVR + surface', SubjectName, tracer));
            [suvrFile, errMsg, surfFile] = pet_process(coregFile, Opts.AtlasName, ...
                Opts.RefROI, Opts.MaskROI, Opts.ApplyMask, Opts.DoProject);
            if ~isempty(errMsg)
                error('pet_process failed for %s / %s: %s', SubjectName, tracer, errMsg);
            end
        end

        o = local_rec(SubjectName, tracer);
        o.Imported = impFile; o.Aggregated = aggFile; o.Coregistered = coregFile;
        o.Suvr = suvrFile; o.Surface = surfFile; o.Skipped = false;
        Out(end+1) = o; %#ok<AGROW>
    end

    panel_protocols('UpdateNode', 'Subject', iSubject);
end


%% ===== helpers =====
function o = local_rec(subj, tracer)
    o = struct('Subject',subj, 'Tracer',tracer, 'Imported','', 'Aggregated','', ...
               'Coregistered','', 'Suvr','', 'Surface','', 'Skipped',false);
end

function MriFile = local_find_t1(sSubject)
% Pick the structural T1: the anatomy that is NOT a volume atlas and NOT a PET/SUVR
% volume. Prefer the subject's default anatomy (iAnatomy) when it qualifies.
    MriFile = '';
    isCand = true(1, numel(sSubject.Anatomy));
    for i = 1:numel(sSubject.Anatomy)
        fn = lower(sSubject.Anatomy(i).FileName);
        cm = lower(sSubject.Anatomy(i).Comment);
        if ~isempty(strfind(fn,'_volatlas')) || ~isempty(strfind(cm,'pet')) ...
                || ~isempty(strfind(cm,'rescaled')) || ~isempty(strfind(cm,'aseg'))
            isCand(i) = false;
        end
    end
    if ~isempty(sSubject.iAnatomy) && sSubject.iAnatomy>=1 && sSubject.iAnatomy<=numel(isCand) && isCand(sSubject.iAnatomy)
        MriFile = sSubject.Anatomy(sSubject.iAnatomy).FileName; return;
    end
    iT1 = find(isCand, 1);
    if ~isempty(iT1), MriFile = sSubject.Anatomy(iT1).FileName; end
end

function petFiles = local_list_pet(BidsPetDir, SubjectName)
% List ses-*/pet/*_pet.nii.gz (skipping macOS AppleDouble '._*' files).
    petFiles = {};
    sesDirs = dir(fullfile(BidsPetDir, SubjectName, 'ses-*'));
    sesDirs = sesDirs([sesDirs.isdir]);
    if isempty(sesDirs)
        sesDirs = struct('name', {''});   % allow pet/ directly under subject
    end
    for s = 1:numel(sesDirs)
        petDir = fullfile(BidsPetDir, SubjectName, sesDirs(s).name, 'pet');
        f = dir(fullfile(petDir, '*_pet.nii.gz'));
        for i = 1:numel(f)
            if strncmp(f(i).name, '._', 2), continue; end   % AppleDouble metadata
            petFiles{end+1} = fullfile(petDir, f(i).name); %#ok<AGROW>
        end
    end
    petFiles = sort(petFiles);
end

function tracer = local_tracer_tag(petFile)
% Extract the BIDS 'trc-<label>' entity; fall back to the file base.
    [~, base] = bst_fileparts(petFile);            % strips .gz
    [~, base] = bst_fileparts(base);               % strips .nii
    tok = regexp(base, 'trc-([A-Za-z0-9]+)', 'tokens', 'once');
    if ~isempty(tok), tracer = tok{1}; else, tracer = base; end
end

function tf = local_already_done(sSubject, tracer)
% A tracer is "done" if some anatomy comment carries both the tracer tag and the
% SUVR rescale tag ('rescaled', stamped by mri_rescale).
    tf = false;
    t = lower(tracer);
    for i = 1:numel(sSubject.Anatomy)
        cm = lower(sSubject.Anatomy(i).Comment);
        if ~isempty(strfind(cm, t)) && ~isempty(strfind(cm, 'rescaled'))
            tf = true; return;
        end
    end
end
