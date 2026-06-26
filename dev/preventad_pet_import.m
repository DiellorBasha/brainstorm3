function Out = preventad_pet_import(BidsPetDir, SubjectName, Opts)
% PREVENTAD_PET_IMPORT: Import PREVENT-AD PET as a 4D (dynamic) volume registered and
% resliced to the subject's T1 - the faithful BASE PET data (no averaging, no SUVR).
%
% For each tracer of the matched (already-imported) subject, this:
%   1. imports the original dynamic 4D PET volume (all frames),
%   2. motion-corrects the frames with SPM realign (no aggregation: frames kept),
%   3. coregisters + reslices ALL frames to the subject's FreeSurfer T1
%      (mri_coregister 4D path: SPM estimate + frame-wise mri_reslice),
%   4. labels the result simply "PET <tracer>" (like the original BIDS file) and
%      DELETES the raw/realigned intermediates (which are not usable in Brainstorm
%      without registration). Realign + coregister run on in-memory structs, so no
%      intermediate DB nodes are created.
%
% The result is a 4D, T1-registered, resliced PET volume per tracer - viewable
% frame-by-frame and a faithful base for any downstream PET analysis (SUVR, PVC,
% surface projection) computed separately. No anatomy is (re)imported.
%
% USAGE:  Out = preventad_pet_import(BidsPetDir, 'sub-MTL0002')
%         Out = preventad_pet_import(BidsPetDir, 'sub-MTL0002', Opts)
%
% INPUTS:
%   BidsPetDir : PREVENT-AD PET BIDS root (default: the SpikeData-2 pet dataset).
%   SubjectName: BIDS subject id, e.g. 'sub-MTL0002' (must already exist in the
%                'preventad' protocol with FreeSurfer anatomy).
%   Opts (all optional; defaults shown):
%       .Aggregation       'ignore'  frame aggregation for mri_realign; 'ignore'
%                                     keeps all frames (4D base). Use 'mean' etc. to
%                                     collapse instead.
%       .CoregMethod       'spm'     mri_coregister method
%       .KeepIntermediates 0         if 1, keep the raw + realigned volumes too
%
% OUTPUT:
%   Out : struct array (one per tracer) with fields .Subject .Tracer .Base .Skipped
%
% SEE ALSO: import_mri, mri_realign, mri_coregister
%
% Author: Diellor Basha, 2026

    Out = struct('Subject',{},'Tracer',{},'Base',{},'Skipped',{});

    % ---- defaults ----
    if (nargin < 1) || isempty(BidsPetDir)
        BidsPetDir = '/Volumes/SpikeData-2/workspace/library/datasets/preventad/pet';
    end
    if (nargin < 2) || isempty(SubjectName)
        error('SubjectName is required (e.g. ''sub-MTL0002'').');
    end
    if (nargin < 3) || isempty(Opts), Opts = struct(); end
    Def = struct('Aggregation','ignore', 'CoregMethod','spm', 'KeepIntermediates',0);
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

    % ---- resolve the EXISTING subject + its T1 ----
    [sSubject, iSubject] = bst_get('Subject', SubjectName);
    if isempty(iSubject) || isempty(sSubject) || isempty(sSubject.Anatomy)
        error(['Subject "%s" not found (or has no anatomy) in protocol "%s". ' ...
               'PET registers to the EXISTING MEG anatomy; import MEG first.'], SubjectName, ProtocolName);
    end
    T1file = local_find_t1(sSubject);
    if isempty(T1file)
        error('Could not identify the structural T1 anatomy for "%s".', SubjectName);
    end
    sT1 = in_mri_bst(T1file);

    % ---- discover PET tracer volumes ----
    petFiles = local_list_pet(BidsPetDir, SubjectName);
    if isempty(petFiles)
        error('No PET volumes found under %s/%s/ses-*/pet/.', BidsPetDir, SubjectName);
    end

    % ---- process each tracer ----
    for k = 1:numel(petFiles)
        petFile  = petFiles{k};
        tracer   = local_tracer_tag(petFile);
        baseName = ['PET ' tracer];

        % idempotency: skip if the 4D base already exists for this tracer
        sSubject = bst_get('Subject', iSubject);
        if any(strcmp({sSubject.Anatomy.Comment}, baseName))
            o = local_rec(SubjectName, tracer); o.Skipped = true;
            Out(end+1) = o; %#ok<AGROW>
            bst_progress('text', sprintf('PET %s / %s: base already exists, skipping', SubjectName, tracer));
            continue;
        end

        % 1) import the raw dynamic 4D volume (DB node, removed below unless kept)
        bst_progress('text', sprintf('PET %s / %s: import (4D)', SubjectName, tracer));
        impFile = import_mri(iSubject, petFile, [], 0, 0, [baseName ' (raw)']);
        sImp = in_mri_bst(impFile);

        % 2) realign frames (keep all frames; in-memory struct -> not saved)
        bst_progress('text', sprintf('PET %s / %s: realign frames', SubjectName, tracer));
        sAgg = mri_realign(sImp, 'spm_realign', 0, Opts.Aggregation);

        % 3) coregister + reslice ALL frames to T1 (in-memory struct -> not saved)
        bst_progress('text', sprintf('PET %s / %s: coregister 4D -> T1', SubjectName, tracer));
        [~, errMsg, ~, sBase] = mri_coregister(sAgg, sT1, Opts.CoregMethod, 1);
        if ~isempty(errMsg) || isempty(sBase)
            local_delete_anat(iSubject, {impFile});   % drop the raw import on failure
            error('Coregistration failed for %s / %s: %s', SubjectName, tracer, errMsg);
        end

        % 4) drop intermediates, save the 4D base labeled like the original tracer
        if ~Opts.KeepIntermediates
            local_delete_anat(iSubject, {impFile});
        end
        % Carry the PET metadata onto the surviving base node. Realign/coregister/
        % reslice preserve the frame count, so the timing in sImp.PET stays valid;
        % the raw import it was read from is deleted above.
        if isfield(sImp, 'PET') && ~isempty(sImp.PET)
            sBase.PET = sImp.PET;
        end
        sBase.Comment = baseName;
        baseFile = local_save_anat(iSubject, T1file, sBase);

        o = local_rec(SubjectName, tracer);
        o.Base = baseFile; o.Skipped = false;
        Out(end+1) = o; %#ok<AGROW>
    end

    panel_protocols('UpdateNode', 'Subject', iSubject);
end


%% ===== helpers =====
function o = local_rec(subj, tracer)
    o = struct('Subject',subj, 'Tracer',tracer, 'Base','', 'Skipped',false);
end

function MriFile = local_find_t1(sSubject)
% Structural T1: anatomy that is not a volume atlas and not a PET/SUVR volume.
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
        sesDirs = struct('name', {''});
    end
    for s = 1:numel(sesDirs)
        petDir = fullfile(BidsPetDir, SubjectName, sesDirs(s).name, 'pet');
        f = dir(fullfile(petDir, '*_pet.nii.gz'));
        for i = 1:numel(f)
            if strncmp(f(i).name, '._', 2), continue; end
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

function local_delete_anat(iSubject, fileList)
% Delete anatomy volume(s) from disk + DB (descending index, fix iAnatomy).
    if isempty(fileList), return; end
    fileList = cellfun(@file_short, fileList, 'UniformOutput', 0);
    sSubject = bst_get('Subject', iSubject);
    file_delete(cellfun(@file_fullpath, fileList, 'UniformOutput', 0), 1);
    [~, iDel] = ismember(fileList, {sSubject.Anatomy.FileName});
    iDel = sort(iDel(iDel > 0), 'descend');
    for k = iDel(:)'
        sSubject.Anatomy(k) = [];
        if ~isempty(sSubject.iAnatomy) && (k < sSubject.iAnatomy)
            sSubject.iAnatomy = sSubject.iAnatomy - 1;
        end
    end
    bst_set('Subject', iSubject, sSubject);
end

function MriFileOut = local_save_anat(iSubject, RefFile, sMri)
% Save an MRI struct as a new anatomy node next to a reference file (e.g. the T1),
% with a unique comment. Mirrors the save pattern of mri_realign/mri_coregister.
    sSubject = bst_get('Subject', iSubject);
    sMri.Comment = file_unique(sMri.Comment, {sSubject.Anatomy.Comment});
    folder = bst_fileparts(file_fullpath(RefFile));
    % The GUI tree types a volume as PET (node 'volpet', with the PET processing
    % options) PURELY from the '_volpet' tag in its FILENAME (node_create_subject);
    % there is no struct field. Keep that tag so the registered 4D base is recognized
    % as PET, not a plain MRI. The Comment stays the clean "PET <tracer>".
    MriFileFull = file_unique(bst_fullfile(folder, ['subjectimage_' file_standardize(sMri.Comment) '_volpet.mat']));
    out_mri_bst(sMri, MriFileFull);
    iAnatomy = length(sSubject.Anatomy) + 1;
    sSubject.Anatomy(iAnatomy) = db_template('Anatomy');
    sSubject.Anatomy(iAnatomy).FileName = file_short(MriFileFull);
    sSubject.Anatomy(iAnatomy).Comment  = sMri.Comment;
    bst_set('Subject', iSubject, sSubject);
    MriFileOut = file_short(MriFileFull);
end
