function preventad_multimodal_subject(BidsDir, OutputDir, SubjectLabel, Module, varargin)
% PREVENTAD_MULTIMODAL_SUBJECT  Build one PREVENT-AD subject's multimodal Brainstorm
%                               protocol: FreeSurfer anatomy + MEG + DWI fibers + PET.
%
% Container worker for the nsp `brainstorm-multimodal` pathway (BST_PIPELINE selects
% it). Standard positional contract:
%   preventad_multimodal_subject(BidsDir, OutputDir, SubjectLabel, Module, ...
%                                'BstDir', D, 'BstDbDir', DB, 'NVertices', N)
% BidsDir / Module / NVertices are accepted for parity but unused.
%
% It starts from the subject's existing per-subject protocol — the fibers protocol
% (anatomy + MEG + sources + fibers) when the subject has DWI, else the base MEG
% protocol — and adds PET in the SAME anatomy, so every modality shares one
% FreeSurfer T1 / cortex and needs no further coregistration:
%
%   per tracer (trc-18FNAV4694 amyloid, trc-18Fflortaucipir tau):
%     preventad_pet_import   4D dynamic PET, realigned + coregistered to the T1 ("PET <trc>")
%     mri_aggregate 'mean'   static volume                                 ("PET <trc>_mean")
%     Gaussian smoothing + pet_suvr (plain cerebellar-cortex mean, no PVC)  ("PET <trc>_suvr")
%       = the VLPP-style method validated against VLPP (global SUVR r=0.99,
%         dev/benchmarks/demo_vlpp_style.m); NSP_PET_METHOD selects the PVC variants
%     mri_interp_vol2tess    SUVR projected onto the cortex (condition "PET")
%     regional Desikan SUVR  -> <Subject>_pet_suvr_<trc>.csv
%
% Brainstorm functions only (import_protocol, import_mri/mri_realign/mri_coregister via
% preventad_pet_import, mri_aggregate, pet_suvr, pet_pvc/pet_gtm, db_add,
% mri_interp_vol2tess, export_protocol). The worker never reads or writes the databank;
% the nsp template stages inputs read-only and collects outputs to staging.
%
% Env inputs (set by the nsp template's apptainer --env):
%   NSP_PROTOCOL_ZIP  exported per-subject protocol .zip to extend    (required)
%   NSP_PET_DIR       PET BIDS root holding <Subject>/ses-*/pet/       (optional; no PET if empty)
%   NSP_PET_METHOD    'vlpp' (default: smooth + plain SUVR, no PVC) | 'mg' | 'gtm'
%   NSP_PET_FWHM      smoothing FWHM in mm for 'vlpp'                  (default 6)
%   NSP_PET_KEEP_4D   '1' keeps the 4D registered base in the protocol (default '0': the
%                     static mean, SUVR volume and cortical map are kept; the 4D base is
%                     ~0.4 GB per tracer and is reproducible from the raw PET)
%   NSP_SESSIONS      'meg=ses-02;pet=ses-01;dwi=ses-FU48A' provenance, recorded as-is
%   NSP_BASE          which protocol was extended ('fibers' | 'meg'), recorded as-is
%
% OUTPUT (to OutputDir):
%   <Subject>_brainstorm.zip            the multimodal protocol
%   <Subject>_multimodal.json           modalities present, sessions, PET reference values
%   <Subject>_pet_suvr_<tracer>.csv     regional (Desikan) + global cortical SUVR
%   <Subject>_multimodal_report.txt     log
%
% Authors: Diellor Basha, 2026 (nsp brainstorm-multimodal pathway)

% ===== Standard container contract =====
p = inputParser;
p.addRequired('BidsDir',      @ischar);
p.addRequired('OutputDir',    @ischar);
p.addRequired('SubjectLabel', @ischar);
p.addRequired('Module',       @ischar);
p.addParameter('BstDir',    '', @ischar);
p.addParameter('BstDbDir',  '', @ischar);
p.addParameter('NVertices', 15000, @isnumeric);   % parity only
p.parse(BidsDir, OutputDir, SubjectLabel, Module, varargin{:});
opts = p.Results;

SubjectLabel = regexprep(opts.SubjectLabel, '^sub-', '');
SubjectName  = ['sub-' SubjectLabel];

% ===== Inputs from environment =====
protocolZip = getenv('NSP_PROTOCOL_ZIP');
petDir      = getenv('NSP_PET_DIR');
petMethod   = lower(getenv('NSP_PET_METHOD'));   if isempty(petMethod), petMethod = 'vlpp'; end
petFwhm     = str2double(getenv('NSP_PET_FWHM')); if isnan(petFwhm),    petFwhm = 6;         end
keep4D      = strcmp(getenv('NSP_PET_KEEP_4D'), '1');
sessionsStr = getenv('NSP_SESSIONS');
baseKind    = getenv('NSP_BASE');
assert(~isempty(protocolZip) && exist(protocolZip, 'file') == 2, ...
    'NSP_PROTOCOL_ZIP not found: %s', protocolZip);
assert(ismember(petMethod, {'vlpp', 'mg', 'gtm'}), 'NSP_PET_METHOD must be vlpp|mg|gtm, got %s', petMethod);

% ===== Resolve BstDir / BstDbDir (mirror import_fibers_subject) =====
if isempty(opts.BstDir)
    if exist('brainstorm', 'file') == 2
        opts.BstDir = fileparts(which('brainstorm'));
    else
        error('BstDir not specified and brainstorm3 not on MATLAB path');
    end
end
if isempty(opts.BstDbDir)
    slurm_tmpdir = getenv('SLURM_TMPDIR');
    if ~isempty(slurm_tmpdir)
        opts.BstDbDir = fullfile(slurm_tmpdir, 'brainstorm_db');
    else
        opts.BstDbDir = fullfile(tempdir, 'brainstorm_db');
    end
end
if exist(opts.BstDbDir, 'dir') ~= 7; mkdir(opts.BstDbDir); end
if exist(opts.OutputDir, 'dir') ~= 7; mkdir(opts.OutputDir); end

reportFile = fullfile(opts.OutputDir, [SubjectName '_multimodal_report.txt']);
diary(reportFile); diary on;
fprintf('=== preventad_multimodal_subject: %s ===\n', SubjectName);
fprintf('protocol : %s (base: %s)\n', protocolZip, baseKind);
fprintf('pet      : %s | method: %s | fwhm: %g | keep 4D: %d\n', petDir, petMethod, petFwhm, keep4D);
fprintf('sessions : %s\n', sessionsStr);

% ===== Init Brainstorm (headless server; mirrors import_fibers_subject) =====
addpath(opts.BstDir);
brainstorm setpath;
bst_user_dir = fullfile(char(java.lang.System.getProperty('user.home')), '.brainstorm');
if exist(bst_user_dir, 'dir') ~= 7; mkdir(bst_user_dir); end
iProtocol             = 0; %#ok<NASGU>
ProtocolsListInfo     = repmat(db_template('ProtocolInfo'), 0);     %#ok<NASGU>
ProtocolsListSubjects = repmat(db_template('ProtocolSubjects'), 0); %#ok<NASGU>
ProtocolsListStudies  = repmat(db_template('ProtocolStudies'), 0);  %#ok<NASGU>
BrainStormDbDir       = opts.BstDbDir; %#ok<NASGU>
DbVersion             = 5.03; %#ok<NASGU>
save(fullfile(bst_user_dir, 'brainstorm.mat'), 'iProtocol', 'ProtocolsListInfo', ...
     'ProtocolsListSubjects', 'ProtocolsListStudies', 'BrainStormDbDir', 'DbVersion');
if ~brainstorm('status'); brainstorm server; end
fprintf('Brainstorm server started. DB = %s\n', opts.BstDbDir);

% pet_pvc / pet_gtm and friends call viewers after saving; headless there is no
% figure to open. Shadow the viewers for this worker only.
if ~bst_get('isGUI')
    local_shadow_viewers(opts.BstDbDir);
end

summary = struct('Subject', SubjectName, 'Base', baseKind, ...
                 'Sessions', local_parse_sessions(sessionsStr), ...
                 'Modalities', struct(), 'PET', struct('Tracer', {}, 'Method', {}, ...
                 'RefValue', {}, 'GlobalCorticalSUVR', {}, 'SurfaceFile', {}, 'RoiCsv', {}));
try
    % ===== Load the per-subject protocol =====
    import_protocol(protocolZip);
    iProtocol = bst_get('iProtocol');
    [sSubject, iSubject] = local_subject(SubjectName);
    fprintf('Loaded protocol #%d, subject #%d: %s\n', iProtocol, iSubject, sSubject.Name);

    % ===== Record what the base protocol already holds =====
    summary.Modalities = local_inventory(iSubject);
    fprintf('Base modalities: %s\n', jsonencode(summary.Modalities));

    % ===== PET =====
    if ~isempty(petDir) && exist(fullfile(petDir, SubjectName), 'dir') == 7
        % SPM (realign/coregister) must be loaded AND its batch system initialised:
        % mri_realign builds cfg_dep batches, which a fresh headless session lacks.
        [isOk, errMsg] = bst_plugin('Load', 'spm12');
        assert(isOk, 'Could not load the spm12 plugin: %s', errMsg);
        spm('defaults', 'PET');
        spm_jobman('initcfg');
        sProtocol = bst_get('ProtocolInfo');
        imp = preventad_pet_import(petDir, SubjectName, struct('ProtocolName', sProtocol.Comment));
        for k = 1:numel(imp)
            summary.PET(end+1) = local_pet_tracer(iSubject, SubjectName, imp(k).Tracer, ...
                                                  petMethod, petFwhm, keep4D, opts.OutputDir); %#ok<AGROW>
        end
        summary.Modalities.PET = {summary.PET.Tracer};
    else
        fprintf('No PET for %s (NSP_PET_DIR=%s)\n', SubjectName, petDir);
        summary.Modalities.PET = {};
    end
    db_save();

    % ===== Export the multimodal protocol =====
    exportZip = fullfile(opts.OutputDir, [SubjectName '_brainstorm.zip']);
    if exist(exportZip, 'file') == 2; delete(exportZip); end
    export_protocol(iProtocol, iSubject, exportZip);
    fprintf('Exported multimodal protocol -> %s\n', exportZip);

    local_write_json(fullfile(opts.OutputDir, [SubjectName '_multimodal.json']), summary);
    fprintf('=== DONE: %s ===\n', SubjectName);
    brainstorm stop;
    diary off;
catch ME
    fprintf(2, 'ERROR in preventad_multimodal_subject(%s): %s\n', SubjectName, ME.message);
    for s = 1:numel(ME.stack)
        fprintf(2, '  at %s (line %d)\n', ME.stack(s).name, ME.stack(s).line);
    end
    try brainstorm stop; catch; end
    diary off;
    rethrow(ME);
end
end


%% ===== PET: one tracer =====
function rec = local_pet_tracer(iSubject, SubjectName, tracer, method, fwhm, keep4D, OutputDir)
    rec = struct('Tracer', tracer, 'Method', method, 'RefValue', NaN, ...
                 'GlobalCorticalSUVR', NaN, 'SurfaceFile', '', 'RoiCsv', '');
    af = @(c) local_anat_file(iSubject, c);
    baseFile = af(['PET ' tracer]);
    assert(~isempty(baseFile), 'PET %s: registered base missing after import', tracer);
    T1File   = local_t1_file(iSubject);
    sAseg    = in_mri_bst(af('ASEG'));

    % static (mean over frames) — kept as its own node
    sMean = mri_aggregate(in_mri_bst(baseFile), 'mean');
    sMean.Comment = ['PET ' tracer '_mean'];
    meanFile = db_add(iSubject, sMean);

    switch method
        case 'vlpp'   % smooth + plain cerebellar-cortex reference (validated vs VLPP)
            sIn = sMean;
            sIn.Cube = imgaussfilt3(double(sMean.Cube(:,:,:,1)), fwhm / 2.355);
            suvrOpts = struct('Erode', 0, 'Robust', 'mean');
        case 'mg'     % Mueller-Gartner PVC, robust eroded reference
            pvcFile = pet_pvc(meanFile, T1File, [], struct());
            sIn = in_mri_bst(pvcFile);  suvrOpts = struct();
        case 'gtm'    % geometric transfer matrix PVC
            pvcFile = pet_gtm(meanFile, [], struct());
            sIn = in_mri_bst(pvcFile);  suvrOpts = struct();
    end
    [sSuvr, info] = pet_suvr(sIn, sAseg, suvrOpts);
    sSuvr.Comment = ['PET ' tracer '_suvr'];
    sSuvr = bst_history('add', sSuvr, 'suvr', sprintf( ...
        'SUVR (%s), ref cerebellar cortex = %.4g', method, info.RefValue));
    suvrFile = db_add(iSubject, sSuvr);
    rec.RefValue = info.RefValue;

    % cortical map on the subject's cortex (condition "PET")
    [surfFile, errProj] = mri_interp_vol2tess(suvrFile, T1File, 'PET', 'SUVR', [0.1 0.8 0.1]);
    if ~isempty(errProj), warning('PET %s: projection failed: %s', tracer, errProj); end
    rec.SurfaceFile = surfFile;

    % regional (Desikan) + global cortical SUVR
    sDK = in_mri_bst(af('Desikan-Killiany'));
    [names, vals] = local_regional(sDK, double(sSuvr.Cube(:,:,:,1)));
    rec.GlobalCorticalSUVR = mean(vals, 'omitnan');
    rec.RoiCsv = fullfile(OutputDir, sprintf('%s_pet_suvr_%s.csv', SubjectName, tracer));
    fid = fopen(rec.RoiCsv, 'w');
    fprintf(fid, 'subject,tracer,method,reference,global_cortical_suvr');
    fprintf(fid, ',%s', names{:}); fprintf(fid, '\n');
    fprintf(fid, '%s,%s,%s,%.6g,%.6g', SubjectName, tracer, method, info.RefValue, rec.GlobalCorticalSUVR);
    fprintf(fid, ',%.6g', vals); fprintf(fid, '\n');
    fclose(fid);
    [~, n, e] = fileparts(rec.RoiCsv); rec.RoiCsv = [n e];

    if ~keep4D
        local_delete_anat(iSubject, baseFile);
    end
    fprintf('PET %s: ref=%.4g global cortical SUVR=%.3f -> %s\n', ...
        tracer, info.RefValue, rec.GlobalCorticalSUVR, surfFile);
end


%% ===== helpers =====
function [sSubject, iSubject] = local_subject(SubjectName)
    [sSubject, iSubject] = bst_get('Subject', SubjectName);
    if isempty(sSubject)   % per-subject exports hold one real subject
        sProt = bst_get('ProtocolSubjects');
        for k = 1:numel(sProt.Subject)
            if ~strcmpi(sProt.Subject(k).Name, bst_get('DirDefaultSubject'))
                iSubject = k; break;
            end
        end
        sSubject = bst_get('Subject', iSubject);
    end
    assert(~isempty(sSubject) && ~isempty(sSubject.iCortex), ...
        'Subject %s has no cortex in the protocol', SubjectName);
end

function M = local_inventory(iSubject)
% What the protocol already holds for this subject (anatomy, MEG, fibers).
    sSubject = bst_get('Subject', iSubject);
    M = struct();
    M.Anatomy = {sSubject.Anatomy.Comment};
    M.Surfaces = {sSubject.Surface.Comment};
    isFib = cellfun(@(c) ~isempty(strfind(lower(c), 'fib')), {sSubject.Surface.FileName});
    M.Fibers = any(isFib);
    [sStudies, ~] = bst_get('StudyWithSubject', sSubject.FileName, 'intra_subject');
    conds = {};
    nRec = 0; nSrc = 0;
    for i = 1:numel(sStudies)
        if ~isempty(sStudies(i).Condition), conds{end+1} = sStudies(i).Condition{1}; end %#ok<AGROW>
        nRec = nRec + numel(sStudies(i).Data);
        nSrc = nSrc + numel(sStudies(i).Result);
    end
    M.MEG = struct('Conditions', {conds}, 'nRecordings', nRec, 'nSourceResults', nSrc);
end

function f = local_anat_file(iSubject, comment)
    sSubject = bst_get('Subject', iSubject);
    i = find(strcmp({sSubject.Anatomy.Comment}, comment), 1);
    if isempty(i), f = ''; else, f = sSubject.Anatomy(i).FileName; end
end

function f = local_t1_file(iSubject)
    sSubject = bst_get('Subject', iSubject);
    f = sSubject.Anatomy(sSubject.iAnatomy).FileName;
end

function local_delete_anat(iSubject, fileName)
    sSubject = bst_get('Subject', iSubject);
    i = find(strcmp({sSubject.Anatomy.FileName}, fileName), 1);
    if isempty(i), return; end
    file_delete(file_fullpath(fileName), 1);
    sSubject.Anatomy(i) = [];
    if sSubject.iAnatomy > i, sSubject.iAnatomy = sSubject.iAnatomy - 1; end
    bst_set('Subject', iSubject, sSubject);
end

function [names, vals] = local_regional(sDK, cube)
% Desikan-Killiany cortical labels (1000-2999), as in preventad_pet_pipeline.
    L = sDK.Labels; v = cell2mat(L(:,1)); nm = L(:,2);
    k = find(v >= 1000 & v < 3000);
    names = nm(k)'; vals = zeros(1, numel(k));
    for i = 1:numel(k), vals(i) = mean(cube(sDK.Cube == v(k(i))), 'omitnan'); end
end

function S = local_parse_sessions(str)
    S = struct();
    for kv = strsplit(str, ';')
        t = strsplit(kv{1}, '=');
        if numel(t) == 2 && ~isempty(t{1}), S.(matlab.lang.makeValidName(t{1})) = t{2}; end
    end
end

function local_write_json(file, s)
    fid = fopen(file, 'w');
    fprintf(fid, '%s\n', jsonencode(s, 'PrettyPrint', true));
    fclose(fid);
end

function local_shadow_viewers(dbDir)
% No-op view_mri / view_surface_data on the path (headless workers only).
    d = fullfile(dbDir, '..', 'nsp_headless_shims');
    if exist(d, 'dir') ~= 7, mkdir(d); end
    for fn = {'view_mri', 'view_surface_data'}
        fid = fopen(fullfile(d, [fn{1} '.m']), 'w');
        fprintf(fid, 'function varargout = %s(varargin)\nvarargout = cell(1, nargout);\nend\n', fn{1});
        fclose(fid);
    end
    addpath(d, '-begin');
end
