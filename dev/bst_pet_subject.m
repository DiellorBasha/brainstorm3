function bst_pet_subject(BidsDir, OutputDir, SubjectLabel, Module, varargin)
% BST_PET_SUBJECT  Build one subject's standalone PET Brainstorm protocol:
%                  FreeSurfer anatomy + per-tracer PET SUVR (volume + cortex).
%
% Container worker for the nsp `brainstorm-pet` pathway (BST_PIPELINE selects it).
% Standard positional contract:
%   bst_pet_subject(BidsDir, OutputDir, SubjectLabel, Module, 'BstDir', D, 'BstDbDir', DB)
% BidsDir / Module are accepted for parity but unused.
%
% Unlike preventad_multimodal_subject it does not extend an existing protocol: it creates
% a per-subject protocol, imports the subject's FreeSurfer recon (import_anatomy_fs,
% icosphere cortex, volume atlases ASEG + Desikan-Killiany), then applies the PET steps per tracer
% (bst_pet_tracer; PVC through the pr/pet-pvc-2 API):
%
%   per tracer (trc-18FNAV4694 amyloid, trc-18Fflortaucipir tau):
%     preventad_pet_import   4D dynamic PET, realigned + coregistered to the T1 ("PET <trc>")
%     mri_aggregate 'mean'   static volume                                 ("PET <trc>_mean")
%     NSP_PET_METHOD:
%       'vlpp'  Gaussian smoothing + pet_suvr (cerebellar-cortex mean, no PVC)
%       'gtm'   pet_pvc Method 'gtm' (Rousset GTM on Desikan-Killiany) + robust pet_suvr
%       'mg'    pet_pvc Method 'mg'  (Muller-Gartner, PETPVE12)        + robust pet_suvr
%     PVC PSF: pet_psf_fwhm = scanner FWHM stored at import from the PET metadata
%              (pet_scanner_fwhm: HRRT 2.5 mm; 6 mm when the scanner is not identified),
%              recorded smoothing added in quadrature                  ("PET <trc>_suvr")
%     mri_interp_vol2tess    SUVR projected onto the cortex (condition "PET")
%     regional Desikan SUVR  -> <Subject>_pet_suvr_<trc>.csv
%
% Env inputs (set by the nsp template's apptainer --env):
%   NSP_FS_DIR        FreeSurfer subject dir (holds mri/, surf/, label/)    (required)
%   NSP_ANAT_METHOD   import_anatomy_fs cortex method (default 'icosphere')
%   NSP_ANAT_SOURCE   the FreeSurfer archive it came from (provenance, recorded as-is)
%   NSP_PET_DIR       PET BIDS root holding <Subject>/ses-*/pet/           (required)
%   NSP_PET_METHOD    'vlpp' (default) | 'mg' | 'gtm'
%   NSP_PET_FWHM      'vlpp': smoothing FWHM in mm (default 6)
%                     'gtm'/'mg': 'auto' or empty = scanner PSF from the PET metadata (default);
%                     a number overrides the scanner FWHM (smoothing is still added)
%   NSP_BST_COMMIT    Brainstorm commit for the provenance (default: read from BstDir/.git)
%   NSP_PET_KEEP_4D   '1' keeps the 4D registered base (default '0')
%   NSP_SESSIONS      'anat=ses-BL00A;pet=ses-01' provenance, recorded as-is
%
% OUTPUT (to OutputDir):
%   <Subject>_brainstorm.zip            the PET protocol
%   <Subject>_pet.json                  tracers, sessions (anat, per-tracer PET file), reference values,
%                                       per-tracer Provenance (method, Brainstorm commit, PSF + source,
%                                       smoothing, reference region + SUVR mode)
%   <Subject>_pet_suvr_<tracer>.csv     regional (Desikan) + global cortical SUVR
%   <Subject>_pet_gtm_<tracer>.csv      'gtm' only: per-region observed + corrected (pet_gtm regTable)
%   <Subject>_pet_report.txt            log
%
% Authors: Diellor Basha, 2026 (nsp brainstorm-pet pathway)

p = inputParser;
p.addRequired('BidsDir',      @ischar);
p.addRequired('OutputDir',    @ischar);
p.addRequired('SubjectLabel', @ischar);
p.addRequired('Module',       @ischar);
p.addParameter('BstDir',    '', @ischar);
p.addParameter('BstDbDir',  '', @ischar);
p.addParameter('NVertices', 15000, @isnumeric);   % parity only (icosphere fixes the cortex)
p.parse(BidsDir, OutputDir, SubjectLabel, Module, varargin{:});
opts = p.Results;

SubjectLabel = regexprep(opts.SubjectLabel, '^sub-', '');
SubjectName  = ['sub-' SubjectLabel];

fsDir       = getenv('NSP_FS_DIR');
anatMethod  = getenv('NSP_ANAT_METHOD');  if isempty(anatMethod), anatMethod = 'icosphere'; end
anatSource  = getenv('NSP_ANAT_SOURCE');
petDir      = getenv('NSP_PET_DIR');
petMethod   = lower(getenv('NSP_PET_METHOD'));   if isempty(petMethod), petMethod = 'vlpp'; end
petFwhm     = str2double(getenv('NSP_PET_FWHM'));   % NaN for '' / 'auto'
if isnan(petFwhm)
    if strcmp(petMethod, 'vlpp'), petFwhm = 6; else, petFwhm = []; end   % [] = scanner PSF from metadata
end
keep4D      = strcmp(getenv('NSP_PET_KEEP_4D'), '1');
sessionsStr = getenv('NSP_SESSIONS');
assert(exist(fullfile(fsDir, 'mri', 'T1.mgz'), 'file') == 2, 'NSP_FS_DIR has no mri/T1.mgz: %s', fsDir);
assert(exist(fullfile(petDir, SubjectName), 'dir') == 7, 'No PET for %s under NSP_PET_DIR=%s', SubjectName, petDir);
assert(ismember(petMethod, {'vlpp', 'mg', 'gtm'}), 'NSP_PET_METHOD must be vlpp|mg|gtm, got %s', petMethod);

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

reportFile = fullfile(opts.OutputDir, [SubjectName '_pet_report.txt']);
diary(reportFile); diary on;
fprintf('=== bst_pet_subject: %s ===\n', SubjectName);
fprintf('anatomy  : %s (%s) from %s\n', fsDir, anatMethod, anatSource);
bstCommit = bst_git_commit(opts.BstDir);
fprintf('pet      : %s | method: %s | fwhm: %s | keep 4D: %d\n', petDir, petMethod, mat2str(petFwhm), keep4D);
fprintf('brainstorm: %s @ %s\n', opts.BstDir, bstCommit);
fprintf('sessions : %s\n', sessionsStr);

% ===== Init Brainstorm (headless server; mirrors preventad_multimodal_subject) =====
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
if ~bst_get('isGUI')
    local_shadow_viewers(opts.BstDbDir);
end

[~, fsName] = fileparts(fsDir);
summary = struct('Subject', SubjectName, ...
                 'Sessions', local_parse_sessions(sessionsStr), ...
                 'Anatomy', struct('Source', anatSource, 'FreeSurferSubject', fsName, 'Method', anatMethod), ...
                 'PET', struct('Tracer', {}, 'Method', {}, 'RefValue', {}, ...
                 'GlobalCorticalSUVR', {}, 'SurfaceFile', {}, 'RoiCsv', {}, 'GtmCsv', {}, ...
                 'Provenance', {}, 'Session', {}, 'File', {}));
try
    % ===== Per-subject protocol + FreeSurfer anatomy =====
    ProtocolName = [SubjectName '_pet'];
    iProtocol = gui_brainstorm('CreateProtocol', ProtocolName, 0, 0, opts.BstDbDir);
    [~, iSubject] = db_add_subject(SubjectName, [], 0, 0);
    errMsg = import_anatomy_fs(iSubject, fsDir, [], 0, [], 0, 1, 0, anatMethod);
    assert(isempty(errMsg), 'import_anatomy_fs: %s', errMsg);
    sSubject = bst_get('Subject', iSubject);
    assert(~isempty(sSubject.iCortex), 'No cortex after the FreeSurfer import of %s', fsDir);
    assert(~isempty(local_anat_file(iSubject, 'ASEG')) && ~isempty(local_anat_file(iSubject, 'Desikan-Killiany')), ...
        'FreeSurfer import lacks the ASEG / Desikan-Killiany volumes (needed for SUVR)');
    fprintf('Protocol #%d, subject #%d: anatomy %s imported\n', iProtocol, iSubject, fsName);

    % ===== PET (as in preventad_multimodal_subject) =====
    [isOk, errMsg] = bst_plugin('Load', 'spm12');
    assert(isOk, 'Could not load the spm12 plugin: %s', errMsg);
    spm('defaults', 'PET');
    spm_jobman('initcfg');
    imp = preventad_pet_import(petDir, SubjectName, struct('ProtocolName', ProtocolName));
    assert(~isempty(imp), 'No PET tracer imported for %s', SubjectName);
    for k = 1:numel(imp)
        rec = bst_pet_tracer(iSubject, SubjectName, imp(k).Tracer, ...
                               petMethod, petFwhm, keep4D, opts.OutputDir, bstCommit);
        [rec.Session, rec.File] = local_pet_source(petDir, SubjectName, imp(k).Tracer);
        summary.PET(end+1) = rec; %#ok<AGROW>
    end
    db_save();

    exportZip = fullfile(opts.OutputDir, [SubjectName '_brainstorm.zip']);
    if exist(exportZip, 'file') == 2; delete(exportZip); end
    export_protocol(iProtocol, iSubject, exportZip);
    fprintf('Exported PET protocol -> %s\n', exportZip);

    local_write_json(fullfile(opts.OutputDir, [SubjectName '_pet.json']), summary);
    fprintf('=== DONE: %s ===\n', SubjectName);
    brainstorm stop;
    diary off;
catch ME
    fprintf(2, 'ERROR in bst_pet_subject(%s): %s\n', SubjectName, ME.message);
    for s = 1:numel(ME.stack)
        fprintf(2, '  at %s (line %d)\n', ME.stack(s).name, ME.stack(s).line);
    end
    try brainstorm stop; catch; end
    diary off;
    rethrow(ME);
end
end


%% ===== PET source file of a tracer (the first, in the order preventad_pet_import takes) =====
function [ses, file] = local_pet_source(petDir, SubjectName, tracer)
    ses = ''; file = '';
    d = dir(fullfile(petDir, SubjectName, 'ses-*', 'pet', ['*trc-' tracer '*_pet.nii*']));
    if isempty(d), return; end
    [~, o] = sort(fullfile({d.folder}, {d.name}));
    d = d(o(1));
    file = d.name;
    t = regexp(d.folder, '(ses-[^/\\]+)', 'tokens', 'once');
    if ~isempty(t), ses = t{1}; end
end


%% ===== helpers =====

function f = local_anat_file(iSubject, comment)
    sSubject = bst_get('Subject', iSubject);
    i = find(strcmp({sSubject.Anatomy.Comment}, comment), 1);
    if isempty(i), f = ''; else, f = sSubject.Anatomy(i).FileName; end
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
