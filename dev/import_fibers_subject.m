function import_fibers_subject(BidsDir, OutputDir, SubjectLabel, Module, varargin)
% IMPORT_FIBERS_SUBJECT  Import DWI tractography (.trk) into an existing per-subject
%                        Brainstorm protocol and build structural connectomes,
%                        using established Brainstorm functions.
%
% Container worker for the nsp `brainstorm-fibers` pathway (BST_PIPELINE selects
% it). Standard positional contract:
%   import_fibers_subject(BidsDir, OutputDir, SubjectLabel, Module, ...
%                         'BstDir', D, 'BstDbDir', DB, 'NVertices', N)
% BidsDir / Module are accepted for parity but unused.
%
% Brainstorm functions used (no reimplemented neuroimaging methods):
%   import_protocol                    load the exported per-subject protocol
%   bst_get('Subject')                 resolve the subject index
%   import_fibers                      trk_read/trk_interp/cs_convert/ComputeColor/save
%   in_mri + mri_coregister('spm')     the DWI's anatomical frame onto the subject MRI (NSP_ACPC_T1)
%   cs_convert                         ACPC world -> that volume's MRI coords -> subject SCS
%   in_tess_bst                        load cortex surface + atlases
%   fibers_helper('AssignToScouts')    assign streamline endpoints to scouts
%   export_protocol                    write the augmented protocol back out
% Region node positions use Brainstorm's own scout-seed convention
% (figure_connect: RowLocs = Vertices([Atlas.Scouts.Seed],:)). The only glue is
% tallying Brainstorm's per-fiber Assignment into an NxN streamline-count matrix
% (Brainstorm has no headless connectome writer).
%
% Env inputs (set by the nsp template's apptainer --env):
%   NSP_PROTOCOL_ZIP  exported protocol .zip to augment          (required)
%   NSP_TRK           TrackVis .trk streamlines                  (required)
%   NSP_ACPC_T1       the T1w that defines the .trk's frame (QSIPrep space-ACPC_desc-preproc_T1w).
%                     When set, the fibres are REGISTERED (below) and NSP_CS is ignored.
%   NSP_CS            legacy, without NSP_ACPC_T1 only: import_fibers' CS {scs|mni|world|mri}.
%                     NOTE import_fibers maps 'world' to the 'mri' conversion and trk_read
%                     returns TrackVis voxmm, so the points are taken as millimetres on the
%                     subject MRI's own grid: correct only when the .trk grid IS that grid.
%   NSP_MIN_NEAR_WHITE  registration gate: fraction of fibre ends within 5 mm of the
%                     subject's white surface below which the subject fails (default 0.6)
%   NSP_NPOINTS       points per streamline after resampling     (default 100)
%   NSP_ATLASES       comma-list of cortical atlases for connectomes
%                     (default 'Desikan-Killiany,Destrieux')
%
% OUTPUT (to OutputDir):
%   <Subject>_brainstorm.zip            augmented protocol (anatomy + fibers)
%   <Subject>_connectome_<atlas>.mat    NxN streamline-count matrix + labels
%   <Subject>_fibers_registration.mat   the registration and its check (NSP_ACPC_T1 path)
%
% REGISTRATION (NSP_ACPC_T1): tractography from QSIPrep/QSIRecon lives in space-ACPC —
% the DWI session's T1w rigidly re-oriented to AC-PC — a frame Brainstorm's subject MRI
% knows nothing of. So the points are read in the .trk's own frame (trk_read_rasmm:
% voxmm + vox_to_ras -> ACPC world RAS mm), the ACPC T1w is loaded and co-registered to
% the subject MRI with Brainstorm's mri_coregister('spm', no reslice) — as the PET
% import does — and each point goes ACPC world -> the ACPC volume's MRI coordinates ->
% the subject's SCS through the registered volume (whose SCS mri_coregister set from the
% subject's). Then a GATE: the fibre ends must lie near the subject's white surface
% (ACT ends at grey/white), or the subject fails rather than save misplaced fibres.
%
% Authors: Diellor Basha, 2026 (nsp brainstorm-fibers pathway)

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

% ===== Fiber inputs from environment =====
protocolZip = getenv('NSP_PROTOCOL_ZIP');
trkFile     = getenv('NSP_TRK');
acpcFile    = getenv('NSP_ACPC_T1');
minNear     = str2double(getenv('NSP_MIN_NEAR_WHITE'));
if isnan(minNear); minNear = 0.6; end
csEnv       = lower(getenv('NSP_CS'));
nPointsEnv  = getenv('NSP_NPOINTS');
atlasEnv    = getenv('NSP_ATLASES');
if isempty(csEnv);      csEnv = 'world'; end
if isempty(nPointsEnv); nPoints = 100; else; nPoints = str2double(nPointsEnv); end
if isempty(atlasEnv);   atlasEnv = 'Desikan-Killiany,Destrieux'; end
atlasList = strtrim(strsplit(atlasEnv, ','));
atlasList = atlasList(~cellfun(@isempty, atlasList));

assert(~isempty(protocolZip) && exist(protocolZip, 'file') == 2, ...
    'NSP_PROTOCOL_ZIP not found: %s', protocolZip);
assert(~isempty(trkFile) && exist(trkFile, 'file') == 2, ...
    'NSP_TRK not found: %s', trkFile);

% ===== Resolve BstDir / BstDbDir (mirror preventad_subject) =====
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

reportFile = fullfile(opts.OutputDir, [SubjectName '_fibers_report.txt']);
diary(reportFile); diary on;
fprintf('=== import_fibers_subject: %s ===\n', SubjectName);
fprintf('protocol : %s\n', protocolZip);
fprintf('trk      : %s\n', trkFile);
if isempty(acpcFile)
    fprintf('CS       : %s (UNREGISTERED legacy path) | nPoints: %d | atlases: %s\n', csEnv, nPoints, strjoin(atlasList, ', '));
else
    assert(exist(acpcFile, 'file') == 2, 'NSP_ACPC_T1 not found: %s', acpcFile);
    csEnv = 'scs (registered: ACPC T1w -> subject MRI, spm)';
    fprintf('ACPC T1w : %s | nPoints: %d | atlases: %s\n', acpcFile, nPoints, strjoin(atlasList, ', '));
end

% ===== Init Brainstorm (headless server) =====
% Pre-seed ~/.brainstorm/brainstorm.mat so `brainstorm server` starts without the
% first-run interactive database-directory prompt (mirrors preventad_subject.m).
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

try
    % ===== Load existing protocol + resolve subject (Brainstorm) =====
    import_protocol(protocolZip);
    iProtocol = bst_get('iProtocol');
    fprintf('Loaded protocol #%d from %s\n', iProtocol, protocolZip);

    [sSubject, iSubject] = bst_get('Subject', SubjectName);
    if isempty(sSubject)
        sProt = bst_get('ProtocolSubjects');
        iSubject = 1;
        for k = 1:numel(sProt.Subject)
            if ~strcmpi(sProt.Subject(k).Name, bst_get('DirDefaultSubject'))
                iSubject = k; break;
            end
        end
        sSubject = bst_get('Subject', iSubject);
    end
    assert(~isempty(sSubject) && ~isempty(sSubject.iCortex), ...
        'Subject %s has no cortex in protocol', SubjectName);
    fprintf('Subject #%d: %s (%d surfaces)\n', iSubject, sSubject.Name, numel(sSubject.Surface));

    % ===== Import fibers =====
    if isempty(acpcFile)
        fprintf('Importing fibers via import_fibers (CS=%s, nPoints=%d)...\n', csEnv, nPoints);
        [iNewFibers, OutputFiles, nFibers] = import_fibers(iSubject, {trkFile}, 'TRK', nPoints, csEnv); %#ok<ASGLU>
    else
        [OutputFiles, nFibers, reg] = local_import_registered(iSubject, sSubject, trkFile, acpcFile, nPoints, minNear, opts.OutputDir, SubjectName);
    end
    if iscell(OutputFiles); fibersFile = OutputFiles{1}; else; fibersFile = OutputFiles; end
    fprintf('Imported %d fibers -> %s\n', nFibers, fibersFile);

    % ===== Structural connectomes (Brainstorm AssignToScouts) =====
    % Region nodes = scout seed vertices (Brainstorm's own connectivity node
    % positions). AssignToScouts does the endpoint->region assignment. The NxN
    % streamline count is a tally of Brainstorm's Assignment output.
    FibMat = load(file_fullpath(fibersFile));
    sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);
    for ia = 1:numel(atlasList)
        atlasName = atlasList{ia};
        iAtlas = find(strcmpi({sCortex.Atlas.Name}, atlasName), 1);
        if isempty(iAtlas)
            warning('Atlas "%s" not on cortex (available: %s)', atlasName, strjoin({sCortex.Atlas.Name}, ', '));
            continue;
        end
        scouts  = sCortex.Atlas(iAtlas).Scouts;
        nScouts = numel(scouts);
        if nScouts < 2; continue; end
        % Brainstorm scout seeds (ensure populated, as Brainstorm does: Seed=Vertices(1))
        seeds = zeros(nScouts, 1);
        for is = 1:nScouts
            if isempty(scouts(is).Seed); seeds(is) = scouts(is).Vertices(1);
            else;                        seeds(is) = scouts(is).Seed; end
        end
        centroids = sCortex.Vertices(seeds, :);      % == figure_connect RowLocs
        labels    = {scouts.Label};
        FibMat = fibers_helper('AssignToScouts', FibMat, sprintf('%s_%s', SubjectName, atlasName), centroids);
        asg = FibMat.Scouts(end).Assignment;         % nFibers x 2 scout indices
        C = zeros(nScouts, nScouts);
        valid = all(asg > 0, 2);
        for f = find(valid)'
            a = asg(f,1); b = asg(f,2);
            C(a,b) = C(a,b) + 1;
            if a ~= b; C(b,a) = C(b,a) + 1; end
        end
        outMat = fullfile(opts.OutputDir, sprintf('%s_connectome_%s.mat', SubjectName, atlasName));
        connectome = struct('Matrix', C, 'Labels', {labels}, 'Atlas', atlasName, ...
                            'Subject', SubjectName, 'nFibers', nFibers, ...
                            'Measure', 'streamline_count', 'CS', csEnv); %#ok<NASGU>
        save(outMat, '-struct', 'connectome');
        fprintf('Connectome[%s]: %dx%d, %d assigned streamlines -> %s\n', ...
            atlasName, nScouts, nScouts, sum(valid), outMat);
    end
    bst_save(file_fullpath(fibersFile), FibMat, 'v7');
    db_save();

    % ===== Re-export the augmented protocol (Brainstorm) =====
    exportZip = fullfile(opts.OutputDir, [SubjectName '_brainstorm.zip']);
    if exist(exportZip, 'file') == 2; delete(exportZip); end
    export_protocol(iProtocol, iSubject, exportZip);
    fprintf('Exported augmented protocol -> %s\n', exportZip);

    fprintf('=== DONE: %s (%d fibers) ===\n', SubjectName, nFibers);
    brainstorm stop;
    diary off;
catch ME
    fprintf(2, 'ERROR in import_fibers_subject(%s): %s\n', SubjectName, ME.message);
    for s = 1:numel(ME.stack)
        fprintf(2, '  at %s (line %d)\n', ME.stack(s).name, ME.stack(s).line);
    end
    try brainstorm stop; catch; end
    diary off;
    rethrow(ME);
end
end


function [OutputFiles, nFibers, reg] = local_import_registered(iSubject, sSubject, trkFile, acpcFile, nPoints, minNear, OutputDir, SubjectName)
% The .trk in its own frame -> co-registered to the subject MRI -> SCS; gated; imported.
fprintf('Reading %s in its own frame (voxmm + vox_to_ras)...\n', trkFile);
[tracks, hdr] = trk_read_rasmm(trkFile);
nFibers = numel(tracks);
fprintf('  %d streamlines, voxel order %s, voxel %s mm\n', nFibers, hdr.voxel_order, mat2str(hdr.voxel_size));
P = trk_interp(tracks, nPoints);                     % [nPoints x 3 x nFibers], TrackVis voxmm
clear tracks;
% voxmm -> voxel (centre-based) -> ACPC world RAS mm (nibabel's trackvis-to-rasmm)
A  = hdr.vox_to_ras;
X  = reshape(permute(P, [2 1 3]), 3, []) ./ hdr.voxel_size(:) - 0.5;
clear P;
Y  = A(1:3,1:3) * X + A(1:3,4);                      % RAS mm
clear X;
pts2D = Y' ./ 1000;                                  % Brainstorm world coordinates, metres
clear Y;

% the ACPC T1w onto the subject MRI (Brainstorm + SPM, no reslice)
sMriRef = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
sAcpc   = in_mri(acpcFile, 'ALL', 0, 0);
[isOk, errMsg] = bst_plugin('Load', 'spm12');
assert(isOk, 'Could not load the spm12 plugin: %s', errMsg);
spm('defaults', 'PET');
spm_jobman('initcfg');
[~, errMsg, ~, sAcpcReg] = mri_coregister(sAcpc, sMriRef, 'spm', 0);
if ~isempty(errMsg) || isempty(sAcpcReg)
    error('Registering the ACPC T1w to the subject MRI failed: %s', errMsg);
end
ptsMri = cs_convert(sAcpc,    'world', 'mri', pts2D);   % the ACPC volume's own MRI coordinates
clear pts2D;
ptsScs = cs_convert(sAcpcReg, 'mri',   'scs', ptsMri);  % -> the subject's SCS, through the registration
clear ptsMri;
if isempty(ptsScs)
    error('cs_convert could not map the fibres to SCS (missing vox2ras or SCS on the volumes)');
end
Points = permute(reshape(ptsScs, nPoints, nFibers, 3), [2 1 3]);   % [nFibers x nPoints x 3]
clear ptsScs;

% THE GATE: ACT tractography ends at the grey/white boundary, so the ends must lie near
% the subject's white surface; a mis-registration puts most of them centimetres away.
iWhite = find(~cellfun(@isempty, regexpi({sSubject.Surface.FileName}, 'cortex_white_low')), 1);
if isempty(iWhite); iWhite = find(~cellfun(@isempty, regexpi({sSubject.Surface.FileName}, 'cortex_white')), 1); end
reg = struct('Method', 'mri_coregister spm (no reslice): ACPC T1w -> subject MRI', 'AcpcT1', acpcFile, ...
             'Trk', trkFile, 'VoxToRas', A, 'nFibers', nFibers, 'nPoints', nPoints);
if isempty(iWhite)
    warning('No white surface in the protocol: the registration gate is skipped');
    reg.Gate = 'skipped (no white surface)';
else
    sWhite = in_tess_bst(sSubject.Surface(iWhite).FileName);
    ends = [reshape(Points(:,1,:), [], 3); reshape(Points(:,end,:), [], 3)];
    [~, d] = bst_nearest(sWhite.Vertices, ends, 1, 0);
    d = d * 1000;
    reg.WhiteSurface = sSubject.Surface(iWhite).FileName;
    reg.EndToWhite_mm = struct('median', median(d), 'p75', local_quantile(d, 0.75), 'p95', local_quantile(d, 0.95), ...
                               'within3mm', mean(d < 3), 'within5mm', mean(d < 5));
    fprintf('  ends -> white surface: median %.2f mm, within 3 mm %.1f%%, within 5 mm %.1f%%\n', ...
            reg.EndToWhite_mm.median, 100*reg.EndToWhite_mm.within3mm, 100*reg.EndToWhite_mm.within5mm);
    if reg.EndToWhite_mm.within5mm < minNear
        save(fullfile(OutputDir, [SubjectName '_fibers_registration.mat']), '-struct', 'reg');
        error('Registration gate failed: only %.1f%% of fibre ends within 5 mm of the white surface (need %.0f%%)', ...
              100*reg.EndToWhite_mm.within5mm, 100*minNear);
    end
    reg.Gate = sprintf('passed (>= %.0f%% of ends within 5 mm of white)', 100*minNear);
end
save(fullfile(OutputDir, [SubjectName '_fibers_registration.mat']), '-struct', 'reg');

% save as a Brainstorm fibres file in SCS and import it as such (no further conversion)
FibMat = db_template('fibersmat');
FibMat.Points  = Points;
FibMat.Header  = hdr;
FibMat.Comment = sprintf('fibers_%dPt_%dFib', nPoints, nFibers);
FibMat = fibers_helper('ComputeColor', FibMat);
tmpFile = fullfile(tempdir, [SubjectName '_fibers_scs.mat']);
save(tmpFile, '-struct', 'FibMat', '-v7.3');
clear FibMat Points;
[~, OutputFiles] = import_fibers(iSubject, {tmpFile}, 'BST', nPoints, 'scs');
delete(tmpFile);
end


function q = local_quantile(x, p)
% The p-quantile of x (nearest rank) — no Statistics Toolbox.
x = sort(x(:));
q = x(max(1, min(numel(x), ceil(p * numel(x)))));
end

