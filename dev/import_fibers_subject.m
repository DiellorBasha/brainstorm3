function import_fibers_subject(BidsDir, OutputDir, SubjectLabel, Module, varargin)
% IMPORT_FIBERS_SUBJECT  Import DWI tractography (.trk) into an existing per-subject
%                        Brainstorm protocol using Brainstorm's own import_fibers.
%
% Container worker for the nsp `brainstorm-fibers` pathway, selected by the
% brainstorm-pipeline entrypoint via BST_PIPELINE=import_fibers_subject and called
% with the standard positional contract:
%
%   import_fibers_subject(BidsDir, OutputDir, SubjectLabel, Module, ...
%                         'BstDir', D, 'BstDbDir', DB, 'NVertices', N)
%
% BidsDir / Module are accepted for contract parity but unused: this worker
% augments an already-computed protocol rather than importing BIDS.
%
% All neuroimaging operations use established Brainstorm functions:
%   import_protocol   - load the exported per-subject protocol into the DB
%   bst_get('Subject')- resolve the subject index
%   import_fibers     - read/interpolate/transform the .trk into the subject
%                       (trk_read -> trk_interp -> cs_convert -> ComputeColor -> save)
%   export_protocol   - write the augmented protocol back out
%
% Fiber-import parameters come from environment variables (set by the nsp
% template's apptainer --env), keeping the entrypoint contract intact:
%   NSP_PROTOCOL_ZIP : exported per-subject protocol .zip to augment   (required)
%   NSP_TRK          : TrackVis .trk streamlines to import             (required)
%   NSP_CS           : coordinate system of the .trk {scs|mni|world|mri} (default 'world')
%   NSP_NPOINTS      : points per streamline after resampling          (default 100)
%
% OUTPUT written to OutputDir:
%   <Subject>_brainstorm.zip   augmented protocol (existing anatomy + imported fibers)
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
csEnv       = lower(getenv('NSP_CS'));
nPointsEnv  = getenv('NSP_NPOINTS');
if isempty(csEnv);      csEnv = 'world'; end
if isempty(nPointsEnv); nPoints = 100; else; nPoints = str2double(nPointsEnv); end

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
fprintf('CS       : %s | nPoints: %d\n', csEnv, nPoints);

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
    % ===== Load the existing exported protocol (Brainstorm) =====
    import_protocol(protocolZip);
    iProtocol = bst_get('iProtocol');
    fprintf('Loaded protocol #%d from %s\n', iProtocol, protocolZip);

    % Resolve subject index (Brainstorm).
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
    assert(~isempty(sSubject), 'Subject %s not found in protocol', SubjectName);
    fprintf('Subject #%d: %s (%d surfaces)\n', iSubject, sSubject.Name, numel(sSubject.Surface));

    % ===== Import fibers (Brainstorm import_fibers) =====
    % import_fibers handles trk_read, trk_interp to nPoints, mm->m, cs_convert to
    % SCS using the subject MRI, ComputeColor, save, and db_add_surface.
    fprintf('Importing fibers via import_fibers (CS=%s, nPoints=%d)...\n', csEnv, nPoints);
    [iNewFibers, OutputFiles, nFibers] = import_fibers(iSubject, {trkFile}, 'TRK', nPoints, csEnv); %#ok<ASGLU>
    fprintf('Imported %d fibers.\n', nFibers);
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
