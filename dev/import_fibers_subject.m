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
%   NSP_CS            .trk coordinate system {scs|mni|world|mri} (default 'world')
%   NSP_NPOINTS       points per streamline after resampling     (default 100)
%   NSP_ATLASES       comma-list of cortical atlases for connectomes
%                     (default 'Desikan-Killiany,Destrieux')
%
% OUTPUT (to OutputDir):
%   <Subject>_brainstorm.zip            augmented protocol (anatomy + fibers)
%   <Subject>_connectome_<atlas>.mat    NxN streamline-count matrix + labels
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
fprintf('CS       : %s | nPoints: %d | atlases: %s\n', csEnv, nPoints, strjoin(atlasList, ', '));

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

    % ===== Import fibers (Brainstorm import_fibers) =====
    fprintf('Importing fibers via import_fibers (CS=%s, nPoints=%d)...\n', csEnv, nPoints);
    [iNewFibers, OutputFiles, nFibers] = import_fibers(iSubject, {trkFile}, 'TRK', nPoints, csEnv); %#ok<ASGLU>
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
