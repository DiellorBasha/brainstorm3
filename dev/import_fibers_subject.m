function import_fibers_subject(BidsDir, OutputDir, SubjectLabel, Module, varargin)
% IMPORT_FIBERS_SUBJECT  Add DWI tractography to an existing Brainstorm protocol
%                        and compute multi-atlas structural connectomes.
%
% Container worker for the nsp `brainstorm-fibers` pathway. It is selected by the
% brainstorm-pipeline entrypoint via BST_PIPELINE=import_fibers_subject and is
% therefore called with the SAME positional contract as every other worker:
%
%   import_fibers_subject(BidsDir, OutputDir, SubjectLabel, Module, ...
%                         'BstDir', D, 'BstDbDir', DB, 'NVertices', N)
%
% BidsDir / Module are accepted for contract parity but not used here: this worker
% augments an already-computed per-subject protocol rather than importing BIDS.
%
% Fiber-specific inputs are passed through the environment (set by the nsp
% template's `apptainer --env`), which keeps the fixed entrypoint contract intact:
%
%   NSP_PROTOCOL_ZIP : exported per-subject protocol .zip to augment   (required)
%   NSP_TRK          : TrackVis .trk streamlines to import             (required)
%   NSP_ATLASES      : comma-list of cortical atlases for connectomes
%                      (default 'Desikan-Killiany,Destrieux')
%   NSP_CS           : coordinate system of the .trk {scs|mni|world}   (default 'world')
%   NSP_NPOINTS      : points per streamline after resampling          (default 100)
%
% OUTPUTS written to OutputDir:
%   <Subject>_brainstorm.zip              augmented protocol (anatomy + fibers)
%   <Subject>_connectome_<atlas>.mat      NxN streamline-count matrix + labels
%   <Subject>_fibers_report.txt           summary + alignment check
%
% Authors: Diellor Basha, 2026 (nsp brainstorm-fibers pathway)

% ===== Parse the standard container contract =====
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

SubjectLabel = regexprep(opts.SubjectLabel, '^sub-', '');   % strip prefix if present
SubjectName  = ['sub-' SubjectLabel];

% ===== Fiber-specific inputs from the environment =====
protocolZip = getenv('NSP_PROTOCOL_ZIP');
trkFile     = getenv('NSP_TRK');
atlasEnv    = getenv('NSP_ATLASES');
csEnv       = lower(getenv('NSP_CS'));
nPointsEnv  = getenv('NSP_NPOINTS');
if isempty(atlasEnv);  atlasEnv = 'Desikan-Killiany,Destrieux'; end
if isempty(csEnv);     csEnv    = 'world';   end
if isempty(nPointsEnv); nPoints = 100; else;  nPoints = str2double(nPointsEnv); end
atlasList = strtrim(strsplit(atlasEnv, ','));
atlasList = atlasList(~cellfun(@isempty, atlasList));

assert(~isempty(protocolZip) && exist(protocolZip, 'file') == 2, ...
    'NSP_PROTOCOL_ZIP not found: %s', protocolZip);
assert(~isempty(trkFile) && exist(trkFile, 'file') == 2, ...
    'NSP_TRK not found: %s', trkFile);
assert(ismember(csEnv, {'scs','mni','world'}), 'NSP_CS must be scs|mni|world (got "%s")', csEnv);

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
fprintf('atlases  : %s\n', strjoin(atlasList, ', '));
fprintf('CS       : %s | nPoints: %d\n', csEnv, nPoints);

% ===== Init Brainstorm (headless server) =====
addpath(opts.BstDir);
brainstorm setpath;
bst_user_dir = fullfile(char(java.lang.System.getProperty('user.home')), '.brainstorm');
if exist(bst_user_dir, 'dir') ~= 7; mkdir(bst_user_dir); end
% Point Brainstorm at the throwaway node-local DB before starting the server.
iProtocol = 1; %#ok<NASGU>
ProtocolsListInfo = []; %#ok<NASGU>
if ~brainstorm('status'); brainstorm server; end
bst_set('BrainstormDbDir', opts.BstDbDir);
fprintf('Brainstorm server started. DB = %s\n', opts.BstDbDir);

try
    % ===== Load the existing exported protocol =====
    import_protocol(protocolZip);               % unzips into DbDir + sets current
    iProtocol = bst_get('iProtocol');
    fprintf('Loaded protocol #%d from %s\n', iProtocol, protocolZip);

    % Resolve the subject inside the imported protocol.
    [sSubject, iSubject] = bst_get('Subject', SubjectName);
    if isempty(sSubject)
        % Fall back to the first non-default subject in a single-subject protocol.
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
        'Subject %s has no cortex surface in the protocol', SubjectName);
    fprintf('Subject #%d: %s (%d surfaces)\n', iSubject, sSubject.Name, numel(sSubject.Surface));

    % ===== Import fibers (.trk) =====
    fprintf('Importing fibers (CS=%s, nPoints=%d)...\n', csEnv, nPoints);
    [iNewFibers, fibersFiles, nFibers] = import_fibers(iSubject, {trkFile}, 'TRK', nPoints, csEnv); %#ok<ASGLU>
    if iscell(fibersFiles); fibersFile = fibersFiles{1}; else; fibersFile = fibersFiles; end
    fprintf('Imported %d fibers -> %s\n', nFibers, fibersFile);

    % ===== Alignment sanity check (fibers bbox vs cortex bbox, SCS metres) =====
    FibMat = load(file_fullpath(fibersFile));
    cortexFile = sSubject.Surface(sSubject.iCortex).FileName;
    sCortex = in_tess_bst(cortexFile);
    fibPts = reshape(FibMat.Points, [], 3);
    fbb = [min(fibPts,[],1); max(fibPts,[],1)];
    cbb = [min(sCortex.Vertices,[],1); max(sCortex.Vertices,[],1)];
    overlap = all(fbb(2,:) >= cbb(1,:)) && all(cbb(2,:) >= fbb(1,:));
    fprintf('Cortex bbox (m): [%.3f %.3f %.3f]..[%.3f %.3f %.3f]\n', cbb(1,:), cbb(2,:));
    fprintf('Fibers bbox (m): [%.3f %.3f %.3f]..[%.3f %.3f %.3f]\n', fbb(1,:), fbb(2,:));
    if ~overlap
        warning(['ALIGNMENT WARNING: fiber and cortex bounding boxes are disjoint. ' ...
                 'The .trk CS (%s) may be wrong for this anatomy.'], csEnv);
    else
        fprintf('Alignment OK: fiber/cortex bounding boxes overlap.\n');
    end

    % ===== Multi-atlas structural connectomes =====
    iAtlasAll = [];
    for ia = 1:numel(sCortex.Atlas)
        for wanted = atlasList
            if strcmpi(sCortex.Atlas(ia).Name, wanted{1})
                iAtlasAll(end+1) = ia; %#ok<AGROW>
            end
        end
    end
    if isempty(iAtlasAll)
        warning('None of the requested atlases (%s) found on cortex; available: %s', ...
            strjoin(atlasList, ', '), strjoin({sCortex.Atlas.Name}, ', '));
    end

    for ia = iAtlasAll
        atlasName = sCortex.Atlas(ia).Name;
        scouts = sCortex.Atlas(ia).Scouts;
        nScouts = numel(scouts);
        if nScouts < 2; continue; end
        centroids = zeros(nScouts, 3);
        labels = cell(nScouts, 1);
        for is = 1:nScouts
            centroids(is,:) = mean(sCortex.Vertices(scouts(is).Vertices, :), 1);
            labels{is} = scouts(is).Label;
        end
        connLabel = sprintf('%s_%s', SubjectName, atlasName);
        FibMat = fibers_helper('AssignToScouts', FibMat, connLabel, centroids);

        % Endpoint assignments -> NxN streamline-count connectome (symmetric).
        asg = FibMat.Scouts(end).Assignment;    % nFibers x 2 scout indices
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
                            'Measure', 'streamline_count', 'CS', csEnv);
        save(outMat, '-struct', 'connectome');
        fprintf('Connectome[%s]: %dx%d, %d assigned streamlines -> %s\n', ...
            atlasName, nScouts, nScouts, sum(valid), outMat);
    end

    % Persist the fiber file with scout assignments back into the protocol.
    bst_save(file_fullpath(fibersFile), FibMat, 'v7');
    db_save();

    % ===== Re-export the augmented protocol =====
    exportZip = fullfile(opts.OutputDir, [SubjectName '_brainstorm.zip']);
    if exist(exportZip, 'file') == 2; delete(exportZip); end
    export_protocol(iProtocol, iSubject, exportZip);
    fprintf('Exported augmented protocol -> %s\n', exportZip);

    fprintf('=== DONE: %s ===\n', SubjectName);
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
