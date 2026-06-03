function info = bench_fixtures(anat, nModesPerHemi)
% BENCH_FIXTURES: Idempotently ensure an anatomy's protocol has cortex eigenmodes
% and a composed eigenmode head model, ready for bst_benchmark_inverse.
%
% USAGE:  info = bench_fixtures(anat, nModesPerHemi)
%   anat          : struct with .protocol, .subject (from bench_config)
%   nModesPerHemi : eigenmodes to compute per hemisphere (e.g. 1000)
% RETURNS info: .iStudy .baseHmFile .ncFile .chFile .eigHmFile
%
% NOTE: The default low-resolution cortex meshes in tutorial protocols are
% non-manifold. bench_fixtures uses repair=1 in process_eigenmodes to fix
% this automatically. When repair changes the vertex count, bench_fixtures
% also recomputes the base surface head model (overlapping spheres) so the
% eigenmode leadfield composition succeeds.

% Select the protocol
iProt = bst_get('Protocol', anat.protocol);
if isempty(iProt); error('bench_fixtures: protocol %s not found.', anat.protocol); end
gui_brainstorm('SetCurrentProtocol', iProt);

% ── Step 1: find a study with noise cov + data + at least one surface HM ──
sStudies = bst_get('ProtocolStudies');
T = find_candidate_study_(sStudies);
if isempty(T)
    error('bench_fixtures: no study with surface HM + noise cov + data in %s.', anat.protocol);
end
s        = sStudies.Study(T.iStudy);
ncFile   = s.NoiseCov(1).FileName;
chFile   = bst_get('ChannelFileForStudy', s.FileName);
dataFile = s.Data(1).FileName;

% ── Step 2: ensure cortex eigenmodes ──────────────────────────────────────
% Use ANY non-eigenmode surface HM to find the cortex surface.
hmAny = in_bst_headmodel(s.HeadModel(T.iAny).FileName, 0, 'SurfaceFile');
surfFile = hmAny.SurfaceFile;

[EigExist, hasEig] = in_tess_eigenmodes(surfFile);
if ~hasEig
    bst_process('CallProcess', 'process_eigenmodes', [], [], ...
        'subjectname', anat.subject, ...
        'surfacetype', 'cortex', ...
        'nmodes',      {nModesPerHemi, '', 0}, ...
        'masstype',    'barycentric', ...
        'removedc',    1, ...
        'repair',      1, ...
        'overwrite',   0);
    [EigExist, hasEig] = in_tess_eigenmodes(surfFile);
    if ~hasEig
        error('bench_fixtures: process_eigenmodes failed to store eigenmodes on %s.', surfFile);
    end
end
nVertEig = size(EigExist.Vectors, 1);

% ── Step 3: pick / ensure a base HM with matching vertex count ────────────
% Scan all non-eigenmode surface HMs in the study for the one whose vertex
% count matches the (possibly repaired) eigenmode surface.
sStudies = bst_get('ProtocolStudies');
s        = sStudies.Study(T.iStudy);
baseHmFile = pick_matching_base_hm_(s, nVertEig);

if isempty(baseHmFile)
    % No matching HM exists yet — recompute with overlapping spheres.
    fprintf(['bench_fixtures: no base surface HM with %d vertices. ' ...
             'Recomputing with overlapping spheres (repair changed vertex count)...\n'], nVertEig);
    bst_process('CallProcess', 'process_headmodel', dataFile, [], ...
        'Comment',     '', ...
        'sourcespace', 1, ...   % Cortex surface
        'meg',         {3, {'<none>','Single sphere','Overlapping spheres','OpenMEEG BEM','DUNEuro FEM'}}, ...
        'eeg',         {1, {'<none>','3-shell sphere','OpenMEEG BEM','DUNEuro FEM'}}, ...
        'ecog',        {1, {'<none>','OpenMEEG BEM','DUNEuro FEM'}}, ...
        'seeg',        {1, {'<none>','OpenMEEG BEM','DUNEuro FEM'}}, ...
        'nirs',        {1, {'<none>','Import from MCXlab'}});
    sStudies   = bst_get('ProtocolStudies');
    s          = sStudies.Study(T.iStudy);
    baseHmFile = pick_matching_base_hm_(s, nVertEig);
    if isempty(baseHmFile)
        error('bench_fixtures: recomputed head model still has wrong vertex count in %s.', anat.protocol);
    end
end

% ── Step 4: ensure eigenmode leadfield ───────────────────────────────────
eigHmFile = find_eig_hm_(s);
if isempty(eigHmFile)
    bst_process('CallProcess', 'process_eigenmode_leadfield', dataFile, [], ...
        'nmodes', {0, '', 0});   % 0 = all available
    sStudies  = bst_get('ProtocolStudies');
    s         = sStudies.Study(T.iStudy);
    eigHmFile = find_eig_hm_(s);
end
if isempty(eigHmFile)
    error('bench_fixtures: failed to create eigenmode head model in %s.', anat.protocol);
end

info = struct('iStudy',     T.iStudy, ...
              'baseHmFile', baseHmFile, ...
              'ncFile',     ncFile, ...
              'chFile',     chFile, ...
              'eigHmFile',  eigHmFile, ...
              'dataFile',   dataFile, ...
              'surfaceFile', surfFile);
end

%% ===== HELPERS =====

function T = find_candidate_study_(sStudies)
% Return T.iStudy / T.iAny for the first study that has >=1 surface HM,
% a noise covariance file, and >=1 data file.  T.iAny is any non-eigenmode
% surface HM index (used only to locate the cortex surface).
T = [];
for iS = 1:numel(sStudies.Study)
    s = sStudies.Study(iS);
    if isempty(s.HeadModel) || isempty(s.NoiseCov) || isempty(s.NoiseCov(1).FileName); continue; end
    if isempty(s.Data); continue; end
    iAny = [];
    for ih = 1:numel(s.HeadModel)
        try hm = in_bst_headmodel(s.HeadModel(ih).FileName, 0); catch; continue; end
        isEig = isfield(hm,'isEigenmode') && ~isempty(hm.isEigenmode) && hm.isEigenmode;
        if ~isEig && strcmpi(hm.HeadModelType,'surface'); iAny = ih; break; end
    end
    if ~isempty(iAny); T = struct('iStudy',iS,'iAny',iAny); return; end
end
end

function baseHmFile = pick_matching_base_hm_(s, nVertEig)
% Return the FileName of the non-eigenmode surface HM whose vertex count
% equals nVertEig, or '' if none found.  Prefers the LAST matching HM
% (most recently computed, in case repair forced a recompute).
baseHmFile = '';
for ih = numel(s.HeadModel):-1:1   % reverse: newest first
    try hm = in_bst_headmodel(s.HeadModel(ih).FileName, 0); catch; continue; end
    isEig = isfield(hm,'isEigenmode') && ~isempty(hm.isEigenmode) && hm.isEigenmode;
    if ~isEig && strcmpi(hm.HeadModelType,'surface')
        nV = size(hm.Gain,2)/3;
        if nV == nVertEig; baseHmFile = s.HeadModel(ih).FileName; return; end
    end
end
end

function eigHmFile = find_eig_hm_(s)
% Return FileName of the first eigenmode head model in study s, or ''.
eigHmFile = '';
for ih = 1:numel(s.HeadModel)
    try hm = in_bst_headmodel(s.HeadModel(ih).FileName, 0); catch; continue; end
    if isfield(hm,'isEigenmode') && ~isempty(hm.isEigenmode) && hm.isEigenmode
        eigHmFile = s.HeadModel(ih).FileName; return;
    end
end
end
