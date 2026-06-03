function info = bench_fixtures(anat, nModesPerHemi)
% BENCH_FIXTURES: Idempotently ensure an anatomy's protocol has cortex eigenmodes
% (>= nModesPerHemi per hemisphere) and a composed eigenmode head model containing
% all those modes, ready for bst_benchmark_inverse.
%
% USAGE:  info = bench_fixtures(anat, nModesPerHemi)
%   anat          : struct with .protocol, .subject (from bench_config)
%   nModesPerHemi : eigenmodes to compute per hemisphere (e.g. 1000)
% RETURNS info: .iStudy .baseHmFile .ncFile .chFile .eigHmFile .dataFile .surfaceFile
%
% NOTE: Default low-res tutorial cortex meshes are non-manifold; repair=1 fixes this
% (and may change vertex count, in which case the base head model is recomputed so the
% leadfield composition stays vertex-consistent). Existing eigenmodes / eigenmode-HM are
% REUSED only when they already contain enough modes; otherwise they are recomputed.

iProt = bst_get('Protocol', anat.protocol);
if isempty(iProt); error('bench_fixtures: protocol %s not found.', anat.protocol); end
gui_brainstorm('SetCurrentProtocol', iProt);

% Step 1: candidate study (surface HM + noise cov + data)
sStudies = bst_get('ProtocolStudies');
T = find_candidate_study_(sStudies);
if isempty(T); error('bench_fixtures: no study with surface HM + noise cov + data in %s.', anat.protocol); end
s        = sStudies.Study(T.iStudy);
ncFile   = s.NoiseCov(1).FileName;
chFile   = bst_get('ChannelFileForStudy', s.FileName);
dataFile = s.Data(1).FileName;

% Step 2: ensure cortex eigenmodes with ENOUGH modes per hemisphere
hmAny = in_bst_headmodel(s.HeadModel(T.iAny).FileName, 0, 'SurfaceFile');
surfFile = hmAny.SurfaceFile;
if ~eig_has_enough_(surfFile, nModesPerHemi)
    bst_process('CallProcess', 'process_eigenmodes', [], [], ...
        'subjectname', anat.subject, 'surfacetype', 'cortex', ...
        'nmodes', {nModesPerHemi, '', 0}, 'masstype', 'barycentric', ...
        'removedc', 1, 'repair', 1, 'overwrite', 1);
    if ~eig_has_enough_(surfFile, nModesPerHemi)
        error('bench_fixtures: process_eigenmodes did not store >= %d modes/hemisphere on %s.', nModesPerHemi, surfFile);
    end
end
[EigExist, ~] = in_tess_eigenmodes(surfFile);
nVertEig  = size(EigExist.Vectors, 1);
nModesEig = size(EigExist.Vectors, 2);

% Step 3: base HM with matching vertex count (recompute if repair changed it)
sStudies = bst_get('ProtocolStudies'); s = sStudies.Study(T.iStudy);
baseHmFile = pick_matching_base_hm_(s, nVertEig);
if isempty(baseHmFile)
    fprintf('bench_fixtures: no base HM with %d vertices; recomputing overlapping spheres...\n', nVertEig);
    bst_process('CallProcess', 'process_headmodel', dataFile, [], ...
        'Comment', '', 'sourcespace', 1, ...
        'meg',  {3, {'<none>','Single sphere','Overlapping spheres','OpenMEEG BEM','DUNEuro FEM'}}, ...
        'eeg',  {1, {'<none>','3-shell sphere','OpenMEEG BEM','DUNEuro FEM'}}, ...
        'ecog', {1, {'<none>','OpenMEEG BEM','DUNEuro FEM'}}, ...
        'seeg', {1, {'<none>','OpenMEEG BEM','DUNEuro FEM'}}, ...
        'nirs', {1, {'<none>','Import from MCXlab'}});
    sStudies = bst_get('ProtocolStudies'); s = sStudies.Study(T.iStudy);
    baseHmFile = pick_matching_base_hm_(s, nVertEig);
    if isempty(baseHmFile); error('bench_fixtures: recomputed HM still has wrong vertex count in %s.', anat.protocol); end
end

% Step 4: ensure a composed eigenmode HM containing ~all stored modes
[eigHmFile, nEigHmModes] = find_best_eig_hm_(s);
if isempty(eigHmFile) || nEigHmModes < nModesEig - 4
    bst_process('CallProcess', 'process_eigenmode_leadfield', dataFile, [], 'nmodes', {0, '', 0});
    sStudies = bst_get('ProtocolStudies'); s = sStudies.Study(T.iStudy);
    [eigHmFile, nEigHmModes] = find_best_eig_hm_(s);
end
if isempty(eigHmFile) || nEigHmModes < nModesEig - 4
    error('bench_fixtures: failed to compose eigenmode HM with enough modes (got %d, need ~%d) in %s.', nEigHmModes, nModesEig, anat.protocol);
end

info = struct('iStudy', T.iStudy, 'baseHmFile', baseHmFile, 'ncFile', ncFile, ...
              'chFile', chFile, 'eigHmFile', eigHmFile, 'dataFile', dataFile, ...
              'surfaceFile', surfFile);
end

%% ===== HELPERS =====
function T = find_candidate_study_(sStudies)
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

function tf = eig_has_enough_(surfFile, nModesPerHemi)
% True if eigenmodes exist AND every connected component has >= nModesPerHemi-2 modes.
tf = false;
[Eig, hasEig] = in_tess_eigenmodes(surfFile);
if ~hasEig || ~isfield(Eig,'Component') || isempty(Eig.Component); return; end
comps  = unique(Eig.Component);
counts = arrayfun(@(c) sum(Eig.Component==c), comps);
tf = all(counts >= nModesPerHemi - 2);
end

function baseHmFile = pick_matching_base_hm_(s, nVertEig)
baseHmFile = '';
for ih = numel(s.HeadModel):-1:1
    try hm = in_bst_headmodel(s.HeadModel(ih).FileName, 0); catch; continue; end
    isEig = isfield(hm,'isEigenmode') && ~isempty(hm.isEigenmode) && hm.isEigenmode;
    if ~isEig && strcmpi(hm.HeadModelType,'surface') && size(hm.Gain,2)/3 == nVertEig
        baseHmFile = s.HeadModel(ih).FileName; return;
    end
end
end

function [eigHmFile, nModes] = find_best_eig_hm_(s)
% Return the eigenmode HM with the MOST modes (columns) and that count.
eigHmFile = ''; nModes = 0;
for ih = 1:numel(s.HeadModel)
    try hm = in_bst_headmodel(s.HeadModel(ih).FileName, 0); catch; continue; end
    if isfield(hm,'isEigenmode') && ~isempty(hm.isEigenmode) && hm.isEigenmode
        nc = size(hm.Gain,2);
        if nc > nModes; nModes = nc; eigHmFile = s.HeadModel(ih).FileName; end
    end
end
end
