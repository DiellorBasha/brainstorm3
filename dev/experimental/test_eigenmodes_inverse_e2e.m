function test_eigenmodes_inverse_e2e
% e2e: drive the eigenmode inverse through BOTH the process wrapper
% (process_eigenmodes_inverse) and the shared core (process_inverse_2018 Compute)
% on a study that has a composed eigenmode head model + noise covariance.
%
% Validates:
%   - the cortex kernel-only results node (Phi*M~) for whitening ON and OFF,
%     and that on/off differ (the whitener has an effect);
%   - the EigenModeKernel scratchpad is stripped from the saved cortex node;
%   - for IMPORTED (non-raw) data, the coefficients matrix node is produced;
%     for RAW data, coefficients are correctly skipped (no matrix node, no error).
% All nodes created here are deleted afterwards so the DB is left unchanged.
% Skips cleanly only if no eigenmode head model + noise cov study exists.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

sStudies = bst_get('ProtocolStudies');
if isempty(sStudies) || ~isfield(sStudies,'Study') || isempty(sStudies.Study)
    disp('SKIP: no protocol loaded.'); return;
end

% Find a study with eigenmode head model + noise cov + a data file.
% Prefer a study that has an IMPORTED (non-raw) recording so the coefficients
% node is exercised; fall back to a raw-only study (cortex path still runs).
iStudyTarget = []; iDataFile = []; isImported = false;
for iS = 1:numel(sStudies.Study)
    s = sStudies.Study(iS);
    if isempty(s.iHeadModel) || s.iHeadModel < 1 || numel(s.HeadModel) < s.iHeadModel; continue; end
    try
        hm = in_bst_headmodel(s.HeadModel(s.iHeadModel).FileName, 0, 'isEigenmode');
    catch; continue; end
    if ~(isfield(hm,'isEigenmode') && ~isempty(hm.isEigenmode) && hm.isEigenmode==1); continue; end
    if ~(isfield(s,'NoiseCov') && ~isempty(s.NoiseCov) && ~isempty(s.NoiseCov(1).FileName)); continue; end
    if isempty(s.Data); continue; end
    % Find an imported (non-raw) recording first; else take any (raw)
    iImp = []; iAny = [];
    for kD = 1:numel(s.Data)
        if isempty(iAny); iAny = kD; end
        if ~isfield(s.Data(kD),'DataType') || ~strcmpi(s.Data(kD).DataType, 'raw')
            iImp = kD; break;
        end
    end
    if ~isempty(iImp)
        iStudyTarget = iS; iDataFile = iImp; isImported = true; break;     % best: imported
    elseif isempty(iStudyTarget)
        iStudyTarget = iS; iDataFile = iAny; isImported = false;            % fallback: raw
    end
end
if isempty(iStudyTarget)
    disp('SKIP: no study with composed eigenmode head model + noise cov + data.'); return;
end
fprintf('Using study %d, data %d, imported=%d\n', iStudyTarget, iDataFile, isImported);

created = {};   % accumulate every node we create, to delete at the end

% ===== Wrapper path: process_eigenmodes_inverse('Run', ...) =====
sInput = struct('iStudy', iStudyTarget, 'iItem', iDataFile, ...
                'FileName', sStudies.Study(iStudyTarget).Data(iDataFile).FileName, 'Comment', 'eigtest');
sProcess = struct('options', struct( ...
    'method',     struct('Value','dspm'), ...
    'prior',      struct('Value','log'), ...
    'snr',        struct('Value',{{3,'',1}}), ...
    'nmodes',     struct('Value',{{0,'',0}}), ...
    'whiten',     struct('Value',1), ...
    'outputtype', struct('Value','both')));
OutFiles = process_eigenmodes_inverse('Run', sProcess, sInput);
created = [created, OutFiles];
assert(~isempty(OutFiles), 'Wrapper must produce output files.');
assert(any(~cellfun(@isempty, regexp(OutFiles, 'results_', 'once'))), 'Wrapper: expected a cortex results node.');
if isImported
    assert(any(~cellfun(@isempty, regexp(OutFiles, 'matrix_', 'once'))), 'Wrapper: expected a coefficient matrix node (imported data).');
end

% ===== Direct Compute() path: whitening ON (default) =====
OPTIONS = struct();
OPTIONS.InverseMethod    = 'eigenmode';
OPTIONS.InverseMeasure   = 'dspm';
OPTIONS.EigenmodePrior   = 'log';
OPTIONS.SnrFixed         = 3;
OPTIONS.SnrMethod        = 'fixed';
OPTIONS.nModes           = 0;
OPTIONS.NoiseMethod      = 'reg';
OPTIONS.NoiseReg         = 0.1;
OPTIONS.SourceOrient     = {'fixed'};
OPTIONS.ComputeKernel    = 1;
OPTIONS.DisplayMessages  = 0;
OPTIONS.DataTypes        = [];
OPTIONS.EigenmodeWhiten  = 1;
OPTIONS.SaveCoefficients = 1;

[outFiles, errMsg] = process_inverse_2018('Compute', iStudyTarget, iDataFile, OPTIONS);
created = [created, outFiles];
assert(isempty(errMsg), sprintf('Compute eigenmode (whiten on) failed: %s', errMsg));
iResNode = find(~cellfun(@isempty, regexp(outFiles, 'results_.*\.mat$', 'once')), 1);
assert(~isempty(iResNode), 'No cortex results node was created.');
ResCheck = in_bst_results(outFiles{iResNode}, 0, 'ImagingKernel', 'nComponents');
assert(~isempty(ResCheck.ImagingKernel), 'Cortex node must carry an ImagingKernel.');
assert(all(isfinite(ResCheck.ImagingKernel(:))), 'Cortex kernel must be finite.');
assert(ResCheck.nComponents == 1, 'Eigenmode cortex node nComponents must be 1.');
% EigenModeKernel scratchpad must not be persisted into the cortex node
Rfull = load(file_fullpath(outFiles{iResNode}));
assert(~isfield(Rfull, 'EigenModeKernel'), 'EigenModeKernel must be stripped from the saved cortex node.');
KWhiten = ResCheck.ImagingKernel;
% Coefficients node only for imported data
if isImported
    assert(any(~cellfun(@isempty, regexp(outFiles, 'matrix_eigencoeffs.*\.mat$', 'once'))), 'Coefficients node expected (imported).');
else
    assert(~any(~cellfun(@isempty, regexp(outFiles, 'matrix_eigencoeffs.*\.mat$', 'once'))), 'Coefficients node must be skipped for raw data.');
end

% ===== Direct Compute() path: whitening OFF (pure sensor-space solve) =====
OPTIONS.EigenmodeWhiten  = 0;
OPTIONS.SaveCoefficients = 0;
[outFiles2, errMsg2] = process_inverse_2018('Compute', iStudyTarget, iDataFile, OPTIONS);
created = [created, outFiles2];
assert(isempty(errMsg2), sprintf('Compute eigenmode (whiten off) failed: %s', errMsg2));
iRes2 = find(~cellfun(@isempty, regexp(outFiles2, 'results_.*\.mat$', 'once')), 1);
assert(~isempty(iRes2), 'Whiten-off path produced no cortex node.');
Res2 = in_bst_results(outFiles2{iRes2}, 0, 'ImagingKernel');
assert(all(isfinite(Res2.ImagingKernel(:))), 'Whiten-off kernel must be finite.');
assert(isequal(size(Res2.ImagingKernel), size(KWhiten)), 'Whiten on/off kernels must share shape.');
assert(max(abs(Res2.ImagingKernel(:) - KWhiten(:))) > 0, 'Whiten on/off should differ (whitener has an effect).');

% ===== Clean up every node we created =====
for k = 1:numel(created)
    f = file_fullpath(created{k});
    if exist(f, 'file'); delete(f); end
end
db_reload_studies(iStudyTarget);

disp('ALL TESTS PASSED');
end
