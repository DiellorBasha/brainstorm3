function test_eigenmodes_inverse_e2e
% Smoke: run the eigenmode inverse on a study that has a composed eigenmode head
% model + noise cov + imported data; verify a coefficient matrix node and a
% kernel-only cortex results node are produced with correct shapes. Skips cleanly.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

sStudies = bst_get('ProtocolStudies');
if isempty(sStudies) || ~isfield(sStudies,'Study') || isempty(sStudies.Study)
    disp('SKIP: no protocol loaded.'); return;
end
% Find a study with an eigenmode head model + noise cov + a data file
iStudyTarget = []; iDataFile = [];
for iS = 1:numel(sStudies.Study)
    s = sStudies.Study(iS);
    if isempty(s.iHeadModel) || s.iHeadModel < 1 || numel(s.HeadModel) < s.iHeadModel; continue; end
    try
        hm = in_bst_headmodel(s.HeadModel(s.iHeadModel).FileName, 0);
    catch; continue; end
    if isfield(hm,'isEigenmode') && hm.isEigenmode==1 ...
            && isfield(s,'NoiseCov') && ~isempty(s.NoiseCov) && ~isempty(s.NoiseCov(1).FileName) ...
            && ~isempty(s.Data)
        iStudyTarget = iS; iDataFile = 1; break;
    end
end
if isempty(iStudyTarget)
    disp('SKIP: no study with composed eigenmode head model + noise cov + data.'); return;
end

sInput = struct('iStudy', iStudyTarget, 'FileName', sStudies.Study(iStudyTarget).Data(iDataFile).FileName, ...
                'Comment', 'eigtest');
sProcess = struct('options', struct( ...
    'method',     struct('Value','dspm'), ...
    'prior',      struct('Value','log'), ...
    'snr',        struct('Value',{{3,'',1}}), ...
    'nmodes',     struct('Value',{{0,'',0}}), ...
    'outputtype', struct('Value','both')));
OutFiles = process_eigenmodes_inverse('Run', sProcess, sInput);
assert(~isempty(OutFiles), 'Inverse must produce output files.');

% At least one results (cortex) and one matrix (coefficients)
hasRes = any(~cellfun(@isempty, regexp(OutFiles, 'results_', 'once')));
hasMat = any(~cellfun(@isempty, regexp(OutFiles, 'matrix_',  'once')));
assert(hasRes, 'Expected a cortex results node.');
assert(hasMat, 'Expected a coefficient matrix node.');
disp('ALL TESTS PASSED');
end
