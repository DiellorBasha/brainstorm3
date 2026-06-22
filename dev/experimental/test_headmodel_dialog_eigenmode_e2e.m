function test_headmodel_dialog_eigenmode_e2e
% Smoke: ComputeHeadModel with the "Cortex surface harmonics" source space
% (passed as sMethod, bypassing the GUI) produces exactly ONE harmonic head
% model node (Gain=[nCh x K], isEigenmode, Comment '... | harmonic') and no
% extra base node. Skips cleanly without a suitable study.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

sStudies = bst_get('ProtocolStudies');
if isempty(sStudies) || ~isfield(sStudies,'Study') || isempty(sStudies.Study)
    disp('SKIP: no protocol loaded.'); return;
end
% Find a study: MEG channel + subject cortex with precomputed eigenmodes
iStudyTarget = [];
for iS = 1:numel(sStudies.Study)
    s = sStudies.Study(iS);
    if isempty(s.Channel) || isempty(s.Channel.FileName); continue; end
    if ~any(strcmpi(s.Channel.DisplayableSensorTypes, 'MEG')); continue; end
    sSubj = bst_get('Subject', s.BrainStormSubject);
    if isempty(sSubj) || isempty(sSubj.iCortex); continue; end
    [~, isEig] = in_tess_eigenmodes(sSubj.Surface(sSubj.iCortex).FileName);
    if isEig; iStudyTarget = iS; break; end
end
if isempty(iStudyTarget)
    disp('SKIP: no study with MEG channel + cortex eigenmodes.'); return;
end

sBefore   = bst_get('Study', iStudyTarget);
nHmBefore = numel(sBefore.HeadModel);

% Build a harmonic sMethod (bypass the GUI)
sMethod = struct('Comment', 'Overlapping spheres | harmonic', ...
    'HeadModelType', 'surface', 'SourceCompression', 'eigenmode', 'nModes', 0, ...
    'MEGMethod', 'os_meg', 'EEGMethod', '', 'ECOGMethod', '', 'SEEGMethod', '', ...
    'NIRSMethod', '', 'SaveFile', 1, 'Interactive', 0);
OutFiles = panel_headmodel('ComputeHeadModel', iStudyTarget, sMethod);
assert(~isempty(OutFiles), 'Harmonic head model must be produced.');

% Exactly one new node, and it is the harmonic one
sAfter = bst_get('Study', iStudyTarget);
assert(numel(sAfter.HeadModel) == nHmBefore + 1, 'Exactly one new head model node expected (no base node).');
HM = in_bst_headmodel(sAfter.HeadModel(sAfter.iHeadModel).FileName, 0);
assert(isfield(HM,'isEigenmode') && HM.isEigenmode==1, 'Active node must be the harmonic head model.');
assert(size(HM.Gain,2) == HM.nModes, 'Gain must be [nCh x nModes].');
assert(~isempty(HM.Eigenvalues), 'Eigenvalues must be stored.');
assert(~isempty(regexp(HM.Comment, '\| harmonic$', 'once')), 'Comment must end with "| harmonic".');
disp('ALL TESTS PASSED');
end
