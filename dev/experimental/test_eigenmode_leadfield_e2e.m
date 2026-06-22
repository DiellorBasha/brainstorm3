function test_eigenmode_leadfield_e2e
% Smoke: build a composed eigenmode head model from a real study and verify it
% is a valid headmodel node with Gain=[nCh x K] and carried metadata.
% Skips cleanly without a suitable protocol.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

sProtocol = bst_get('ProtocolStudies');
if isempty(sProtocol) || ~isfield(sProtocol,'Study') || isempty(sProtocol.Study)
    disp('SKIP: no protocol loaded.'); return;
end
iStudyTarget = [];
for iS = 1:numel(sProtocol.Study)
    s = sProtocol.Study(iS);
    if isfield(s,'iHeadModel') && ~isempty(s.iHeadModel) && s.iHeadModel >= 1 ...
            && numel(s.HeadModel) >= s.iHeadModel
        try
            hm = in_bst_headmodel(s.HeadModel(s.iHeadModel).FileName, 0, 'HeadModelType','SurfaceFile');
            if strcmpi(hm.HeadModelType,'surface')
                [~, isEig] = in_tess_eigenmodes(hm.SurfaceFile);
                if isEig; iStudyTarget = iS; break; end
            end
        catch; end
    end
end
if isempty(iStudyTarget)
    disp('SKIP: no study with surface head model + eigenmodes.'); return;
end

OutFiles = process_eigenmode_leadfield('Run', ...
    struct('options', struct('nmodes', struct('Value',{{0,'',0}}))), ...
    struct('iStudy', iStudyTarget, 'FileName', '', 'Comment', 'test'));
% Process returns the (pass-through) input files; the head model is a DB side-effect.
sStudyAfter = bst_get('Study', iStudyTarget);
HM = in_bst_headmodel(sStudyAfter.HeadModel(sStudyAfter.iHeadModel).FileName, 0);
assert(isfield(HM,'isEigenmode') && HM.isEigenmode==1, ...
    'Active head model after the process must be an eigenmode leadfield.');
assert(size(HM.Gain,2) == HM.nModes, 'Gain columns must equal nModes.');
assert(~isempty(HM.Eigenvalues), 'Eigenvalues must be stored.');
disp('ALL TESTS PASSED');
end
