function test_eigenmode_canonical_consistency_e2e
% No-regression safety net: canonical-default reconstruction equals explicit-ModeIndices,
% and the canonical-path leadfield Gain equals the existing composed HM.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end
if isempty(bst_get('Protocol','TutorialAuditory')); disp('SKIP: TutorialAuditory missing.'); return; end
iProt=bst_get('Protocol','TutorialAuditory'); gui_brainstorm('SetCurrentProtocol',iProt);
sStudies=bst_get('ProtocolStudies'); T=[];
for iS=1:numel(sStudies.Study)
    s=sStudies.Study(iS); if isempty(s.HeadModel)||isempty(s.NoiseCov)||isempty(s.NoiseCov(1).FileName); continue; end
    iBase=[]; iEig=[];
    for ih=1:numel(s.HeadModel)
        try hm=in_bst_headmodel(s.HeadModel(ih).FileName,0); catch; continue; end
        isE=isfield(hm,'isEigenmode')&&hm.isEigenmode;
        if isE && isempty(iEig); iEig=ih; elseif ~isE && strcmpi(hm.HeadModelType,'surface') && isempty(iBase); iBase=ih; end
    end
    if ~isempty(iBase)&&~isempty(iEig); T=struct('iS',iS,'iBase',iBase,'iEig',iEig); break; end
end
if isempty(T); disp('SKIP: need base + eigenmode head models.'); return; end
s=sStudies.Study(T.iS);
oldEigHM=in_bst_headmodel(s.HeadModel(T.iEig).FileName,0);

% Derive GoodChannel from Gain rows that have a valid forward model (non-NaN).
GainTmp = double(oldEigHM.Gain);
GoodChannel = ~any(isnan(GainTmp), 2);

% Reconstruct two ways: canonical default vs explicit ModeIndices (== canonical Order).
[Inv,err]=bst_inverse_eigenmodes(s.HeadModel(T.iEig).FileName, s.NoiseCov(1).FileName, ...
    bst_get('ChannelFileForStudy',s.FileName), GoodChannel, 'Method','mne','Prior','log','SNR',3);
assert(isempty(err), 'inverse failed: %s', err);
Kidx = bst_eigenmode_reconstruct(oldEigHM.SurfaceFile, Inv.ImagingKernel, oldEigHM.ModeIndices);
Kdef = bst_eigenmode_reconstruct(oldEigHM.SurfaceFile, Inv.ImagingKernel);   % canonical default
rel = norm(Kidx - Kdef,'fro')/norm(Kidx,'fro');
fprintf('reconstruct canonical-default vs explicit rel.diff = %.3e\n', rel);
assert(rel < 1e-9, 'canonical-default reconstruction must equal explicit ModeIndices');

% Leadfield consistency (current leadfield already sorts globally, so should already match).
% Restrict to valid (non-NaN) rows — auxiliary channels have no forward model.
baseHM=in_bst_headmodel(s.HeadModel(T.iBase).FileName,0);
Eig=in_tess_eigenmodes(baseHM.SurfaceFile);
newEigHM=bst_eigenmode_leadfield(baseHM, Eig);
nC=min(size(newEigHM.Gain,2), size(oldEigHM.Gain,2));
oldG=double(oldEigHM.Gain(:,1:nC)); newG=double(newEigHM.Gain(:,1:nC));
validRows = ~any(isnan(oldG),2) & ~any(isnan(newG),2);
relL=norm(newG(validRows,:)-oldG(validRows,:),'fro')/norm(oldG(validRows,:),'fro');
fprintf('leadfield Gain rel.diff = %.3e\n', relL);
assert(relL < 1e-9, 'canonical leadfield Gain diverged from existing composed HM');
disp('ALL TESTS PASSED');
end
