function test_eigenmodes_transform_canonical_e2e
% The unregularized eigenmode transform must select modes via the canonical Order, so
% K<1000 spans BOTH hemispheres (previously Vectors(:,1:K) was hemisphere-1 only). We
% verify the canonical selection and exercise the transform's kernel math directly
% (bst_eigenmodes_transform), without the heavier CallProcess pipeline.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end
if isempty(bst_get('Protocol','TutorialAuditory')); disp('SKIP: TutorialAuditory missing.'); return; end
iProt=bst_get('Protocol','TutorialAuditory'); gui_brainstorm('SetCurrentProtocol',iProt);
sSubj=bst_get('Subject','Subject01'); ctxFile=sSubj.Surface(sSubj.iCortex).FileName;
Eig=in_tess_eigenmodes(ctxFile);

% Canonical selection (what the fixed transform uses) vs the old naive first-K slice.
K = 600;
assert(numel(unique(Eig.Component(1:K)))==1, 'precondition: naive first-K is single-hemisphere');
sel = Eig.Order(1:K);
assert(numel(unique(Eig.Component(sel)))==2, 'canonical first-K must span both hemispheres');
assert(max(Eig.Values(sel)) < max(Eig.Values(1:K)), 'canonical selects strictly lower frequencies');

% Exercise the transform kernel math on the canonical Phi (mirrors the process internals,
% which now use Order(1:K)). Use the study's base surface leadfield.
sStudies=bst_get('ProtocolStudies'); s=sStudies.Study(6);
ih=[];
for k=1:numel(s.HeadModel)
    hm=in_bst_headmodel(s.HeadModel(k).FileName,0);
    if (~isfield(hm,'isEigenmode')||~hm.isEigenmode) && strcmpi(hm.HeadModelType,'surface'); ih=k; break; end
end
assert(~isempty(ih), 'no base surface head model found');
HM=in_bst_headmodel(s.HeadModel(ih).FileName,1);   % [nCh x nVert] constrained
iGood=all(isfinite(double(HM.Gain)),2);
Phi=double(Eig.Vectors(:, sel));
[Kernel, Info]=bst_eigenmodes_transform(double(HM.Gain(iGood,:)), Phi);
assert(size(Kernel,1)==K && all(isfinite(Kernel(:))), 'transform kernel must be finite [K x nCh]');
fprintf('canonical transform: K=%d, rank=%d, both hemispheres, maxLam %.4g < naive %.4g\n', ...
    K, Info.Rank, max(Eig.Values(sel)), max(Eig.Values(1:K)));
disp('ALL TESTS PASSED');
end
