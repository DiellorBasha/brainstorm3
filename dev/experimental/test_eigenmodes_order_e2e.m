function test_eigenmodes_order_e2e
% tess_eigenmodes writes a global-sort .Order; in_tess_eigenmodes backfills it. On the
% 2-hemisphere cortex, the first-K canonical modes span BOTH components.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end
if isempty(bst_get('Protocol','TutorialAuditory')); disp('SKIP: TutorialAuditory missing.'); return; end
iProt=bst_get('Protocol','TutorialAuditory'); gui_brainstorm('SetCurrentProtocol',iProt);
sSubj=bst_get('Subject','Subject01'); ctxFile=sSubj.Surface(sSubj.iCortex).FileName;
Surf=in_tess_bst(ctxFile);

Eig = tess_eigenmodes(Surf.Vertices, Surf.Faces, 'nModes', 50);   % ~50/hemi, 2 components
assert(isfield(Eig,'Order') && numel(Eig.Order)==size(Eig.Vectors,2), 'Order not written');
assert(issorted(Eig.Values(Eig.Order)), 'Order does not sort Values globally');
assert(isequal(sort(Eig.Order(:)), (1:numel(Eig.Order))'), 'Order not a permutation');
K = 40; comp = Eig.Component(Eig.Order(1:K));
assert(numel(unique(comp))==2, 'first-K canonical must include both hemispheres');

% Backfill: strip Order, save to a temp surface, reload via in_tess_eigenmodes.
TessMat = in_tess_bst(ctxFile, 0);
EigStrip = rmfield(Eig,'Order'); TessMat.Eigenmodes = EigStrip;
tmpFile = fullfile(tempdir, sprintf('tess_order_test_%d.mat', feature('getpid')));
bst_save(tmpFile, TessMat, 'v7');
EigBack = in_tess_eigenmodes(tmpFile);
assert(isfield(EigBack,'Order') && ~isempty(EigBack.Order), 'Order not backfilled');
assert(issorted(EigBack.Values(EigBack.Order)), 'backfilled Order wrong');
delete(tmpFile);
disp('ALL TESTS PASSED');
end
