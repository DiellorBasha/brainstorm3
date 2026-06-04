function test_eigenmodes_ensure_e2e
% bst_eigenmodes_ensure returns canonical eigenmodes when present (idempotent, no
% recompute), with a valid canonical Order.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end
if isempty(bst_get('Protocol','TutorialAuditory')); disp('SKIP: TutorialAuditory missing.'); return; end
iProt=bst_get('Protocol','TutorialAuditory'); gui_brainstorm('SetCurrentProtocol',iProt);
sSubj=bst_get('Subject','Subject01'); ctxFile=sSubj.Surface(sSubj.iCortex).FileName;

t0=tic; Eig=bst_eigenmodes_ensure(ctxFile); el=toc(t0);
assert(~isempty(Eig) && isfield(Eig,'Order'), 'ensure did not return canonical eigenmodes');
assert(issorted(Eig.Values(Eig.Order)), 'returned Order invalid');
assert(el < 30, 'ensure must be fast when eigenmodes already exist (no recompute)');
Eig2=bst_eigenmodes_ensure(ctxFile);
assert(size(Eig2.Vectors,2)==size(Eig.Vectors,2), 'ensure not idempotent');
fprintf('ensure: %d modes returned in %.2fs (existing, no recompute)\n', size(Eig.Vectors,2), el);
disp('ALL TESTS PASSED');
end
