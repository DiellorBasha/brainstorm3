function results = test_import_bids_ico(varargin)
% TEST_IMPORT_BIDS_ICO: Confirm process_import_bids DEFAULTS produce an ico5 manifold cortex.
% Runs the real importer on OMEGA sub-0002 (local) with DEFAULT options (no method set ->
% GetDescription default icosphere/ico5) and checks the cortex.
%
% The protocol is deliberately KEPT after the run (Gate 0 GUI inspection). It lives only in
% the harness' isolated Brainstorm user dir; the user's main database is never touched.
%
% USAGE: results = test_import_bids_ico()
%        results = test_import_bids_ico('BidsDir', '/path/to/omega-tutorial')

Def.BidsDir      = '/Volumes/SpikeData-2/workspace/library/datasets/omega-tutorial';
Def.ProtocolName = 'omega-tutorial-cortical-flow';
Def.SubjTag      = 'sub-0002';
OPT = Def;
for i = 1:2:numel(varargin), OPT.(varargin{i}) = varargin{i+1}; end
assert(exist(OPT.BidsDir, 'dir') == 7, 'BIDS dataset folder not found: %s', OPT.BidsDir);

if ~brainstorm('status')
    brainstorm nogui
end

% Delete first: makes reruns idempotent and cleans up any protocol left behind by a
% previous run that aborted mid-way.
gui_brainstorm('DeleteProtocol', OPT.ProtocolName);
% A run that errors out never reaches "brainstorm stop", so the protocol registration
% may never have been persisted to brainstorm.mat while its folders WERE created on
% disk. DeleteProtocol no-ops on such an orphan (bst_get('Protocol',...) is empty) and
% CreateProtocol then refuses because the anat/data folders are not empty. Remove the
% orphaned folder tree directly so reruns stay idempotent.
OrphanDir = bst_fullfile(bst_get('BrainstormDbDir'), OPT.ProtocolName);
if (exist(OrphanDir, 'dir') == 7)
    fprintf('Removing orphaned protocol folder: %s\n', OrphanDir);
    rmdir(OrphanDir, 's');
end
gui_brainstorm('CreateProtocol', OPT.ProtocolName, 0, 0);
bst_report('Start');

% Force-reload the process registry from disk so any stale cached entry for
% process_import_bids (e.g. from a session started before this branch added the
% downsamplemethod/icolevel options) is replaced before CallProcess runs. Without
% this, sProcess.options.downsamplemethod can be missing and Run errors.
panel_process_select('ParseProcessFolder', 1);

% Run the real importer with DEFAULT downsampling (do NOT set downsamplemethod/icolevel —
% this exercises the GetDescription defaults: icosphere / ico5 -> 20484 vertices).
bst_process('CallProcess', 'process_import_bids', [], [], ...
    'bidsdir',      {OPT.BidsDir, 'BIDS'}, ...
    'selectsubj',   OPT.SubjTag, ...
    'channelalign', 0);

% Find the imported subject's cortex using the proven idiom from test_import_fs_ico:
%   bst_get('Subject', iSubject) -> sSubject.iCortex -> sSubject.Surface(iCortex).FileName
ProtocolSubjects = bst_get('ProtocolSubjects');
% The DeleteProtocol above guarantees this protocol holds only the just-imported
% subject, so the first subject with a cortex surface is sub-0002's.
iCortexSubj = [];
for iS = 1:numel(ProtocolSubjects.Subject)
    if ~isempty(ProtocolSubjects.Subject(iS).iCortex) && ...
            ~isempty(ProtocolSubjects.Subject(iS).Surface)
        iCortexSubj = iS;
        break;
    end
end
assert(~isempty(iCortexSubj), 'No subject with a cortex surface was imported.');

% Re-fetch via the single-subject form to get a fully populated struct
sSubject   = bst_get('Subject', iCortexSubj);
CortexFile = sSubject.Surface(sSubject.iCortex).FileName;
TessMat    = in_tess_bst(CortexFile);

nVertices = size(TessMat.Vertices, 1);
nFaces    = size(TessMat.Faces, 1);

% Connected components (two hemispheres => 2)
VertConn  = tess_vertconn(TessMat.Vertices, TessMat.Faces);
nComp     = max(conncomp(graph(VertConn)));

% Per-hemisphere 2-manifold check: split on the Structures atlas (tess_hemisplit),
% never conncomp; tess_repair is the harness-local manifold oracle.
isManifold = true;
try
    [rH, lH] = tess_hemisplit(TessMat);
    for hemi = {rH, lH}
        iV  = hemi{1};
        iF  = all(ismember(TessMat.Faces, iV), 2);
        lut = zeros(nVertices, 1); lut(iV) = 1:numel(iV);
        [~, ~, isM] = tess_repair(TessMat.Vertices(iV,:), lut(TessMat.Faces(iF,:)));
        isManifold = isManifold && logical(isM);
    end
catch ME
    fprintf('Manifold check errored: %s\n', ME.message);
    isManifold = false;
end

results = struct( ...
    'subject',    sSubject.Name, ...
    'cortexFile', CortexFile, ...
    'nVertices',  nVertices, ...
    'nFaces',     nFaces, ...
    'nComp',      nComp, ...
    'isManifold', logical(isManifold));
results.pass = (nVertices == 20484) && (nComp == 2) && isequal(logical(isManifold), true);

fprintf('Cortex: %s\n  vertices=%d (expect 20484), faces=%d, components=%d (expect 2), manifold=%d (expect 1)\n', ...
    sSubject.Name, nVertices, nFaces, nComp, isManifold);

% Save the report so diagnostics are left behind (non-fatal).
try
    bst_report('Save', []);
catch
    % non-fatal
end

% NOTE: the protocol is intentionally NOT deleted here — it is kept for Gate 0 inspection.

if results.pass
    fprintf('ALL TESTS PASSED: test_import_bids_ico\n');
else
    error('test_import_bids_ico FAILED: cortex did not match ico5/2-component/manifold (vertices=%d, components=%d, manifold=%d).', ...
        nVertices, nComp, isManifold);
end
end
