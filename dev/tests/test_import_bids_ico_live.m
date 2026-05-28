function results = test_import_bids_ico_live(varargin)
% TEST_IMPORT_BIDS_ICO_LIVE: Confirm process_import_bids DEFAULTS produce an ico5 manifold cortex.
% Runs the real importer on OMEGA sub-0002 (local) with DEFAULT options (no method set ->
% GetDescription default icosphere/ico5), checks the cortex, then deletes the disposable protocol.
%
% USAGE: results = test_import_bids_ico_live()
%        results = test_import_bids_ico_live('BidsDir', '/path/to/omega-tutorial')

Def.BidsDir      = '/Users/diellorbasha/workspace/library/datasets/omega-tutorial';
Def.ProtocolName = 'BidsIcoDefaultCheck';
Def.SubjTag      = 'sub-0002';
OPT = Def;
for i = 1:2:numel(varargin), OPT.(varargin{i}) = varargin{i+1}; end
assert(exist(OPT.BidsDir, 'dir') == 7, 'BIDS dataset folder not found: %s', OPT.BidsDir);

if ~brainstorm('status')
    brainstorm nogui
end

% Delete first: makes reruns idempotent and cleans up any protocol left behind by a
% previous run that aborted mid-way (an assertion failure skips the end-of-run delete).
gui_brainstorm('DeleteProtocol', OPT.ProtocolName);
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

% Find the imported subject's cortex using the proven idiom from test_omega_icosphere_sourcemap:
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

% 2-manifold check
try
    [~, ~, isManifold] = tess_manifold(TessMat.Vertices, TessMat.Faces);
catch
    isManifold = NaN;
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

% Save the report before teardown so a failure leaves diagnostics behind (non-fatal).
try
    bst_report('Save', []);
catch
    % non-fatal
end

% Tear down the disposable protocol
gui_brainstorm('DeleteProtocol', OPT.ProtocolName);

if results.pass
    fprintf('ALL TESTS PASSED: test_import_bids_ico_live\n');
else
    error('test_import_bids_ico_live FAILED: cortex did not match ico5/2-component/manifold (vertices=%d, components=%d, manifold=%d).', ...
        nVertices, nComp, isManifold);
end
end
