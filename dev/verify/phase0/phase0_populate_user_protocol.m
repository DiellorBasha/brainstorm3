% PHASE0_POPULATE_USER_PROTOCOL: BIDS-import OMEGA sub-0002 (clean-branch code,
% default icosphere/ico5 options) into the USER'S existing protocol
% omega-tutorial-cortical-flow on SpikeData-2, loaded ALONE into this isolated
% Brainstorm instance (the rest of the user's brainstorm_db is never touched).
BidsDir      = '/Volumes/SpikeData-2/workspace/library/datasets/omega-tutorial';
ProtocolDir  = '/Volumes/SpikeData-2/workspace/library/datasets/brainstorm_db/omega-tutorial-cortical-flow';
ProtocolName = 'omega-tutorial-cortical-flow';
SubjTag      = 'sub-0002';

assert(exist(BidsDir, 'dir') == 7, 'BIDS dir not found: %s', BidsDir);
assert(exist(fullfile(ProtocolDir, 'data'), 'dir') == 7, 'Protocol data dir not found: %s', ProtocolDir);

% --- Register ONLY this protocol; NEVER reuse a same-named entry that points
%     elsewhere (the bug that shadowed the user's protocol with the isolated
%     Gate-0 copy). Run this script in a FRESH user dir (BST_USERDIR_CLEAN). ---
iProtocol = bst_get('Protocol', ProtocolName);
if isempty(iProtocol)
    sProtocol = db_template('ProtocolInfo');
    sProtocol.Comment  = ProtocolName;
    sProtocol.SUBJECTS = fullfile(ProtocolDir, 'anat');
    sProtocol.STUDIES  = fullfile(ProtocolDir, 'data');
    iProtocol = db_edit_protocol('load', sProtocol);
    assert(iProtocol > 0, 'Failed to load protocol from %s', ProtocolDir);
end
gui_brainstorm('SetCurrentProtocol', iProtocol);
% HARD GUARD: the registered protocol must point at the intended folders
ProtocolInfo = bst_get('ProtocolInfo');
assert(strcmp(ProtocolInfo.SUBJECTS, fullfile(ProtocolDir, 'anat')), ...
    'Registered protocol points at %s, expected %s — refusing to import.', ...
    ProtocolInfo.SUBJECTS, fullfile(ProtocolDir, 'anat'));
db_reload_database(iProtocol);
bst_report('Start');

% Refresh process registry so the new downsamplemethod/icolevel options exist
panel_process_select('ParseProcessFolder', 1);

% --- The real importer, DEFAULT downsampling options (icosphere / ico5) ---
bst_process('CallProcess', 'process_import_bids', [], [], ...
    'bidsdir',      {BidsDir, 'BIDS'}, ...
    'selectsubj',   SubjTag, ...
    'channelalign', 0);

% --- Verify: ico5 cortex, per-hemisphere manifold ---
ProtocolSubjects = bst_get('ProtocolSubjects');
iSubj = [];
for iS = 1:numel(ProtocolSubjects.Subject)
    if ~isempty(strfind(ProtocolSubjects.Subject(iS).Name, SubjTag))
        iSubj = iS; break;
    end
end
assert(~isempty(iSubj), 'Imported subject %s not found', SubjTag);
sSubject = bst_get('Subject', iSubj);
assert(~isempty(sSubject.iCortex), 'No cortex surface for %s', SubjTag);
CortexFile = sSubject.Surface(sSubject.iCortex).FileName;
Cortex = in_tess_bst(CortexFile);
nVert = size(Cortex.Vertices, 1);
assert(nVert == 20484, 'Expected ico5 cortex (20484), got %d', nVert);
[rH, lH] = tess_hemisplit(Cortex);
for hemi = {rH, lH}
    iV = hemi{1};
    iF = all(ismember(Cortex.Faces, iV), 2);
    lut = zeros(nVert, 1); lut(iV) = 1:numel(iV);
    [~, ~, isM] = tess_repair(Cortex.Vertices(iV,:), lut(Cortex.Faces(iF,:)));
    assert(isM, 'Hemisphere fails 2-manifold check');
end

% --- Count linked recordings (raw links) ---
sStudies = bst_get('StudyWithSubject', sSubject.FileName);
nRaw = 0; rawNames = {};
for iSt = 1:numel(sStudies)
    for iD = 1:numel(sStudies(iSt).Data)
        if strcmpi(sStudies(iSt).Data(iD).DataType, 'raw')
            nRaw = nRaw + 1;
            rawNames{end+1} = sStudies(iSt).Data(iD).Comment; %#ok<AGROW>
        end
    end
end

% Persist the database index for this protocol
db_save(1);
fprintf('PHASE0 USER-PROTOCOL IMPORT PASSED: cortex=%d verts, manifold OK, raw recordings linked=%d\n', nVert, nRaw);
for k = 1:numel(rawNames), fprintf('  raw: %s\n', rawNames{k}); end
