function tests = test_tess_manifold
% Tests for tess_manifold: assemble the nxr geometry backbone as a manifold_ DB node.
tests = functiontests(localfunctions);
end

% ---------------------------------------------------------------------------
function SurfaceFile = local_cortex()
    if ~brainstorm('status'); brainstorm nogui; end
    SurfaceFile = local_find_cortex(20484);
end

function s = local_find_cortex(nVertTarget)
    % Single source of truth for the canonical test cortex (never a hand-built mesh).
    s = bst_canonical_cortex(nVertTarget);
end

% Clean up a manifold node created during a test:
%   - delete the .mat file from disk
%   - remove its entry from the subject's Surface.Manifold list
%   - write back via bst_set + db_save + db_reload_database('current')
function local_delete_manifold(manifestFileName)
    if isempty(manifestFileName)
        return
    end
    ManifoldFileFull = file_fullpath(manifestFileName);
    % 1. Delete file from disk
    if exist(ManifoldFileFull, 'file')
        file_delete(ManifoldFileFull, 1);
    end
    % 2. Remove DB entry
    ProtocolSubjects = bst_get('ProtocolSubjects');
    % Search default subject
    changed = false;
    if ~isempty(ProtocolSubjects.DefaultSubject)
        for iSurf = 1:numel(ProtocolSubjects.DefaultSubject.Surface)
            if isfield(ProtocolSubjects.DefaultSubject.Surface(iSurf),'Manifold') && ...
               ~isempty(ProtocolSubjects.DefaultSubject.Surface(iSurf).Manifold)
                entries = ProtocolSubjects.DefaultSubject.Surface(iSurf).Manifold;
                keep = ~strcmp(file_short(manifestFileName), {entries.FileName});
                ProtocolSubjects.DefaultSubject.Surface(iSurf).Manifold = entries(keep);
                changed = true;
            end
        end
    end
    % Search all subjects
    for iSubj = 1:numel(ProtocolSubjects.Subject)
        for iSurf = 1:numel(ProtocolSubjects.Subject(iSubj).Surface)
            if isfield(ProtocolSubjects.Subject(iSubj).Surface(iSurf),'Manifold') && ...
               ~isempty(ProtocolSubjects.Subject(iSubj).Surface(iSurf).Manifold)
                entries = ProtocolSubjects.Subject(iSubj).Surface(iSurf).Manifold;
                keep = ~strcmp(file_short(manifestFileName), {entries.FileName});
                ProtocolSubjects.Subject(iSubj).Surface(iSurf).Manifold = entries(keep);
                changed = true;
            end
        end
    end
    if changed
        bst_set('ProtocolSubjects', ProtocolSubjects);
        db_save();
        db_reload_database('current');
    end
end

% ---------------------------------------------------------------------------
function test_nosave_returns_struct_1x2(tc)
% tess_manifold with NoSave=1, ForceRecompute=1 returns a valid ManifoldMat
% with 1x2 Topology/Embedded/Intrinsic/Extrinsic/Gauge arrays.
    SurfaceFile = local_cortex();

    ManifoldMat = tess_manifold(SurfaceFile, 'NoSave', 1, 'ForceRecompute', 1);

    verifyTrue(tc, isstruct(ManifoldMat), 'ManifoldMat should be a struct');
    for f = {'Topology','Embedded','Intrinsic','Extrinsic','Gauge'}
        verifyTrue(tc, isfield(ManifoldMat, f{1}), sprintf('Missing field: %s', f{1}));
        verifyEqual(tc, size(ManifoldMat.(f{1})), [1 2], ...
            sprintf('%s should be 1x2', f{1}));
    end

    % Hemisphere tags
    verifyEqual(tc, ManifoldMat.Embedded(1).Hemisphere, 'L');
    verifyEqual(tc, ManifoldMat.Embedded(2).Hemisphere, 'R');

    % Provenance
    verifyEqual(tc, ManifoldMat.Provenance.Backend, 'nxr');
    verifyEqual(tc, ManifoldMat.Provenance.Gauge, 'trivial');

    % Gauge type
    verifyEqual(tc, lower(ManifoldMat.Gauge(1).type), 'trivial');
end

% ---------------------------------------------------------------------------
function test_embedded_grid_orthonormal(tc)
% Embedded(hh).vertex.grid should be an nVh x 3 complex matrix whose real
% and imaginary parts are unit-length and mutually orthogonal per vertex.
    SurfaceFile = local_cortex();

    ManifoldMat = tess_manifold(SurfaceFile, 'NoSave', 1, 'ForceRecompute', 1);

    tol = 1e-4;
    for hh = 1:2
        grid = ManifoldMat.Embedded(hh).vertex.grid;  % nVh x 3 complex
        e1   = real(grid);
        e2   = imag(grid);
        % e1 unit-norm
        n1 = sqrt(sum(e1.^2, 2));
        verifyLessThan(tc, max(abs(n1 - 1)), tol, ...
            sprintf('Hemi %d: e1 not unit-length (max dev = %g)', hh, max(abs(n1-1))));
        % e2 unit-norm
        n2 = sqrt(sum(e2.^2, 2));
        verifyLessThan(tc, max(abs(n2 - 1)), tol, ...
            sprintf('Hemi %d: e2 not unit-length (max dev = %g)', hh, max(abs(n2-1))));
        % e1 perp e2
        dot12 = abs(sum(e1 .* e2, 2));
        verifyLessThan(tc, max(dot12), tol, ...
            sprintf('Hemi %d: e1 and e2 not orthogonal (max dot = %g)', hh, max(dot12)));
    end
end

% ---------------------------------------------------------------------------
function test_gauss_bonnet_index_sum(tc)
% For each hemisphere (topological sphere), the sum of singularity indices
% in the trivial gauge should equal 2 (Poincare-Hopf / Gauss-Bonnet).
    SurfaceFile = local_cortex();

    ManifoldMat = tess_manifold(SurfaceFile, 'NoSave', 1, 'ForceRecompute', 1);

    for hh = 1:2
        G = ManifoldMat.Gauge(hh);
        verifyTrue(tc, isfield(G,'singularity') && isfield(G.singularity,'indices'), ...
            sprintf('Hemi %d: Gauge.singularity.indices missing', hh));
        idx_sum = sum(G.singularity.indices);
        verifyEqual(tc, idx_sum, 2, ...
            sprintf('Hemi %d: Gauss-Bonnet index sum = %g (expected 2)', hh, idx_sum));
    end
end

% ---------------------------------------------------------------------------
function test_save_creates_manifold_node(tc)
% tess_manifold with default save=true (ForceRecompute=1) should:
%   1. Create a manifold_ node registered in the DB.
%   2. The cortex's .Manifold list becomes non-empty.
%   3. bst_get('ManifoldFile', <thatFile>) resolves it.
% Cleanup removes the node.
    SurfaceFile = local_cortex();

    % Ensure we start clean (remove any pre-existing manifold entries for this surface)
    [sSubjPre, iSubjPre] = bst_get('SurfaceFile', SurfaceFile);
    iSurfPre = find(strcmpi({sSubjPre.Surface.FileName}, file_short(SurfaceFile)), 1);
    preExisting = {};
    if isfield(sSubjPre.Surface(iSurfPre),'Manifold')
        for k = 1:numel(sSubjPre.Surface(iSurfPre).Manifold)
            preExisting{end+1} = sSubjPre.Surface(iSurfPre).Manifold(k).FileName; %#ok<AGROW>
        end
    end

    % Run with default save
    ManifoldMat = tess_manifold(SurfaceFile, 'ForceRecompute', 1);

    % Capture the newly created file for cleanup
    createdFile = '';
    cleanup = onCleanup(@() local_delete_manifold(createdFile));

    % Query the DB for the surface's Manifold list
    [sSubj, iSubj] = bst_get('SurfaceFile', SurfaceFile);
    iSurf = find(strcmpi({sSubj.Surface.FileName}, file_short(SurfaceFile)), 1);
    verifyFalse(tc, isempty(iSurf), 'Surface not found in subject after tess_manifold');

    mList = sSubj.Surface(iSurf).Manifold;
    verifyFalse(tc, isempty(mList), 'Surface.Manifold should be non-empty after tess_manifold');

    % Find the newly added entry (not in preExisting)
    newEntries = {};
    for k = 1:numel(mList)
        if ~any(strcmp(mList(k).FileName, preExisting))
            newEntries{end+1} = mList(k).FileName; %#ok<AGROW>
        end
    end
    verifyFalse(tc, isempty(newEntries), 'No new manifold entry found after tess_manifold');
    createdFile = newEntries{1};

    % Update the cleanup handle with the actual file name
    % (re-assign since onCleanup captured the variable by value — use a nested func workaround)
    cleanup2 = onCleanup(@() local_delete_manifold(createdFile));  %#ok<NASGU>
    clear cleanup;

    % 2. Verify bst_get('ManifoldFile',...) resolves
    [sRes, iRes, iSurfRes, iMRes] = bst_get('ManifoldFile', createdFile);
    verifyFalse(tc, isempty(sRes), 'bst_get(ManifoldFile) should resolve the new node');
    verifyEqual(tc, iRes, iSubj, 'Resolved subject index mismatch');
    verifyEqual(tc, iSurfRes, iSurf, 'Resolved surface index mismatch');
    verifyGreaterThan(tc, iMRes, 0, 'Resolved manifold index should be > 0');

    % 3. Verify the file exists on disk
    verifyTrue(tc, logical(exist(file_fullpath(createdFile),'file')), ...
        'manifold_*.mat should exist on disk');

    % 4. Verify the saved content has 1x2 arrays
    SavedMat = load(file_fullpath(createdFile));
    for f = {'Topology','Embedded','Intrinsic','Extrinsic','Gauge'}
        verifyTrue(tc, isfield(SavedMat, f{1}), ...
            sprintf('Saved manifoldmat missing field: %s', f{1}));
        verifyEqual(tc, size(SavedMat.(f{1})), [1 2], ...
            sprintf('Saved %s should be 1x2', f{1}));
    end
    % ParentSurface is set
    verifyFalse(tc, isempty(SavedMat.ParentSurface), ...
        'Saved manifoldmat.ParentSurface should not be empty');
    % Cleanup happens via onCleanup(cleanup2)
end

% ---------------------------------------------------------------------------
function test_skip_when_already_registered(tc)
% Without ForceRecompute, if a manifold node is already registered,
% tess_manifold should return an empty struct (skip signal).
    SurfaceFile = local_cortex();

    % First, create a node
    tess_manifold(SurfaceFile, 'ForceRecompute', 1);

    % Capture the new node for cleanup
    [sSubj, ~] = bst_get('SurfaceFile', SurfaceFile);
    iSurf = find(strcmpi({sSubj.Surface.FileName}, file_short(SurfaceFile)), 1);
    mList = sSubj.Surface(iSurf).Manifold;
    createdFile = '';
    if ~isempty(mList)
        createdFile = mList(end).FileName;
    end
    cleanup = onCleanup(@() local_delete_manifold(createdFile));  %#ok<NASGU>

    % Now call without ForceRecompute — should return empty (cached)
    result = tess_manifold(SurfaceFile);
    verifyTrue(tc, isempty(fieldnames(result)) || isempty(result), ...
        'Without ForceRecompute, tess_manifold should return empty struct when node exists');
end
