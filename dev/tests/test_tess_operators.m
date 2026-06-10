function tests = test_tess_operators
% Tests for tess_operators: per-hemisphere {operator,mass} as an operator_ DB node.
tests = functiontests(localfunctions);
end

% ---------------------------------------------------------------------------
function setupOnce(tc)
    if ~brainstorm('status'); brainstorm nogui; end
    [isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
    assert(isOk, 'nxr-compute required for tess_operators tests: %s', errMsg);
    tc.TestData.SurfaceFile = local_find_cortex(20484);
    % Back up the parent surface file so any accidental mutation is restored.
    src = file_fullpath(tc.TestData.SurfaceFile);
    bak = [src '.test_operators.bak'];
    copyfile(src, bak);
    tc.TestData.SurfaceBackup = bak;
    tc.TestData.SurfaceSrc    = src;
    tc.TestData.cleanup = onCleanup(@() local_restore_surface(src, bak));
end

function s = local_find_cortex(nVertTarget)
    % Single source of truth for the canonical test cortex (never a hand-built mesh).
    s = bst_canonical_cortex(nVertTarget);
end

function local_restore_surface(src, bak)
    if exist(bak, 'file')
        copyfile(bak, src);
        delete(bak);
    end
end

% Clean up an operator node created during a test:
%   - delete the .mat file from disk
%   - remove its entry from the subject's Surface.Operator list
%   - write back via bst_set + db_save + db_reload_database('current')
function local_delete_operator(opFileName)
    if isempty(opFileName)
        return
    end
    OpFileFull = file_fullpath(opFileName);
    if exist(OpFileFull, 'file')
        file_delete(OpFileFull, 1);
    end
    ProtocolSubjects = bst_get('ProtocolSubjects');
    changed = false;
    if ~isempty(ProtocolSubjects.DefaultSubject)
        [ProtocolSubjects.DefaultSubject, c] = local_strip(ProtocolSubjects.DefaultSubject, opFileName);
        changed = changed || c;
    end
    for iSubj = 1:numel(ProtocolSubjects.Subject)
        [ProtocolSubjects.Subject(iSubj), c] = local_strip(ProtocolSubjects.Subject(iSubj), opFileName);
        changed = changed || c;
    end
    if changed
        bst_set('ProtocolSubjects', ProtocolSubjects);
        db_save();
        db_reload_database('current');
    end
end

function [sSubject, changed] = local_strip(sSubject, opFileName)
    changed = false;
    for iSurf = 1:numel(sSubject.Surface)
        if isfield(sSubject.Surface(iSurf),'Operator') && ...
           ~isempty(sSubject.Surface(iSurf).Operator)
            entries = sSubject.Surface(iSurf).Operator;
            keep = ~strcmp(file_short(opFileName), {entries.FileName});
            if ~all(keep)
                sSubject.Surface(iSurf).Operator = entries(keep);
                changed = true;
            end
        end
    end
end

% ---------------------------------------------------------------------------
function test_lbo_shapes(tc)
% Laplace-Beltrami: 1x2 Operator/Mass, real symmetric [nVh x nVh].
    SurfaceFile = tc.TestData.SurfaceFile;
    OM = tess_operators(SurfaceFile, 'Laplace-Beltrami', 'NoSave',1, 'ForceRecompute',1);

    verifyEqual(tc, OM.Variant, 'Laplace-Beltrami');
    verifyEqual(tc, size(OM.Operator), [1 2], 'Operator should be 1x2');
    verifyEqual(tc, size(OM.Mass),     [1 2], 'Mass should be 1x2');

    for hh = 1:2
        A   = OM.Operator{hh};
        B   = OM.Mass{hh};
        nVh = numel(OM.GlobalVertices{hh});
        verifyEqual(tc, size(A), [nVh nVh], sprintf('Hemi %d: A not nVh x nVh', hh));
        verifyEqual(tc, size(B), [nVh nVh], sprintf('Hemi %d: B not nVh x nVh', hh));
        verifyTrue(tc, isreal(A), sprintf('Hemi %d: LBO A should be real', hh));
        symErr = full(max(max(abs(A - A.'))));
        verifyLessThan(tc, symErr, 1e-9, sprintf('Hemi %d: LBO A not symmetric (%g)', hh, symErr));
    end
end

% ---------------------------------------------------------------------------
function test_connection_shapes(tc)
% Connection Laplacian: 1x2 Operator/Mass, Hermitian complex [nVh x nVh].
    SurfaceFile = tc.TestData.SurfaceFile;
    OM = tess_operators(SurfaceFile, 'Connection Laplacian', 'NoSave',1, 'ForceRecompute',1);

    verifyEqual(tc, OM.Variant, 'Connection Laplacian');
    verifyEqual(tc, size(OM.Operator), [1 2], 'Operator should be 1x2');
    verifyEqual(tc, size(OM.Mass),     [1 2], 'Mass should be 1x2');

    for hh = 1:2
        A   = OM.Operator{hh};
        nVh = numel(OM.GlobalVertices{hh});
        verifyEqual(tc, size(A), [nVh nVh], sprintf('Hemi %d: A not nVh x nVh', hh));
        hermErr = full(max(max(abs(A - A'))));
        verifyLessThan(tc, hermErr, 1e-9, sprintf('Hemi %d: connection A not Hermitian (%g)', hh, hermErr));
    end
end

% ---------------------------------------------------------------------------
function test_dirac_shapes(tc)
% Dirac: 1x2 Operator/Mass, A is [4nVh x 4nVh], B == kron(Mg, I4).
    SurfaceFile = tc.TestData.SurfaceFile;
    OM = tess_operators(SurfaceFile, 'Dirac', 'Tau',0.5, 'NoSave',1, 'ForceRecompute',1);

    verifyEqual(tc, OM.Variant, 'Dirac');
    verifyEqual(tc, size(OM.Operator), [1 2], 'Operator should be 1x2');
    verifyEqual(tc, size(OM.Mass),     [1 2], 'Mass should be 1x2');
    verifyEqual(tc, OM.Provenance.Tau, 0.5, 'Provenance.Tau should be 0.5');

    for hh = 1:2
        A   = OM.Operator{hh};
        B   = OM.Mass{hh};
        nVh = numel(OM.GlobalVertices{hh});
        verifyEqual(tc, size(A), [4*nVh 4*nVh], sprintf('Hemi %d: Dirac A not 4nVh x 4nVh', hh));
        verifyEqual(tc, size(B), [4*nVh 4*nVh], sprintf('Hemi %d: Dirac B not 4nVh x 4nVh', hh));
        % B should be block-diagonal kron(Mg, I4): reconstruct Mg from B and compare.
        Mg = B(1:4:end, 1:4:end);
        verifyEqual(tc, full(max(max(abs(B - kron(Mg, speye(4)))))), 0, ...
            sprintf('Hemi %d: Dirac B not kron(Mg, I4)', hh));
    end
end

% ---------------------------------------------------------------------------
function test_eigensolve_smoke(tc)
% Per-variant smoke eigensolve of the (A,B) pencil returns finite values.
    SurfaceFile = tc.TestData.SurfaceFile;

    % LBO: smallest eig ~ 0 (constant mode)
    OM = tess_operators(SurfaceFile, 'Laplace-Beltrami', 'NoSave',1, 'ForceRecompute',1);
    dA = eigs(OM.Operator{1}, OM.Mass{1}, 6, 'smallestabs');
    verifyTrue(tc, all(isfinite(dA)), 'LBO eigs not finite');

    % Connection: Hermitian pencil, eigs finite
    OM = tess_operators(SurfaceFile, 'Connection Laplacian', 'NoSave',1, 'ForceRecompute',1);
    dC = eigs(OM.Operator{1}, OM.Mass{1}, 6, 'smallestabs');
    verifyTrue(tc, all(isfinite(dC)), 'Connection eigs not finite');

    % Dirac: ~4-fold near-zero clusters
    OM = tess_operators(SurfaceFile, 'Dirac', 'Tau',0.5, 'NoSave',1, 'ForceRecompute',1);
    opts = struct('tol',1e-6, 'maxit',1000, 'disp',0);
    dL = eigs(OM.Operator{1}, OM.Mass{1}, 8, 'smallestabs', opts);
    verifyTrue(tc, all(isfinite(real(dL))), 'Dirac eigs not finite');
end

% ---------------------------------------------------------------------------
function test_save_creates_operator_node(tc)
% Default save=true should create an operator_ node registered in the DB,
% make the cortex .Operator list non-empty, and bst_get('OperatorFile',...)
% should resolve it.  Cleanup removes the node.
    SurfaceFile = tc.TestData.SurfaceFile;

    % Snapshot pre-existing operator entries for this surface.
    [sSubjPre, ~] = bst_get('SurfaceFile', SurfaceFile);
    iSurfPre = find(strcmpi({sSubjPre.Surface.FileName}, file_short(SurfaceFile)), 1);
    preExisting = {};
    if isfield(sSubjPre.Surface(iSurfPre),'Operator')
        for k = 1:numel(sSubjPre.Surface(iSurfPre).Operator)
            preExisting{end+1} = sSubjPre.Surface(iSurfPre).Operator(k).FileName; %#ok<AGROW>
        end
    end

    tess_operators(SurfaceFile, 'Dirac', 'Tau',0.5, 'ForceRecompute',1);

    [sSubj, iSubj] = bst_get('SurfaceFile', SurfaceFile);
    iSurf = find(strcmpi({sSubj.Surface.FileName}, file_short(SurfaceFile)), 1);
    opList = sSubj.Surface(iSurf).Operator;
    verifyFalse(tc, isempty(opList), 'Surface.Operator should be non-empty after save');

    newEntries = {};
    for k = 1:numel(opList)
        if ~any(strcmp(opList(k).FileName, preExisting))
            newEntries{end+1} = opList(k).FileName; %#ok<AGROW>
        end
    end
    verifyFalse(tc, isempty(newEntries), 'No new operator entry found');
    createdFile = newEntries{1};

    % Register cleanup only after createdFile is known (anonymous functions
    % capture by value at creation time, so an empty placeholder would never
    % remove the node).
    cleanup = onCleanup(@() local_delete_operator(createdFile));  %#ok<NASGU>

    [sRes, iRes, iSurfRes, iOpRes] = bst_get('OperatorFile', createdFile);
    verifyFalse(tc, isempty(sRes), 'bst_get(OperatorFile) should resolve the new node');
    verifyEqual(tc, iRes, iSubj, 'Resolved subject index mismatch');
    verifyEqual(tc, iSurfRes, iSurf, 'Resolved surface index mismatch');
    verifyGreaterThan(tc, iOpRes, 0, 'Resolved operator index should be > 0');

    verifyTrue(tc, logical(exist(file_fullpath(createdFile),'file')), ...
        'operator_*.mat should exist on disk');

    SavedMat = load(file_fullpath(createdFile));
    verifyEqual(tc, SavedMat.Variant, 'Dirac');
    verifyEqual(tc, size(SavedMat.Operator), [1 2], 'Saved Operator should be 1x2');
    verifyFalse(tc, isempty(SavedMat.ParentSurface), 'Saved ParentSurface should not be empty');
end
