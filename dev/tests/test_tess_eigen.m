function tests = test_tess_eigen
% Tests for tess_eigen: per-hemisphere eigenbasis as an eigen_ DB node, with
% find-or-create of the underlying operator_ node.
tests = functiontests(localfunctions);
end

% ---------------------------------------------------------------------------
function setupOnce(tc)
    if ~brainstorm('status'); brainstorm nogui; end
    [isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
    assert(isOk, 'nxr-compute required for tess_eigen tests: %s', errMsg);
    tc.TestData.SurfaceFile = local_find_cortex(20484);
    % Back up the parent surface file so any accidental mutation is restored.
    src = file_fullpath(tc.TestData.SurfaceFile);
    bak = [src '.test_eigen.bak'];
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

% ---------------------------------------------------------------------------
% Cleanup helpers (mirror test_tess_operators): delete a node .mat, strip its
% entry from the subject's Surface child lists, write back, reload.

function local_delete_operator(opFileName)
    local_delete_node(opFileName, 'Operator');
end

function local_delete_eigen(eigFileName)
    local_delete_node(eigFileName, 'Eigen');
end

function local_delete_node(fileName, listField)
    if isempty(fileName)
        return
    end
    FileFull = file_fullpath(fileName);
    if exist(FileFull, 'file')
        file_delete(FileFull, 1);
    end
    ProtocolSubjects = bst_get('ProtocolSubjects');
    changed = false;
    if ~isempty(ProtocolSubjects.DefaultSubject)
        [ProtocolSubjects.DefaultSubject, c] = local_strip(ProtocolSubjects.DefaultSubject, fileName, listField);
        changed = changed || c;
    end
    for iSubj = 1:numel(ProtocolSubjects.Subject)
        [ProtocolSubjects.Subject(iSubj), c] = local_strip(ProtocolSubjects.Subject(iSubj), fileName, listField);
        changed = changed || c;
    end
    if changed
        bst_set('ProtocolSubjects', ProtocolSubjects);
        db_save();
        db_reload_database('current');
    end
end

function [sSubject, changed] = local_strip(sSubject, fileName, listField)
    changed = false;
    for iSurf = 1:numel(sSubject.Surface)
        if isfield(sSubject.Surface(iSurf), listField) && ...
           ~isempty(sSubject.Surface(iSurf).(listField))
            entries = sSubject.Surface(iSurf).(listField);
            keep = ~strcmp(file_short(fileName), {entries.FileName});
            if ~all(keep)
                sSubject.Surface(iSurf).(listField) = entries(keep);
                changed = true;
            end
        end
    end
end

% Snapshot the FileNames in a surface child list (Operator/Eigen) for the cortex.
function names = local_list_names(SurfaceFile, listField)
    names = {};
    [sSubj, ~] = bst_get('SurfaceFile', SurfaceFile);
    iSurf = find(strcmpi({sSubj.Surface.FileName}, file_short(SurfaceFile)), 1);
    if isfield(sSubj.Surface(iSurf), listField)
        L = sSubj.Surface(iSurf).(listField);
        for k = 1:numel(L)
            names{end+1} = L(k).FileName; %#ok<AGROW>
        end
    end
end

% Count operator entries of a given Variant under the cortex (loads files).
function [n, files] = local_count_operator_variant(SurfaceFile, Variant)
    n = 0; files = {};
    [sSubj, ~] = bst_get('SurfaceFile', SurfaceFile);
    iSurf = find(strcmpi({sSubj.Surface.FileName}, file_short(SurfaceFile)), 1);
    if ~isfield(sSubj.Surface(iSurf), 'Operator')
        return;
    end
    ops = sSubj.Surface(iSurf).Operator;
    for k = 1:numel(ops)
        try
            S = load(file_fullpath(ops(k).FileName), 'Variant');
            if isfield(S,'Variant') && strcmpi(S.Variant, Variant)
                n = n + 1;
                files{end+1} = ops(k).FileName; %#ok<AGROW>
            end
        catch
        end
    end
end

% ---------------------------------------------------------------------------
function test_lbo_shapes_and_ortho(tc)
% Laplace-Beltrami: 1x2 Phi/Lambda, Phi{hh} is [nVh x K], B-orthonormal.
    SurfaceFile = tc.TestData.SurfaceFile;
    K = 12;
    EM = tess_eigen(SurfaceFile, 'Laplace-Beltrami', 'K',K, 'NoSave',1, 'ForceRecompute',1);
    cleanup = onCleanup(@() local_cleanup_operators(SurfaceFile, 'Laplace-Beltrami')); %#ok<NASGU>

    verifyEqual(tc, EM.Variant, 'Laplace-Beltrami');
    verifyEqual(tc, size(EM.Phi),    [1 2], 'Phi should be 1x2');
    verifyEqual(tc, size(EM.Lambda), [1 2], 'Lambda should be 1x2');

    Op = local_load_operator(EM.OperatorFile);
    for hh = 1:2
        nVh = numel(EM.GlobalVertices{hh});
        verifyEqual(tc, size(EM.Phi{hh}),    [nVh K], sprintf('Hemi %d: Phi not nVh x K', hh));
        verifyEqual(tc, size(EM.Lambda{hh}), [K 1],   sprintf('Hemi %d: Lambda not K x 1', hh));
        B = Op.Mass{hh};
        res = norm(full(EM.Phi{hh}' * (B * EM.Phi{hh}) - eye(K)));
        fprintf('[test] LBO hemi %d B-orthonormality residual = %.3e\n', hh, res);
        verifyLessThan(tc, res, 1e-6, sprintf('Hemi %d: LBO Phi not B-orthonormal (%g)', hh, res));
        lam = EM.Lambda{hh};
        verifyTrue(tc, all(diff(lam) >= -1e-9), sprintf('Hemi %d: LBO Lambda not ascending', hh));
    end
end

% ---------------------------------------------------------------------------
function test_connection_shapes_and_ortho(tc)
% Connection Laplacian: 1x2 Phi/Lambda, Phi{hh} is [nVh x K], B-orthonormal.
    SurfaceFile = tc.TestData.SurfaceFile;
    K = 12;
    EM = tess_eigen(SurfaceFile, 'Connection Laplacian', 'K',K, 'NoSave',1, 'ForceRecompute',1);
    cleanup = onCleanup(@() local_cleanup_operators(SurfaceFile, 'Connection Laplacian')); %#ok<NASGU>

    verifyEqual(tc, EM.Variant, 'Connection Laplacian');
    verifyEqual(tc, EM.Provenance.Ortho, 'Rayleigh-Ritz', 'Connection Ortho provenance');
    verifyEqual(tc, size(EM.Phi),    [1 2], 'Phi should be 1x2');
    verifyEqual(tc, size(EM.Lambda), [1 2], 'Lambda should be 1x2');

    Op = local_load_operator(EM.OperatorFile);
    for hh = 1:2
        nVh = numel(EM.GlobalVertices{hh});
        verifyEqual(tc, size(EM.Phi{hh}),    [nVh K], sprintf('Hemi %d: Phi not nVh x K', hh));
        verifyEqual(tc, size(EM.Lambda{hh}), [K 1],   sprintf('Hemi %d: Lambda not K x 1', hh));
        B = Op.Mass{hh};
        res = norm(full(EM.Phi{hh}' * (B * EM.Phi{hh}) - eye(K)));
        fprintf('[test] Connection hemi %d B-orthonormality residual = %.3e\n', hh, res);
        verifyLessThan(tc, res, 1e-6, sprintf('Hemi %d: Connection Phi not B-orthonormal (%g)', hh, res));
        lam = EM.Lambda{hh};
        verifyTrue(tc, all(diff(lam) >= -1e-9), sprintf('Hemi %d: Connection Lambda not ascending', hh));
        verifyLessThan(tc, max(abs(imag(lam))), 1e-6, sprintf('Hemi %d: Connection Lambda not real-ish', hh));
    end
end

% ---------------------------------------------------------------------------
function test_dirac_shapes_and_ortho(tc)
% Dirac: 1x2 Phi/Lambda, Phi{hh} is [4nVh x K], B-orthonormal (Rayleigh-Ritz),
% eigenvalues form near-4-fold clusters.
    SurfaceFile = tc.TestData.SurfaceFile;
    K = 12;
    EM = tess_eigen(SurfaceFile, 'Dirac', 'K',K, 'Tau',0.5, 'NoSave',1, 'ForceRecompute',1);
    cleanup = onCleanup(@() local_cleanup_operators(SurfaceFile, 'Dirac')); %#ok<NASGU>

    verifyEqual(tc, EM.Variant, 'Dirac');
    verifyEqual(tc, EM.Provenance.Ortho, 'Rayleigh-Ritz', 'Dirac Ortho provenance');
    verifyEqual(tc, size(EM.Phi),    [1 2], 'Phi should be 1x2');
    verifyEqual(tc, size(EM.Lambda), [1 2], 'Lambda should be 1x2');

    Op = local_load_operator(EM.OperatorFile);
    for hh = 1:2
        nVh = numel(EM.GlobalVertices{hh});
        verifyEqual(tc, size(EM.Phi{hh}),    [4*nVh K], sprintf('Hemi %d: Dirac Phi not 4nVh x K', hh));
        verifyEqual(tc, size(EM.Lambda{hh}), [K 1],     sprintf('Hemi %d: Dirac Lambda not K x 1', hh));
        B = Op.Mass{hh};
        res = norm(full(EM.Phi{hh}' * (B * EM.Phi{hh}) - eye(K)));
        fprintf('[test] Dirac hemi %d B-orthonormality residual = %.3e\n', hh, res);
        verifyLessThan(tc, res, 1e-6, sprintf('Hemi %d: Dirac Phi not B-orthonormal (%g)', hh, res));
        lam = EM.Lambda{hh};
        verifyTrue(tc, all(diff(lam) >= -1e-9), sprintf('Hemi %d: Dirac Lambda not ascending', hh));
        % Loose multiplet check: the within-quadruplet spread should be much
        % smaller than the gap between consecutive quadruplets (first 2 groups).
        if K >= 8
            g1 = lam(1:4); g2 = lam(5:8);
            spread1 = max(g1) - min(g1);
            gap = min(g2) - max(g1);
            fprintf('[test] Dirac hemi %d group1 spread=%.3e, gap to group2=%.3e\n', hh, spread1, gap);
            verifyGreaterThanOrEqual(tc, gap, -1e-9, sprintf('Hemi %d: groups not ordered', hh));
        end
    end
end

% ---------------------------------------------------------------------------
function test_find_or_create_operator(tc)
% Starting from a surface with NO operator node of the variant, the first
% tess_eigen call must CREATE the operator node; a second call must REUSE it.
    SurfaceFile = tc.TestData.SurfaceFile;
    Variant = 'Laplace-Beltrami';
    K = 12;

    % Make sure no operator of this variant exists to begin with.
    [~, preFiles] = local_count_operator_variant(SurfaceFile, Variant);
    for k = 1:numel(preFiles)
        local_delete_operator(preFiles{k});
    end
    [n0, ~] = local_count_operator_variant(SurfaceFile, Variant);
    verifyEqual(tc, n0, 0, 'Precondition: no operator of this variant should exist');

    % First call: must create exactly one operator node of the variant.
    EM1 = tess_eigen(SurfaceFile, Variant, 'K',K, 'NoSave',1, 'ForceRecompute',1);
    [n1, files1] = local_count_operator_variant(SurfaceFile, Variant);
    createdOp = files1{1};
    cleanup = onCleanup(@() local_delete_operator(createdOp)); %#ok<NASGU>
    verifyEqual(tc, n1, 1, 'First tess_eigen should create one operator node');
    verifyEqual(tc, file_short(EM1.OperatorFile), file_short(createdOp), ...
        'EigenMat.OperatorFile should reference the created operator');

    % Second call: must REUSE the same operator node (no new operator entry).
    EM2 = tess_eigen(SurfaceFile, Variant, 'K',K, 'NoSave',1, 'ForceRecompute',1);
    [n2, files2] = local_count_operator_variant(SurfaceFile, Variant);
    verifyEqual(tc, n2, 1, 'Second tess_eigen should REUSE, not create a new operator node');
    verifyEqual(tc, file_short(EM2.OperatorFile), file_short(createdOp), ...
        'Second call should reuse the same operator file');

    fprintf('[test] find-or-create: created 1 operator, reused on 2nd call (n1=%d, n2=%d)\n', n1, n2);
end

% ---------------------------------------------------------------------------
function test_save_creates_eigen_node(tc)
% Default save=true should create an eigen_ node registered in the DB; both the
% eigen file and the referenced operator file must resolve. Cleanup removes both.
    SurfaceFile = tc.TestData.SurfaceFile;
    K = 12;

    preEig = local_list_names(SurfaceFile, 'Eigen');

    EM = tess_eigen(SurfaceFile, 'Dirac', 'K',K, 'Tau',0.5);

    postEig = local_list_names(SurfaceFile, 'Eigen');
    verifyFalse(tc, isempty(postEig), 'Surface.Eigen should be non-empty after save');
    newEig = setdiff(postEig, preEig);
    verifyFalse(tc, isempty(newEig), 'No new eigen entry found');
    createdEig = newEig{1};

    % Register cleanup once both file names are known.
    createdOp = file_short(EM.OperatorFile);
    cleanup = onCleanup(@() local_cleanup_save(createdEig, createdOp)); %#ok<NASGU>

    % eigen file resolves via bst_get
    [sEig, ~, ~, iEig] = bst_get('EigenFile', createdEig);
    verifyFalse(tc, isempty(sEig), 'bst_get(EigenFile) should resolve the new node');
    verifyGreaterThan(tc, iEig, 0, 'Resolved eigen index should be > 0');
    verifyTrue(tc, logical(exist(file_fullpath(createdEig),'file')), 'eigen_*.mat should exist on disk');

    % operator file resolves via bst_get
    [sOp, ~, ~, iOp] = bst_get('OperatorFile', EM.OperatorFile);
    verifyFalse(tc, isempty(sOp), 'bst_get(OperatorFile) should resolve EigenMat.OperatorFile');
    verifyGreaterThan(tc, iOp, 0, 'Resolved operator index should be > 0');

    SavedMat = load(file_fullpath(createdEig));
    verifyEqual(tc, SavedMat.Variant, 'Dirac');
    verifyEqual(tc, size(SavedMat.Phi), [1 2], 'Saved Phi should be 1x2');
    verifyEqual(tc, SavedMat.K, K, 'Saved K mismatch');
    verifyFalse(tc, isempty(SavedMat.ParentSurface), 'Saved ParentSurface should not be empty');
end

% ---------------------------------------------------------------------------
% local utilities

function Op = local_load_operator(OperatorFile)
    [sOp, ~, ~, ~] = bst_get('OperatorFile', OperatorFile);
    assert(~isempty(sOp), 'Operator file did not resolve: %s', OperatorFile);
    Op = load(file_fullpath(OperatorFile));
end

% Delete all operator nodes of a given variant under the surface.
function local_cleanup_operators(SurfaceFile, Variant)
    [~, files] = local_count_operator_variant(SurfaceFile, Variant);
    for k = 1:numel(files)
        local_delete_operator(files{k});
    end
end

% Cleanup for the save test: delete eigen first, then its operator.
function local_cleanup_save(eigFile, opFile)
    local_delete_eigen(eigFile);
    local_delete_operator(opFile);
end
