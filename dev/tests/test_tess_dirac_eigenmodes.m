function tests = test_tess_dirac_eigenmodes
% Property tests for tess_dirac_eigenmodes (per-hemisphere Dirac eigenbasis).
tests = functiontests(localfunctions);
end

function SurfaceFile = local_cortex()
    if ~brainstorm('status'); brainstorm nogui; end
    SurfaceFile = local_find_cortex(20484);
end

function s = local_find_cortex(nVertTarget)
    s = '';
    P = bst_get('ProtocolSubjects');
    subj = P.Subject;
    if isfield(P,'DefaultSubject') && ~isempty(P.DefaultSubject)
        subj = [P.DefaultSubject, subj];
    end
    for k = 1:numel(subj)
        surfs = subj(k).Surface;
        for i = 1:numel(surfs)
            if strcmpi(surfs(i).SurfaceType, 'Cortex')
                m = load(file_fullpath(surfs(i).FileName), 'Vertices');
                if size(m.Vertices,1) == nVertTarget
                    s = surfs(i).FileName; return;
                end
            end
        end
    end
    error('No %d-vertex cortex found in the loaded protocol.', nVertTarget);
end

function test_store_1x2_shapes(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    K = 40;
    DE = tess_dirac_eigenmodes(SurfaceFile, 'K', K, 'Tau', 0.5, 'ForceRecompute', 1);
    verifyEqual(tc, size(DE), [1 2]);
    T = in_tess_bst(SurfaceFile, 0);
    verifyTrue(tc, isfield(T,'DiracEigen') && isequal(size(T.DiracEigen),[1 2]));
    for hh = 1:2
        nVh = numel(T.DiracEigen(hh).GlobalVertices);
        verifyEqual(tc, size(T.DiracEigen(hh).Vectors), [4*nVh, K]);
        verifyEqual(tc, numel(T.DiracEigen(hh).Values), K);
        verifyEqual(tc, T.DiracEigen(hh).Tau, 0.5);
        verifyEqual(tc, T.DiracEigen(hh).Provenance.Backend, 'nxr');
    end
    verifyEqual(tc, T.DiracEigen(1).Hemisphere, 'L');
    verifyEqual(tc, T.DiracEigen(2).Hemisphere, 'R');
end

function test_B_orthonormal(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    K = 40;
    tess_dirac_eigenmodes(SurfaceFile, 'K', K, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    for hh = 1:2
        Phi = T.DiracEigen(hh).Vectors;
        B   = kron(T.DiracEigen(hh).Mass, speye(4));
        G   = Phi' * (B * Phi);
        verifyLessThan(tc, full(max(abs(G - speye(K)), [], 'all')), 1e-7);
    end
end

function test_values_ascending_nonneg_multiplets(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    K = 40;
    tess_dirac_eigenmodes(SurfaceFile, 'K', K, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    for hh = 1:2
        lam = T.DiracEigen(hh).Values;
        verifyGreaterThanOrEqual(tc, min(lam), -1e-9);
        verifyTrue(tc, issorted(lam));
        % 4-fold quaternionic multiplets: find the first non-null cluster
        % The null space may be 0, 3, or 4 dimensional; skip zeros and verify
        % the next 4 non-trivial eigenvalues form a tight cluster (relative spread < 1e-3).
        iFirst = find(lam > 1e-6, 1);
        if ~isempty(iFirst) && iFirst + 3 <= K
            cluster4 = lam(iFirst:iFirst+3);
            verifyLessThan(tc, cluster4(4) - cluster4(1), 1e-3 * (1 + abs(cluster4(4))));
        end
    end
end

function test_eigenpair_residual(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    K = 40;
    tess_dirac_eigenmodes(SurfaceFile, 'K', K, 'Tau', 0.5, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    [isOk] = bst_plugin('Install', 'nxr-compute'); verifyTrue(tc, logical(isOk));
    Vtx = double(T.Vertices); Fcs = double(T.Faces); nVtot = size(Vtx,1);
    for hh = 1:2
        vH = T.DiracEigen(hh).GlobalVertices;
        isV = false(nVtot,1); isV(vH) = true;
        fMask = all(isV(Fcs), 2);
        map = zeros(nVtot,1); map(vH) = 1:numel(vH);
        Vloc = Vtx(vH,:); Floc = map(Fcs(fMask,:));
        h = nxr_compute('create', Vloc, Floc);
        L = nxr_compute('operators', h, 'dirac', T.DiracEigen(hh).Tau);
        nxr_compute('destroy', h);
        B  = kron(T.DiracEigen(hh).Mass, speye(4));
        Phi = T.DiracEigen(hh).Vectors; mu = T.DiracEigen(hh).Values;
        Res = L*Phi - B*Phi*spdiags(mu, 0, numel(mu), numel(mu));
        relResid = max(sqrt(sum(Res.^2, 1))) / normest(L);
        verifyLessThan(tc, relResid, 1e-6);   % every stored column is a genuine eigenpair
    end
end

function test_cache_return_no_recompute(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    K = 40;
    tess_dirac_eigenmodes(SurfaceFile, 'K', K, 'Tau', 0.5, 'ForceRecompute', 1);
    TF = load(TessFile);
    TF.DiracEigen(1).Provenance.ComputeDate = 'SENTINEL';
    bst_save(TessFile, TF, 'v7');

    tess_dirac_eigenmodes(SurfaceFile, 'K', K, 'Tau', 0.5);   % must cache-return
    T = in_tess_bst(SurfaceFile, 0);
    verifyEqual(tc, T.DiracEigen(1).Provenance.ComputeDate, 'SENTINEL');
end

function test_tau_mismatch_recomputes(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    K = 40;
    tess_dirac_eigenmodes(SurfaceFile, 'K', K, 'Tau', 0.5, 'ForceRecompute', 1);
    TF = load(TessFile);
    TF.DiracEigen(1).Provenance.ComputeDate = 'SENTINEL';
    bst_save(TessFile, TF, 'v7');

    tess_dirac_eigenmodes(SurfaceFile, 'K', K, 'Tau', 0.75);  % different Tau -> recompute
    T = in_tess_bst(SurfaceFile, 0);
    verifyNotEqual(tc, T.DiracEigen(1).Provenance.ComputeDate, 'SENTINEL');
    verifyEqual(tc, T.DiracEigen(1).Tau, 0.75);
end

function test_nosave_returns_without_writing(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    TF = load(TessFile);
    if isfield(TF,'DiracEigen'), TF = rmfield(TF,'DiracEigen'); end
    bst_save(TessFile, TF, 'v7');

    DE = tess_dirac_eigenmodes(SurfaceFile, 'K', 40, 'NoSave', 1, 'ForceRecompute', 1);
    verifyEqual(tc, size(DE), [1 2]);
    T = in_tess_bst(SurfaceFile, 0);
    verifyFalse(tc, isfield(T,'DiracEigen'));
end
