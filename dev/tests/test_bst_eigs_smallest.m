function test_bst_eigs_smallest()
% TEST_BST_EIGS_SMALLEST: bst_eigs_smallest on singular symmetric AND Hermitian pencils.
    test_real_case();
    test_hermitian_case();
    test_real_operator_case();
end

function test_real_case()
% Real symmetric singular pencil: A = Neumann path-graph Laplacian (constant null
% vector, smallest eig 0); B = 1e-6*I (tiny mass mimics meters-scale conditioning).
    n = 60; k = 6;
    e = ones(n,1);
    A = spdiags([-e, 2*e, -e], -1:1, n, n);
    A(1,1) = 1; A(n,n) = 1;
    B = 1e-6 * speye(n);
    opts = struct('tol',1e-10,'disp',0,'maxit',1000);

    lam_ref = sort(eig(full(A), full(B)));
    lam_ref = lam_ref(1:k);

    lastwarn('');
    [V, D] = bst_eigs_smallest(A, B, k, opts);
    [wmsg, wid] = lastwarn;
    assert(isempty(wid), 'real: new path emitted a warning [%s] %s', wid, wmsg);

    lam = sort(real(diag(D)));
    scale = max(1, abs(lam_ref(end)));
    assert(max(abs(lam - lam_ref)) < 1e-6 * scale, 'real: spectrum mismatch %.3e', max(abs(lam - lam_ref)));
    assert(abs(lam(1)) < 1e-7 * scale, 'real: kernel mode not captured %.3e', lam(1));

    nrm = sqrt(real(diag(V' * (B * V)))); nrm(nrm < eps) = eps;
    Vn = V * diag(1 ./ nrm);
    Lam = diag(real(diag(D)));
    res = norm(A*Vn - B*Vn*Lam, 'fro') / (norm(full(A)) * sqrt(k));
    assert(res < 1e-6, 'real: eigenpair residual too large %.3e', res);
    G = Vn' * (B * Vn);
    assert(max(max(abs(G - eye(k)))) < 1e-8, 'real: not B-orthonormal %.3e', max(max(abs(G-eye(k)))));

    disp('PASS: test_bst_eigs_smallest (real symmetric case)');
end

function test_hermitian_case()
% Complex Hermitian singular pencil: a path "magnetic" Laplacian with a fixed edge
% phase. On a path (no cycle) the flux gauges away -> unitarily equivalent to the real
% path Laplacian -> still singular (one zero mode). Exercises the complex eigs path and
% the real-eigenvalue enforcement.
    n = 40; k = 5; theta = 0.3;
    A = diag(2*ones(n,1)) ...
        + diag(-exp(1i*theta)*ones(n-1,1), 1) ...
        + diag(-exp(-1i*theta)*ones(n-1,1), -1);
    A(1,1) = 1; A(n,n) = 1;
    A = sparse(A);
    B = 1e-6 * speye(n);
    opts = struct('tol',1e-10,'disp',0,'maxit',1000);

    lam_ref = sort(real(eig(full(A), full(B))));
    lam_ref = lam_ref(1:k);

    lastwarn('');
    [V, D] = bst_eigs_smallest(A, B, k, opts);
    [wmsg, wid] = lastwarn;
    assert(isempty(wid), 'hermitian: new path emitted a warning [%s] %s', wid, wmsg);

    lam = diag(D);
    assert(isreal(lam), 'hermitian: eigenvalues not enforced real (max|imag|=%.3e)', max(abs(imag(lam))));
    scale = max(1, abs(lam_ref(end)));
    assert(max(abs(lam - lam_ref)) < 1e-6 * scale, 'hermitian: spectrum mismatch %.3e', max(abs(lam - lam_ref)));
    assert(abs(lam(1)) < 1e-7 * scale, 'hermitian: kernel mode not captured %.3e', lam(1));

    nrm = sqrt(real(diag(V' * (B * V)))); nrm(nrm < eps) = eps;
    Vn = V * diag(1 ./ nrm);
    G = Vn' * (B * Vn);
    assert(max(max(abs(G - eye(k)))) < 1e-8, 'hermitian: not B-orthonormal %.3e', max(max(abs(G-eye(k)))));

    disp('PASS: test_bst_eigs_smallest (complex Hermitian case)');
end

function test_real_operator_case()
% Reuse a stored Dirac operator node (if present) to confirm on REAL data: legacy
% 'smallestabs' warns, bst_eigs_smallest does not, and the K smallest eigenvalues agree.
% Skips cleanly if no protocol / Dirac operator node is available.
    PI = bst_get('ProtocolInfo');
    if isempty(PI) || ~isfield(PI,'SUBJECTS') || isempty(PI.SUBJECTS)
        disp('SKIP: test_bst_eigs_smallest (real operator) -- no protocol loaded'); return;
    end
    opFiles = dir(fullfile(PI.SUBJECTS, '**', 'operator_*.mat'));
    f = '';
    for i = 1:numel(opFiles)
        p = fullfile(opFiles(i).folder, opFiles(i).name);
        S = load(p, 'Variant');
        if isfield(S,'Variant') && strcmpi(S.Variant,'Dirac'); f = p; break; end
    end
    if isempty(f)
        disp('SKIP: test_bst_eigs_smallest (real operator) -- no Dirac operator node found'); return;
    end
    Op = load(f); A = Op.Operator{1}; B = Op.Mass{1};
    k = 12; opts = struct('tol',1e-6,'maxit',1000,'disp',0);

    lastwarn('');
    lam_old = sort(real(eigs(A, B, k, 'smallestabs', opts)));
    [~, wid_old] = lastwarn;

    lastwarn('');
    [~, D] = bst_eigs_smallest(A, B, k, opts);
    [wmsg_new, wid_new] = lastwarn;
    lam_new = sort(real(diag(D)));

    assert(isempty(wid_new), 'real operator: new path emitted a warning [%s] %s', wid_new, wmsg_new);
    sc = max(abs(lam_old));
    assert(max(abs(lam_old - lam_new)) < 1e-3 * sc, ...
        'real operator: spectrum mismatch %.3e (rel %.1e)', max(abs(lam_old-lam_new)), max(abs(lam_old-lam_new))/sc);
    fprintf('PASS: test_bst_eigs_smallest (real operator; legacy warned id=%s)\n', wid_old);
end
