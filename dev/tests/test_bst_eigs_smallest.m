function test_bst_eigs_smallest()
% TEST_BST_EIGS_SMALLEST: bst_eigs_smallest on singular symmetric AND Hermitian pencils.
    test_real_case();
    test_hermitian_case();
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
    assert(~strcmp(wid, 'MATLAB:nearlySingularMatrix'), 'real: nearly-singular warning: %s', wmsg);

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
    assert(~strcmp(wid, 'MATLAB:nearlySingularMatrix'), 'hermitian: nearly-singular warning: %s', wmsg);

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
