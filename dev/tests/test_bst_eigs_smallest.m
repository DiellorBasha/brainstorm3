function test_bst_eigs_smallest()
% TEST_BST_EIGS_SMALLEST: synthetic singular symmetric generalized pencil.
% A = Neumann path-graph Laplacian (singular: constant null vector, smallest eig 0).
% B = 1e-6*I  (tiny mass, mimics meters-scale conditioning: generalized eig ~1e6*eig(A)).
    n = 60; k = 6;
    e = ones(n,1);
    A = spdiags([-e, 2*e, -e], -1:1, n, n);
    A(1,1) = 1; A(n,n) = 1;                 % Neumann ends -> row sums 0 -> constant null
    B = 1e-6 * speye(n);
    opts = struct('tol',1e-10,'disp',0,'maxit',1000);

    % Dense reference for the k smallest generalized eigenvalues
    lam_ref = sort(eig(full(A), full(B)));
    lam_ref = lam_ref(1:k);

    % Solver under test: must not raise the nearly-singular warning
    lastwarn('');
    [V, D] = bst_eigs_smallest(A, B, k, opts);
    [wmsg, wid] = lastwarn;
    assert(~strcmp(wid, 'MATLAB:nearlySingularMatrix'), ...
        'bst_eigs_smallest raised the nearly-singular warning: %s', wmsg);

    lam = sort(real(diag(D)));
    % (1) spectrum matches the dense reference
    scale = max(1, abs(lam_ref(end)));
    assert(max(abs(lam - lam_ref)) < 1e-6 * scale, ...
        'spectrum mismatch: max diff %.3e', max(abs(lam - lam_ref)));
    % (2) the zero (kernel) mode is captured
    assert(abs(lam(1)) < 1e-7 * scale, 'zero/kernel mode not captured: %.3e', lam(1));
    % (3) B-orthonormalize the returned vectors, then check residual + orthonormality
    nrm = sqrt(real(diag(V' * (B * V))));
    nrm(nrm < eps) = eps;
    Vn = V * diag(1 ./ nrm);
    Lam = diag(real(diag(D)));
    res = norm(A*Vn - B*Vn*Lam, 'fro') / (normest(A) * sqrt(k));
    assert(res < 1e-6, 'eigenpair residual too large: %.3e', res);
    G = Vn' * (B * Vn);
    assert(max(max(abs(G - eye(k)))) < 1e-8, 'not B-orthonormal: %.3e', max(max(abs(G-eye(k)))));

    disp('PASS: test_bst_eigs_smallest');
end
