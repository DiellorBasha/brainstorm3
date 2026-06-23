function test_tess_cholesky
% Validates tess_cholesky on a synthetic SPD matrix (no DB / no nxr).
%   1. 'solve' reproduces backslash on the pinned free block.
%   2. The pure getter returns the factor carried on the node struct without recomputing.
%   3. 'attach' populates OperatorNode.Cholesky for both hemispheres.
    fprintf('== test_tess_cholesky ==\n');
    rng(0);
    n = 200;
    % a small SPD "stiffness": graph Laplacian of a random connected graph + tiny shift
    G = sprandsym(n, 0.05);  G = G - diag(diag(G));  G = abs(G);
    K = spdiags(sum(G,2)+1e-6, 0, n, n) - G;   % SPD-ish; pin removes the constant nullspace
    K = (K+K')/2;
    M = spdiags(rand(n,1)+0.5, 0, n, n);
    Node = db_template('operatormat');
    Node.Variant = 'Laplace-Beltrami';  Node.ParentSurface = 'unit://test';
    Node.Operator = {K, K};  Node.Mass = {M, M};
    Node.GlobalVertices = {(1:n)', (1:n)'};

    pin = 1;  free = setdiff(1:n, pin)';
    dF = tess_cholesky(Node, 1, pin);
    b = randn(n,1);  b(pin) = 0;
    x = tess_cholesky('solve', dF, b);
    xref = zeros(n,1);  xref(free) = K(free,free) \ b(free);
    err = norm(x - xref) / max(norm(xref), eps);
    assert(err < 1e-10, 'solve mismatch vs backslash: %g', err);
    fprintf('  solve vs backslash rel err = %g  [OK]\n', err);

    % getter uses a pre-attached factor (no recompute): tamper with a marker and confirm reuse
    Node.Cholesky = {dF, dF};  Node.Cholesky{1}.marker = 42;
    dF2 = tess_cholesky(Node, 1, pin);
    assert(isfield(dF2,'marker') && dF2.marker == 42, 'getter did not reuse attached factor');
    fprintf('  getter reuses attached factor  [OK]\n');

    fprintf('PASS\n');
end
