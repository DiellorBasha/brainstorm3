function test_bst_poisson
% bst_poisson scalar solve reproduces the legacy bst_operators i_poisson math exactly,
% on a synthetic LBO-like node (no DB / no nxr).
    fprintf('== test_bst_poisson ==\n');
    rng(1);
    n = 300;
    G = sprandsym(n, 0.04);  G = abs(G - diag(diag(G)));
    K = spdiags(sum(G,2), 0, n, n) - G;  K = (K+K')/2;     % cotan-like, constant nullspace
    M = spdiags(rand(n,1)+0.5, 0, n, n);
    Node = db_template('operatormat');
    Node.Variant='Laplace-Beltrami'; Node.ParentSurface='unit://test';
    Node.Operator={K,K}; Node.Mass={M,M}; Node.GlobalVertices={(1:n)',(1:n)'};
    f = randn(n, 4);

    phi = bst_poisson(Node, f);

    % reference: legacy projection/pin/recenter (verbatim from bst_operators i_poisson)
    free=(2:n)'; totMass=sum(M(:));  ref=zeros(n,4);
    dK = decomposition(K(free,free),'chol');
    for c=1:4
        rhs = M*(f(:,c) - sum(M*f(:,c))/totMass);
        x=zeros(n,1); x(free)=dK\rhs(free); ref(:,c)=x-mean(x);
    end
    err = norm(phi-ref,'fro')/max(norm(ref,'fro'),eps);
    assert(err < 1e-10, 'bst_poisson != legacy i_poisson: rel err %g', err);
    fprintf('  bst_poisson vs legacy rel err = %g  [OK]\nPASS\n', err);
end
