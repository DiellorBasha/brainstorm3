function tests = test_lift_connectome_dirac
tests = functiontests(localfunctions);
end

function test_lift_shape_orthonormal(t)
    n = 6; K = 4;  rng(0);
    Ms = speye(n);                                  % identity mass -> Phis just needs orthonormal cols
    [Q,~] = qr(randn(n,K),0);  Phis = Q;  Lams = (1:K)';
    [Phiq, Lamq, Mq] = bst_lift_connectome_dirac(Phis, Lams, Ms);
    verifyEqual(t, size(Phiq), [4*n 3*K]);
    verifyEqual(t, size(Mq), [4*n 4*n]);
    verifyEqual(t, numel(Lamq), 3*K);
    % orthonormal in Mq
    G = Phiq' * Mq * Phiq;
    verifyLessThan(t, max(abs(G(:) - reshape(eye(3*K),[],1))), 1e-10);
    % each lifted column has w=0 and lives in exactly one imag axis
    for c = 1:size(Phiq,2)
        col = reshape(Phiq(:,c), 4, n);            % [w;x;y;z] x n
        verifyEqual(t, col(1,:), zeros(1,n), 'AbsTol', 0);      % w = 0
        nzAx = find(any(col(2:4,:) ~= 0, 2));
        verifyEqual(t, numel(nzAx), 1);            % exactly one of x/y/z active
    end
    % Lamq contains each Lams value 3 times
    verifyEqual(t, sort(Lamq), sort(repelem(Lams,3)));
end

function test_kron_mass(t)
    n=3; Ms=sprand(n,n,0.5); Ms=(Ms+Ms')/2+n*speye(n);
    [~,~,Mq]=bst_lift_connectome_dirac(eye(n),(1:n)',Ms);
    verifyLessThan(t, max(abs(full(Mq)-kron(full(Ms),eye(4))),[],'all'), 1e-12);
end
