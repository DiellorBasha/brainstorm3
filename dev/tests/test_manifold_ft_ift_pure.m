function test_manifold_ft_ift_pure
% Pure manifold Fourier pair: ft = Phi'*M*U, ift = Phi*C. Round-trip onto the basis
% subspace recovers the field; ift handles a matrix C (kernel); dimension errors raise.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);

nV = 20; K = 8;
[Q, ~] = qr(randn(nV, K), 0);     % Phi: orthonormal columns
Phi = Q; M = speye(nV);            % M = I -> Phi M-orthonormal
u = randn(nV, 3);

c = manifold_ft(Phi, M, u);
assert(isequal(size(c), [K 3]), 'manifold_ft output shape');
uproj = Phi * (Phi' * u);          % projection of u onto span(Phi)
urec  = manifold_ift(Phi, c);
assert(norm(urec - uproj, 'fro')/norm(uproj,'fro') < 1e-10, 'round-trip onto subspace failed');

% ift handles a matrix C (mode-space kernel -> vertex kernel)
Kk = manifold_ift(Phi, randn(K, 5));
assert(isequal(size(Kk), [nV 5]), 'manifold_ift matrix shape');

% dimension errors
ok=false; try; manifold_ft(Phi, M, randn(nV+1,2)); catch; ok=true; end; assert(ok, 'ft dim check');
ok=false; try; manifold_ift(Phi, randn(K+1,2)); catch; ok=true; end; assert(ok, 'ift dim check');
disp('ALL TESTS PASSED');
end
