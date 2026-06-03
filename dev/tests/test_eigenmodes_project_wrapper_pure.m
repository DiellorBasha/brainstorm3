function test_eigenmodes_project_wrapper_pure
% project's coefficients equal manifold_ft; its reconstruct output equals manifold_ift;
% ModeRange selects over the canonical Order.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);

nV=12; K=6; [Q,~]=qr(randn(nV,K),0);
Eig=struct('Vectors',Q,'Values',[10;30;50;20;40;60],'Component',[1;1;1;2;2;2], ...
           'CompRank',[1;2;3;1;2;3]);
[~,Eig.Order]=sort(Eig.Values,'ascend');   % [1 4 2 5 3 6]
M=speye(nV); u=randn(nV,4);

C = bst_eigenmodes_project(Eig, u, M);
assert(isequal(C, manifold_ft(Q, M, u)), 'project coeffs must equal manifold_ft');

% Reconstruct over canonical ranks 1..3 (lowest 3 eigenvalues: 10,20,30 -> stored [1 4 2])
[~, R] = bst_eigenmodes_project(Eig, u, M, 'ModeRange', [1 3]);
iSel = Eig.Order(1:3);
Rexp = manifold_ift(Q(:,iSel), C(iSel,:));
assert(norm(R - Rexp,'fro') < 1e-12, 'ranged reconstruct must use canonical order');
disp('ALL TESTS PASSED');
end
