function tests = test_apply_cache_equiv
tests = functiontests(localfunctions);
end
function test_cached_equals_analysis(tc)
% Controller runs this with a surface on the path. Build a small LB eigenbasis, a random
% scalar field, and verify C-then-gain == bst_eigenfilter('Analysis').
% bst_eigen('Axes') needs the time-axis fields (as the panel's i_atom_axes passes); the equivalence
% only uses the spatial basis (Phi/Lambda/Mass), so a short placeholder window is fine.
ax = bst_eigen('Axes', struct('SurfaceFile', getenv('BST_TEST_SURF'), 'Variant','Laplace-Beltrami', ...
               'nModes',40, 'TimeWindow',[0 0.04], 'SampleRate',100));
nV = 0; for h=1:numel(ax.GlobalVertices), nV = max(nV, max(ax.GlobalVertices{h}(:))); end
F  = randn(nV, 5);
kp = struct('t', 0.02, 'lmax', max(ax.Lambda{1}(:)));
g  = bst_eigfilter_kernel('heat', kp);
% cached path
Ffilt = zeros(nV,5);
for h=1:numel(ax.Phi)
    gv = ax.GlobalVertices{h}(:); C = manifold_ft(ax.Phi{h}, ax.Mass{h}, F(gv,:));
    hg = g(ax.Lambda{h}(:)); Ffilt(gv,:) = manifold_ift(ax.Phi{h}, hg(:).*C);
end
% direct Analysis
EigenMat = struct('Phi',{ax.Phi},'Lambda',{ax.Lambda},'Variant','Laplace-Beltrami','GlobalVertices',{ax.GlobalVertices});
OperatorMat = struct('Mass',{ax.Mass});
[Fana,~,err] = bst_eigenfilter('Analysis', F, EigenMat, OperatorMat, 'heat', kp);
verifyEqual(tc, err, 0);
verifyLessThan(tc, max(abs(Ffilt(:)-Fana(:))), 1e-9);
end
