function tests = test_atom_realise_v3
tests = functiontests(localfunctions);
end
function test_dirac_v3_shape_and_dir(tc)
    surf = getenv('BST_TEST_SURF'); if isempty(surf), surf='sub-MTL0002/tess_cortex_pial_low.mat'; end
    ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Dirac','nModes',80,'TimeWindow',[0 .04],'SampleRate',100));
    seed = ax.GlobalVertices{1}(1);  kp = struct('lmax', max(ax.Lambda{1}(:)));
    % call the pure realise-core directly (exposed for tests): [W,gv,V3] = i_atom_realise_core(ax,k,kp,seed,dir)
    [~, gv, V3] = panel_bst_dynamics('i_atom_realise_core', ax, 'diffusion', kp, seed, [0 0 1]); %#ok
    nV = 0; for h=1:numel(ax.GlobalVertices), nV=max(nV,max(ax.GlobalVertices{h}(:))); end
    verifyEqual(tc, size(V3), [nV 3]);
    [~,ipk] = max(sqrt(sum(V3.^2,2)));
    verifyGreaterThan(tc, abs(V3(ipk,3)), max(abs(V3(ipk,1)), abs(V3(ipk,2))), 'peak dipole ~ +Z');
end
function test_scalar_v3_empty(tc)
    surf = getenv('BST_TEST_SURF'); if isempty(surf), surf='sub-MTL0002/tess_cortex_pial_low.mat'; end
    ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Laplace-Beltrami','nModes',40,'TimeWindow',[0 .04],'SampleRate',100));
    seed = ax.GlobalVertices{1}(1);  kp = struct('lmax', max(ax.Lambda{1}(:)));
    [~,~,V3] = panel_bst_dynamics('i_atom_realise_core', ax, 'diffusion', kp, seed, 1);
    verifyEmpty(tc, V3);
end
