function tests = test_atom_default_dir
tests = functiontests(localfunctions);
end

function test_dirac_normal(tc)
    surf = getenv('BST_TEST_SURF'); if isempty(surf), surf='sub-MTL0002/tess_cortex_pial_low.mat'; end
    ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Dirac','nModes',6,'TimeWindow',[0 .04],'SampleRate',100));
    seed = ax.GlobalVertices{1}(10);
    d = panel_bst_dynamics('i_atom_default_dir', ax, seed);
    verifyEqual(tc, numel(d), 3);  verifyEqual(tc, norm(d), 1, 'AbsTol', 1e-6);
    S = in_tess_bst(surf, 0);
    verifyGreaterThan(tc, abs(dot(d(:), S.VertNormals(seed,:)')), 0.99, 'must be the seed surface normal');
end

function test_scalar_amplitude(tc)
    surf = getenv('BST_TEST_SURF'); if isempty(surf), surf='sub-MTL0002/tess_cortex_pial_low.mat'; end
    ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Laplace-Beltrami','nModes',6,'TimeWindow',[0 .04],'SampleRate',100));
    d = panel_bst_dynamics('i_atom_default_dir', ax, ax.GlobalVertices{1}(1));
    verifyEqual(tc, d, 1);
end
