function tests = test_atom_fiber
tests = functiontests(localfunctions);
end
function test_fieldtype_map(tc)
    surf = getenv('BST_TEST_SURF'); if isempty(surf), surf='sub-MTL0002/tess_cortex_pial_low.mat'; end
    map = {'Laplace-Beltrami',1,'scalar'; 'Connection Laplacian',2,'tangent'; 'Dirac',4,'quaternion'};
    for i=1:size(map,1)
        ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant',map{i,1},'nModes',6,'TimeWindow',[0 .04],'SampleRate',100));
        [C,kind] = bst_eigenfilter('Fiber', ax);
        verifyEqual(tc, C, map{i,2}, sprintf('%s C', map{i,1}));
        verifyEqual(tc, kind, map{i,3}, sprintf('%s kind', map{i,1}));
    end
end
function test_layout_fallback(tc)
    surf = getenv('BST_TEST_SURF'); if isempty(surf), surf='sub-MTL0002/tess_cortex_pial_low.mat'; end
    ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Dirac','nModes',6,'TimeWindow',[0 .04],'SampleRate',100));
    ax.Operator = [];                          % strip Registry -> must fall back to Phi layout
    [C,kind] = bst_eigenfilter('Fiber', ax);
    verifyEqual(tc, C, 4); verifyEqual(tc, kind, 'quaternion');
end
