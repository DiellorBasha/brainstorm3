function tests = test_atom_seed
tests = functiontests(localfunctions);
end
function test_dirac_seed_locates_at_vertex(tc)
    surf = getenv('BST_TEST_SURF'); if isempty(surf), surf='sub-MTL0002/tess_cortex_pial_low.mat'; end
    ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Dirac','nModes',80,'TimeWindow',[0 .04],'SampleRate',100));
    gv = ax.GlobalVertices{1};  seed = gv(round(numel(gv)/2));       % a vertex that is NOT local index 1
    kp = struct('lmax', max(ax.Lambda{1}(:)));
    [W, gvb] = bst_eigenfilter('Atom', ax, 'diffusion', kp, seed, [0;0;1]);  % +Z dipole
    n = numel(gvb);
    imag3 = [W(2:4:end,1) W(3:4:end,1) W(4:4:end,1)];  nrm = sqrt(sum(imag3.^2,2));
    [~,ipk] = max(nrm);
    Vxyz = getfield(in_tess_bst(surf,0),'Vertices'); %#ok
    dmm = norm(Vxyz(gvb(ipk),:) - Vxyz(seed,:))*1000;
    verifyLessThan(tc, dmm, 15, 'Dirac impulse must peak within 15mm of the seed (mis-seed regression)');
end
function test_dirac_seed_direction(tc)
    surf = getenv('BST_TEST_SURF'); if isempty(surf), surf='sub-MTL0002/tess_cortex_pial_low.mat'; end
    ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Dirac','nModes',80,'TimeWindow',[0 .04],'SampleRate',100));
    gv = ax.GlobalVertices{1};  seed = gv(1);  kp = struct('lmax', max(ax.Lambda{1}(:)));
    Wx = bst_eigenfilter('Atom', ax, 'diffusion', kp, seed, [1;0;0]);
    Wz = bst_eigenfilter('Atom', ax, 'diffusion', kp, seed, [0;0;1]);
    verifyGreaterThan(tc, norm(Wx - Wz), 1e-9, 'different seed directions must give different fields');
end
function test_scalar_unchanged(tc)
    surf = getenv('BST_TEST_SURF'); if isempty(surf), surf='sub-MTL0002/tess_cortex_pial_low.mat'; end
    ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Laplace-Beltrami','nModes',40,'TimeWindow',[0 .04],'SampleRate',100));
    seed = ax.GlobalVertices{1}(1);  kp = struct('lmax', max(ax.Lambda{1}(:)));
    [W, gvb] = bst_eigenfilter('Atom', ax, 'diffusion', kp, seed);      % no seedDir -> default amplitude 1
    verifyEqual(tc, size(W,1), numel(gvb));                             % scalar: one row per vertex
    verifyGreaterThan(tc, W(1), max(W)*0.5);                            % peak at the seed (local idx 1)
end
