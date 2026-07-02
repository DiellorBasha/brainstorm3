function test_realise_core_delegates
% i_atom_realise_core must now be a thin passthrough of bst_eigenfilter('Atom'): identical [W,gv,V3].
    surf = getenv('BST_TEST_SURF');
    if isempty(surf), fprintf('SKIP test_realise_core_delegates: set BST_TEST_SURF to a cortex surface file.\n'); return; end
    kp = struct('lmax', []);

    for variant = {'Laplace-Beltrami','Dirac'}
        ax   = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant',variant{1},'nModes',50,'TimeWindow',[0 0.5],'SampleRate',100));
        seed = ax.GlobalVertices{1}(1);
        dir  = bst_atom_default_dir(ax, seed);
        [Wc, gvc, V3c]   = panel_bst_dynamics('i_atom_realise_core', ax, 'heat', kp, seed, dir);
        [Wa, gva, V3a]   = bst_eigenfilter('Atom', ax, 'heat', kp, seed, dir);
        assert(isequal(Wc,Wa) && isequal(gvc,gva) && isequal(V3c,V3a), ...
            'realise_core must equal Atom for %s', variant{1});
    end
    disp('test_realise_core_delegates PASSED');
end
