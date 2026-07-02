function test_atom_default_dir
% bst_atom_default_dir: scalar -> 1; quaternion -> unit seed normal (byte-equiv to the former panel local).
    surf = getenv('BST_TEST_SURF');
    if isempty(surf), fprintf('SKIP test_atom_default_dir: set BST_TEST_SURF to a cortex surface file.\n'); return; end

    axS = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Laplace-Beltrami','nModes',40,'TimeWindow',[0 0.5],'SampleRate',100));
    seedS = axS.GlobalVertices{1}(1);
    assert(isequal(bst_atom_default_dir(axS, seedS), 1), 'scalar default must be 1');

    axQ = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Dirac','nModes',40,'TimeWindow',[0 0.5],'SampleRate',100));
    seedQ = axQ.GlobalVertices{1}(1);
    d = bst_atom_default_dir(axQ, seedQ);
    assert(isequal(size(d),[1 3]) && abs(norm(d)-1) < 1e-9, 'quaternion default must be a unit row 3-vector');
    S = in_tess_bst(axQ.SurfaceFile, 0);
    nref = S.VertNormals(seedQ,:);  nref = nref / norm(nref);
    assert(max(abs(d - nref)) < 1e-9, 'quaternion default must equal the unit seed normal');
    disp('test_atom_default_dir PASSED');
end
