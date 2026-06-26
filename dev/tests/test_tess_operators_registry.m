function test_tess_operators_registry()
% Build operators NoSave on a real cortex and check registry population.
% Picks the first cortex surface with a usable Structures atlas in the
% current protocol; errors with guidance if none is loaded.

    SurfaceFile = bst_canonical_cortex(20484);   % ico5 canonical test cortex (has VertNormals)

    Op = tess_operators(SurfaceFile, 'Laplace-Beltrami', 'NoSave', true);
    assert(isfield(Op,'Registry') && ~isempty(Op.Registry), 'Registry not populated');
    assert(strcmp(Op.Registry.Primary.id, 'laplaceBeltrami'), ...
        'wrong primary id: %s', Op.Registry.Primary.id);

    Od = tess_operators(SurfaceFile, 'Dirac', 'NoSave', true, 'Tau', 0.5);
    assert(strcmp(Od.Registry.Primary.id, 'relativeDirac'), ...
        'wrong Dirac primary id: %s', Od.Registry.Primary.id);
    cids = {Od.Registry.Components.id};
    assert(any(strcmp(cids,'intrinsicDirac')) && any(strcmp(cids,'extrinsicDirac')), ...
        'Dirac components missing from Registry');

    fprintf('test_tess_operators_registry: ALL PASS\n');
end

