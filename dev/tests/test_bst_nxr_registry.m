function test_bst_nxr_registry()
% Unit tests for the registry helper. The operator()/field() calls hit the
% live MEX (requires nxr v0.2.0 loaded); the map accessors are pure.

    % --- pure map accessors ---
    assert(strcmp(bst_nxr_registry('idForVariant','Laplace-Beltrami'),     'laplaceBeltrami'));
    assert(strcmp(bst_nxr_registry('idForVariant','Connection Laplacian'), 'leviCivitaConnectionLaplacian'));
    assert(strcmp(bst_nxr_registry('idForVariant','Dirac'),                'relativeDirac'));
    assert(strcmp(bst_nxr_registry('idForVariant','Dirac-Face'),           'relativeFaceDirac'));
    assert(strcmp(bst_nxr_registry('idForVariant','Hodge-Face'),           'faceLaplacianGreenGauss'));
    assert(strcmp(bst_nxr_registry('idForVariant','Covariant'),            'flatCovariantLaplacian'));
    assert(isempty(bst_nxr_registry('idForVariant','Bogus')));

    comps = bst_nxr_registry('componentsForVariant','Dirac');
    assert(iscell(comps) && any(strcmp(comps,'intrinsicDirac')) && any(strcmp(comps,'extrinsicDirac')));
    assert(isempty(bst_nxr_registry('componentsForVariant','Bogus')));

    % --- live registry passthrough ---
    oi = bst_nxr_registry('operator','laplaceBeltrami');
    assert(isstruct(oi) && strcmp(oi.id,'laplaceBeltrami') && isfield(oi,'bundle') && isfield(oi,'role'));
    assert(isempty(bst_nxr_registry('operator','definitelyNotAnId')));   % graceful []

    fprintf('test_bst_nxr_registry: ALL PASS\n');
end
