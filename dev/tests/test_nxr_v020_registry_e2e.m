function test_nxr_v020_registry_e2e()
% Build every Variant NoSave on a real cortex and assert the registry primary
% id matches the expected mapping. Face/Covariant variants require the L/R
% Structures atlas (present on FreeSurfer-imported cortex).

    SurfaceFile = bst_canonical_cortex(20484);   % ico5 canonical test cortex (has VertNormals)

    cases = { ...
        'Laplace-Beltrami',     'laplaceBeltrami'; ...
        'Connection Laplacian', 'leviCivitaConnectionLaplacian'; ...
        'Dirac',                'relativeDirac'; ...
        'Dirac-Face',           'relativeFaceDirac'; ...
        'Hodge-Face',           'faceLaplacianGreenGauss'; ...
        'Covariant',            'flatCovariantLaplacian'};

    for i = 1:size(cases,1)
        V = cases{i,1}; want = cases{i,2};
        Op = tess_operators(SurfaceFile, V, 'NoSave', true, 'Tau', 0.5);
        assert(~isempty(Op.Registry), '%s: Registry empty', V);
        got = Op.Registry.Primary.id;
        assert(strcmp(got, want), '%s: primary id %s, expected %s', V, got, want);
        fprintf('  %-22s -> %s  OK\n', V, got);
    end
    fprintf('test_nxr_v020_registry_e2e: ALL PASS\n');
end

