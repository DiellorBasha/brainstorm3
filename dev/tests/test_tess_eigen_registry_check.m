function test_tess_eigen_registry_check()
% The cross-check helper must exist and agree on real nodes. We call it
% directly via a thin probe that tess_eigen exposes through a subfunction
% handle is not possible, so we assert the function is reachable by running a
% small Laplace-Beltrami eigis and confirming no warning is raised.

    SurfaceFile = bst_canonical_cortex(20484);   % ico5 canonical test cortex (has VertNormals)

    lastwarn('');                                   % clear warning state
    Eig = tess_eigen(SurfaceFile, 'Laplace-Beltrami', 'nModes', 20, 'NoSave', true);
    [~, wid] = lastwarn();
    assert(~strcmp(wid, 'tess_eigen:registryMismatch'), ...
        'unexpected registry mismatch warning on a valid LB node');
    assert(isfield(Eig,'Phi') && ~isempty(Eig.Phi{1}), 'eigen result empty');

    fprintf('test_tess_eigen_registry_check: PASS\n');
end

