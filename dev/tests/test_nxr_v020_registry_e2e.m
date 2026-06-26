function test_nxr_v020_registry_e2e()
% Build every Variant NoSave on a real cortex and assert the registry primary
% id matches the expected mapping. Face/Covariant variants require the L/R
% Structures atlas (present on FreeSurfer-imported cortex).

    SurfaceFile = local_pick_cortex();
    assert(~isempty(SurfaceFile), 'No cortex surface in the current protocol.');

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

function SurfaceFile = local_pick_cortex()
    SurfaceFile = '';
    ProtocolSubjects = bst_get('ProtocolSubjects');
    if isempty(ProtocolSubjects), return; end
    allSubj = [ProtocolSubjects.Subject];
    for i = 1:numel(allSubj)
        S = allSubj(i);
        if ~isfield(S,'Surface') || isempty(S.Surface), continue; end
        k = find(strcmpi({S.Surface.SurfaceType}, 'Cortex'), 1);
        if ~isempty(k), SurfaceFile = S.Surface(k).FileName; return; end
    end
end
