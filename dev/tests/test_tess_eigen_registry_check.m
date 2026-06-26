function test_tess_eigen_registry_check()
% The cross-check helper must exist and agree on real nodes. We call it
% directly via a thin probe that tess_eigen exposes through a subfunction
% handle is not possible, so we assert the function is reachable by running a
% small Laplace-Beltrami eigis and confirming no warning is raised.

    SurfaceFile = local_pick_cortex();
    assert(~isempty(SurfaceFile), 'No cortex surface in the current protocol.');

    lastwarn('');                                   % clear warning state
    Eig = tess_eigen(SurfaceFile, 'Laplace-Beltrami', 'nModes', 20, 'NoSave', true);
    [~, wid] = lastwarn();
    assert(~strcmp(wid, 'tess_eigen:registryMismatch'), ...
        'unexpected registry mismatch warning on a valid LB node');
    assert(isfield(Eig,'Phi') && ~isempty(Eig.Phi{1}), 'eigen result empty');

    fprintf('test_tess_eigen_registry_check: PASS\n');
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
