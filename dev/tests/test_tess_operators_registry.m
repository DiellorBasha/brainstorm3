function test_tess_operators_registry()
% Build operators NoSave on a real cortex and check registry population.
% Picks the first cortex surface with a usable Structures atlas in the
% current protocol; errors with guidance if none is loaded.

    SurfaceFile = local_pick_cortex();
    assert(~isempty(SurfaceFile), ...
        'No cortex surface found in the current protocol — load a subject first.');

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
