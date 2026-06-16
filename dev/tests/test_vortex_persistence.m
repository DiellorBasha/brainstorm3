function test_vortex_persistence()
% Pure unit tests for bst_vortex_persistence on hand-computed tiny graphs.
% Author: Diellor Basha, 2026
    nFail = 0;

    % Chain 1-2-3-4-5, field = [5 2 4 1 3].
    % +field: global max v1 (Inf); v3 dies at saddle v2 -> pers 4-2=2;
    %         v5 dies at saddle v4 -> pers 3-1=2.
    % -field: global max v4 (Inf, min of field); v2 dies at saddle v3 -> pers 2.
    nb = {2; [1;3]; [2;4]; [3;5]; 4};
    f  = [5;2;4;1;3];
    C  = bst_vortex_persistence(f, [], 'Neighbors', nb);

    nFail = nFail + chk('5 features total', numel(C.vertex) == 5);
    nFail = nFail + chk('v1 is +global',  any(C.vertex==1 & C.chirality==1 & isinf(C.persistence)));
    nFail = nFail + chk('v4 is -global',  any(C.vertex==4 & C.chirality==-1 & isinf(C.persistence)));
    nFail = nFail + chk('v3 +core pers=2', any(C.vertex==3 & C.chirality==1  & abs(C.persistence-2)<1e-12));
    nFail = nFail + chk('v5 +core pers=2', any(C.vertex==5 & C.chirality==1  & abs(C.persistence-2)<1e-12));
    nFail = nFail + chk('v2 -core pers=2', any(C.vertex==2 & C.chirality==-1 & abs(C.persistence-2)<1e-12));
    nFail = nFail + chk('sorted by persistence desc', issorted(C.persistence,'descend'));

    % MinPersistence=3 keeps only the two Inf globals.
    C2 = bst_vortex_persistence(f, [], 'Neighbors', nb, 'MinPersistence', 3);
    nFail = nFail + chk('MinPersistence keeps only globals', numel(C2.vertex)==2 && all(isinf(C2.persistence)));

    % Flat field -> one global per sign, no spurious finite-persistence cores.
    Cf = bst_vortex_persistence(zeros(5,1), [], 'Neighbors', nb);
    nFail = nFail + chk('flat field -> 2 globals only', numel(Cf.vertex)==2 && all(isinf(Cf.persistence)));

    % Faces path builds 1-ring: single triangle, peak at the highest vertex.
    Ct = bst_vortex_persistence([3;1;2], [1 2 3]);
    nFail = nFail + chk('triangle global is vertex 1', any(Ct.vertex==1 & isinf(Ct.persistence) & Ct.chirality==1));

    fprintf('\n==== test_vortex_persistence: %d failed ====\n', nFail);
    if nFail > 0, error('test_vortex_persistence FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
