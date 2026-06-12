function test_eigenmode_vector_field()
% TEST_EIGENMODE_VECTOR_FIELD  Headless regression for the Dirac eigenmode
% vector-field reconstruction (view_eigenmodes pure core). Requires Brainstorm
% on path so view_eigenmodes('ReconstructModeField', ...) dispatches.
% Authors: Diellor Basha, 2026
    nPass = 0; nFail = 0;

    % Synthetic 2-hemi EigenMat: nVh=2 per hemi (Phi{hh} is [4*2 x K]), K=3.
    % Global vertices: hemi L -> [1;2], hemi R -> [3;4]; nVert=4.
    % Per mode k, quaternion [w i j k] per vertex:
    %   L v1 = [0 k 0 0] (x=k), L v2 = [0 0 k 0] (y=k),
    %   R v3 = [0 0 0 k] (z=k), R v4 = [0 k k k] (k,k,k).
    EM.Variant = 'Dirac';
    EM.GlobalVertices = {[1;2], [3;4]};
    PhiL = zeros(8,3); PhiR = zeros(8,3);
    for k = 1:3
        PhiL(1:4,k) = [0;k;0;0];
        PhiL(5:8,k) = [0;0;k;0];
        PhiR(1:4,k) = [0;0;0;k];
        PhiR(5:8,k) = [0;k;k;k];
    end
    EM.Phi = {PhiL, PhiR};
    EM.Lambda = {(1:3)', (1:3)'};

    V3 = view_eigenmodes('ReconstructModeField', EM, 2, 4);
    [nPass,nFail] = chk('v1 = (2,0,0)', isequal(V3(1,:),[2 0 0]), nPass,nFail);
    [nPass,nFail] = chk('v2 = (0,2,0)', isequal(V3(2,:),[0 2 0]), nPass,nFail);
    [nPass,nFail] = chk('v3 = (0,0,2)', isequal(V3(3,:),[0 0 2]), nPass,nFail);
    [nPass,nFail] = chk('v4 = (2,2,2)', isequal(V3(4,:),[2 2 2]), nPass,nFail);

    % off-support vertices stay zero (nVert larger than mapped indices)
    V3b = view_eigenmodes('ReconstructModeField', EM, 1, 6);
    [nPass,nFail] = chk('off-support vertices zero', ...
        isequal(V3b(5,:),[0 0 0]) && isequal(V3b(6,:),[0 0 0]), nPass,nFail);

    % w slot is dropped: a large w must not change the vector part
    EM2 = EM; EM2.Phi{1}(1,2) = 99;   % w of L v1, mode 2
    V3c = view_eigenmodes('ReconstructModeField', EM2, 2, 4);
    [nPass,nFail] = chk('w slot ignored', isequal(V3c(1,:),[2 0 0]), nPass,nFail);

    % out-of-range mode errors
    err = false; try, view_eigenmodes('ReconstructModeField', EM, 9, 4); catch, err = true; end
    [nPass,nFail] = chk('out-of-range mode errors', err, nPass,nFail);

    fprintf('\n==== test_eigenmode_vector_field: %d passed, %d failed ====\n', nPass, nFail);
    if nFail > 0, error('test_eigenmode_vector_field: %d test(s) FAILED.', nFail); end
end

function [p,f] = chk(name, c, p, f)
    if c, fprintf('  PASS  %s\n', name); p=p+1; else, fprintf('  FAIL  %s\n', name); f=f+1; end
end
