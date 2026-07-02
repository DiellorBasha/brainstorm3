function tests = test_quat_imag
tests = functiontests(localfunctions);
end
function test_extracts_imag(t)
    % Q = [4n x T] interleaved [w,x,y,z]; V = [3n x T] [x,y,z]
    n = 5; T = 3;  rng(0);
    w = randn(n,T); x = randn(n,T); y = randn(n,T); z = randn(n,T);
    Q = zeros(4*n,T);
    Q(1:4:end,:)=w; Q(2:4:end,:)=x; Q(3:4:end,:)=y; Q(4:4:end,:)=z;
    V = manifold_quat_imag(Q);
    verifyEqual(t, size(V), [3*n T]);
    verifyEqual(t, V(1:3:end,:), x, 'AbsTol', 0);   % x rows
    verifyEqual(t, V(2:3:end,:), y, 'AbsTol', 0);   % y rows
    verifyEqual(t, V(3:3:end,:), z, 'AbsTol', 0);   % z rows (w dropped)
end
function test_bad_size_errors(t)
    verifyError(t, @() manifold_quat_imag(zeros(7,2)), 'manifold_quat_imag:size');
end
