function test_time_derivative()
% Unit tests for bst_time_derivative (backward finite differences).
% Author: Diellor Basha, 2026
    nFail = 0;
    a=[1;2]; b=[4;6]; c=[9;12]; dt=0.5;
    nFail = nFail + chk('order 0 = current frame',     isequal(bst_time_derivative([a b c], dt, 0), c));
    nFail = nFail + chk('order 1 = (b-a)/dt',          isequal(bst_time_derivative([a b],   dt, 1), (b-a)/dt));
    nFail = nFail + chk('order 2 = (c-2b+a)/dt^2',     isequal(bst_time_derivative([a b c], dt, 2), (c-2*b+a)/dt^2));
    nFail = nFail + chk('uses newest cols (extra ok)', isequal(bst_time_derivative([a b c], dt, 1), (c-b)/dt));
    ok = false; try, bst_time_derivative(a, dt, 1); catch, ok = true; end
    nFail = nFail + chk('too few frames errors', ok);
    fprintf('\n==== test_time_derivative: %d failed ====\n', nFail);
    if nFail > 0, error('test_time_derivative FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
