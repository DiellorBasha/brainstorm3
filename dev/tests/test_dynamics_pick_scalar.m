function test_dynamics_pick_scalar
    Ht = struct('Div',[1;2;3], 'Curl',[4;5;6], 'Phi',[7;8;9], 'Psi',[10;11;12]);
    assert(isequal(view_dynamics('PickScalar', Ht, 'Divergence'), [1;2;3]), 'Divergence->Div');
    assert(isequal(view_dynamics('PickScalar', Ht, 'Curl'),       [4;5;6]), 'Curl->Curl');
    assert(isequal(view_dynamics('PickScalar', Ht, 'Potential'),  [7;8;9]), 'Potential->Phi');
    assert(isequal(view_dynamics('PickScalar', Ht, 'Stream'),     [10;11;12]), 'Stream->Psi');
    ok = false; try, view_dynamics('PickScalar', Ht, 'Bogus'); catch, ok = true; end
    assert(ok, 'unknown operator must error');
    fprintf('PASS test_dynamics_pick_scalar\n');
end
