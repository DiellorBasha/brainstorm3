function test_operatormat_registry_field()
    t = db_template('operatormat');
    assert(isfield(t, 'Registry'), 'operatormat missing Registry field');
    assert(isempty(t.Registry), 'operatormat.Registry default should be []');
    fprintf('test_operatormat_registry_field: PASS\n');
end
