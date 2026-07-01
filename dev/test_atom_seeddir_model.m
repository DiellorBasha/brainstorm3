function tests = test_atom_seeddir_model
tests = functiontests(localfunctions);
end
function test_template_has_seeddir(tc)
    G = db_template('atomgroup');
    verifyTrue(tc, isfield(G, 'SeedDir'));
    verifyEmpty(tc, G.SeedDir);
end
function test_roundtrip(tc)
    T = bst_dynamics('New', 'test');
    G = bst_dynamics('NewGroup', 'a'); G.vertices = 5; G.Operator = 'Dirac'; G.SeedDir = [0 0 1];
    T = bst_dynamics('AddGroup', T, G);
    f = [tempname '.mat'];  bst_dynamics('Save', f, T);  T2 = bst_dynamics('Load', f);  delete(f);
    verifyEqual(tc, T2.Groups(1).SeedDir, [0 0 1]);
end
