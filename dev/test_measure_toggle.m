function tests = test_measure_toggle
tests = functiontests(localfunctions);
end
function test_measure_default_and_set(t)
    verifyEqual(t, panel_bst_dynamics('i_measure_default'), 'amplitude');
    verifyTrue(t, ismember('dspm', panel_bst_dynamics('i_measure_options')));
    verifyTrue(t, ismember('amplitude', panel_bst_dynamics('i_measure_options')));
end
