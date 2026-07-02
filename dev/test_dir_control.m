function tests = test_dir_control
tests = functiontests(localfunctions);
end
function test_control_spec(tc)
    s = panel_bst_dynamics('i_dir_control_spec', 'scalar');     verifyFalse(tc, s.show);
    a = panel_bst_dynamics('i_dir_control_spec', 'tangent');    verifyTrue(tc, a.show); verifyEqual(tc, a.type, 'angle');
    q = panel_bst_dynamics('i_dir_control_spec', 'quaternion'); verifyTrue(tc, q.show); verifyEqual(tc, q.type, 'preset');
    verifyTrue(tc, ismember('Normal', q.presets));
    verifyTrue(tc, ismember('Pick-on-surface', q.presets));
end
function test_preset_to_dir(tc)
    verifyEqual(tc, panel_bst_dynamics('i_preset_dir', '+Z'), [0 0 1]);
    verifyEqual(tc, panel_bst_dynamics('i_preset_dir', '+X'), [1 0 0]);
end
