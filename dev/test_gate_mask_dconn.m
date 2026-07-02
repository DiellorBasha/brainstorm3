function tests = test_gate_mask_dconn
tests = functiontests(localfunctions);
end
function test_masks(t)
    % Order: {Laplace-Beltrami, LB-Connectome, Connection Laplacian, Dirac, Dirac-Connectome}
    verifyEqual(t, panel_bst_dynamics('i_gate_mask', 1), logical([1 1 1 0 0]));  % constrained
    verifyEqual(t, panel_bst_dynamics('i_gate_mask', 3), logical([1 1 0 1 1]));  % unconstrained (Dirac + Dirac-Connectome)
    verifyEqual(t, panel_bst_dynamics('i_gate_mask', []), logical([1 1 1 1 1]));  % unknown -> permissive
end
