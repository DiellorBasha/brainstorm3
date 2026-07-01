function tests = test_operator_gate
tests = functiontests(localfunctions);
end
function test_scalar_source_disables_vector_ops(tc)
% opVariants order: {'Laplace-Beltrami','LB-Connectome','Connection Laplacian','Dirac'}
m = panel_bst_dynamics('i_gate_mask', 1);          % scalar source
verifyEqual(tc, m, logical([1 1 0 0]));            % LB/LB-Conn enabled; Connection/Dirac disabled
end
function test_vector_source_enables_all(tc)
m = panel_bst_dynamics('i_gate_mask', 3);          % unconstrained vector source
verifyEqual(tc, m, logical([1 1 1 1]));
end
function test_unknown_defaults_permissive(tc)
m = panel_bst_dynamics('i_gate_mask', []);         % unknown -> permissive (all enabled)
verifyEqual(tc, m, logical([1 1 1 1]));
end
