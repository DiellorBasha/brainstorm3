function tests = test_operator_gate
tests = functiontests(localfunctions);
end
function test_constrained_scalar_source(tc)
% opVariants order: {'Laplace-Beltrami','LB-Connectome','Connection Laplacian','Dirac'}
% Constrained (scalar) supports the two scalar bases + the tangent operator (gradient of the scalar
% map), but NOT Dirac (no vectors).
m = panel_bst_dynamics('i_gate_mask', 1);
verifyEqual(tc, m, logical([1 1 1 0]));
end
function test_unconstrained_vector_source(tc)
% Unconstrained (vector) supports Dirac (native) + the two scalar bases (on the norm |J|), but NOT
% the tangent-only Connection Laplacian.
m = panel_bst_dynamics('i_gate_mask', 3);
verifyEqual(tc, m, logical([1 1 0 1]));
end
function test_unknown_defaults_permissive(tc)
m = panel_bst_dynamics('i_gate_mask', []);         % unknown -> permissive (all enabled)
verifyEqual(tc, m, logical([1 1 1 1]));
end
