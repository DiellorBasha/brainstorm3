% test_dynamics_operator - per-atom operator + launch-derived default
kp = struct('lmax',40,'tau',0.3,'vals',[400 0 0]);
G  = panel_bst_dynamics('i_default_atom', 'diffusion', kp, 13, 's.mat', 'atom1', 'Dirac');
assert(strcmp(G.Operator,'Dirac'), 'atom carries its operator');
s = panel_bst_dynamics('i_atom_detail', G);  assert(contains(s,'Dirac'), 'detail shows operator');
% operator survives AddGroup (must be a template field, not a dropped extra)
T = bst_dynamics('AddGroup', bst_dynamics('New','t'), G);
assert(strcmp(T.Groups(1).Operator,'Dirac'), 'operator persists through AddGroup');
% launch-derived default from the source-result comment (test seam: st.srcComment)
stD = struct('srcComment','MN: shared dirac kernel');
assert(strcmp(panel_bst_dynamics('i_launch_operator', stD), 'Dirac'), 'dirac comment -> Dirac');
stL = struct('srcComment','MN: 2018 (Constr) 2018');
assert(strcmp(panel_bst_dynamics('i_launch_operator', stL), 'Laplace-Beltrami'), 'else -> LBO');
% default operator when omitted = Laplace-Beltrami (5-arg back-compat)
G2 = panel_bst_dynamics('i_default_atom', 'diffusion', kp, 7, 's.mat', 'atom2');
assert(strcmp(G2.Operator,'Laplace-Beltrami'), 'omitted operator defaults to LBO');
disp('OK');
