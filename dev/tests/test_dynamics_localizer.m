% test_dynamics_localizer - Get/Set round-trip on a group, via bst_dynamics (moved from bst_atom)
G = bst_dynamics('NewGroup','t');
loc = bst_dynamics('NewLoc','time');  loc.center = 0.4; loc.extent = 0.1;
G = bst_dynamics('Set', G, 'time', 1, loc);
out = bst_dynamics('Get', G, 'time', 1);
assert(abs(out.center - 0.4) < 1e-12 && abs(out.extent - 0.1) < 1e-12, 'time loc round-trip');
assert(strcmp(out.state,'window'), 'extent>0 => window');
M = bst_dynamics('AxisMeta');
assert(any(strcmp({M.name},'source')) && any(strcmp({M.name},'scale')), 'axis metadata');
disp('OK');
