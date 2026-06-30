% test_atomgroup_generator - atom group carries the kernel generator, round-trips through Save/Load
G = bst_dynamics('NewGroup','t');
assert(isfield(G,'KernelName') && isfield(G,'KernelParams') && isfield(G,'Threshold'), 'generator fields present');
G.KernelName = 'gabor';  G.KernelParams = struct('f0',10,'k0',209,'sf',2);  G.Threshold = 0.5;
T = bst_dynamics('New','gen'); T.SurfaceFile = 'x'; T = bst_dynamics('AddGroup', T, G);
tmp = [tempname '.mat'];  bst_dynamics('Save', tmp, T);  T2 = bst_dynamics('Load', tmp);  delete(tmp);
g2 = T2.Groups(1);
assert(strcmp(g2.KernelName,'gabor') && g2.KernelParams.f0==10 && g2.Threshold==0.5, 'generator round-trips');
disp('OK');
