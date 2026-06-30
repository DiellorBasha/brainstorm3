% test_dynamics_filter_atom - a created atom is a FILTER (generator set, threshold/scout/event unset)
kp = struct('lmax',40,'tau',0.3);
G  = panel_bst_dynamics('i_default_atom', 'diffusion', kp, 13, 'surf.mat', 'atom1');
assert(strcmp(G.label,'atom1') && strcmp(G.KernelName,'diffusion'), 'label + kernel');
assert(isequal(G.KernelParams, kp) && isequal(G.vertices,13), 'generator carried');
assert(isempty(G.Threshold) && isempty(G.region) && isempty(G.times), 'NOT thresholded (filter, not marker)');
assert(strcmp(G.SurfaceFile,'surf.mat'), 'surface provenance');
s = panel_bst_dynamics('i_atom_detail', G);
assert(ischar(s) && contains(s,'diffusion') && contains(s,'13'), 'detail shows kernel + seed');
assert(contains(s,'tau'), 'detail shows a param (not lmax)');
disp('OK');
