% test_atom_from_kernel - realise+threshold an atom into a populated group; threshold monotonicity
nV=60; K=20; nT=64; Fs=100;
[Q,~]=qr(reshape(cos(1:(nV*K)),nV,K),0); Phi=Q; Lam=(linspace(0,5,K)').^2; M=speye(nV);
ax=struct('nT',nT,'NFFT',nT,'Fs',Fs,'SurfaceFile','synthetic'); ax.Phi={Phi}; ax.Lambda={Lam}; ax.Mass={M};
ax.GlobalVertices={(1:nV)'}; ax.tlag=(0:nT-1)/Fs; ax.omega=(0:nT-1)*(Fs/nT);
kp = struct('lmax',max(Lam),'tau',0.3);
G = bst_dynamics('AtomFromKernel', ax, 'diffusion', kp, 13, 0.5);
assert(G.vertices==13 && strcmp(G.KernelName,'diffusion') && G.Threshold==0.5, 'generator fields');
assert(~isempty(G.region{1}) && numel(G.times)==2 && G.times(2)>=G.times(1), 'Scout + Event populated');
% threshold monotonicity: higher threshold -> Scout is a subset
G9 = bst_dynamics('AtomFromKernel', ax, 'diffusion', kp, 13, 0.9);
assert(all(ismember(G9.region{1}, G.region{1})), 'higher threshold -> subset Scout');
disp('OK');
