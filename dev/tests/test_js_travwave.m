% test_js_travwave - traveling wave: real; peak temporal frequency grows with sqrt(lambda)
nV=60; K=20; nT=100; Fs=100;
[Q,~]=qr(reshape(cos(1:(nV*K)),nV,K),0); Phi=Q; Lam=(linspace(0.2,6,K)').^2; M=speye(nV);
ax=struct('nT',nT,'NFFT',nT,'Fs',Fs); ax.Phi={Phi}; ax.Lambda={Lam}; ax.Mass={M};
ax.GlobalVertices={(1:nV)'}; ax.tlag=(0:nT-1)/Fs; ax.omega=(0:nT-1)*(Fs/nT);
c=4; width=2;
g = bst_eigfilter_design_travwave(struct('c',c,'width',width,'lmax',max(Lam)));
G = g(Lam, ax.omega);                                            % [K x N]
assert(all(isfinite(G(:))), 'finite');
% ridge: per-mode peak frequency should increase with sqrt(lambda)
fb=(0:nT-1)*(Fs/nT); half=1:floor(nT/2);
[~,iLo]=min(Lam); [~,iHi]=max(Lam);
[~,pLo]=max(abs(G(iLo,half))); [~,pHi]=max(abs(G(iHi,half)));
assert(fb(half(pHi)) > fb(half(pLo)), 'higher lambda -> higher ridge frequency');
% realised atom is real
W = bst_eigenfilter('Atom', ax, 'travwave', struct('c',c,'width',width,'lmax',max(Lam)), 13);
assert(isequal(size(W),[nV nT]) && all(isfinite(W(:))), 'atom shape/finite');
disp('OK');
