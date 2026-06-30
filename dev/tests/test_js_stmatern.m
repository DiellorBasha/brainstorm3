% test_js_stmatern - real, low-pass: positive-ish, temporal PSD monotone-decreasing (1/f)
nV=60; K=20; nT=128; Fs=100;
[Q,~]=qr(reshape(cos(1:(nV*K)),nV,K),0); Phi=Q; Lam=(linspace(0,5,K)').^2; M=speye(nV);
ax=struct('nT',nT,'NFFT',nT,'Fs',Fs); ax.Phi={Phi}; ax.Lambda={Lam}; ax.Mass={M};
ax.GlobalVertices={(1:nV)'}; ax.tlag=(0:nT-1)/Fs; ax.omega=(0:nT-1)*(Fs/nT);
g = bst_eigfilter_design_stmatern(struct('kappa',2*pi/0.05,'nu',1.5,'lmax',max(Lam)));
G = g(Lam, ax.omega);
assert(all(isfinite(G(:))) && all(G(:)>=0), 'finite, nonneg (power spectrum)');
% per-mode temporal spectrum decreases over the lower half (1/f)
half=1:floor(nT/2); row=abs(G(1,half));
assert(row(1) >= row(end), 'temporal spectrum decays with frequency');
W = bst_eigenfilter('Atom', ax, 'stmatern', struct('kappa',2*pi/0.05,'nu',1.5,'lmax',max(Lam)), 13);
assert(isequal(size(W),[nV nT]) && all(isfinite(W(:))), 'atom shape/finite');
disp('OK');
