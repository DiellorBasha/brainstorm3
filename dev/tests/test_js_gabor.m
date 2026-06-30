% test_js_gabor - gabor js atom: real, and temporal spectrum peaks at f0
nV=60; K=20; nT=100; Fs=100;
[Q,~]=qr(reshape(cos(1:(nV*K)),nV,K),0); Phi=Q; Lam=(linspace(0,5,K)').^2; M=speye(nV);
ax=struct('nT',nT,'NFFT',nT,'Fs',Fs); ax.Phi={Phi}; ax.Lambda={Lam}; ax.Mass={M};
ax.GlobalVertices={(1:nV)'}; ax.tlag=(0:nT-1)/Fs; ax.omega=(0:nT-1)*(Fs/nT);
f0=12; k0=2*pi/0.03; sf=2;                                   % 12 Hz, ~30 mm scale
[W,gv]=bst_eigenfilter('Atom', ax, 'gabor', struct('f0',f0,'k0',k0,'sf',sf,'lmax',max(Lam)), 13);
assert(isequal(size(W),[nV nT]), 'shape');
assert(all(isfinite(W(:))), 'finite');
% temporal spectrum of the most active vertex peaks near f0
[~,iv]=max(sum(W.^2,2)); P=abs(fft(W(iv,:))).^2; fb=(0:nT-1)*(Fs/nT);
half=2:floor(nT/2); [~,ip]=max(P(half)); fpk=fb(half(ip));
assert(abs(fpk-f0) <= Fs/nT*1.5, sprintf('peak %.1f near f0=%.1f', fpk, f0));
disp('OK');
