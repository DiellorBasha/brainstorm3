% test_js_resonator - real; oscillates at f0; higher Q -> longer-lasting envelope
nV=60; K=20; nT=200; Fs=100;
[Qr,~]=qr(reshape(cos(1:(nV*K)),nV,K),0); Phi=Qr; Lam=(linspace(0,5,K)').^2; M=speye(nV);
ax=struct('nT',nT,'NFFT',nT,'Fs',Fs); ax.Phi={Phi}; ax.Lambda={Lam}; ax.Mass={M};
ax.GlobalVertices={(1:nV)'}; ax.tlag=(0:nT-1)/Fs; ax.omega=(0:nT-1)*(Fs/nT);
f0=10;
Wlo = bst_eigenfilter('Atom', ax, 'resonator', struct('f0',f0,'Q',3, 'lmax',max(Lam)), 13);
Whi = bst_eigenfilter('Atom', ax, 'resonator', struct('f0',f0,'Q',12,'lmax',max(Lam)), 13);
assert(isequal(size(Wlo),[nV nT]) && all(isfinite(Wlo(:))), 'shape/finite');
[~,iv]=max(sum(Whi.^2,2)); fb=(0:nT-1)*(Fs/nT); half=2:floor(nT/2);
P=abs(fft(Whi(iv,:))).^2;
[~,ip]=max(P(half)); assert(abs(fb(half(ip))-f0) <= Fs/nT*2, 'oscillates at f0');
% envelope persistence: late-window energy fraction larger for high Q
lateE = @(W) sum(sum(W(:,round(0.6*nT):end).^2)) / max(sum(W(:).^2),eps);
assert(lateE(Whi) > lateE(Wlo), 'higher Q -> more late-window energy');
disp('OK');
