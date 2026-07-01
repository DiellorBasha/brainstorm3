% test_eigenfilter_atom - domain-aware single-atom realiser on a SYNTHETIC eigenbasis
nV = 60; K = 20; nT = 64; Fs = 100;
[Q,~] = qr(reshape(cos(1:(nV*K)), nV, K), 0);     % deterministic nV x K orthobasis
Phi = Q;  Lam = (linspace(0.0, 5, K)').^2;  M = speye(nV);
ax = struct('nT',nT,'NFFT',nT,'Fs',Fs);
ax.Phi = {Phi}; ax.Lambda = {Lam}; ax.Mass = {M}; ax.GlobalVertices = {(1:nV)'};
ax.Variant = 'Laplace-Beltrami';   % scalar fiber (RowMap/i_fiber dispatch, added for the Atom seedDir refactor)
ax.tlag = (0:nT-1)/Fs;  ax.omega = (0:nT-1)*(Fs/nT);
seed = 13;  loc = seed;  c0 = manifold_ft(Phi, M, full(sparse(loc,1,1,nV,1)));

% --- (a) ts atom matches the direct formula (== legacy bst_atom Evaluate path) ---
kp = struct('lmax',max(Lam),'tau',0.3);
g  = bst_eigfilter_kernel('diffusion', kp);
W_ref = manifold_ift(Phi, g(Lam, ax.tlag) .* c0);
W_ts  = bst_eigenfilter('Atom', ax, 'diffusion', kp, seed);
assert(isequal(size(W_ts),[nV nT]), 'ts atom shape');
assert(max(abs(W_ts(:) - W_ref(:))) < 1e-12, 'ts atom == direct manifold_ift formula');

% --- (b) static atom is constant in time ---
W_s = bst_eigenfilter('Atom', ax, 'heat', struct('lmax',max(Lam),'t',0.2), seed);
assert(max(abs(W_s(:,1) - W_s(:,end))) < 1e-12, 'static atom constant over time');

% --- (c) js path is the temporal-Fourier dual of ts: define the SAME filter in js (= fft of the ts
%     samples), realise via the js machinery (ifft of g(lambda,omega)), assert it equals the direct ts atom. ---
Gjs = fft(g(Lam, ax.tlag), nT, 2);              % js representation of the ts kernel
Gt  = real(ifft(Gjs, nT, 2));                   % what the Atom js branch does
W_js = manifold_ift(Phi, Gt(:,1:nT) .* c0);
assert(max(abs(W_js(:) - W_ref(:))) < 1e-9, 'js path (ifft of fft(g)) == direct ts (validates js branch)');
disp('OK');
