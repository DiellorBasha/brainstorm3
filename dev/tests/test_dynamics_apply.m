% test_dynamics_apply - the atom filter-apply core (real source field -> Analysis -> filtered field)
nV=60; nT=40;
[Q,~]=qr(reshape(cos(1:(nV*nV)),nV,nV),0); Phi=Q; Lam=(linspace(0,5,nV)').^2; M=speye(nV);   % complete orthonormal basis
ax = struct('Phi',{{Phi}}, 'Lambda',{{Lam}}, 'Mass',{{M}}, 'GlobalVertices',{{(1:nV)'}});
F  = randn(nV, nT);
% low-pass (heat) must reduce energy and keep shape
Ff = panel_bst_dynamics('i_atom_filter_field', F, ax, 'Laplace-Beltrami', 'heat', struct('lmax',max(Lam),'t',0.2));
assert(isequal(size(Ff), size(F)), 'filtered field keeps shape');
assert(all(isfinite(Ff(:))), 'finite');
assert(norm(Ff,'fro') < norm(F,'fro'), 'low-pass reduces energy');
% flat (all-pass) ~ identity
Fa = panel_bst_dynamics('i_atom_filter_field', F, ax, 'Laplace-Beltrami', 'flat', struct('lmax',max(Lam)));
assert(norm(Fa - F,'fro')/norm(F,'fro') < 1e-6, 'flat kernel is ~identity');
% cursor window: a 4-sample window from a synthetic time vector + cursor seam
[iWin] = panel_bst_dynamics('i_cursor_window_test', (0:99)/100, 0.50, 0.04);  % Fs=100, cursor=0.50s, 0.04s window
assert(numel(iWin)==4 && iWin(1)==51 && iWin(end)==54, sprintf('window samples [%s]', num2str(iWin([1 end]))));
disp('OK');
