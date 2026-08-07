% TEST_TESS_LAPLACIAN: analytic right triangle + null space + oracle parity.
V = [0 0 0; 1 0 0; 0 1 0]; F = [1 2 3];
A = tess_laplacian(V, F);
% right angle at v1 (cot 0), 45 deg at v2,v3 (cot 1):
% w(edge v2v3) = cot(v1)/2 = 0 ; w(v1v3) = cot(v2)/2 = 0.5 ; w(v1v2) = cot(v3)/2 = 0.5
Aexp = [1 -0.5 -0.5; -0.5 0.5 0; -0.5 0 0.5];
assert(norm(full(A) - Aexp, 'fro') < 1e-14, 'right-triangle cotan weights wrong');
% --- sphere checks ---
[Vs, Fs] = tess_sphere(2562);
As = tess_laplacian(Vs, Fs);
assert(norm(As - As', 'fro') < 1e-14, 'stiffness must be symmetric');
assert(max(abs(As * ones(size(As,1),1))) < 1e-12, 'constants must be in the null space');
e = eigs(As, 3, 'smallestreal');
assert(min(e) > -1e-10, 'stiffness must be positive semidefinite');
% --- parity vs nxr oracle on the real cortex ---
S = load(fullfile(fileparts(mfilename('fullpath')), 'oracle_lbo_sub0002.mat'));
T = load(S.meta.SurfaceFile);
Vtx = double(T.Vertices); Fcs = double(T.Faces); nVtot = size(Vtx,1);
for hh = 1:2
    v = S.vH{hh};
    isV = false(nVtot,1); isV(v) = true;
    mapV = zeros(nVtot,1); mapV(v) = 1:numel(v);
    Ah = tess_laplacian(Vtx(v,:), mapV(Fcs(all(isV(Fcs),2),:)));
    relErr = norm(Ah - S.A{hh}, 'fro') / norm(S.A{hh}, 'fro');
    fprintf('hemi %d stiffness parity rel err: %g\n', hh, relErr);
    assert(relErr <= 1e-12, 'stiffness parity failed (hemi %d): %g', hh, relErr);
end
disp('test_tess_laplacian PASSED');
