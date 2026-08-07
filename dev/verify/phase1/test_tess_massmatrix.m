% TEST_TESS_MASSMATRIX: analytic single-triangle + sphere-area + oracle parity.
% --- analytic: unit right triangle, Area = 0.5 ---
V = [0 0 0; 1 0 0; 0 1 0]; F = [1 2 3];
B = tess_massmatrix(V, F);
Bexp = 0.5 * [1/6 1/12 1/12; 1/12 1/6 1/12; 1/12 1/12 1/6];
assert(norm(full(B) - Bexp, 'fro') < 1e-15, 'single-triangle Galerkin entries wrong');
% --- sphere: total mass = surface area ---
[Vs, Fs] = tess_sphere(2562);          % unit sphere approximation
Bs = tess_massmatrix(Vs, Fs);
v1 = Vs(Fs(:,1),:); v2 = Vs(Fs(:,2),:); v3 = Vs(Fs(:,3),:);
totalArea = sum(0.5 * sqrt(sum(cross(v2-v1, v3-v1, 2).^2, 2)));
assert(abs(full(sum(Bs(:))) - totalArea) < 1e-12 * totalArea, 'total mass ~= mesh area');
assert(norm(Bs - Bs', 'fro') == 0, 'mass must be exactly symmetric');
% --- parity vs nxr oracle on the real cortex ---
S = load(fullfile(fileparts(mfilename('fullpath')), 'oracle_lbo_sub0002.mat'));
T = load(S.meta.SurfaceFile);
Vtx = double(T.Vertices); Fcs = double(T.Faces); nVtot = size(Vtx,1);
for hh = 1:2
    v = S.vH{hh};
    isV = false(nVtot,1); isV(v) = true;
    mapV = zeros(nVtot,1); mapV(v) = 1:numel(v);
    Bh = tess_massmatrix(Vtx(v,:), mapV(Fcs(all(isV(Fcs),2),:)));
    relErr = norm(Bh - S.B{hh}, 'fro') / norm(S.B{hh}, 'fro');
    fprintf('hemi %d mass parity rel err: %g\n', hh, relErr);
    assert(relErr <= 1e-12, 'mass parity failed (hemi %d): %g', hh, relErr);
end
disp('test_tess_massmatrix PASSED');
