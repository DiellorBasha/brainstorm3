function tests = test_jtvatoms_localize
tests = functiontests(localfunctions);
end
function test_one_atom_per_band_at_peak(tc)
% Synthetic W [nV x nT x M]: band m peaks at vertex 10*m, time m.
nV = 60; nT = 4; M = 3;  W = zeros(nV,nT,M);
for m=1:M, W(10*m, m, m) = 1; end
ax = struct('SurfaceFile','', 'TimeFile','', 'Time', (1:nT), 'tlag',(1:nT));
T = bst_eigenwavelet('JTVAtoms', W, ax, 0.5);
verifyEqual(tc, numel(T.Groups), M);
for m=1:M, verifyEqual(tc, T.Groups(m).vertices, 10*m); end
end
