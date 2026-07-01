function tests = test_jtvatoms_localize
tests = functiontests(localfunctions);
end
function test_one_atom_per_band_at_peak(tc)
% Controller runs this with BST_TEST_SURF = a real cortex on the path (JTVAtoms reads vertex
% positions from the surface; the seed INDICES come from W). Skips if the surface is not set.
surf = getenv('BST_TEST_SURF');
assumeTrue(tc, ~isempty(surf));
nV = 60; nT = 4; M = 3;  W = zeros(nV,nT,M);
for m=1:M, W(10*m, m, m) = 1; end                    % band m peaks at vertex 10*m, time m
ax = struct('SurfaceFile', surf, 'TimeFile','', 'Time', (1:nT), 'tlag',(1:nT));
T = bst_eigenwavelet('JTVAtoms', W, ax, 0.5);
verifyEqual(tc, numel(T.Groups), M);                 % one atom per band
for m=1:M, verifyEqual(tc, T.Groups(m).vertices, 10*m); end   % at the expected peak vertex
end
