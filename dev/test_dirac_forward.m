function tests = test_dirac_forward
tests = functiontests(localfunctions);
end
function test_allpass_roundtrip(tc)
% Toy: 1 hemi, a FULL quaternion basis (Phi = eye(4nV), 4nV modes, Mass = I) so
% project+reconstruct is an exact identity. With all-pass g, i_dirac_forward must round-trip
% the embedded source EXACTLY (Jfilt == J), and Dfilt has the right shape. The full-basis
% round-trip is a real check (not the tautological "Dfilt==Leig*cfilt"); physics is live-verified.
nV = 4;  K = 4*nV;  nT = 3;  nCh = 5;
ax = struct('Phi',{{eye(4*nV)}}, 'Mass',{{speye(4*nV)}}, 'Lambda',{{(1:K)'}}, 'GlobalVertices',{{(1:nV)'}});
J = randn(3*nV, nT);                       % source 3-vector field
Leig = randn(nCh, K);                      % eigenbasis leadfield (1 hemi -> K cols)
g = @(l) ones(size(l));                    % all-pass
[Dfilt, Jfilt, cfilt] = panel_bst_dynamics('i_dirac_forward', ax, Leig, J, g);
verifyEqual(tc, size(Dfilt), [nCh nT]);
verifyEqual(tc, size(cfilt), [K nT]);
verifyEqual(tc, size(Jfilt), [3*nV nT]);
verifyLessThan(tc, max(abs(Jfilt(:) - J(:))), 1e-12);      % all-pass full-basis round-trip is exact
end
function test_zero_gain_zeros(tc)
% g==0 -> cfilt=0 -> Jfilt=0, Dfilt=0 (the filter removes everything).
nV = 4;  K = 4*nV;  nT = 2;  nCh = 5;
ax = struct('Phi',{{eye(4*nV)}}, 'Mass',{{speye(4*nV)}}, 'Lambda',{{(1:K)'}}, 'GlobalVertices',{{(1:nV)'}});
[Dfilt, Jfilt] = panel_bst_dynamics('i_dirac_forward', ax, randn(nCh,K), randn(3*nV,nT), @(l) zeros(size(l)));
verifyLessThan(tc, max(abs(Jfilt(:))), 1e-12);
verifyLessThan(tc, max(abs(Dfilt(:))), 1e-12);
end
