function test_eigenmodes_transform_pure
% Verify the unregularized SVD pseudoinverse transform: exact recovery in the
% noise-free case, and the left/right-inverse invariants in both regimes.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

rng(7);

% ----- Overdetermined regime: K < nch (well-determined least squares) -----
nCh = 40; nVert = 200; K = 25;
Gain    = randn(nCh, nVert);
Phi     = orth(randn(nVert, K));        % nVert x K, full column rank
L_tilde = Gain * Phi;

c0 = randn(K, 4);                        % known coefficients
D  = L_tilde * c0;                       % noise-free data
[Kernel, Info] = bst_eigenmodes_transform(Gain, Phi);

assert(isequal(size(Kernel), [K, nCh]), 'Kernel must be [K x nch].');
Theta = Kernel * D;
assert(max(abs(Theta(:) - c0(:))) < 1e-8, 'Did not recover coefficients (overdetermined).');
assert(norm(Kernel * L_tilde - eye(K), 'fro') < 1e-8, 'Left-inverse invariant Kernel*L=I_K failed.');
assert(Info.Rank == K, 'Rank should equal K when K<nch and full rank.');
assert(isfinite(Info.ConditionNumber) && Info.ConditionNumber >= 1, 'Bad condition number.');
assert(numel(Info.SingularValues) == K, 'SingularValues length mismatch.');

% Custom Tol floors more singular values -> lower effective rank.
sAll   = svd(L_tilde);                       % descending
TolMid = (sAll(9) + sAll(10)) / 2;           % keep exactly the first 9 modes
[~, InfoTol] = bst_eigenmodes_transform(Gain, Phi, 'Tol', TolMid);
assert(InfoTol.Rank == 9, 'Custom Tol did not reduce rank as expected (got %d).', InfoTol.Rank);

% Tol above the largest singular value -> rank 0, condition number Inf.
[~, InfoZero] = bst_eigenmodes_transform(Gain, Phi, 'Tol', sAll(1) * 2);
assert(InfoZero.Rank == 0 && isinf(InfoZero.ConditionNumber), 'Rank-0 path: expected Inf condition number.');

% Invalid Tol must error.
threw = false;
try
    bst_eigenmodes_transform(Gain, Phi, 'Tol', -1);
catch
    threw = true;
end
assert(threw, 'Negative Tol should error.');

% ----- Underdetermined regime: K > nch (min-norm right inverse) -----
nCh2 = 20; nVert2 = 200; K2 = 35;
Gain2    = randn(nCh2, nVert2);
Phi2     = orth(randn(nVert2, K2));
L_tilde2 = Gain2 * Phi2;
[Kernel2, Info2] = bst_eigenmodes_transform(Gain2, Phi2);

assert(isequal(size(Kernel2), [K2, nCh2]), 'Kernel must be [K x nch] (underdetermined).');
assert(norm(L_tilde2 * Kernel2 - eye(nCh2), 'fro') < 1e-8, 'Right-inverse invariant L*Kernel=I_nch failed.');
assert(Info2.Rank == nCh2, 'Rank should equal nch in the underdetermined regime.');

fprintf('ALL TESTS PASSED: test_eigenmodes_transform_pure\n');
end
