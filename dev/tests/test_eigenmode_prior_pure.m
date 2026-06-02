function test_eigenmode_prior_pure
% Verify the spectral prior R from eigenvalues:
%   - 'flat'  -> all ones
%   - 'power' -> normalized lambda^(-alpha), DC handled
%   - 'log'   -> positive, decreasing in lambda, ratio-preserving, max==1
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

lambdas = [0; 1; 4; 9; 16; 25];     % includes DC=0; ascending
K = 5;

% flat
Rf = bst_eigenmode_prior(lambdas, K, 'flat', 0);
assert(isequal(size(Rf), [K 1]), 'R must be [K x 1].');
assert(all(Rf == 1), 'flat prior must be all ones.');

% power (alpha=1): proportional to lambda^-1 after DC swap, normalized to max 1
Rp = bst_eigenmode_prior(lambdas, K, 'power', 1);
assert(all(Rp > 0), 'power prior must be positive.');
assert(abs(max(Rp) - 1) < 1e-12, 'prior must be normalized to max 1.');
assert(all(diff(Rp) <= 1e-12), 'power prior must be non-increasing in lambda.');

% log (2026): positive, decreasing, normalized; smoother modes favored
Rl = bst_eigenmode_prior(lambdas, K, 'log', 0);
assert(all(Rl > 0), 'log prior must be positive (lambda normalized into (0,1)).');
assert(abs(max(Rl) - 1) < 1e-12, 'log prior must be normalized to max 1.');
assert(all(diff(Rl) <= 1e-9), 'log prior must be non-increasing in lambda.');

% ratio-preservation invariant: scaling all eigenvalues by c only shifts log-space,
% so the resulting log prior is unchanged after max-normalization.
Rl_scaled = bst_eigenmode_prior(lambdas * 7.3, K, 'log', 0);
assert(max(abs(Rl - Rl_scaled)) < 1e-9, 'log prior must be invariant to uniform eigenvalue scaling.');

disp('ALL TESTS PASSED');
end
