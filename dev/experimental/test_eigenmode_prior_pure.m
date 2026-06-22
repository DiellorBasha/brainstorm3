function test_eigenmode_prior_pure
% Verify the spectral prior R from eigenvalues:
%   - 'flat'  -> all ones
%   - 'power' -> normalized lambda^(-alpha), DC handled
%   - 'log'   -> positive, decreasing in lambda, gentle high-mode rolloff, max==1
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

% DC mode (lambda=0) is swapped to lambda(2), so modes 1 and 2 get equal prior
assert(abs(Rp(1) - Rp(2)) < 1e-12, 'DC swap: mode-1 and mode-2 prior should be equal.');

% log (2026): positive, decreasing, normalized; smoother modes favored
Rl = bst_eigenmode_prior(lambdas, K, 'log', 0);
assert(all(Rl > 0), 'log prior must be positive (lambda normalized into (0,1)).');
assert(abs(max(Rl) - 1) < 1e-12, 'log prior must be normalized to max 1.');
assert(all(diff(Rl) <= 1e-9), 'log prior must be non-increasing in lambda.');

% GBF-exact prior is INTENTIONALLY scale-dependent (millimetre scale is physical):
% the key fix is that high modes are NOT over-suppressed (gentle rolloff).
fprintf('log prior Rl = [%s]\n', sprintf('%.4f ', Rl));
assert(Rl(end) > 0.5, 'log prior must keep substantial variance on the edge mode (gentle rolloff).');

% Analysis-only kernels must be rejected as priors BY THEIR ADMISSIBILITY FLAG
% (not merely as "unknown") -- assert the specific error identifier.
eid = '';
try, bst_eigenmode_prior(lambdas, K, 'mexhat', 0); catch e, eid = e.identifier; end
assert(strcmp(eid, 'bst_eigenmode_prior:Inadmissible'), 'mexhat must be rejected as inadmissible.');
eid = '';
try, bst_eigenmode_prior(lambdas, K, 'dog', 0); catch e, eid = e.identifier; end
assert(strcmp(eid, 'bst_eigenmode_prior:Inadmissible'), 'dog must be rejected as inadmissible.');
% A registered-but-unwired admissible kernel is a distinct, clear error
eid = '';
try, bst_eigenmode_prior(lambdas, K, 'matern', 0); catch e, eid = e.identifier; end
assert(strcmp(eid, 'bst_eigenmode_prior:UnsupportedPrior'), 'matern must report unsupported-prior.');
% A truly unknown name is UnknownPrior
eid = '';
try, bst_eigenmode_prior(lambdas, K, 'no_such', 0); catch e, eid = e.identifier; end
assert(strcmp(eid, 'bst_eigenmode_prior:UnknownPrior'), 'unknown name must report unknown-prior.');

disp('ALL TESTS PASSED');
end
