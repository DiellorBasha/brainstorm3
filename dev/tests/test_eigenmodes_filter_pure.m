function test_eigenmodes_filter_pure
% Verify the pure spectral filter (lowpass + heat kernel limits).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

rng(1);
n = 40; k = 10;
d = 0.5 + rand(n, 1);
M = spdiags(d, 0, n, n);
A = randn(n, k);
G = A' * (M * A);
R = chol((G + G') / 2);
Phi = A / R;

lambdas = (0:k-1)'.^2;   % eigenvalues; lambda(1)=0 is the DC mode
Em = struct('Vectors', Phi, 'Values', lambdas, 'nModes', k, 'MassType', 'barycentric');

c0 = randn(k, 1);
Data = Phi * c0;

% Lowpass cutoff 4 keeps modes 1..4 only.
LP = bst_eigenmodes_filter(Em, Data, M, 'lowpass', 'CutoffMode', 4);
assert(max(abs(LP - Phi(:, 1:4) * c0(1:4))) < 1e-9, 'lowpass mismatch.');

% Heat kernel with t -> 0 is the identity.
H0 = bst_eigenmodes_filter(Em, Data, M, 'heat', 'DiffusionTime', 1e-12);
assert(max(abs(H0 - Data)) < 1e-6, 'heat t->0 is not identity.');

% Heat kernel with large t isolates the DC mode (lambda=0 -> gain 1; rest -> 0).
HT = bst_eigenmodes_filter(Em, Data, M, 'heat', 'DiffusionTime', 1e6);
assert(max(abs(HT - Phi(:, 1) * c0(1))) < 1e-6, 'heat large t did not isolate DC.');

fprintf('ALL TESTS PASSED: test_eigenmodes_filter_pure\n');
end
