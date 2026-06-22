function test_eigenmodes_project_pure
% Verify the pure projection recovers M-orthonormal coefficients exactly.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

rng(0);
n = 50; k = 8;
d = 0.5 + rand(n, 1);
M = spdiags(d, 0, n, n);

% Build an M-orthonormal basis: Phi = A * inv(chol(A'MA)) -> Phi'M Phi = I.
A = randn(n, k);
G = A' * (M * A);
R = chol((G + G') / 2);
Phi = A / R;
assert(norm(Phi' * (M * Phi) - eye(k), 'fro') < 1e-9, 'setup: Phi not M-orthonormal');

Em = struct('Vectors', Phi, 'Values', (1:k)', 'nModes', k, 'MassType', 'barycentric');

% Known coefficients -> data -> recover exactly (proves M weighting is applied).
c0 = randn(k, 3);
Data = Phi * c0;
[C, Recon] = bst_eigenmodes_project(Em, Data, M);
assert(max(abs(C(:) - c0(:))) < 1e-9, 'Projection did not recover coefficients (M weighting bug?).');
assert(max(abs(Recon(:) - Data(:))) < 1e-9, 'Full reconstruction not exact.');

% Parseval in the M-norm.
uMnorm2 = sum(sum(Data .* (M * Data)));
assert(abs(uMnorm2 - sum(C(:).^2)) < 1e-7, 'Parseval identity failed.');

% ModeRange reconstruction uses only the selected modes.
[~, ReconLP] = bst_eigenmodes_project(Em, Data, M, 'ModeRange', [1 3]);
expected = Phi(:, 1:3) * C(1:3, :);
assert(max(abs(ReconLP(:) - expected(:))) < 1e-9, 'ModeRange reconstruction wrong.');

fprintf('ALL TESTS PASSED: test_eigenmodes_project_pure\n');
end
