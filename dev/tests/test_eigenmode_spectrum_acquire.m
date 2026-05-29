function test_eigenmode_spectrum_acquire
% Verify CollapseProject: constrained passthrough + unconstrained rms-collapse,
% both projected onto an M-orthonormal eigenbasis.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end
addpath(fullfile(repoRoot, 'toolbox', 'gui'));

rng(0);
nV = 40; k = 6;
% Diagonal mass matrix and an M-orthonormal basis Phi (Phi' M Phi = I).
d = 0.5 + rand(nV, 1);
M = spdiags(d, 0, nV, nV);
A = randn(nV, k);
G = A' * (M * A);
R = chol((G + G') / 2);
Phi = A / R;
assert(norm(Phi' * (M * Phi) - eye(k), 'fro') < 1e-9, 'setup: Phi not M-orthonormal');
Em = struct('Vectors', Phi, 'Values', (1:k)', 'nModes', k, 'MassType', 'barycentric', ...
            'Component', [ones(3,1); 2*ones(3,1)], 'CompRank', [1;2;3;1;2;3], 'nComponents', 2);

% ---- Constrained (nComponents=1): a pure eigenmode field -> unit coeff at that mode ----
S = Phi(:, 3) * [1, 0, 5];                    % [nV x 3], column 1 is mode 3, scaled across time
Theta = view_eigenmode_spectrum('CollapseProject', Em, S, 1, M);
assert(isequal(size(Theta), [k 3]), 'constrained Theta wrong shape');
assert(abs(Theta(3,1) - 1) < 1e-9, 'pure mode-3 field must give unit coeff at mode 3');
assert(max(abs(Theta([1 2 4 5 6], 1))) < 1e-9, 'other modes must be ~0');
assert(max(abs(Theta(:,2))) < 1e-9, 'zero field column must give zero coeffs');
% Constrained path must equal a direct projection.
ThetaDirect = bst_eigenmodes_project(Em, S, M);
assert(max(abs(Theta(:) - ThetaDirect(:))) < 1e-9, 'constrained path must equal direct projection');

% ---- Unconstrained (nComponents=3): rms-collapse then project ----
A3 = randn(3*nV, 4);                           % [3nV x nTime]
expectedS = bst_source_orient([], 3, [], A3, 'rms');   % [nV x nTime]
expectedTheta = bst_eigenmodes_project(Em, expectedS, M);
ThetaU = view_eigenmode_spectrum('CollapseProject', Em, A3, 3, M);
assert(isequal(size(ThetaU), [k 4]), 'unconstrained Theta wrong shape');
assert(max(abs(ThetaU(:) - expectedTheta(:))) < 1e-9, 'unconstrained path must equal collapse-then-project');

fprintf('ALL TESTS PASSED: test_eigenmode_spectrum_acquire\n');
end
