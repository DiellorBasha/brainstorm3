function test_bst_eigenmodes_modekernel_pure
% Pure test: mode kernel equals the explicit projection, plus capping/no-kernel.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end
addpath(fullfile(repoRoot, 'toolbox', 'math'));

rng(7);
nV = 40; nModesAll = 12; nChan = 25; nT = 9;
Phi = randn(nV, nModesAll);
M   = spdiags(rand(nV,1) + 0.5, 0, nV, nV);   % SPD diagonal mass matrix
Eig = struct('Vectors', Phi);
K   = randn(nV, nChan);                         % imaging kernel [nV x nChan]
F   = randn(nChan, nT);                         % sensor data

% ---- Kernel case: ModeKernel*F == Phi(:,1:nModes)' * M * (K*F) ----
nModes = 7;
MK = bst_eigenmodes_modekernel(Eig, M, K, nModes);
assert(isequal(size(MK), [nModes, nChan]), 'ModeKernel size wrong');
D = MK*F - Phi(:,1:nModes)' * (M * (K * F));
assert(max(abs(D(:))) < 1e-9, 'kernel-folded projection mismatch');

% ---- No-kernel case: returns the projector Phi(:,1:nModes)' * M ----
P = bst_eigenmodes_modekernel(Eig, M, [], nModes);
assert(isequal(size(P), [nModes, nV]), 'projector size wrong');
Dp = P - Phi(:,1:nModes)' * M;
assert(max(abs(Dp(:))) < 1e-9, 'projector mismatch');

% ---- Capping: requesting more than available uses all; first-K consistency ----
MKall = bst_eigenmodes_modekernel(Eig, M, K, 999);
assert(size(MKall,1) == nModesAll, 'cap to available modes failed');
Dc = MKall(1:nModes,:) - MK;
assert(max(abs(Dc(:))) < 1e-12, 'first-K rows must match');

fprintf('ALL TESTS PASSED: test_bst_eigenmodes_modekernel_pure\n');
end
