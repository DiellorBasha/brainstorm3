function test_eigenmodes_harmonic_pure
% Verify the harmonic eigenmode kernel M = pinv(iW*L*Phi)*iW:
%   - reconstructs the whitened compressed system (Phi*M is the min-norm-in-mode-space inverse)
%   - is rank-safe when L*Phi is rank-deficient (no blow-up)
%   - equals bst_eigenmodes_transform(iW*L,Phi)*iW
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% Synthetic problem: 8 channels, 12 vertices, 5 modes (deterministic, no rng).
% L uses sin(i*j*pi/13) so rows are linearly independent (full rank over K=5).
% Phi uses DCT-II basis (normalized) so L*Phi is guaranteed full column rank.
nCh = 8; nV = 12; K = 5;
L = zeros(nCh, nV);
for i = 1:nCh
    for j = 1:nV
        L(i,j) = sin(i * j * pi / 13);
    end
end
Phi = zeros(nV, K);
for k = 1:K
    for n = 1:nV
        Phi(n,k) = cos(pi/nV * (n-0.5) * (k-1));
    end
end
Phi = Phi ./ vecnorm(Phi);                         % normalize columns
iW  = diag(linspace(1, 2, nCh));                  % deterministic whitener [nCh x nCh]

M = bst_eigenmodes_harmonic(L, Phi, iW);          % [K x nCh]
assert(isequal(size(M), [K, nCh]), 'M must be [K x nCh].');

% Equals the explicit construction via the transform
[Kt, ~] = bst_eigenmodes_transform(iW*L, Phi);    % pinv(iW*L*Phi) [K x nCh]
Mref = Kt * iW;
assert(max(abs(M(:) - Mref(:))) < 1e-9, 'M must equal pinv(iW*L*Phi)*iW.');

% M maps RAW sensor data to coefficients, so M*(L*Phi) ~ I_K (full column rank):
%   M*(L*Phi) = pinv(iW*L*Phi)*iW*(L*Phi) = pinv(iW*L*Phi)*(iW*L*Phi) = I_K
ID = M * (L * Phi);                                % should be ~ I_K
assert(max(abs(ID(:) - reshape(eye(K),[],1))) < 1e-6, 'M*(L*Phi) must be ~ I_K for full-rank case.');

% Rank-safety: duplicate two mode columns so L*Phi is rank-deficient -> no Inf/NaN
Phi2 = Phi; Phi2(:,5) = Phi2(:,1);
M2 = bst_eigenmodes_harmonic(L, Phi2, iW);
assert(all(isfinite(M2(:))), 'Rank-deficient input must not produce Inf/NaN (rank-safe pinv).');

disp('ALL TESTS PASSED');
end
