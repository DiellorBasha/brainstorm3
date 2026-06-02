function test_inverse_eigenmodes_pure
% Verify the mode-space solver math (no DB), via the 'SolvePure' entry:
%   Kernel = SolvePure(L_tilde, lambdas, iW, Proj, Method, Prior, Alpha, SNR, Unreg)
%   - shapes [K x nCh]; finite
%   - projector folding: Kernel annihilates the projected-out direction
%   - harmonic limit: flat prior + Unreg reproduces pinv(iW*L_tilde)*iW
%   - dSPM rows are unit noise-normalized in whitened space
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

nCh = 8; K = 5;
L_tilde = zeros(nCh, K);
for i = 1:nCh
    for j = 1:K
        L_tilde(i,j) = sin(i * j * pi / 13);
    end
end
lambdas = (0:K-1)';                       % ascending, DC=0
iW   = diag(linspace(1, 2, nCh));         % whitener
Proj = eye(nCh);                          % no projector for base cases

% --- harmonic limit: flat + Unreg == pinv(iW*L)*iW ---
Kh = bst_inverse_eigenmodes('SolvePure', L_tilde, lambdas, iW, Proj, 'mne', 'flat', 0, 3, true);
Mref = pinv(iW * L_tilde) * iW;
assert(isequal(size(Kh), [K nCh]), 'Kernel must be [K x nCh].');
assert(max(abs(Kh(:) - Mref(:))) < 1e-8, 'flat+Unreg must equal pinv(iW*L)*iW.');

% --- regularized MNE with log prior: finite, correct shape ---
Km = bst_inverse_eigenmodes('SolvePure', L_tilde, lambdas, iW, Proj, 'mne', 'log', 0, 3, false);
assert(all(isfinite(Km(:))), 'Regularized MNE kernel must be finite.');

% --- projector folding: a rank-1 projector removes channel direction p ---
p = zeros(nCh,1); p(3) = 1;               % project out channel 3 subspace
Proj1 = eye(nCh) - (p*p')/(p'*p);
Kp = bst_inverse_eigenmodes('SolvePure', L_tilde, lambdas, iW, Proj1, 'mne', 'flat', 0, 3, false);
assert(max(abs(Kp * p)) < 1e-9, 'Kernel must annihilate the projected-out direction.');

% --- dSPM: rows unit-normalized by noise std in whitened space ---
Kd = bst_inverse_eigenmodes('SolvePure', L_tilde, lambdas, iW, Proj, 'dspm', 'flat', 0, 3, false);
Kd_white = Kd / iW;                        % undo final de-whitening
rn = sqrt(sum(Kd_white.^2, 2));
assert(max(abs(rn - 1)) < 1e-6, 'dSPM rows must have unit noise norm in whitened space.');

disp('ALL TESTS PASSED');
end
