function R = bst_eigenmode_prior(lambdas, K, priorType, alpha)
% BST_EIGENMODE_PRIOR: Diagonal source-covariance prior R from LBO eigenvalues.
%
% USAGE:  R = bst_eigenmode_prior(lambdas, K, priorType, alpha)
%
% DESCRIPTION:
%     Returns the diagonal of the source-covariance prior R [K x 1] in eigenmode
%     space. R plays the role of the source covariance in the standard inverse
%     J = R*L'*(L*R*L' + lambda*C)^-1*d, replacing depth weighting. Larger R_k =
%     more prior variance for mode k.
%
%     priorType:
%       'flat'  : R = ones(K,1)                         (no spectral prior)
%       'power' : R proportional to lambda_k^(-alpha)   (legacy 1/f-like)
%       'log'   : R = -log(lambda_mm), GBF millimetre-scale eigenvalues (gentle high-mode rolloff)
%
%     For 'log', eigenvalues are normalized into (0,1) by lambda_ref (the first
%     discarded eigenvalue, i.e. lambda(K+1), else lambda(K)*(1+eps)). R_k =
%     log(lambda_ref/lambda_k) is positive and decreasing in lambda (smoother
%     modes get more prior variance). Scaling all eigenvalues by a constant c
%     only shifts log-space uniformly; after max-normalization R is unchanged
%     (scale invariance). The DC mode (lambda~0) is swapped to lambda(2).
%     R is normalized so max(R) = 1; absolute scale is absorbed by the global
%     regularizer in the inverse.
%
% Authors: Diellor Basha, 2026

lambdas = double(lambdas(:));
nAvail  = numel(lambdas);
K = max(1, min(K, nAvail));

% Need a usable spectrum: a single zero eigenvalue cannot be normalized/swapped.
if nAvail < 1 || (nAvail == 1 && lambdas(1) <= 0)
    error('bst_eigenmode_prior:InvalidInput', ...
        'lambdas must contain at least one positive eigenvalue.');
end

% DC handling: replace a (near-)zero leading eigenvalue with the next one
lam = lambdas;
if lam(1) <= max(lam) * 1e-12 && nAvail >= 2
    lam(1) = lam(2);
end

switch lower(priorType)
    case 'flat'
        R = ones(K, 1);
        return;

    case 'power'
        lamK = lam(1:K);
        lamK = max(lamK, max(lamK) * 1e-12);
        R = lamK .^ (-alpha);

    case 'log'
        % GBF 2026 prior on RAW millimetre-scale eigenvalues. Brainstorm surfaces
        % are in metres and the cotan stiffness is scale-invariant (only the mass
        % matrix scales with area), so lambda_mm = lambda_m * 1e-6 reproduces GBF's
        % millimetre mesh exactly (eigenvectors unchanged). With lambda_mm in (0,1),
        % R = -log(lambda_mm) > 0 and decreases gently with mode index (a large
        % additive offset, ~ -log(1e-6) ≈ 13.8, so high modes are softened, not
        % annihilated). The global regulariser in the inverse absorbs absolute scale.
        M2MM2 = 1e-6;                                  % m^-2 -> mm^-2 (coords scale 1e3)
        lamMM = lam(1:K) * M2MM2;
        lamMM = min(max(lamMM, eps), 1 - 1e-12);       % keep strictly in (0,1)
        R = -log(lamMM);                               % > 0, gently decreasing

    otherwise
        error('bst_eigenmode_prior:UnknownPrior', 'Unknown priorType: %s', priorType);
end

% Normalize to max 1
R = R / max(R);
end
