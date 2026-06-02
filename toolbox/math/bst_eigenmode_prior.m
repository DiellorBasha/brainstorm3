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
%       'log'   : R proportional to log(lambda_ref/lambda_k)      (2026 GBF log prior)
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
        % Reference scale = first discarded eigenvalue (else just above max used)
        if nAvail >= K+1
            lamRef = lambdas(K+1);
        else
            lamRef = lam(K) * (1 + 1e-6);
        end
        lamRef = max(lamRef, lam(K) * (1 + 1e-12));   % ensure strictly > lam(K)
        lamTilde = lam(1:K) / lamRef;                  % in (0,1)
        lamTilde = min(lamTilde, 1 - 1e-12);
        R = -log(lamTilde);                            % = log(lamRef/lam) > 0, decreasing in lambda

    otherwise
        error('bst_eigenmode_prior:UnknownPrior', 'Unknown priorType: %s', priorType);
end

% Normalize to max 1
R = R / max(R);
end
