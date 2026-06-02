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
%     The raw spectral shape is delegated to the eigfilter kernel library
%     (bst_eigfilter_kernel / bst_eigfilter_evaluate: 'flat'/'power'/'log'); this
%     function keeps the prior's own conventions (DC-mode swap, the 'log' millimetre
%     rescaling, max(R)=1 normalization, and an admissibility check).
%
%     For 'log', eigenvalues are taken on GBF's millimetre scale: Brainstorm
%     surfaces are in metres and the cotan stiffness is scale-invariant (only the
%     mass matrix scales with area), so lambda_mm = lambda_m * 1e-6 reproduces a
%     millimetre mesh (eigenvectors unchanged). With lambda_mm in (0,1),
%     R = -log(lambda_mm) > 0 and decreases gently with mode index (a large
%     additive offset ~ -log(1e-6) ~ 13.8 softens, rather than annihilates, the
%     high modes). The DC mode (lambda~0) is swapped to lambda(2). R is then
%     normalized so max(R) = 1; the inverse's global regulariser absorbs absolute
%     scale. Unlike a reference-normalized prior, this is intentionally
%     scale-dependent (the millimetre scale is physical).
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
        g = bst_eigfilter_kernel('flat');
        R = bst_eigfilter_evaluate(g, lam(1:K));      % ones
    case 'power'
        lamK = lam(1:K);
        lamK = max(lamK, max(lamK) * 1e-12);
        g = bst_eigfilter_kernel('power', struct('alpha', alpha));
        R = bst_eigfilter_evaluate(g, lamK);
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
        g = bst_eigfilter_kernel('log');
        R = bst_eigfilter_evaluate(g, lamMM);
    otherwise
        % Not one of the wired prior shapes. Use the kernel registry to give a
        % precise reason: a registered analysis-only kernel (e.g. mexhat/dog) is
        % rejected by its admissibility flag; a registered admissible kernel is
        % simply not wired here yet (only flat/power/log carry the prior-specific
        % scaling and parameterization); an unregistered name is unknown.
        if any(strcmp(bst_eigfilter_kernel('list'), lower(priorType)))
            minfo = bst_eigfilter_kernel('info', lower(priorType));
            if isfield(minfo,'priorAdmissible') && ~minfo.priorAdmissible
                error('bst_eigenmode_prior:Inadmissible', ...
                    'Kernel ''%s'' is not admissible as a prior (analysis-only; zero/negative spectral density).', priorType);
            end
            error('bst_eigenmode_prior:UnsupportedPrior', ...
                'Kernel ''%s'' is registered but not wired as a prior shape (use flat, power, or log).', priorType);
        end
        error('bst_eigenmode_prior:UnknownPrior', 'Unknown priorType: %s', priorType);
end

% Defensive backstop: a prior is a covariance, so it must be finite and >= 0.
if any(~isfinite(R)) || any(R < 0)
    error('bst_eigenmode_prior:Inadmissible', ...
        'Prior ''%s'' produced an invalid (non-finite or negative) covariance.', priorType);
end

% Normalize to max 1 (flat is already ones; this leaves it unchanged)
R = R / max(R);
end
