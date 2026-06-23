function test_bst_eigs_smallest()
% TEST_BST_EIGS_SMALLEST: the smallest-mode singular-pencil eigensolver, exercised
% end-to-end through tess_eigen (its sole caller and host).
%
% bst_eigs_smallest exists to return the K smallest generalized eigenpairs of a
% (near-)singular symmetric/Hermitian pencil (A,B) WITHOUT triggering MATLAB's RCOND
% "matrix close to singular" warning that the naive eigs(A,B,k,'smallestabs') shift-
% invert raises on the singular A. It is now a LOCAL function of tess_eigen (not a
% standalone callable), so this test drives it through the public entry point: it loads
% a real operator pencil, confirms the legacy 'smallestabs' path warns on it (the
% premise), then runs tess_eigen and asserts the solve
%   (1) emits NO warning,
%   (2) returns a real, ascending spectrum that agrees with the legacy eigenvalues,
%   (3) returns B-orthonormal modes.
%
% Environment-dependent: needs a loaded protocol whose Subject 1 has a cortex surface.
% SKIPs cleanly otherwise.
%
% Author: Diellor Basha, 2026
    PI = bst_get('ProtocolInfo');
    if isempty(PI) || ~isfield(PI,'SUBJECTS') || isempty(PI.SUBJECTS)
        disp('SKIP: test_bst_eigs_smallest -- no protocol loaded'); return;
    end
    sSubject = bst_get('Subject', 1);
    if isempty(sSubject) || isempty(sSubject.Surface)
        disp('SKIP: test_bst_eigs_smallest -- Subject 1 has no surfaces'); return;
    end
    % Resolve a cortex surface (default cortex if tagged, else first Cortex-type surface).
    iCx = sSubject.iCortex;
    if isempty(iCx)
        iCx = find(strcmpi({sSubject.Surface.SurfaceType}, 'Cortex'), 1);
    end
    if isempty(iCx)
        disp('SKIP: test_bst_eigs_smallest -- no cortex surface on Subject 1'); return;
    end
    SurfaceFile = sSubject.Surface(iCx).FileName;

    % Laplace-Beltrami: a real-symmetric singular pencil (constant null mode), the
    % cleanest case that bst_eigs_smallest is built to solve without warning.
    K = 12;
    Op = tess_operators(SurfaceFile, 'Laplace-Beltrami', 'NoSave', true);
    A = Op.Operator{1};  B = Op.Mass{1};      % left hemisphere pencil

    % --- Premise: the legacy 'smallestabs' shift-invert warns on this singular A. ---
    %     (Environment-robust: if it does NOT warn here, the no-warn contrast cannot be
    %      exercised, but the correctness + spectrum-agreement assertions below still run.)
    opts = struct('tol', 1e-6, 'maxit', 1000, 'disp', 0);
    lastwarn('');
    lam_legacy = sort(real(eigs(A, B, K, 'smallestabs', opts)));
    [~, wid_legacy] = lastwarn;
    legacyWarned = ~isempty(wid_legacy);
    if ~legacyWarned
        disp('NOTE: legacy ''smallestabs'' did not warn in this environment; testing correctness only.');
    end

    % --- tess_eigen runs bst_eigs_smallest internally; it must NOT warn. ---
    lastwarn('');
    Eig = tess_eigen(SurfaceFile, 'Laplace-Beltrami', 'K', K, 'NoSave', true);
    [wmsg_new, wid_new] = lastwarn;
    assert(isempty(wid_new), 'tess_eigen solve emitted a warning [%s] %s', wid_new, wmsg_new);
    if legacyWarned
        fprintf('PASS: legacy ''smallestabs'' warned (id=%s), tess_eigen solve was clean\n', wid_legacy);
    end

    lam_new = Eig.Lambda{1}(:);

    % --- (2) real, ascending spectrum, agreeing with the legacy eigenvalues ---
    assert(isreal(lam_new), 'eigenvalues not real (max|imag|=%.3e)', max(abs(imag(lam_new))));
    assert(all(diff(lam_new) >= -1e-8), 'eigenvalues not ascending');
    sc = max(1, max(abs(lam_legacy)));
    dspec = max(abs(lam_legacy - lam_new(1:K)));
    assert(dspec < 1e-3 * sc, 'spectrum mismatch %.3e (rel %.1e)', dspec, dspec / sc);
    fprintf('PASS: spectrum agrees with legacy eigs (max abs diff %.2e, rel %.1e)\n', dspec, dspec / sc);

    % --- (3) B-orthonormal modes ---
    Phi = Eig.Phi{1}(:, 1:K);
    G = Phi' * (B * Phi);
    ortho = norm(G - eye(K), 'fro') / sqrt(K);
    assert(ortho < 1e-8, 'modes not B-orthonormal (||Phi''BPhi - I||/sqrt(K) = %.3e)', ortho);
    fprintf('PASS: modes B-orthonormal (||Phi''BPhi - I||/sqrt(K) = %.2e)\n', ortho);

    disp('PASS: test_bst_eigs_smallest (end-to-end via tess_eigen)');
end
