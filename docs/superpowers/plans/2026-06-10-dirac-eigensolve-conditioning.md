# Dirac Eigensolve Conditioning Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the `sigma=0` `'smallestabs'` shift-invert at `tess_eigen.m:217` (which factorizes the singular operator and warns `RCOND≈1e-18`) with a symmetric, scale-aware negative-shift solve, shared across the LBO / Connection / Dirac variants.

**Architecture:** Extract a standalone numerical helper `bst_eigs_smallest(A,B,k,opts)` that (1) symmetrizes/Hermitian-izes the pencil to force `eigs` onto the Lanczos path, and (2) selects the k smallest generalized eigenpairs via a small negative sigma shift `σ=−1e-7·λ_max`, where `λ_max` comes from a factorization-free `'largestabs'` estimate. `tess_eigen` calls it in place of the bare `eigs`. The returned `(V,D)` have identical meaning to today's call, so all downstream normalization/Rayleigh–Ritz and the stored basis are unchanged.

**Tech Stack:** MATLAB (Brainstorm toolbox), `eigs`, MATLAB MCP for running tests (`mcp__plugin_brainstorm-dev_MATLAB__run_matlab_file` / `evaluate_matlab_code`), `checkcode` for lint.

---

## Spec → plan note

This plan refines the approved spec in one way: the helper is a **standalone function** `toolbox/math/bst_eigs_smallest.m` rather than a file-local subfunction of `tess_eigen.m`. The logic is identical to the spec's `local_eigs_smallest`; the standalone form is unit-testable with a synthetic singular pencil (no nxr-compute, no surface). Everything else matches the spec exactly. Issue B remains out of scope.

## File structure

- **Create** `toolbox/math/bst_eigs_smallest.m` — the robust smallest-eigenpair solver (one responsibility: symmetrize + scale-aware shifted `eigs`).
- **Create** `dev/tests/test_bst_eigs_smallest.m` — portable synthetic unit test (singular symmetric pencil with a tiny mass mimicking meters scale) + an optional real-operator check.
- **Modify** `toolbox/anatomy/tess_eigen.m:217` — call `bst_eigs_smallest` instead of the bare `eigs`.

---

## Task 1: Standalone solver `bst_eigs_smallest` (TDD)

**Files:**
- Create: `toolbox/math/bst_eigs_smallest.m`
- Test: `dev/tests/test_bst_eigs_smallest.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_bst_eigs_smallest.m`:

```matlab
function test_bst_eigs_smallest()
% TEST_BST_EIGS_SMALLEST: synthetic singular symmetric generalized pencil.
% A = Neumann path-graph Laplacian (singular: constant null vector, smallest eig 0).
% B = 1e-6*I  (tiny mass, mimics meters-scale conditioning: generalized eig ~1e6*eig(A)).
    n = 60; k = 6;
    e = ones(n,1);
    A = spdiags([-e, 2*e, -e], -1:1, n, n);
    A(1,1) = 1; A(n,n) = 1;                 % Neumann ends -> row sums 0 -> constant null
    B = 1e-6 * speye(n);
    opts = struct('tol',1e-10,'disp',0,'maxit',1000);

    % Dense reference for the k smallest generalized eigenvalues
    lam_ref = sort(eig(full(A), full(B)));
    lam_ref = lam_ref(1:k);

    % Solver under test: must not raise the nearly-singular warning
    lastwarn('');
    [V, D] = bst_eigs_smallest(A, B, k, opts);
    [wmsg, wid] = lastwarn;
    assert(~strcmp(wid, 'MATLAB:nearlySingularMatrix'), ...
        'bst_eigs_smallest raised the nearly-singular warning: %s', wmsg);

    lam = sort(real(diag(D)));
    % (1) spectrum matches the dense reference
    scale = max(1, abs(lam_ref(end)));
    assert(max(abs(lam - lam_ref)) < 1e-6 * scale, ...
        'spectrum mismatch: max diff %.3e', max(abs(lam - lam_ref)));
    % (2) the zero (kernel) mode is captured
    assert(abs(lam(1)) < 1e-7 * scale, 'zero/kernel mode not captured: %.3e', lam(1));
    % (3) B-orthonormalize the returned vectors, then check residual + orthonormality
    nrm = sqrt(real(diag(V' * (B * V))));
    nrm(nrm < eps) = eps;
    Vn = V * diag(1 ./ nrm);
    Lam = diag(real(diag(D)));
    res = norm(A*Vn - B*Vn*Lam, 'fro') / (normest(A) * sqrt(k));
    assert(res < 1e-6, 'eigenpair residual too large: %.3e', res);
    G = Vn' * (B * Vn);
    assert(max(max(abs(G - eye(k)))) < 1e-8, 'not B-orthonormal: %.3e', max(max(abs(G-eye(k)))));

    disp('PASS: test_bst_eigs_smallest');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run (MATLAB MCP `run_matlab_file` on `dev/tests/test_bst_eigs_smallest.m`):
```matlab
run('dev/tests/test_bst_eigs_smallest.m')
```
Expected: FAIL — `Undefined function 'bst_eigs_smallest'`.

- [ ] **Step 3: Write the minimal implementation**

Create `toolbox/math/bst_eigs_smallest.m`. Prepend the standard Brainstorm GPL license header block (copy verbatim from `toolbox/anatomy/tess_eigen.m` lines 52–68, adjusting the function name in the leading comment). Function body:

```matlab
function [V, D] = bst_eigs_smallest(A, B, k, opts)
% BST_EIGS_SMALLEST: k smallest generalized eigenpairs of a (near-)singular,
% symmetric/Hermitian pencil (A, B) with B SPD, avoiding the sigma=0 shift-invert
% of the singular A that triggers MATLAB's RCOND "close to singular" warning.
%
% USAGE:  [V, D] = bst_eigs_smallest(A, B, k, opts)
%
% INPUTS:
%    - A    : [n x n] symmetric (real) or Hermitian (complex) operator, may be singular
%    - B    : [n x n] symmetric positive-definite mass matrix
%    - k    : number of smallest-magnitude generalized eigenpairs to return
%    - opts : struct passed through to eigs (e.g. tol, maxit, disp)
% OUTPUTS:
%    - V, D : eigenvectors / diagonal eigenvalues, same meaning as eigs(A,B,k,'smallestabs')
%
% Forces the symmetric/Hermitian Lanczos path (real spectrum, faster, clean
% degenerate multiplets) and selects the smallest modes via a small negative
% sigma shift so (A - sigma*B) = A + |sigma|*B is SPD/well-conditioned.

    if nargin < 4 || isempty(opts)
        opts = struct('tol', 1e-6, 'maxit', 1000, 'disp', 0);
    end
    % Symmetrize / Hermitian-ize: no-op to ~1e-16 for already-symmetric operators,
    % but makes issymmetric()/ishermitian() true so eigs uses Lanczos, not Arnoldi.
    A = (A + A') / 2;
    B = (B + B') / 2;
    % Factorization-free spectrum-scale estimate: 'largestabs' factorizes only the
    % SPD, well-conditioned mass B (never the singular A), so it cannot warn.
    lmax = abs(eigs(A, B, 1, 'largestabs', opts));
    if ~isfinite(lmax) || (lmax <= 0)
        [V, D] = eigs(A, B, k, 'smallestabs', opts);   % degenerate scale: legacy path
        return;
    end
    % Small negative shift below the (PSD) spectrum bottom. sigma < 0 can never
    % coincide with an eigenvalue, and 'nearest sigma' still returns the k smallest
    % (including the near-zero kernel).
    try
        [V, D] = eigs(A, B, k, -1e-7 * lmax, opts);
    catch
        try
            [V, D] = eigs(A, B, k, -1e-4 * lmax, opts);   % larger lift
        catch
            [V, D] = eigs(A, B, k, 'smallestabs', opts);  % legacy fallback
        end
    end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run:
```matlab
run('dev/tests/test_bst_eigs_smallest.m')
```
Expected: prints `PASS: test_bst_eigs_smallest`, no error, no `MATLAB:nearlySingularMatrix` warning.

- [ ] **Step 5: Lint the new function**

Run (MATLAB MCP `evaluate_matlab_code`):
```matlab
msg = checkcode('toolbox/math/bst_eigs_smallest.m','-id','-string'); if isempty(msg); disp('(clean)'); else; disp(msg); end
```
Expected: `(clean)` (or only pre-existing Brainstorm-idiom messages — there should be none for a new file).

- [ ] **Step 6: Commit**

```bash
git add toolbox/math/bst_eigs_smallest.m dev/tests/test_bst_eigs_smallest.m
git commit -m "feat(eigen): bst_eigs_smallest — symmetric scale-aware shifted eigs

Standalone robust solver for the k smallest generalized eigenpairs of a
near-singular symmetric/Hermitian pencil. Symmetrizes to force the Lanczos
path and uses sigma=-1e-7*lambda_max to avoid factorizing the singular A.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Wire `bst_eigs_smallest` into `tess_eigen`

**Files:**
- Modify: `toolbox/anatomy/tess_eigen.m:216-217`

- [ ] **Step 1: Write the failing integration check**

Add this function to the end of `dev/tests/test_bst_eigs_smallest.m` and a call to it at the top-level (after the synthetic test). It reuses a real stored Dirac operator node if present (no nxr/surface needed — it only loads matrices), and asserts the legacy path warns while the new path does not, with matching spectra:

```matlab
function test_bst_eigs_smallest_real_operator()
% Reuse a stored Dirac operator node (if any) to confirm: legacy 'smallestabs'
% warns, bst_eigs_smallest does not, and the K smallest eigenvalues agree.
    PI = bst_get('ProtocolInfo');
    if isempty(PI) || ~isfield(PI,'SUBJECTS') || isempty(PI.SUBJECTS)
        disp('SKIP: no protocol loaded'); return;
    end
    opFiles = dir(fullfile(PI.SUBJECTS, '**', 'operator_*.mat'));
    f = '';
    for i = 1:numel(opFiles)
        p = fullfile(opFiles(i).folder, opFiles(i).name);
        S = load(p, 'Variant');
        if isfield(S,'Variant') && strcmpi(S.Variant,'Dirac'); f = p; break; end
    end
    if isempty(f); disp('SKIP: no Dirac operator node found'); return; end

    Op = load(f); A = Op.Operator{1}; B = Op.Mass{1};
    k = 12; opts = struct('tol',1e-6,'maxit',1000,'disp',0);

    lastwarn('');
    lam_old = sort(real(eigs(A, B, k, 'smallestabs', opts)));
    [~, wid_old] = lastwarn;

    lastwarn('');
    [~, D] = bst_eigs_smallest(A, B, k, opts);
    [~, wid_new] = lastwarn;
    lam_new = sort(real(diag(D)));

    assert(~strcmp(wid_new, 'MATLAB:nearlySingularMatrix'), 'new path warned');
    % Compare the well-conditioned (non-kernel) modes as a set; kernel modes ~0.
    sc = max(abs(lam_old));
    assert(max(abs(lam_old - lam_new)) < 1e-3 * sc, ...
        'real-operator spectrum mismatch: %.3e (rel %.1e)', ...
        max(abs(lam_old-lam_new)), max(abs(lam_old-lam_new))/sc);
    fprintf('PASS: test_bst_eigs_smallest_real_operator (legacy warned id=%s)\n', wid_old);
end
```

Update the top of the file so both run:

```matlab
function test_bst_eigs_smallest()
    run_synthetic();
    test_bst_eigs_smallest_real_operator();
end
function run_synthetic()
    % ... (the synthetic body from Task 1 Step 1 moves here unchanged) ...
end
```

- [ ] **Step 2: Run to verify the integration check passes against the legacy-wired tess_eigen path is not yet relevant — run the real-operator check now**

Run:
```matlab
run('dev/tests/test_bst_eigs_smallest.m')
```
Expected: `PASS: test_bst_eigs_smallest` and either `PASS: test_bst_eigs_smallest_real_operator (legacy warned id=MATLAB:nearlySingularMatrix)` or a `SKIP:` line. This proves `bst_eigs_smallest` removes the warning on the real operator before we wire it in.

- [ ] **Step 3: Modify `tess_eigen.m` to call the helper**

In `toolbox/anatomy/tess_eigen.m`, replace lines 216–217:

```matlab
        opts = struct('tol', 1e-6, 'maxit', 1000, 'disp', 0);
        [V, D] = eigs(A, B, nRequest, 'smallestabs', opts);
```

with:

```matlab
        opts = struct('tol', 1e-6, 'maxit', 1000, 'disp', 0);
        % Robust smallest-eigenpair solve: A has a near-zero kernel (constant /
        % constant-quaternion modes), so eigs(...,'smallestabs') would shift-invert
        % at sigma=0 and factorize the singular A (RCOND ~ 1e-18 warning). The helper
        % symmetrizes A,B (forces the Lanczos path) and uses a small negative sigma.
        [V, D] = bst_eigs_smallest(A, B, nRequest, opts);
```

- [ ] **Step 4: Lint `tess_eigen.m`**

Run:
```matlab
msg = checkcode('toolbox/anatomy/tess_eigen.m','-id','-string'); if isempty(msg); disp('(clean)'); else; disp(msg); end
```
Expected: `(clean)` (no new messages relative to before the edit).

- [ ] **Step 5: Re-run the test file (now exercises the wired helper path indirectly)**

Run:
```matlab
run('dev/tests/test_bst_eigs_smallest.m')
```
Expected: both PASS (or real-operator SKIP if no node), no `nearlySingularMatrix` warning.

- [ ] **Step 6: Commit**

```bash
git add toolbox/anatomy/tess_eigen.m dev/tests/test_bst_eigs_smallest.m
git commit -m "fix(eigen): route tess_eigen eigensolve through bst_eigs_smallest

Replaces the shared eigs(A,B,K,'smallestabs') at tess_eigen line 217 (which
factorizes the singular operator at sigma=0 and warns RCOND~1e-18) with the
symmetric scale-aware shifted solve, for the LBO/Connection/Dirac variants.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-review

**Spec coverage:**
- Symmetrize + scale-aware negative shift → Task 1 Step 3 (`bst_eigs_smallest`). ✓
- Shared across all 3 variants → Task 2 (single shared call site). ✓
- Factorization-free `λ_max` via `'largestabs'` → Task 1 Step 3. ✓
- Error/convergence fallback (retry larger lift → legacy) → Task 1 Step 3. ✓
- Invariant: returned `(V,D)` same meaning, downstream unchanged → Task 2 Step 3 (only the call swapped). ✓
- Tests: no-warning (1), eigenpair residual (2), B-orthonormality (3), spectrum match vs legacy (4), real/≥0 — residual+spectrum assertions in Task 1 Step 1; spectrum-match-vs-legacy in Task 2 Step 1. ✓
- Out of scope: issue B, process migration, dev-bench updates — untouched. ✓

**Placeholder scan:** no TBD/TODO; all code blocks complete; the only "copy from sibling" instruction (license header) names exact source lines. ✓

**Type/name consistency:** `bst_eigs_smallest(A,B,k,opts)` signature identical across the helper, the synthetic test, the real-operator test, and the `tess_eigen` call site. ✓
