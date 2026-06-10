# Design: Robust smallest-eigenpair solve in `tess_eigen` (eigensolve conditioning)

**Date:** 2026-06-10
**Branch:** `feat/bst-dirac-leadfield`
**Status:** Approved design — ready for implementation plan
**Scope:** Issue **A** only (numerical conditioning of the eigensolve). Issue **B** (τ block-imbalance / operator definition) is explicitly **deferred to a later phase**.

## Problem

Computing the Dirac eigenbasis via `tess_eigen(SurfaceFile, 'Dirac')` emits:

```
Warning: First input matrix is close to singular or badly scaled
(RCOND = 1.242899e-18) and results may be inaccurate. Consider
specifying a small nonzero numeric sigma value instead of 'smallestabs'...
 In eigs>AminusSigmaBSolve / eigs>WarnIfIllConditioned
 In tess_eigen (line 217)
```

The same `eigs(A, B, nRequest, 'smallestabs', opts)` at `toolbox/anatomy/tess_eigen.m:217` is shared by all three operator variants (Laplace–Beltrami, Connection Laplacian, Dirac), all of which have a near-zero kernel and therefore hit the same issue.

## Root cause (investigated and confirmed)

`eigs(..., 'smallestabs')` performs shift-invert at **sigma = 0**, which factorizes **A** itself (`A − 0·B = A`). The operators are genuinely singular:

- The relative-Dirac operator has a ~4-dimensional constant-quaternion kernel (LBO: 1-dim constant null; Connection: its own kernel).
- Measured on the real operator node (TutorialAuditory/Subject01, left hemi, 40968×40968): `condest(A) = 4.1e20` (rcond ≈ 2.4e-21), matching the warned RCOND.
- A small nonzero shift (`sigma = −1e3`) removes the warning **and runs faster** (0.5 s vs 1.4 s).

Two compounding contributors map onto the warning's two clauses:

| Clause | Contributor | Evidence |
|---|---|---|
| "close to singular" | genuine kernel → sigma-0 factorizes a singular A | condest(A)=4.1e20; ~4 zero eigenvalues |
| "badly scaled" | meters units → mass ~ area ~1e-6, spectrum 0 → 2.3e13 | λ_max = 2.3e13 |

**Not a bug in the operator.** The nxr-compute C++ assembly (`src/facets.cpp diracFamily_`, `src/dirac_operator.cpp` `E = DᵀW_F D`) and the CSC→MATLAB marshaling were verified correct: the operator is real-symmetric to 1e-16, PSD, real spectrum, with genuine quaternion coupling (not scalar·I₄). One actionable detail: A is symmetric only to ~1e-16, **not bit-exact**, so `issymmetric(A)=0` and `eigs` falls back to the general (Arnoldi) path instead of symmetric Lanczos — slower, complex iterates, worse degenerate-multiplet recovery.

## Decisions

- **Approach 1 (chosen): symmetrize + scale-aware negative shift.** Selected over scalar-pencil-normalization (Approach 2) because it achieves the same unit-robustness (σ tracks `λ_max`) without a normalize/un-normalize round-trip that risks the B-orthonormality invariant.
- **Kernel deflation rejected:** the near-zero modes are the smoothest, most leadfield-relevant modes — they must remain in the basis.
- **Mesh-rescaling to mm rejected for Dirac:** the operator mixes a scale-invariant intrinsic block (`cotanL`, s⁰) with a scale-dependent extrinsic block (`E ~ s⁻²`); geometric rescaling changes the intrinsic/extrinsic balance → *different* eigenvectors (not merely rescaled), corrupting the basis. (GBFs can rescale only because its LBO is purely intrinsic.)
- **Shared fix** across LBO / Connection / Dirac (single code path), per scope decision.

## Design

### Location & structure

Single change in `toolbox/anatomy/tess_eigen.m`: replace the bare call at line 217

```matlab
[V, D] = eigs(A, B, nRequest, 'smallestabs', opts);
```

with a call to a new file-local helper `local_eigs_smallest(A, B, nRequest, opts)`. All downstream code (`sort`, `real`, the LBO column normalization, `local_ritz_basis`) is unchanged.

### The helper

```matlab
function [V, D] = local_eigs_smallest(A, B, k, opts)
% Robust smallest-magnitude generalized eigenpairs for a (near-)singular,
% symmetric/Hermitian operator A against an SPD mass B. Avoids the sigma=0
% shift-invert of the singular A that triggers MATLAB's RCOND warning.
    % Force the symmetric/Hermitian Lanczos path: real spectrum, faster, clean
    % degenerate multiplets. No-op to ~1e-16 (A,B already symmetric) — only strips
    % last-bit noise that makes issymmetric()/ishermitian() return false.
    A = (A + A')/2;
    B = (B + B')/2;
    % Factorization-free spectrum-scale estimate. 'largestabs' factorizes only the
    % SPD, well-conditioned mass B (never the singular A), so it cannot warn.
    lmax = abs(eigs(A, B, 1, 'largestabs', opts));
    if ~isfinite(lmax) || lmax <= 0
        % Degenerate scale estimate: fall back to legacy behavior.
        [V, D] = eigs(A, B, k, 'smallestabs', opts);
        return;
    end
    % Small negative shift below the PSD spectrum bottom: (A - sigma*B) = A + |sigma|*B
    % is SPD/well-conditioned, yet 'nearest sigma' still returns the k smallest modes
    % (including the near-zero kernel we keep). sigma<0 on a PSD spectrum can never
    % coincide with an eigenvalue, so there is no shift-collision risk.
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

### Error handling / fallback

- Guard `lmax` (finite, > 0); otherwise fall back to legacy `'smallestabs'`.
- If the shifted solve fails to converge, retry once with a larger lift (`−1e-4·lmax`), then fall back to legacy. A numerically pathological surface degrades to today's behavior rather than erroring.

### Invariant preserved

`local_eigs_smallest` returns `(V, D)` with identical meaning to the old call. Downstream code still:
- sorts ascending and takes `real` of eigenvalues,
- B-orthonormalizes eigenvectors against the **true** meters-scale mass `B` (LBO column normalization; Rayleigh–Ritz for Connection/Dirac).

So the **stored basis is mathematically the one already intended** (true generalized eigenvalues; eigenvectors B-orthonormal vs the true mass). Issue B is untouched. Only the computation changes: no warning, guaranteed-real spectrum, faster, cleaner multiplets.

## Testing (verification)

Per-variant check on a canonical test surface (LBO, Connection, Dirac):

1. **No warning** fires across the solve (`lastwarn('')` before; assert empty after).
2. **Eigenpair residual** `‖Aφ − λBφ‖ / (‖A‖·‖φ‖) < 1e-6` for every stored column.
3. **B-orthonormality** `‖ΦᵀBΦ − I‖∞ < 1e-8`.
4. **Spectrum match vs legacy path:** the K eigenvalues agree with the old `'smallestabs'` result within tolerance (degenerate multiplets compared as sets), confirming downstream `bst_dirac` output is unchanged.
5. **Spectrum real and ≥ 0** after clamping.

Per standing scope, the dev-test wiring lands with the later test pass; this spec records the assertions so the implementation plan includes them.

## Out of scope (deferred)

- **Issue B** — τ block-imbalance / operator-definition normalization (`L(τ)` mixes dimensionally-incoherent blocks; at meters scale τ=0.5 is ≈99.999% extrinsic). Modeling change that *redefines the basis*; needs its own benchmark before committing 400-mode bases.
- `process_dirac_eigenmode_leadfield.m` migration; the numerical `bst_dirac` Reconstruct equivalence check; `dev/` test/bench updates.

## References

- `toolbox/anatomy/tess_eigen.m` (line 217, shared eigs call)
- `toolbox/anatomy/tess_operators.m` (operator/mass assembly via nxr-compute)
- nxr-compute `src/facets.cpp` (`diracFamily_`), `src/dirac_operator.cpp` (`extrinsicBlock`)
- Memory: `dirac-eigenmode-leadfield` (eigensolve investigation; issue B)
