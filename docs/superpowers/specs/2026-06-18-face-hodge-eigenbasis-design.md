# Face Hodge vector eigenbasis ('Hodge-Face') — Design

**Date:** 2026-06-18
**Status:** design (approved approach; detailed schema below)
**Repo:** brainstorm3 (`development`). No nxr changes.
**Supersedes:** the face Dirac eigenbasis as the inverse basis (proven non-localizing — `face-dirac-inverse` memory). Keeps everything else (face leadfield, bst_dirac/bst_inverse_dirac face path).

---

## 1. Motivation

The face Dirac eigenbasis fails (43 mm mislocalization) because the face Dirac square is built from **wide roots** (`D̃ [4V×4F]`) → non-smooth modes that miss ~half the observable subspace. The fix: a **full-rank, smooth** face vector eigenbasis from the scalar face Laplacian `lapFace` (built from the **tall** `gradFace` root `[3F×F]` → Weyl spectrum), lifted to vectors via the Hodge/Helmholtz construction.

## 2. Construction (the math)

Per hemisphere submesh (closed), all operators from the already-built nxr `gradFace`/`lapFace` + geometry:
- Scalar face Laplacian eigenproblem: `lapFace ψ_k = λ_k M_F ψ_k`, `ψ_k ∈ ℝ^F`, `M_F = diag(faceArea)`, `ψ` `M_F`-orthonormal. Smallest `K` via `bst_eigs_smallest` (well-conditioned, simple spectrum — like LBO). Drop the constant mode (λ≈0).
- Hodge lift to two smooth vector families per scalar mode (`gradFace` = `G̃ [3F×F]`, `SkewG = n_f×G̃`):
  - irrotational `û_k = G̃ ψ_k / √λ_k` `[3F]`
  - solenoidal  `ŵ_k = SkewG ψ_k / √λ_k` `[3F]`
  - `‖û_k‖²_{W_F} = ψ_kᵀ lapFace ψ_k = λ_k` (since `lapFace = G̃ᵀW_F G̃`), so `/√λ_k` ⇒ unit `W_F`-norm; `⟨û_j,û_k⟩_{W_F}=δ_jk` and likewise for `ŵ` (rotation is a `W_F`-isometry). The two families are **individually** `W_F`-orthonormal; the cross term `⟨û,ŵ⟩` is small but nonzero (the discrete skew-grad isn't exactly div-free), so a final Gram-Cholesky `W_F`-orthonormalizes the stacked `[3F × 2K]` set.
- Embed as pure-imaginary quaternions `Phi [4F × 2K]`: rows 2,3,4 of each face block = the vector mode, row 1 (w) = 0. Then `Phiᵀ W_F⁽⁴⁾ Phi = I` (W_F⁽⁴⁾ = `kron(faceArea,I₄)`), exactly what `bst_dirac` requires.

## 3. Pipeline integration (minimal changes — reuses the face path)

- **`tess_operators 'Hodge-Face'`**: a node carrying `Operator{hh}=lapFace [F×F]` (record), `Mass{hh}=W_F⁽⁴⁾ [4F×4F]` (what `bst_dirac` reads as `B`), `GlobalFaces{hh}`. (Co-located with the existing face-variant plumbing.)
- **`tess_eigen 'Hodge-Face'`**: self-contained — recompute `gradFace`/`lapFace`/`M_F`/normals from nxr on the hemi submesh, scalar-eigensolve, Hodge-lift, Gram-Cholesky orthonormalize, store `Phi [4F×2K]` (pseudo-quaternion), `Lambda` (the `[λ_k;λ_k]` paired or the vector-mode Rayleigh quotients), `Variant='Hodge-Face'`, `GlobalFaces`, `K=2K_scalar`. References the `'Hodge-Face'` operator node for `Mass`.
- **`bst_dirac`**: the face branch currently hardcodes `Variant='Dirac-Face'`. Generalize: `Variant = HeadModel.FaceBasis` if set, else default `'Hodge-Face'` (the working basis). Everything else (quaternion embed, `Psiᵀ B Phi`, reconstruct) is unchanged — the `w=0` modes pass through identically.
- **`bst_inverse_dirac`**: **no change** (already domain-agnostic; consumes whatever `bst_dirac` returns).

## 4. Correctness contract (tests) + decisive gate

- **Phase A (`tess_eigen 'Hodge-Face'`)**: `Phi [4F×2K]`; `w`-rows are zero; **`Phiᵀ W_F⁽⁴⁾ Phi = I` to <1e-6** (orthonormal — the embed gate); each mode tangent (imag part ⟂ face normal). Scalar λ_k ascending, Weyl-like (ratio ≈ vertex LBO, NOT 845).
- **Phase B (`bst_dirac` Hodge face)**: Transform/Reconstruct match an independent re-derivation (rel err ~0); **mode-forward observability rank ≫ the Dirac-Face one** — expect it to approach the vertex ceiling like the vertex Dirac does (the smoothness test).
- **Phase C — THE GATE (`bench`)**: the face **Hodge** inverse localizes — peak within a few mm of the vertex Dirac inverse AND the plain face wMNE (vs the 43 mm Dirac-Face failure); focal, not diffuse; source-power corr high. This is the decisive empirical test the whole exercise is for.

## 5. Notes / risks

- This is a **Hodge/Helmholtz** vector basis, not a Dirac basis — no extrinsic (Gauss-map) quaternion coupling. That coupling is exactly what made the face Dirac square wide-rooted and unusable; the smooth scalar-Laplacian lift is the price for a working face vector basis. (The vertex pipeline keeps the true Dirac basis; the face pipeline uses the Hodge basis. Documented asymmetry, justified by the wide-root obstruction.)
- 2K vector modes from K scalar modes: to match the vertex inverse's ~400 modes/hemi, use K_scalar≈200 (→ 400 vector modes).
- If the cross-family non-orthogonality is large, the Gram-Cholesky may drop near-dependent columns (rank-reveal) — keep the survivors; log the count.

## 6. Self-review

- **Placeholders:** none — the lift is closed-form from existing operators; the gate is the runnable localization benchmark.
- **Reuse:** bst_dirac face branch + bst_inverse_dirac unchanged (pseudo-quaternion embedding); only the eigenbasis construction is new.
- **Decision honored:** face-native (lapFace/gradFace), full-rank, smooth; no vertex borrowing, no new nxr.
- **Decisive gate:** unlike Dirac-Face (which passed structural gates but failed localization), Phase C localization is the explicit accept/reject criterion.
