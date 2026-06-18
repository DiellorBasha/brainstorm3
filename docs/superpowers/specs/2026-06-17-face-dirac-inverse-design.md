# Face-Dirac eigenbasis inverse — Design

**Date:** 2026-06-17
**Status:** design (awaiting user review)
**Repo:** brainstorm3 (`development`); no nxr changes (MATLAB assembly per decision)

---

## 1. Goal & motivation

A full-3D **face-based** minimum-norm inverse expressed in the **face Dirac eigenbasis** — the exact dual of the validated vertex `bst_inverse_dirac`. Sources are unconstrained dipoles at face centroids (the leadfield from the previous phase); the source side is expressed in curvature-aware face Dirac eigenmodes. This completes end-to-end parity between the vertex and face pipelines.

The vertex inverse chain is: `tess_operators`(Dirac operator `A=(1−τ)D²_int+τE` `[4nV×4nV]`, mass `B`) → `tess_eigen`(eigenmodes `Phi [4nV×K]`) → `bst_dirac`(Transform leadfield→mode-forward; Reconstruct modes→vertices) → `bst_inverse_dirac`(whitened MNE in mode space). The face pipeline mirrors every stage on the **face** domain.

**Decisions (locked):** eigenbasis = `(1−τ)Ẽ_int + τẼ_ext` mirroring the vertex Dirac (τ=0.5 default); `Ẽ_int = D̃_intᵀ W_V D̃_int` assembled in MATLAB in `tess_operators` (parity with `local_dirac_intrinsic_sq`).

## 2. Architecture — four phases

### Phase 1 — `tess_operators` `'Dirac-Face'` variant
Per hemisphere submesh (closed, via `tess_hemisplit`):
- `D̃_int = nxr_compute('operators', h, 'diracFaceIntrinsicD')` `[4V×4F]`; `D̃_ext = diracFaceD` `[4V×4F]`.
- vertex dual-area mass `W_V = kron(diag(vertexDualAreas), I₄)` `[4V×4V]`; face-area mass `W_F = kron(diag(faceArea), I₄)` `[4F×4F]`.
- intrinsic face Dirac²: `Ẽ_int = D̃_intᵀ W_V D̃_int` `[4F×4F]` (MATLAB; new helper `local_dirac_face_intrinsic_sq`).
- extrinsic face Dirac²: `Ẽ_ext = nxr_compute('operators', h, 'diracFace', 1)` — i.e. `extrinsicBlockFace` `[4F×4F]`.
- **co-normalize** each block to unit largest generalized eigenvalue vs `W_F` (reuse the vertex Dirac co-normalization logic so τ is dimensionless), then `Ã = (1−τ)·Ẽ_int^ + τ·Ẽ_ext^` `[4F×4F]`, `B̃ = W_F`.
- store `Operator{hh}=Ã`, `Mass{hh}=B̃`, and a per-hemi **`GlobalFaces{hh}`** scatter map (NEW field on operatormat — faces are the domain now; keep `GlobalVertices` for reference). `Variant='Dirac-Face'`, `diracScales`, `Tau` in Provenance.

### Phase 2 — `tess_eigen` `'Dirac-Face'` variant
- Recognize `OperatorName='Dirac-Face'` → `Variant='Dirac-Face'`, treated like Dirac (quaternion multiplets → Rayleigh-Ritz `local_ritz_basis`).
- Load the face pencil `(Ã [4F×4F], B̃ [4F×4F])` via `local_find_operator`/`tess_operators`.
- `bst_eigs_smallest(Ã, B̃, K+over-fetch)` → Rayleigh-Ritz → `Phi [4nF×K]` (B̃-orthonormal), `Lambda [K]`.
- store `eigen_` node: `Phi{1×2} [4nFh×K]`, `Lambda{1×2}`, `Variant='Dirac-Face'`, `GlobalFaces{1×2}`, `K`, `Tau`, Provenance. `db_add_eigen`.

### Phase 3 — `bst_dirac` face branch
Detect the face domain from the head model (`HeadModel.isFaceBased==1`) and fetch the `'Dirac-Face'` eigen node (Tau/K match).
- **Transform** `CompHM = bst_dirac(faceHM, 'nModes',K, 'Tau',τ)`: per hemi, embed the face leadfield `[nCh×3nFh]` as pure-imaginary quaternions `Psi [4nFh×nCh]` (rows 2,3,4 of each 4-block = x,y,z), then `L̃_h = Psiᵀ B̃_h Phi_h` `[nCh×K]`; stack `[L̃_L L̃_R]` `[nCh×2K]`. `CompHM.isFaceBased=1`, `isDiracEigenmode=1`, `ModeHemisphere`, `HemiGlobalFaces`.
- **Reconstruct** `J = bst_dirac(CompHM,'Reconstruct',c)`: per hemi `R_h = Phi_h c_hᵀ` `[4nFh×m]`, take rows 2,3,4 → `J[:, faces]` `[m×3nF]`.
- Implementation: add a domain switch inside `bst_dirac` keyed on `isFaceBased` (modes/mass/scatter come from faces); the quaternion embed/extract math is identical to the vertex path.

### Phase 4 — `bst_inverse_dirac` face path
`bst_inverse_dirac` already routes Transform (line 105) and Reconstruct (line 243) through `bst_dirac`. Make it pass a face head model through unchanged:
- accept `HeadModel.isFaceBased==1` (face leadfield `[nCh×3nF]`); `bst_dirac` then returns face modes, and STAGE 4 `Wres = bst_dirac(CompHM,'Reconstruct',VL')'` returns `[3nF×r]`, so the kernel is `[3nF×nCh]` (the correct face-source kernel — fixes the silent vertex-fallback observed in the leadfield benchmark).
- the dSPM/sLORETA per-source (3-component) normalization is domain-agnostic (operates on kernel rows) → works on faces unchanged. `nVert` → `nFace` naming only.

## 3. Correctness contract (tests, per phase)

**Phase 1** (`test_dirac_face_operator_bst.m`): `Ã`/`B̃` are `[4F×4F]`, symmetric, `Ã` PSD; co-normalization gives both blocks unit top-eigenvalue; `B̃ == W_F`; `GlobalFaces` covers all faces.

**Phase 2** (extend `test_tess_eigen` or new): `Phi` is `[4nFh×K]`; **B̃-orthonormal** `‖Phiᵀ B̃ Phi − I‖<1e-6`; `Lambda` ascending ≥ −ε; the lowest modes are smooth (low Dirichlet energy); multiplet structure resolved (no NaN/duplicate-collapse from Rayleigh-Ritz).

**Phase 3** (`test_bst_dirac_face.m`): Transform shapes `[nCh×2K]`; **Reconstruct∘Transform** of a leadfield ≈ its `W_F`-orthogonal projection onto the mode span (residual ↓ as K↑); a field built from the modes round-trips to machine precision; energy/`B̃`-norm preserved in mode coefficients.

**Phase 4** (`test_face_dirac_inverse.m`): face kernel is `[3nF×nCh]` (NOT `3nV`); on the real alpha frame the face Dirac-inverse source map co-localizes with the vertex Dirac-inverse (peak within a few mm) and with the plain face wMNE already validated (corr high); observability rank matches.

## 4. Validation benchmark

`dev/benchmarks/bench_face_dirac_inverse.m`: build the face `'Dirac-Face'` eigenmodes (K≈400), run `bst_inverse_dirac` with the unconstrained face head model, compare to (a) the vertex Dirac inverse and (b) the plain face wMNE (from `bench_face_leadfield`). Report peak separation (mm), source-power correlation, observability rank, mode-truncation HarmFrac; side-by-side PNG. This is the end-to-end parity demonstration.

## 5. Risks

- **Conditioning:** the intrinsic face Dirac (barycentric-dual) and the face representation are rougher than vertices (prior finding); the face Dirac spectrum may have larger degenerate multiplets / rank deficiency. Mitigation: the Rayleigh-Ritz path (already used for vertex Dirac) handles this; Phase-2 B̃-orthonormality test is the gate. If conditioning fails, fall back to intrinsic-only `Ẽ_int` (the decision's runner-up) — an operator-assembly swap, no pipeline change.
- **Schema:** adding `GlobalFaces` to operatormat/eigenmat is a small `db_template` addition; face nodes set it, vertex nodes leave it empty (back-compatible).

## 6. Out of scope

EEG/OpenMEEG; GUI/process wrappers; retiring the scalar-LBO `bst_face_eigenmode_leadfield`; depth weighting beyond what the vertex `bst_inverse_dirac` already has.

## 7. Self-review

- **Placeholders:** none — each phase has concrete operators (existing nxr `diracFaceIntrinsicD`/`diracFace`, existing `bst_eigs_smallest`/`local_ritz_basis`) and a runnable gate.
- **Parity:** every stage mirrors a named vertex stage; the quaternion embed/extract and the MNE stages are reused verbatim (domain-agnostic).
- **Decisions honored:** τ-mix intrinsic+extrinsic; `Ẽ_int` in MATLAB; extend existing functions (no new nxr).
- **Scope:** 4 phases, MATLAB-only, each independently testable; fallback path stated.
