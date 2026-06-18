# Face-native dual-mesh gradient `gradFace` + face Laplacian `lapFace` — Design

**Date:** 2026-06-17
**Status:** design (awaiting user review)
**Repos:** `nxr-compute` (new C++ operators + bindings), `brainstorm3` (rewrite `bst_dirac_helmholtz_face`)
**Supersedes the open piece in:** `2026-06-17-face-domain-helmholtz-design.md` (the "dual reconstruction operator, still to be derived")

---

## 1. Motivation & decisive evidence

The face-domain Helmholtz needs to decompose a **per-face** 3D vector field `J` (full unconstrained dipoles at face centroids) into irrotational + solenoidal + harmonic parts, mirroring the validated **vertex** pipeline (`bst_dirac_helmholtz`).

The committed INTERMEDIATE (`bst_dirac_helmholtz_face`, commit `a6328fd6`) extracts curl/div on *vertices* via the intrinsic face Dirac `D̃_int`, then solves Poisson on the *vertex* cotan `K` and reconstructs with the *vertex* FEM gradient `G`. These three pieces are not mutually consistent → a planted pure-skew-gradient field recovers at corr 0.997 but leaves ~10% spurious HarmFrac (a Dirac-curl-vs-FEM-gradient discretization gap).

A decisive experiment (2026-06-17, validated to a Pythagoras identity at 1e-35) **ruled out** the tempting "no-new-operator" shortcut — projecting the per-face field directly onto `range(D̃*)` (the true adjoint `D̃* = W_F⁻¹ D̃ᵀ W_V`). Because `D̃_int` is `[4V×4F]` with `F≈2V`, `range(D̃*)` is only ~half of `ℍ^F`; a **random** imaginary face field shows **52% HarmFrac**, when a genus-0 closed hemisphere's true harmonic space is constants-only (~0%). That projection conflates `kernel(D̃)` with "harmonic" — physically wrong.

**Conclusion:** the only correct face-native Helmholtz uses a scalar Poisson on a **face Laplacian `K̃` with a constants-only kernel**, plus a **dual-mesh face gradient `G̃`** consistent with it by construction (`K̃ = G̃ᵀ W_F G̃`) — the exact dual of the vertex pipeline's `K = GᵀW_F G`.

## 2. The two new nxr operators (everything else is trivial MATLAB)

Define, per closed hemisphere (genus-0; reuse the existing closed-mesh-v1 guard):

- **`gradFace`** — dual-mesh gradient of a per-face scalar. Maps `ψ̃ ∈ ℝ^F` → a per-face ambient 3-vector field, returned as a sparse `[3F × F]` matrix (rows `3f+{0,1,2}` = x,y,z of face `f`). Each output vector lies in face `f`'s tangent plane. **Barycentric dual** construction (dual vertices at face barycenters `C_f`, matching `diracFaceIntrinsicD`'s centroid convention): the per-face gradient is assembled from barycentric-dual edge-neighbor differences (`ψ̃_g − ψ̃_f` across each of `f`'s three edges), so `gradFace · 1 = 0` (annihilates constants).

- **`lapFace`** — the face Laplacian `K̃ = gradFaceᵀ · W_F · gradFace`, `[F × F]`, where `W_F = diag(face areas)` lifted to 3 components per face. Built in C++ from the **same** `gradFace` triplets (single source of truth, exactly as `extrinsicBlock` builds `E` from `matrix`'s triplets — they never drift). Symmetric PSD; `kernel = constants` (genus-0).

**Why only these two:** rotation by the face normal `n_f` preserves the area-weighted inner product, so the **skew-gradient** `SkewG = n_f × gradFace` satisfies `SkewGᵀ W_F SkewG = gradFaceᵀ W_F gradFace = K̃`. Hence the stream-function Poisson and the potential Poisson use the **same** `K̃`. `SkewG`, the divergence map `Div = gradFaceᵀ W_F`, and the curl map `Curl = SkewGᵀ W_F` are all one-line MATLAB products of `gradFace` + per-face normals — no extra C++.

## 3. Correctness contract (the operators' acceptance tests)

These are the properties the implementation must satisfy; they double as the TDD tests.

| Property | Check |
|---|---|
| Constant precision | `gradFace · 1 = 0` (machine precision) |
| Tangency | each output 3-vector ⟂ its face normal (machine precision) |
| Single source of truth | `lapFace == gradFaceᵀ W_F gradFace` (machine precision) |
| Symmetry + PSD | `lapFace == lapFaceᵀ`; eigenvalues ≥ −ε |
| Harmonic space | `dim ker(lapFace) == 1` (constants only), on a closed hemisphere |
| Closed-mesh guard | throws on an open boundary (reuse existing guard message) |
| **Round-trip (the headline)** | plant `J = SkewG ψ0`; solve `K̃ ψ̃ = SkewGᵀ W_F J`; then `‖J − SkewG ψ̃‖²_{W_F} / ‖J‖²_{W_F} → 0` (HarmFrac→0) and `corr(ψ̃, ψ0) → 1` |

Linear precision (reproducing the gradient of a globally linear field) is **not** required — the barycentric-dual gradient trades it for robustness/no-negative-weights; the round-trip exactness is guaranteed by the `G̃ᵀ` adjoint regardless.

## 4. Brainstorm integration (`bst_dirac_helmholtz_face` rewrite)

`Prepare` (per hemisphere): pull `G̃ = gradFace` and `K̃ = lapFace` from a fresh nxr handle on the hemisphere submesh (same `nxr_safe_create` pattern already in the file). Build `W_F` (face areas), per-face normals `n_f`, `SkewG = n_f × G̃`, and cache the pinned-`K̃` Cholesky factor. Drop the vertex cotan `K`, the vertex FEM gradient `Gx/Gy/Gz`, and `D̃_int` from the Helmholtz path (D̃_int may stay only if we still want its vertex-domain vorticity as a *diagnostic* scalar — decide during implementation; default: drop it, all on faces).

`Frame(Op, Jf)`:
1. `divS = G̃ᵀ W_F Jf`; `curlS = SkewGᵀ W_F Jf` (per-face scalars).
2. Poisson on `K̃` (mean-zero, pinned): `φ̃ = K̃⁺ divS`, `ψ̃ = K̃⁺ curlS`.
3. `Virr = G̃ φ̃`; `Vsol = SkewG ψ̃ = n_f × G̃ ψ̃`; `Vharm = Jf − Virr − Vsol` (all `[F×3]`).
4. `HarmFrac = ‖Vharm‖²_{W_F} / ‖Jf‖²_{W_F}`.
5. Diagnostics now live on **faces**: vorticity `ω_f = curlS` (or `K̃ ψ̃`), divergence `div_f = divS`. Cores/sources: persistence on the per-face potentials `ψ̃`/`φ̃` over the **dual** adjacency (faces sharing an edge) — extend `bst_vortex_persistence` to accept a face-neighbor graph (it already takes a generic `'Neighbors'` list).

Return struct keeps the same field names where possible; scalars (`Curl`, `Div`, `Psi`, `Phi`) move from `[nV×1]` to `[nF×1]`. Update `test_dirac_helmholtz_face` gates to the strict round-trip (`HarmFrac < 0.02`, `corr > 0.99`).

## 5. Phasing

- **Phase C (nxr-compute, C++):** implement `gradFace` + `lapFace` (`src/`, `OperatorId`, cache, facet accessor), MEX + WASM dispatch, C++ test (`test/`) + MATLAB binding test, all properties in §3. Build, install the fresh mex into the plugin folder (with a dated backup, per the stale-binary trap), commit on a feature branch. **Do not** cut a `v*` release tag.
- **Phase D (brainstorm3, MATLAB):** rewrite `bst_dirac_helmholtz_face` per §4; tighten `test_dirac_helmholtz_face`; re-run `bench_dirac_face_helmholtz` (vertex vs face, now apples-to-apples); commit on `development`.

## 6. Out of scope (later)

Face leadfield + face-Dirac eigenbasis inverse (these reuse `gradFace`/`lapFace` and the rewritten Helmholtz). Circumcentric-dual or least-squares variants (only if barycentric conditioning proves inadequate on real cortex — the round-trip exactness does not depend on the choice).

## 7. Self-review

- **Placeholders:** none — the one genuinely-derived piece (the exact barycentric `gradFace` assembly) is bounded by the §3 property tests, which are concrete and runnable.
- **Consistency:** `K̃ = G̃ᵀW_FG̃` is enforced by building `lapFace` from `gradFace`'s triplets; the shared-`K̃` claim for both Poissons rests on `n×` being an area-weighted isometry (stated in §2).
- **Scope:** two C++ operators + one MATLAB rewrite — single coherent deliverable, two phases. Leadfield explicitly deferred.
- **Ambiguity:** whether to retain `D̃_int` as a diagnostic is flagged as an implementation-time decision with a stated default (drop it).
