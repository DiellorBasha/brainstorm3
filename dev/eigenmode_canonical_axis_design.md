# Canonical Eigenmode Axis — Design (Foundation A)

**Date:** 2026-06-03 (rev 2)
**Author:** Diellor Basha (with Claude)
**Status:** Design — approved in brainstorming; revised for the pure `manifold_ft`/`manifold_ift` pair; pending spec review before the implementation plan
**Scope:** Foundation refactor only. The eigenmode time-series ↔ cortical-activation viewer (Feature B) is a separate, later cycle that depends on this.

## Motivation

A cortical surface file has **one** vertex axis (Vertices/faces) shared by all clients without divergence. The Laplace–Beltrami eigenmodes are the **Fourier space of that same surface**, so they should likewise be a single canonical axis per surface, shared by every consumer.

Today they are not. The stored eigenmodes are **grouped by hemisphere** (cols 1–1000 = hemisphere-1, 1001–2000 = hemisphere-2; sorted only *within* each, not globally):

```
total modes = 2000, components = 2, issorted(Values) = 0
component-2 starts at stored index 1001   (so first 1000 stored modes are all hemisphere-1)
```

So consumers that take a subset re-invent mode selection differently:
- `bst_eigenmode_leadfield` re-sorts globally by eigenvalue and records its own `ModeIndices` on the composed head model.
- `process_eigenmodes_transform` naively takes `Vectors(:,1:K)` — for K<1000 that is **entirely hemisphere-1** (a real bias bug).
- `bst_eigenmode_reconstruct` takes `ModeIndices` if given, else falls back to first-K.

These are different selections of the same modes (the same bug class fixed in the benchmark's `bst_benchmark_inverse`).

## Principles

1. **One canonical, globally eigenvalue-sorted ordering** of a surface's eigenmodes, stored with the eigenmode structure (`Eigenmodes.Order`). It is the single index; every consumer selects through it (`Vectors(:, Order(1:K))`).
2. **The analytic transform is pure math, named for what it is — a manifold Fourier transform.** Forward `manifold_ft(Φ,M,U) = Φ'·M·U` (vertex → mode coefficients); inverse `manifold_ift(Φ,C) = Φ·C` (coefficients/kernel → vertex). Pure functions in `toolbox/math/`, independent of files, sensors, or selection.
3. **Source-mapping machinery is separate** from this analytic core: the eigenmode leadfield (`L·Φ`) + inverse + reconstruction build *on top of* the same `Φ` and the same canonical `Order`.

## Design

### 1. Canonical mode index — `Eigenmodes.Order`

The `Eigenmodes` struct in the surface file gains:
- **`Eigenmodes.Order`** : `[nModes×1]` permutation sorting the stored modes by **global eigenvalue ascending** (interleaving hemispheres).

It is *informative*: with `.Component`/`.CompRank`, canonical rank `k` is fully described — stored column `Order(k)`, eigenvalue `Values(Order(k))`, hemisphere `Component(Order(k))`, within-hemisphere rank `CompRank(Order(k))`. "First K canonical modes" = `Vectors(:, Order(1:K))` = whole-brain lowest spatial frequencies, never hemisphere-biased.

- `tess_eigenmodes` writes `.Order` at computation time (`sort(Values,'ascend')`).
- `in_tess_eigenmodes` **backfills** `.Order` on read for legacy files (same pattern it uses for `.Component`/`.CompRank`) — transparent, no migration.

**No dedicated accessor.** Selection is the inline idiom `Vectors(:, Order(1:K))`. The stored `Order` *is* the single shared index; a wrapper function would add a layer without adding correctness.

### 2. Pure analytic transforms — `manifold_ft` / `manifold_ift`

New, in `toolbox/math/`:
- **`manifold_ft(Phi, M, U)` → `C = Phi' * (M * U)`** — forward manifold Fourier transform, `[K×nTime]`. `Phi` is the (already-selected) eigenvector matrix; `M` the mass matrix (basis is M-orthonormal, so the forward needs `M`).
- **`manifold_ift(Phi, C)` → `U = Phi * C`** — inverse manifold Fourier transform, `[nV×nTime]` or `[nV×nCh]`. Works identically for coefficient vectors and mode-space kernels. No `M` needed.

These are pure (no I/O, no selection, no validation beyond dimension checks). They are the analytic machinery in named form.

### 3. `project` / `reconstruct` become thin wrappers (callers untouched)

`bst_eigenmodes_project` and `bst_eigenmode_reconstruct` have **12 callers** between them, including the production inverse `process_inverse_2018` and the pushed `bst_benchmark_inverse`. Rather than retire them (high churn/risk to the validated Eigen-MNE path), they become thin wrappers over the pure pair, holding only the non-math concerns:

- **`bst_eigenmodes_project`** = validate + `manifold_ft(Vectors, M, Data)`; optional reconstruct output via `manifold_ift`, with any `ModeRange` selecting over the canonical `Order`.
- **`bst_eigenmode_reconstruct`** = load surface eigenmodes (if given a file) + select canonical `Order` + `manifold_ift(Phi(:,Order(1:K)), ModeKernel)`. The hemisphere-bias fix is applied here, inside the wrapper.

Full retirement of these two wrappers (callers migrating directly to `manifold_ft`/`manifold_ift`) is an **optional later cleanup**, out of scope for this foundation.

### 4. Source-mapping consumers re-grounded on `Order`

- **`bst_eigenmode_leadfield`** — stop computing its own global sort + `ModeIndices`. Compose `L̃ = L · Vectors(:, Order)` over **all** modes (truncation moves to the inverse); copy `Order` into the composed head model as `ModeIndices` (derived, not recomputed).
- **`bst_inverse_eigenmodes`** — unchanged math; its `nModes` cap now means "first-K of the canonical axis," correct for free (the composed Gain is canonical-ordered).
- **`process_eigenmodes_transform`** — replace `Vectors(:,1:K)` with `Vectors(:, Order(1:K))` (the bias fix). Its reconstruction stays self-consistent.

**Viewers** (`view_eigenmodes`, `panel_eigenmodes`, `view_eigenmode_spectrum`) carry **no** first-K slice — they pair by `CompRank` tags — so they need no change for correctness here. The mode-k↔canonical-rank UX alignment is part of Feature B.

### 5. Leadfield contract + ensure-if-empty

- **Leadfield uses all canonical modes — no independent truncation.** K-capping moves to `bst_inverse_eigenmodes` (first-K of canonical).
- **`bst_eigenmodes_ensure(SurfaceFile)`** — new helper: return the surface's canonical eigenmodes if present; else compute the default (**1000 modes/hemisphere**, barycentric mass, remove-DC), store, return. The default does **not** repair: a non-manifold surface errors with a clear "remesh to ico" message (repair changes vertex count and breaks surface↔leadfield↔eigenmode consistency; ico5 avoids it). This is the principled home for what `bench_fixtures` does ad hoc.

### 6. Backward compatibility — consistency-preserving by construction

The *old* `bst_eigenmode_leadfield` already sorted globally by eigenvalue; it merely stored that index on the head model instead of the surface. Therefore:
- Legacy surface eigenmodes without `.Order` → backfilled on read. Transparent.
- Legacy composed head models' `ModeIndices` **already equal** the canonical `Order(1:K)`, so existing eigenmode kernels — including the validated `Eigen-MNE` in TutorialAuditory — **stay correct**. No migration.
- The only genuinely buggy behavior (`Vectors(:,1:K)` in the transform) is the one behavior that changes.

The refactor is **consistency-preserving for all source-mapping data** and **bug-fixing for the analytic transform**.

## Testing

**Pure transforms:** `manifold_ft`/`manifold_ift` round-trip on an M-orthonormal basis (`manifold_ift(Φ, manifold_ft(Φ,M,u)) == u` within tol); dimension errors raised; `manifold_ift` works for both a vector and a matrix `C`.

**Canonical index (pure + e2e):** `Values(Order)` globally ascending; `Order` a valid permutation. Discriminating e2e: on the 2-hemisphere cortex, `Order(1:600)` spans **both** components (old first-K = one hemisphere). Backfill: a legacy struct lacking `.Order` gets it filled on read.

**Wrapper equivalence (pure):** `bst_eigenmodes_project` equals `manifold_ft` (+ ranged `manifold_ift`); `bst_eigenmode_reconstruct` with explicit `ModeIndices` equals with canonical default.

**Consistency / no-regression (e2e, safety net):** on TutorialAuditory — new canonical-path leadfield Gain **equals** the existing composed eigenmode-HM Gain; the stored, validated Eigen-MNE kernel reconstructs **identically**.

**Bug-fix (behavior that should change):** `process_eigenmodes_transform` at K<1000 selects whole-brain modes (assert both components present); end-to-end transform runs and is finite.

**ensure-if-empty (e2e):** returns canonical fast when present (idempotent); errors clearly (no silent repair) on a non-manifold surface.

**Regression sweep:** re-run existing eigenmode tests (`test_eigenmode_leadfield_*`, `test_inverse_eigenmodes_*`, `test_eigenmode_reconstruct_*`, `test_eigenmodes_project_pure`, `test_kernel_comparison`) green.

All follow the `dev/tests` pattern (function, `addpath`, brainstorm, `assert`, `ALL TESTS PASSED`); pure where possible, e2e for protocol-dependent equivalence checks.

## Out of scope

- **Feature B** — the eigenmode time-series ↔ cortical-activation viewer.
- **Full retirement** of `bst_eigenmodes_project` / `bst_eigenmode_reconstruct` (caller migration to the pure pair) — optional later cleanup.
- Re-sorting the stored `Vectors` on disk (we use the `Order` index instead).
- Changing inverse math, priors, or the eigenfilter library.
- Viewer mode-k↔canonical-rank UX (Feature B).
