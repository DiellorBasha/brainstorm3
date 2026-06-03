# Canonical Eigenmode Axis — Design (Foundation A)

**Date:** 2026-06-03
**Author:** Diellor Basha (with Claude)
**Status:** Design — approved in brainstorming; pending spec review before the implementation plan
**Scope:** Foundation refactor only. The eigenmode time-series ↔ cortical-activation viewer (Feature B) is a separate, later cycle that depends on this.

## Motivation

A cortical surface file has **one** vertex axis (Vertices/faces) shared by all clients without divergence. The Laplace–Beltrami eigenmodes are the **Fourier space of that same surface**, so they should likewise be a single canonical axis per surface, shared by every consumer.

Today they are not. The surface's stored eigenmodes are **grouped by hemisphere** (component 1's modes, then component 2's), sorted only *within* each component, not globally:

```
total modes = 2000, components = 2
issorted(Values) = 0
component-2 starts at stored index 1001   (so first 1000 stored modes are all hemisphere-1)
```

Because of this, consumers that take a subset re-invent mode selection differently:
- `bst_eigenmode_leadfield` re-sorts globally by eigenvalue and records its own `ModeIndices` on the composed head model.
- `process_eigenmodes_transform` / `bst_eigenmodes_transform` naively take `Vectors(:,1:K)` — which for K<1000 is **entirely hemisphere-1** (a real bias bug).
- `bst_eigenmode_reconstruct` takes `ModeIndices` if given, else falls back to first-K.

These are different selections of the same modes. (This is the same class of bug found and fixed in the benchmark's `bst_benchmark_inverse`.)

## Principle

There is **one canonical, globally eigenvalue-sorted ordering of a surface's eigenmodes**, stored with the eigenmode structure in the surface file, and **every consumer selects modes through one shared accessor**. The mode-selection index is not retired — it is relocated to its canonical home and made informative.

This cleanly separates two machineries that both rest on the same axis:
- **Source-mapping machinery:** eigenmode leadfield (sensor↔eigenmode) + inverse + reconstruction (eigenmode↔vertex).
- **Analytic machinery (separate):** forward eigenmode transform (eigenspectrum: any vertex pattern → coefficients, `Φ'M·u`) and inverse eigenmode transform (coefficients → vertices, `Φ·c`).

## Design

### 1. Canonical mode index — `Eigenmodes.Order`

The `Eigenmodes` struct in the surface file gains a stored field:

- **`Eigenmodes.Order`** : `[nModes×1]` permutation that orders the stored modes by **global eigenvalue ascending** (interleaving hemispheres).

It is *informative*: with the existing `.Component` and `.CompRank` tags, every canonical mode is fully described. For canonical rank `k`:
- stored column = `Order(k)`
- spatial frequency = `Values(Order(k))`
- hemisphere = `Component(Order(k))`, within-hemisphere rank = `CompRank(Order(k))`

"Give me the first K canonical modes" → `Vectors(:, Order(1:K))` = whole-brain lowest-spatial-frequency modes, never hemisphere-biased.

**Who fills it:**
- `tess_eigenmodes` writes `.Order` at computation time (one global `sort(Values, 'ascend')`).
- `in_tess_eigenmodes` **backfills** `.Order` on read if a legacy file lacks it (same pattern it already uses for `.Component`/`.CompRank`), so old protocols work transparently with no migration step.

### 2. Single shared accessor — `bst_eigenmodes_canonical`

```
[Phi, lambdas, prov] = bst_eigenmodes_canonical(Eig, K)
```
Returns `Vectors(:, Order(1:K))`, the matching eigenvalues `Values(Order(1:K))`, and provenance (`Component`/`CompRank` for the selected modes). `K` omitted/empty ⇒ all modes. This is the **one** place that turns "I want K modes" into eigenvector columns. Every consumer goes through it.

### 3. Consumers re-grounded on the accessor

Wherever code does `Vectors(:,1:K)` or computes its own sort, it calls the accessor instead.

**Source-mapping machinery:**
- `bst_eigenmode_leadfield` — stop computing its own global sort + `ModeIndices`. Compose `L̃ = L · Vectors(:, Order)` over **all** modes via the accessor; copy `Eig.Order` into the composed head model for self-containment (derived, never recomputed).
- `bst_eigenmode_reconstruct` — select via the accessor (`Φ(:, Order(1:K)) · ModeKernel`); the index it consumes *is* the canonical `Order` (from the surface, or the head model's copy). The first-K fallback is removed.
- `bst_inverse_eigenmodes` — unchanged math; its `nModes` cap now means "first-K of the canonical axis," whole-brain-correct for free.

**Analytic machinery (the bug-carriers today):**
- `process_eigenmodes_transform` / `bst_eigenmodes_transform` — replace naive `Vectors(:,1:K)` with the accessor. **This is where the hemisphere-bias bug lived; it is fixed here.** Reconstruction becomes `Φ(:, Order(1:K)) · Θ`.
- `bst_eigenmodes_project` — when it truncates to a mode range, that range is over the canonical `Order`.

**Viewers** (`view_eigenmodes`, `panel_eigenmodes`, `view_eigenmode_spectrum`) — wherever the lever indexes "mode k" or truncates a band, `k` maps to canonical rank via `Order`. The pair-by-`CompRank` grid is unaffected (tags travel with the modes).

Net: exactly one selection path, used by source-mapping and analytic consumers alike.

### 4. Leadfield contract + ensure-if-empty

- **Leadfield uses all canonical modes — no independent truncation.** `bst_eigenmode_leadfield` composes `L̃ = L · Vectors(:, Order)` over all stored modes. K-capping (e.g. the benchmark K-sweep) moves to the **inverse** (`bst_inverse_eigenmodes` `nModes` = first-K of canonical). The forward node is purely "the leadfield expressed in the surface's Fourier basis."
- **`bst_eigenmodes_ensure(SurfaceFile)`** — new helper: return the surface's canonical eigenmodes if present; otherwise compute the default (**1000 modes/hemisphere**, barycentric mass, remove-DC), store them, and return. The leadfield and any source-mapping consumer call this. It is the principled home for what `bench_fixtures` does ad hoc.
  - The default does **not** repair. On a non-manifold surface it errors with a clear "remesh to ico / repair manually" message, because repair changes vertex count and breaks surface↔leadfield↔eigenmode consistency (the benchmark lesson; ico5 avoids it).

### 5. Backward compatibility — consistency-preserving by construction

The *old* `bst_eigenmode_leadfield` already sorted globally by eigenvalue; it merely stored that index in the wrong place (the composed head model instead of the surface). Therefore:
- Legacy surface eigenmodes without `.Order` → backfilled on read. Transparent.
- Legacy composed head models' `ModeIndices` **already equal** the new canonical `Order(1:K)` (both are the global eigenvalue sort), so existing eigenmode kernels — including the validated `Eigen-MNE` in TutorialAuditory — **stay correct**. No migration, no recompute.
- The only genuinely buggy behavior (hemisphere-biased `Vectors(:,1:K)` in the analytic transform) is the one behavior that changes.

The refactor is **consistency-preserving for all source-mapping data** and **bug-fixing for the analytic transform**.

## Testing

**Canonical index (pure):**
- `Values(Order)` globally ascending; `Order` a valid permutation of `1:nModes`.
- Discriminating test: on a 2-hemisphere surface, `Order(1:600)` spans **both** components in a balanced split (old first-K was 600/0).
- Provenance: `Component/CompRank/Values` indexed by `Order(k)` are mutually consistent.

**Backfill + accessor (pure):** a legacy struct lacking `.Order` gets it filled on read, matching a fresh compute; `bst_eigenmodes_canonical(Eig,K)` returns the right columns/eigenvalues/provenance.

**Consistency / no-regression (e2e, safety net):** on TutorialAuditory —
- new canonical-path leadfield Gain **equals** the existing composed eigenmode-HM Gain;
- Eigen-MNE reconstructed through the accessor **equals** the stored, validated Eigen-MNE kernel.

**Bug-fix (behavior that should change):** `process_eigenmodes_transform` at K<1000 selects whole-brain modes (assert both components present); transform→reconstruct round-trip still recovers the field.

**ensure-if-empty (e2e):** `bst_eigenmodes_ensure` computes 1000/hemi + stores + returns canonical on a bare surface; idempotent on re-call; errors clearly (no silent repair) on a non-manifold surface.

**Regression sweep:** re-run existing eigenmode tests (`test_eigenmode_leadfield_*`, `test_inverse_eigenmodes_*`, `test_eigenmode_reconstruct_*`, `test_kernel_comparison`) green.

All follow the `dev/tests` pattern (function, `addpath`, brainstorm, `assert`, `ALL TESTS PASSED`); pure where possible, e2e for protocol-dependent equivalence checks.

## Out of scope

- **Feature B** — the eigenmode time-series ↔ cortical-activation linked viewer (complex amplitude/phase, eigenfilters, time-frequency integration). Separate spec/plan once this foundation lands.
- Re-sorting the stored `Vectors` array on disk (we use the stored `Order` index instead, leaving `Vectors` stable).
- Reorganizing the analytic transform/project into a new module beyond adopting the canonical accessor.
- Changing the inverse math, priors, or the eigenfilter library.
