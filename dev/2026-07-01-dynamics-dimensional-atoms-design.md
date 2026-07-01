# Dimensionality-aware atom construction — design

**Date:** 2026-07-01
**Status:** design (approved; spec written for review)
**Part of:** the Dynamics-portal program. Builds on the operator-default work
(`atom-operator-applicability`, commit dc09c040/802a3bd8).
**Depends on:** A/B/C/D (shipped); the D Dirac-basis-alignment guard (2c76a34d).

---

## 1. Motivation

An atom is a **kernel `g(λ)` on an operator's eigen-coefficients** — a diagonal in eigen-coefficient
space. The kernel is dimension-agnostic; it only weights the eigenvalue axis. But the **realised atom**
(drop a delta, propagate it through `g(λ)`, reconstruct) lives in the operator's per-vertex value space —
its *fiber* — and that fiber grows with the operator:

| Operator | `field_type` | fiber | impulse needs |
|---|---|---|---|
| Laplace-Beltrami | `scalar` | ℝ (1 comp) | vertex [+ amplitude] |
| Connection Laplacian | `tangent`(*) | ℝ² tangent plane (2 comp) | vertex + **tangent direction** (S¹) |
| Dirac | `quaternion` | ℍ (4 comp; data in imag 3-vector) | vertex + **3-D direction** (S²) |

(*) exact registry string for the tangent operator confirmed headless during implementation; the map
below and the `size(Phi,1)/nV` fallback are correct regardless.

**The bug this fixes.** `bst_eigenfilter('Atom', ax, kernel, kp, seedVert)` seeds the delta with
`sparse(loc, 1, 1, size(Phi,1), 1)` where `loc` is the seed's vertex index. For a scalar `Phi [n×K]`
that is a correct spatial delta. For a multi-component `Phi [C·n×K]` (interleaved `[c1..cC]` per vertex)
row `loc` only lands on the seed vertex when `loc==1`; for any other seed it hits vertex `⌈loc/C⌉`,
component `((loc-1) mod C)` — a mislocated, wrong-component seed. Verified live on a Dirac atom: with the
default seed (vertex 1) the delta fell on the **real (w) slot** of vertex 1 and diffused; for any other
seed it would be wrong. There is also **no direction argument at all**, so vector operators cannot express
the orientation their fiber requires.

**The quaternion fact (verified, informs the contract).** The Dirac eigenvectors are full quaternions —
the real (w) slot carries a robust ~25% of each mode's energy, uniformly across mode count (measured:
mean w-fraction 0.250 at K = 60/150/300). So a filtered atom is a **full quaternion**, not pure-imaginary;
reconstructing pure-imaginary data leaves a few-percent, non-vanishing real companion. The conformal
coupling is genuinely baked into the basis (the reason Dirac was chosen). Consequence for the pipeline:
**carry the full quaternion through project → filter → reconstruct; extract the imaginary 3-vector only at
the end** as the physical current (for display + the sensor forward). This is already what `i_dirac_forward`
does; the realiser and the display must do the same.

## 2. Decisions (from brainstorming)

- **Dimensionality source of truth = `ax.Operator.Registry.Primary.field_type`** (`scalar`/`tangent`/
  `quaternion` → C = 1/2/4), read via the same pattern as `tess_eigen local_registry_consistency`.
  Fallback `C = size(Phi,1)/nV` when Registry is absent (pre-registry binary).
- **Default direction** (fallback when the caller omits one): **surface normal** for Dirac (a normal/radial
  dipole — anatomically natural for cortex); **frame e1** (angle 0) for the Connection Laplacian (the
  normal has no tangent part). Both geometry-derived from stored operator/surface data.
- **v1 scope includes the GUI direction picker** (interactive orientation), not default-only.

## 3. Architecture

### 3.1 Dimensionality resolver (the one place that knows the layout)

`[C, kind, embed, decode] = i_operator_fiber(ax)`:
- `C` — components/vertex; `kind` ∈ {`scalar`,`tangent`,`quaternion`}.
- Reads `ax.Operator.Registry.Primary.field_type`; maps to C; falls back to `size(ax.Phi{1},1)/numel(ax.GlobalVertices{1})`.
- `embed(seed, dir, amp)` → the operator's DOF column (length `C·n_block`), placing the impulse at the
  **seed vertex's block** `(loc-1)*C + (1:C)`:
  - scalar: `U((loc-1)*1+1) = amp`.
  - tangent: `U((loc-1)*2+[1 2]) = amp*[cos θ, sin θ]`, θ = `dir` interpreted in the operator's `Frame`.
  - quaternion: `U((loc-1)*4+[2 3 4]) = amp*dhat` (imaginary slots), `w = 0`.
- `decode(field)` → `[nV×3]` ambient vectors for display:
  - scalar: `[]` (no vectors).
  - tangent: lift `(a,b)` per vertex through the stored `Frame` → `a*e1 + b*e2` (a 3-vector).
  - quaternion: the imaginary slots `field((v-1)*4+[2 3 4])`.

Lives in `bst_eigenfilter` (a local, exposed as a verb for tests). One implementation; `bst_eigenwavelet`
atom paths call it.

### 3.2 Realiser signature

`[W, gv] = bst_eigenfilter('Atom', ax, kernelName, kernelParams, seedVert, seedDir)`:
- Backward-compatible: `seedDir` optional. Scalar operator ignores it. Vector operator with `seedDir`
  omitted uses the **default** (§2), never a silent mis-seed.
- Body: `c0 = manifold_ft(Phi, M, embed(seed, dir, amp))`, then the existing domain-aware `g(λ)` path
  (static / ts / js) unchanged. Returns the full field `W` (full quaternion for Dirac).

### 3.3 Default-direction sourcing

- Dirac normal: from the displayed surface's per-vertex normals (`in_tess_bst(ax.SurfaceFile).VertNormals`;
  ico5 canonical carries them), fallback compute from faces. Embedded as the unit imag 3-vector.
- Connection Laplacian e1: angle θ = 0 in the operator's stored `Frame` (`OperatorMat.Frame`, the
  `(e1,e2,normal)` the eigenmodes decode in) — no extra geometry needed.

### 3.4 Atom data model

`bst_dynamics('NewGroup')` gains `G.SeedDir` (default `[]`): `[]` for scalar, scalar angle (rad) for
tangent, unit `[nx ny nz]` for Dirac. Persisted with the atom; round-trips through Save/Load.

### 3.5 GUI direction picker (panel Atom section)

Adapts to the selected operator's `kind`:
- scalar → the picker is hidden.
- tangent → an **angle** control (0–360°, interpreted in the Frame; default 0 = e1).
- quaternion → a **preset dropdown** `{Normal (default), +X, +Y, +Z, Pick-on-surface}`. `Pick-on-surface`
  arms a one-shot cortex vertex pick (reuse the existing `figure_3d` pick) and sets
  `dir = unit(target_xyz − seed_xyz)`.
Editing the direction re-realises + re-previews the selected atom (routes through `i_atom_preview`).

### 3.6 Vector display

The realised atom draws as **quivers over its norm colormap** (the user's original ask), reusing the
native `QuiverVectorOverride` mechanism (`figure_3d PlotSourceVectors` — an explicit `[nVert×3]` field over
an independent scalar colormap). `view_dynamics('SetAtomField', …)` gains a `V3` argument = `decode(W)` on
the full surface (zeros off-support); when non-empty it sets `QuiverVectorOverride` + `SetShowSourceVectors
(…,1)`; when empty (scalar atom) it clears them. `ClearAtomField` clears the override + disables vectors.
`i_paintable_scalar(W)` remains the norm colormap. Physical field = the imaginary/decoded 3-vector; the
real companion is not drawn (its magnitude is a legitimate auxiliary scalar, out of scope).

## 4. Components & files

- `toolbox/eigen/bst_eigenfilter.m` — `i_operator_fiber` resolver (+ `Fiber` verb for tests); `Atom` gains
  `seedDir`, embeds via the resolver (fixes the mis-seed).
- `toolbox/eigen/bst_eigenwavelet.m` — atom/`JTVAtoms` seeding delegates to the resolver (one impl).
- `toolbox/anatomy/bst_dynamics.m` — `G.SeedDir` field + Save/Load carry-through.
- `toolbox/gui/panel_bst_dynamics.m` — direction picker in the Atom section; `i_atom_realise`/
  `OnCreateAtom`/`i_atom_preview_impulse` pass `SeedDir` (default per operator); `i_atom_realise` returns
  the decoded `V3` alongside the scalar norm.
- `toolbox/gui/view_dynamics.m` — `SetAtomField`(+V3) / `ClearAtomField` quiver override; find/enable on the
  Design figure.
- `toolbox/gui/figure_3d.m` — reuse `QuiverVectorOverride` (no change expected; confirm the Design figure's
  `DataSource.Type=='Source'` so `PlotSourceVectors` runs).
- **Reuse:** `manifold_ft/ift`, `bst_nxr_registry`/`Op.Registry.Primary.field_type`, `OperatorMat.Frame`,
  `i_paintable_scalar`, `QuiverVectorOverride`/`SetShowSourceVectors`, the existing cortex vertex pick.

## 5. Scope & edge cases

- **In:** scalar/tangent/Dirac realise with a dimensionality-correct, direction-carrying impulse; default
  direction; the GUI picker; quiver+norm display; `SeedDir` persistence.
- **Out (v1):** face-domain operators (`Dirac-Face`, `Hodge-Face`); the conformal/real companion as a
  displayed field; a saved realised-atom file; per-atom direction in the Apply (real-data) path (Apply
  filters the real source, which already carries its own orientation — the direction concerns the *impulse*
  preview only).
- **Edge:** pre-registry binary → `field_type` absent → `size(Phi,1)/nV` fallback. Degenerate/umbilic
  normal → normal still defined (VertNormals); a zero/NaN normal falls back to `+Z`. Connection Laplacian
  is gated OFF for unconstrained sources and Dirac OFF for scalar (per `i_gate_mask`), so the tangent
  picker only appears for constrained-source sessions and the Dirac picker for unconstrained — consistent.

## 6. Testing

- **Headless (`matlab -batch`, controller-run):**
  - resolver: `field_type` → C for scalar/tangent/quaternion; fallback `size(Phi,1)/nV` when Registry `[]`;
    confirm the exact tangent `field_type` string here.
  - `embed`: a unit impulse at seed vertex `p` (p ≠ 1) places all its energy in block `(p-1)*C+(1:C)` and
    zero elsewhere (the mis-seed regression test).
  - `decode`: quaternion → imag slots; tangent → `a*e1+b*e2` matches the Frame; scalar → `[]`.
  - `Atom` with `seedDir`: for Dirac, the realised field's imaginary part near the seed aligns with the
    requested direction; default (no dir) aligns with the surface normal.
- **Live (MCP/headless):** Create Dirac atom → quivers along the normal at the seed + norm colormap; pick a
  direction → quivers reorient; a constrained-source session with the Connection Laplacian → angle control
  rotates the tangent impulse; scalar operator → no picker, no quivers (unchanged).

## 7. Risks / notes

- **Exact tangent `field_type` string** — confirm headless; the Variant→field_type map + the layout
  fallback make the resolver correct either way.
- **Connection Laplacian embedding/decoding** — the `(a,b)`↔Frame convention must match how the eigenmodes
  encode tangent vectors; test `decode(embed(·))` round-trips a known tangent vector.
- **Live-session instability** (Apple-silicon `GlobalData` drops) — prefer `matlab -batch` for the headless
  tests (per `matlab-mcp-wedged-fallback`); do live GUI checks in short, self-contained passes.
- **Backward compatibility** — `Atom`'s new `seedDir` is trailing/optional; existing scalar callers and the
  B/C cached-projection Apply path are unaffected (the change is in the impulse seed, not the kernel).
