# SP1 — Metadata-driven eigen stack + Dirac-Connectome as a factory operator (design)

**Date:** 2026-07-02 · **Status:** design (approved in brainstorm; pending spec review)
**Context:** foundation for extracting the analytic (data-free) atom designer out of the Dynamics panel
(concerns doc `dev/2026-07-02-dynamics-panel-concerns.md` C3; north-star N1). "Round out the designer"
(SP2) rides on this. This spec is **SP1 only** — the eigen-stack foundation. SP2 (the designer) is a
separate spec.

---

## 1. Problem

`bst_eigen` routes eigenmode ↔ data-type by the operator's **variant name string**, not by the operator's
properties. `switch EigenMat.Variant` (and sibling name-matches) appear at **four sites** in
`toolbox/eigen/bst_eigen.m`:

| Site | Role |
|---|---|
| `:398–463` main `switch EigenMat.Variant` | embed source map `F` → operator's native row layout `U_h` |
| `:220` `strcmp(Variant,'Connection Laplacian')` | complex handling in the filter path |
| `:238` `any(strcmp(Variant,{'Dirac','Dirac-Face','Hodge-Face'}))` | project quaternion result → 3-vector |
| `:541 / :590` same name-list | set `nComponents` (3 vs 1) on the output results file |

Each case only actually needs the **field kind** (component count + embedding) and the **index domain**
(vertex vs face). The name is a proxy for those, and it is a **strict allowlist with `otherwise → error`**,
so operators that are *layout-identical* to supported ones are rejected:

- `LB-Connectome` is scalar (same layout as `Laplace-Beltrami`) → falls to `otherwise` →
  `"unsupported eigen variant 'LB-Connectome'"` (observed live).
- `Dirac-Connectome` is quaternion (same layout as `Dirac`) → cannot route through `bst_eigen` at all →
  the panel builds it with an **inline lift** (`bst_lift_connectome_dirac` in `panel_bst_dynamics`
  `i_atom_axes`), off the factory path, with **no operator node and no metadata**.

Meanwhile the metadata **already exists but is unused here and incomplete**: `bst_nxr_registry` maps each
variant → a registry id whose `field_type` (`real`/`complex`/`quaternion`) comes from the nxr binary, and
`tess_operators` stamps it onto `OperatorMat.Registry.Primary` (`tess_operators.m:462`). But `bst_eigen`
ignores it, and two operator families never get stamped: `LB-Connectome`/`Connectome Laplacian` (built in
the connectome early-return branch `tess_operators.m:167–173`, before the stamping) and `Dirac-Connectome`
(no node at all).

## 2. Architecture / layering (the frame)

- **`tess_operators` = the operator factory.** Find-or-create for *every* operator, including **composed**
  ones (e.g. Dirac-Connectome, built by composing find-or-create'd sub-operators). Stamps registry metadata.
- **`tess_eigen` = the eigen-file factory.** Find-or-create an eigenbasis from an operator; **structure-aware**
  for lifted/composed operators (decompose the small base, lift — never the large lifted operator).
- **`bst_eigen` = the metadata-driven orchestrator.** Routes eigenmode ↔ data-type by operator metadata;
  **no `switch Variant`, no UX defaults.**
- **Applications own visualization defaults** (e.g. a default seed direction) — never the library (SP2).

## 3. The metadata model — a two-axis `(field_type, domain)` field spec

Routing is driven by two metadata axes read from the operator, not by the name:

- **`field_type` ∈ {`real`, `complex`, `quaternion`}** → component count, source `nComponents`, and embedding.
- **`domain` ∈ {`vertex`, `face`}** → which index set (`GlobalVertices` vs `GlobalFaces`). **Default `vertex`**
  if metadata does not supply it (vertex is Brainstorm's universal domain; this is a completeness fallback,
  not a UX default).

Derived field spec (single source of truth; replaces every name-switch):

| `field_type` | internal width / vertex (or face) | source-map `nComponents` | embedding of `F` → `U_h` |
|---|---|---|---|
| `real` | 1 | 1 | identity: `U_h = F(idx,:)` |
| `complex` | 1 (complex) | 1 (complex) | tangent in operator frame (requires complex `F`; see §7) |
| `quaternion` | 4 | 3 | imag-slot embed `[0; x; y; z]` per element |

`domain=vertex` indexes `GlobalVertices`; `domain=face` indexes `GlobalFaces` (same embedding, face count).
So `Dirac` = `(quaternion, vertex)`, `Dirac-Face`/`Hodge-Face` = `(quaternion, face)`,
`Laplace-Beltrami`/`LB-Connectome` = `(real, vertex)`, `Dirac-Connectome` = `(quaternion, vertex)`,
`Connection Laplacian` = `(complex, vertex)`.

## 4. Components

### A. Complete the operator metadata
Every built operator exposes a non-empty `Registry.Primary.field_type` (+ `domain`):
- Add registry ids for `LB-Connectome` and `Connectome Laplacian` in `bst_nxr_registry`'s `local_map`
  (`field_type='real'`, `domain='vertex'`), and ensure the connectome early-return branch in
  `tess_operators` (`:167`) reaches the registry-stamping (move/duplicate the stamp before the return).
- `Dirac-Connectome` gets `field_type='quaternion'` via component B.
- Verify whether the nxr `fieldInfo` already carries a vertex/face **support/domain**; if yes, read it; if
  not, add `domain` to the registry descriptor (small, additive). Absent → default `vertex`.

### B. Dirac-Connectome as a first-class factory operator
- **`tess_operators(surf,'Dirac-Connectome')`** — find-or-create: locate/build the **LB-Connectome**
  operator, then produce a Dirac-Connectome operator node = the quaternion lift `L_conn ⊗ I₃`, stamped
  `Registry.Primary.field_type='quaternion'`, `domain='vertex'`, `Registry.Components=[LB-Connectome]`.
- **`tess_eigen(surf,'Dirac-Connectome')`** — **structure-aware**: find-or-create the LB-Connectome **eigen
  file**, then lift its modes (relocate `bst_lift_connectome_dirac`: `Φ_q`, `Λ_q=repelem(Λ,3)`,
  `M_q=kron(M,I₄)`; whole-brain single block). Decompose the small base + lift — preserves the
  no-3×-eigensolve win.
- **Retire** the inline lift in `panel_bst_dynamics` `i_atom_axes`; it calls `bst_eigen('Axes','Dirac-Connectome')`
  like every other variant.

### C. `bst_eigen` routes by metadata, not name
- Add a shared helper, e.g. **`bst_eigen('FieldSpec', ax)`** → `struct('field_type',…, 'domain',…, 'width',…,
  'nComponents',…, 'embed',@…)`, reading `ax.Operator.Registry.Primary.field_type` (+ `domain`, default
  `vertex`). Promote `i_fiber`'s registry read into this one helper so `bst_eigen` and `bst_eigenfilter`
  share it.
- Replace all four name-sites (`:398` U-construction, `:220` complex, `:238` project-back, `:541/:590`
  `nComponents`) with `FieldSpec`-driven dispatch. The `otherwise → error` becomes "operator has no/unknown
  `field_type`" — a metadata error, not a name-allowlist rejection.
- Keep `i_fiber`'s Phi-layout inference **only** as a guarded fallback for genuinely pre-registry binaries.

### D. Interfaces / isolation
`(field_type, domain)` is the single source of truth for eigenmode ↔ data-type. `tess_operators` owns
construction, `tess_eigen` owns eigenbases (structure-aware for lifts), `bst_eigen` owns routing. No
operator carries a UX default; no application logic leaks into the library. A new operator becomes routable
by declaring its `field_type` (+ `domain`) — no `bst_eigen` edit.

### E. Testing
- **Unit — factory equivalence:** `tess_operators`/`tess_eigen` `Dirac-Connectome` eigenbasis == today's
  panel inline lift (`Φ`, `Λ`, mass) to ~1e-12.
- **Unit — routing regression:** `bst_eigen` `FieldSpec`-driven `U_h` / `nComponents` == the old
  `switch Variant` output for **every** existing variant (`Laplace-Beltrami`, `Dirac`, `Connection
  Laplacian`, `Dirac-Face`, `Hodge-Face`).
- **Unit — metadata completeness:** every operator built by `tess_operators` exposes a non-empty
  `field_type`; `domain` absent ⇒ resolves to `vertex`.
- **Integration:** the Dynamics panel (now via `bst_eigen('Axes','Dirac-Connectome')`) still produces the
  fiber-spread atom (contra-hemi spread), Apply/Analyze/Helmholtz unchanged.

## 5. Sequencing (within SP1)
1. Metadata helper + `FieldSpec` (C) with the pre-registry fallback — no behavior change yet.
2. Complete metadata (A) — LB-Connectome/Connectome-Laplacian stamped `real`.
3. Swap `bst_eigen`'s four name-sites to `FieldSpec` (C) — regression-tested against the old switch.
4. Dirac-Connectome factory operator + structure-aware eigen (B); retire the panel lift.
Each step is independently testable; 1–3 are pure de-duplication (no new capability), 4 relocates existing math.

## 6. SP2 (next spec, sketch)
Designer reads `field_type` to render scalar/tangent/quaternion generically; operator dropdown expands to the
4 working operators; the seed-direction **default** is added app-side in the designer (never the library);
then the Dynamics design half is excised.

## 7. Risks / open items
- **Connection Laplacian is complex + frame-not-persisted.** Its embedding needs the per-vertex tangent
  frame, not yet stored (`bst_eigen.m:436`). SP1 keeps it in the metadata model (`field_type='complex'`) and
  routes it generically, but the frame-persistence gap is **out of scope** (it's why Connection Laplacian is
  excluded from SP2's designer). Note the existing inconsistency: `:541` treats it as `nComp=1` while `:590`
  treats it as vector — the `FieldSpec` model must resolve this one way (proposed: `complex` → `nComponents=1`
  complex field) and the regression test will pin it.
- **`fieldInfo` domain availability** — verify the nxr descriptor carries vertex/face; add `domain` if absent.
- **Quaternion width-4 vs source nComponents-3** — `field_type='quaternion'` means internal width 4 but source
  map `nComponents=3`; the `FieldSpec` must carry both (`width` for `U_h`, `nComponents` for the results file).
