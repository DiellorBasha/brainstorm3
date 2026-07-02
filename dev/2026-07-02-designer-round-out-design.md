# SP2a — Round out the atom designer to the full operator range (design)

**Date:** 2026-07-02 · **Status:** design (approved in brainstorm; pending spec review)
**Context:** part of SP2 (concerns C3 — one design surface; the analytic, data-free atom design lives in the
designer, not the Dynamics panel). SP2 = **round out the designer (SP2a, this spec)** then **excise the
design half from Dynamics + embed the shared designer panel (SP2b, next spec)**. SP2a is independently
shippable — the standalone designer gains all operators immediately. Builds on **SP1** (metadata-driven
`bst_eigen`: `bst_eigen('Axes', variant)` now routes any operator by `field_type`, and `Dirac-Connectome`
is a factory operator).

---

## 1. Problem

`view_atom_designer` (+ its docked `panel_atom_designer`) is **scalar-only**: the operator selector offers
just Geometric (`Laplace-Beltrami`) / Connectomic (`LB-Connectome`); `i_eval_atom` returns a scalar field
(`bst_eigenfilter('Atom')` → `[W, gv]`, painted as `W(gv,:)`); there is no vector (quaternion) decode and no
seed-direction control. The richer design logic (Dirac/Dirac-Connectome operators, the quaternion/tangent →
V3 decode, the seed-direction control) exists only in `panel_bst_dynamics` (`i_atom_realise_core`,
`i_atom_default_dir`, `jDirCombo`), duplicated and diverged. SP1 made routing metadata-driven, so the
designer can now support the full range — it just needs the UI + the shared realise.

## 2. The Atom / Fiber taxonomy (naming reorganization)

**An atom is one supercategory; its scalar/vector nature is its *fiber*, a descriptor — never a separate
verb.** There is NO `AtomVec`. Two questions characterize an atom: **(1) `Atom`** (realise it) and **(2)
`Fiber`** (which fiber — `scalar` / `complex-tangent` / `quaternion`, driven by `field_type`).

- **`bst_eigenfilter('Fiber', ax)`** — RETAINED. Reports the ax's fiber (`[C, kind]`). The subcategory descriptor.
- **`bst_eigenfilter('Atom', ax, kernel, kp, seed, seedDir)`** — the ONE realiser. Extended to return the
  atom **already decoded for its fiber**: `[W, gv, V3, isSigned]` where `W` is the raw modal field, `V3` is
  the decoded physical 3-vector (`[]` for the scalar fiber; imag-3-vec for quaternion; `a·e1+b·e2` for
  tangent), `isSigned`/`kind` per fiber. It consults `Fiber` internally. Scalar and vector atoms come out of
  the **same** verb. Backward-compatible: existing `[W, gv]` callers are unaffected (extra outputs only).

The quaternion/tangent decode currently stranded in `panel_bst_dynamics` `i_atom_realise_core` (`:679–697`)
**moves into `bst_eigenfilter('Atom')`** — it is part of what realising an atom *means*. `i_atom_realise_core`
is then updated to **consume `Atom`'s `V3`** (thin; no longer decodes itself). Other `Atom` callers that only
read `[W, gv]` are unaffected (extra outputs only). SP2b later removes `i_atom_realise_core` entirely.

## 3. Round out the designer (data-free, all operators)

1. **Operator selector** (`panel_atom_designer`): replace the 2 toggles with the 4 SP1-routable operators —
   Geometric (`Laplace-Beltrami`), Connectomic (`LB-Connectome`), **Dirac**, **Dirac (connectome)**. On change,
   `view_atom_designer` `i_build_basis(variant)` (already `bst_eigen('Axes', variant)`) rebuilds the axes —
   now works for the vector operators via SP1. (Connection Laplacian stays OUT — its tangent frame is
   unpersisted/incomplete.)
2. **Fiber-driven realise + render** (`view_atom_designer` `i_eval_atom`): call `bst_eigenfilter('Atom')`,
   take its `V3`. For the scalar fiber → keep the existing density/peak scalar paint. For the quaternion
   fiber → paint the per-vertex **magnitude** (`source` colormap) + **quivers** from `V3` (reuse the panel's
   quiver display idiom via `figure_3d('SetShowSourceVectors')` + a vector override). `Fiber` decides which.
3. **Seed-direction control** (`panel_atom_designer`): add the direction control (Dirac quaternion presets
   Normal/+X/+Y/+Z/Pick-on-surface; hidden for the scalar fiber), mirroring the panel's `jDirCombo`. The
   **default direction is app-side** (never the library): relocate `panel_bst_dynamics`'s local
   `i_atom_default_dir` into a **standalone shared GUI helper `toolbox/gui/bst_atom_default_dir.m`**
   (`dir = bst_atom_default_dir(ax, seed)`) that both the designer and (SP2b) Dynamics call; update
   `panel_bst_dynamics` to call the shared one (local copy deleted — a small SP2a-scoped change to avoid
   duplication, byte-equivalent). `Pick-on-surface` resolves `unit(target − seed)` via the existing `WaveletDesignerPick`
   hook (which the designer already uses for seed-picking — no geodesic Region tool, cf. C6).
4. **Scale/physical calibration** unchanged (`view_atom_designer` already calibrates mm via `2π/√λ`); the
   vector operators reuse the same kernel/param sliders (`bst_eigfilter_controls`).

## 4. Components

| File | Change |
|---|---|
| `toolbox/eigen/bst_eigenfilter.m` | extend `Atom` to return `[W, gv, V3, isSigned]` fiber-decoded (fold in the quaternion/tangent decode; consult `Fiber`). `Fiber` unchanged. |
| `toolbox/gui/view_atom_designer.m` | `i_eval_atom` uses `Atom`'s `V3`; render magnitude+quivers for the quaternion fiber, scalar paint for the scalar fiber; pass `seedDir`. |
| `toolbox/gui/panel_atom_designer.m` | operator selector → 4 operators; add the seed-direction control (hidden for scalar fiber); `CurrentOperator` returns the 4. |
| `toolbox/gui/bst_atom_default_dir.m` (new) | app-side default-direction helper (relocated from `panel_bst_dynamics`'s local `i_atom_default_dir`); both the designer and (SP2b) Dynamics call it. |
| `toolbox/gui/panel_bst_dynamics.m` | (1) delete the local `i_atom_default_dir`; call `bst_atom_default_dir`. (2) `i_atom_realise_core` consumes `Atom`'s `V3` instead of decoding itself. Both byte-equivalent; the only SP2a touches to this file. |

## 5. Interfaces / isolation
- `bst_eigenfilter('Atom')` is the single realiser; `Fiber` the fiber descriptor. The designer never
  branches on operator NAME — it branches on `Fiber` (scalar vs vector), so a new operator with a known
  fiber renders with no designer edit.
- The library takes an **explicit `seedDir`, no default**; the app supplies the default (`i_atom_default_dir`).
- Data-free: SP2a touches nothing about the inverse kernel / sensor data (that is Dynamics' concern, SP2b).

## 6. Testing
- **Unit (headless):** `bst_eigenfilter('Atom')` fiber-decoded output (`[W, gv, V3]`) == the panel's existing
  `i_atom_realise_core` output for the scalar fiber AND the Dirac quaternion fiber (byte-equivalence, the
  decode simply relocated). Scalar `V3=[]`; quaternion `V3` = imag-3-vec.
- **Live:** open `view_atom_designer(surf)`; for each of the 4 operators, drop a seed → the impulse renders
  (scalar density, or vector magnitude+quivers for Dirac/Dirac-Connectome); the direction control reorients
  the Dirac atom (quiver `|dot|≈1` to the chosen preset) and is hidden for scalar operators; Pick-on-surface
  reorients to `unit(target−seed)`.

## 7. Out of scope (SP2b, next spec)
Excising `panel_bst_dynamics`'s design half; embedding `panel_atom_designer` as a shared docked component in
the Dynamics session bound to the session's filterbank + cortex figure; sharing the current-atom spec. In
SP2a `panel_bst_dynamics`'s design half stays functional — its ONLY SP2a change is the `i_atom_default_dir`
relocation (§4); its `i_atom_realise_core` keeps working via the backward-compatible `Atom` outputs (SP2b removes it).

## 8. Risks
- The quaternion-decode relocation into `bst_eigenfilter('Atom')` must be byte-equivalent to the panel's
  current decode (pinned by the unit test) so `panel_bst_dynamics` (unchanged in SP2a) keeps behaving.
- The designer's vector render (quivers on its managed working-results overlay) must not fight the native
  source-vector display — reuse the panel's established quiver idiom.
