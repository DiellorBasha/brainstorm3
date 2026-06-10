# Dirac Eigenmode Leadfield — Save + Visualize (Forward Change-of-Basis) — Design

**Date:** 2026-06-10
**Status:** Design (pending review)
**Repo:** `brainstorm3` (MATLAB). Builds on `bst_dirac_eigenmode_leadfield` (merged).

## Goal

Expose the **Dirac eigenmode leadfield** — the forward operator re-expressed in the
Dirac eigenbasis — as a first-class, inspectable artifact, and let the user
**visualize it sensor-by-sensor** the way `view_leadfield_vectors` shows a regular
leadfield. This is purely a **change of basis of the forward solution** (vertex
deltas → Dirac eigenmode source patterns); it answers "given activation of a
single Dirac eigenmode, what sensor pattern is observed," and conversely lets you
pick a sensor and see its leadfield as represented by the K Dirac modes. **No
inverse / source-from-data step** is involved.

## Background

`bst_dirac_eigenmode_leadfield(HeadModel, DiracEigen)` already composes the
unconstrained leadfield into the Dirac eigenbasis, returning `CompHM` with
`Gain [nCh × 2K]` (channels × eigenmodes — Brainstorm's native orientation),
`Eigenvalues`, `ModeHemisphere [2K×1]`, `HemiGlobalVertices {L,R}`,
`isDiracEigenmode=1`, `SurfaceFile`. This design adds:
1. a **process** that runs that composer and **saves** the result as a
   `headmodel_dirac_eigenmode_*.mat` DB node;
2. an **extension to the leadfield viewers** to display it.

## Decisions (resolved during brainstorming)

- **Storage orientation:** keep `Gain [nSensors × nEigenmodes]` — that *is*
  Brainstorm's `HeadModel.Gain` convention (channels rows, sources/modes columns);
  the inverse and viewers already expect it. No transpose.
- **Primary view:** **sensor-selected** → cortical eigenmode-reconstructed
  leadfield (vectors + sensitivity), mirroring `view_leadfield_vectors`. The
  per-eigenmode ("select a mode → sensor topography") view is deferred.
- **Save** the composed model as `headmodel_dirac_eigenmode_*.mat` via a process
  (mirroring `process_eigenmode_leadfield`).
- **Viewer approach A:** extend `view_leadfield_vectors` / `view_leadfield_sensitivity`
  with an `isDiracEigenmode` branch (localized, max reuse, regular leadfields
  untouched).
- **Forward-only:** the `Φ_D·(coeffs)` map here acts on a *leadfield row*, not on
  data — it is the leadfield band-limited to the K Dirac modes, not an inverse.

## The reconstruction (shared, pure)

For channel `c`, per hemisphere `hh` (`Φ_D,ₕ = DiracEigen(hh).Vectors(:,1:Kₕ)`,
`Kₕ = CompHM.nModes/2`, hemi columns `g_h = CompHM.Gain(c, ModeHemisphere==hh)'`):

```
psi_h = Φ_D,h * g_h            % [4Vh × 1] real quaternion field
J_h(:,1) = psi_h(2:4:end)      % x  (drop the w rows 1:4:end)
J_h(:,2) = psi_h(3:4:end)      % y
J_h(:,3) = psi_h(4:4:end)      % z
J(HemiGlobalVertices{hh}, :) = J_h     % scatter to [nVtot × 3]
```

`J [nVtot×3]` is sensor `c`'s leadfield represented by the K Dirac modes (→ the
raw `view_leadfield_vectors` field as K→all). This is extracted into a small,
**unit-testable pure helper** `bst_dirac_eigenmode_field` so both viewers and any
later consumer share one implementation:

```
J = bst_dirac_eigenmode_field(DiracEigen, gainRow, ModeHemisphere, HemiGlobalVertices, nVtot)
```

---

## Phase 1 — `process_dirac_eigenmode_leadfield` (save the artifact)

New process, sibling to `process_eigenmode_leadfield`:
1. Resolve the **base unconstrained** surface head model in the study
   (`in_bst_headmodel`, ApplyOrient=0); guard it is surface + unconstrained;
   refuse to chain on an already-eigenmode model.
2. Get `DiracEigen` from `HeadModel.SurfaceFile` (`in_tess_bst`); if absent,
   compute via `tess_dirac_eigenmodes` (options `Tau` default 0.5, `K` default
   0=all-available → use the stored `nModes`; process exposes `nModes`/`Tau`).
3. `CompHM = bst_dirac_eigenmode_leadfield(HeadModel, DiracEigen, 'nModes', K)`.
4. `OutputFile = bst_process('GetNewFilename', StudyDir, 'headmodel_dirac_eigenmode')`;
   `bst_save(OutputFile, CompHM, 'v7')`; register the DB node with `CompHM.Comment`
   (exactly the `process_eigenmode_leadfield` save/register pattern).

**Options:** `nModes` (0=all), `Tau` (0.5).

**Tests:** a process smoke run on a synthetic study is heavy; instead unit-test
the composition path is already covered (`bst_dirac_eigenmode_leadfield`). Phase 1
adds a focused test that, given a synthetic base head model struct + synthetic
`DiracEigen`, the process's compute helper returns a `CompHM` with
`isDiracEigenmode`, `Gain [nCh×2K]`, and a `headmodel_dirac_eigenmode` filename is
produced (save path verified via a temp dir / `onCleanup`).

---

## Phase 2 — viewer extension (Approach A) + the pure helper

### 2a. `bst_dirac_eigenmode_field` (pure helper)
The reconstruction above as a standalone function. **Unit-tested** synthetically:
with `Mass=I`, orthonormal `Φ_D`, and a gain row built as `Φ_Dᵀ·Ψ` for a known
embedded field `Ψ`, the returned `J` equals `imag(P·Ψ)` (the B-orthogonal
projection's imaginary part); and `J` has the right shape `[nVtot×3]` with the
correct hemisphere scatter. (Pure, fast, no GUI.)

### 2b. Extend `view_leadfield_vectors`
At load, detect `isDiracEigenmode` on the head model. When set:
- load `DiracEigen` from `SurfaceFile`; set `GridLoc` = cortex `Vertices`
  (`[nVtot×3]`);
- replace the per-channel data fetch (`LeadField = LF_finale{iLF}(iChannel,:)`
  + reshape) with `J = bst_dirac_eigenmode_field(DiracEigen, Gain(iChannel,:), …)`
  → feed `J [nVtot×3]` to the existing `quiver3` render.
All scrolling / modality / sizing / colormap logic is reused unchanged. The
branch only activates for `isDiracEigenmode` models — regular leadfields are
byte-for-byte unaffected.

### 2c. Extend `view_leadfield_sensitivity`
Same `isDiracEigenmode` branch: render the **magnitude** `‖J(v)‖` of the same
per-channel reconstructed field on the cortex/MRI, reusing the existing display
modes.

**Tests:** the pure helper (2a) carries the testable math. The viewers (2b/2c)
are GUI wiring verified **visually** (launch on the saved Dirac head model, scroll
sensors, confirm vectors/sensitivity render and converge toward the raw leadfield
as K grows). Document the manual verification steps; do not attempt headless GUI
assertions.

---

## Edge cases / notes

- **K consistency:** the viewer uses `Kₕ = CompHM.nModes/2` columns of
  `DiracEigen.Vectors` (the first `Kₕ`, matching how the composer truncated);
  `ModeHemisphere`/`HemiGlobalVertices` come from the saved `CompHM`, so the
  viewer needs `DiracEigen` only for `Φ_D` (`Vectors`).
- **Modality referencing:** EEG/SEEG average/ref subtraction in the viewers
  operates on the *reconstructed* per-channel field consistently (subtract the
  reconstructed reference field), matching the raw-leadfield behavior.
- **Stale DiracEigen:** if `SurfaceFile` lacks `DiracEigen`, the viewer errors
  clearly ("run Compute Dirac eigenmode leadfield first").

## Out of scope

- **Inverse / source-from-data reconstruction** (mode kernel × data → currents).
  Deferred (this is forward-only).
- **Per-eigenmode view** (select a mode → sensor topography). Deferred.
- **Changes to regular-leadfield behavior** of the two viewers.

## Decomposition

- **Plan 1 (Phase 1):** `process_dirac_eigenmode_leadfield` + save/register + test.
- **Plan 2 (Phase 2):** `bst_dirac_eigenmode_field` (pure helper + unit test) →
  extend `view_leadfield_vectors` then `view_leadfield_sensitivity` (visual
  verification).
