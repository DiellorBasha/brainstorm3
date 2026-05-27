# Consolidation: `tess_check_manifold` + `tess_fix_manifold` → `tess_manifold`

Folds the two manifold utility functions developed for the icosphere downsampling work
(see `dev/ico-downsize.md`) into a single check-and-repair entry point. Pure consolidation —
the validation and repair **algorithms are copied verbatim**; no behavioral change.

## API

```matlab
[Vertices, Faces, isManifold, report] = tess_manifold(Vertices, Faces, ...)
```

- **Default (`'Repair', 0`)** — validation only; `Vertices`/`Faces` returned **unchanged**.
  A pure check is `[~,~,ok] = tess_manifold(V,F)`.
- **`'Repair', 1`** — validate first, repair **only if defective** (no-op on a clean mesh),
  re-validate, return the fixed mesh. `isManifold` reflects the **final** (post-repair) state.

Options (union of both former functions, names preserved): `Repair` (new, default 0),
`RequireClosed`, `RequireConnected`, `Verbose`, `AreaTol`, `DupVertexTol`, `MaxIter`.

The `report` top level always describes the **returned** mesh (same fields as the former
`tess_check_manifold`). When repair runs, `report.repair` is added:
`.performed`, `.nFixed`, `.iterations`, `.removedFaces`, `.removedVertices`,
`.nConsistencyFlips`, `.nOutwardFlips`, `.nFlippedFaces`, `.validationBefore`
(the pre-repair validation report). `report.repair` is absent when `Repair=0`.

## Internal structure (one file, `toolbox/anatomy/tess_manifold.m`)

- Main: parse opts → `doCheck` → if `Repair && ~ok` → `doRepair` → re-`doCheck` → assemble report.
- `doCheck` — the 10 validation checks (silent; the former `tess_check_manifold` body).
- `doRepair` — non-manifold face removal, orphan cleanup, BFS orientation + outward-normal pass
  (the former `tess_fix_manifold` body, minus its trailing self-validation — the main function
  re-validates).
- Shared helpers merged/de-duplicated: `computeFaceNormalsAndAreas` (was two area helpers),
  `fixOrientation`, `fillDefaults`, `printSummary`.

## Callers updated (all 3 dependents)

| File | Change |
|---|---|
| `toolbox/anatomy/tess_eigenmodes.m` | `FixMesh` block: replaced the check + conditional fix two-step with one `tess_manifold(V,F,'Repair',1,'Verbose',Verbose)` call. Help + SEE-ALSO updated. |
| `dev/tests/test_omega_icosphere_sourcemap.m` | `[~,~,chk.isManifold] = tess_manifold(...)`. |
| `toolbox/anatomy/tess_laplacian.m` | SEE-ALSO comment → `tess_manifold`. |

`tess_check_manifold.m` and `tess_fix_manifold.m` deleted (`git rm`). No remaining `.m`
references to the old names (the only mention is the historical note in the new file's docstring).

## Verification (live, MATLAB R2023b, dev Brainstorm, protocol `TutorialOmega_IcoTest`)

- `checkcode`: clean on the new file and all 3 edited callers (only the two pre-existing
  `isempty(find)` / `find` style hints carried over verbatim from the originals).
- Old names removed: `exist('tess_check_manifold')` = `exist('tess_fix_manifold')` = 0.
- **Icosphere ico5 cortex** (`sub-0002-test`, 20484 V): CHECK → `isManifold=1`, mesh unchanged,
  no `.repair` field. `Repair=1` → **no-op** (`performed=0`, mesh unchanged).
- **Reducepatch cortex** (`sub-0002-test-rp`, 15002 V): CHECK → `isManifold=0`, mesh unchanged.
  `Repair=1` → repaired (11 faces removed, orientation fixed) → final `isManifold=1`
  (14996 V, 29974 F); `validationBefore.ok=0`.
- **`tess_eigenmodes` FixMesh path** end-to-end on the reducepatch cortex: repair produced a
  valid 2-manifold, then 10 eigenmodes solved cleanly (2 DC modes removed, post-normalization
  orthogonality error 9e-16, eigenvalues finite/ascending, L/M assembled).

Behavior is identical to the former two functions.
