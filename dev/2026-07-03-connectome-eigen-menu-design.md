# Connectome eigenmodes: GUI menu + Subject01 eigen files

**Date:** 2026-07-03
**Branch:** development

## Goal

1. Expose the two connectome operators (`LB-Connectome`, `Dirac-Connectome`) — already
   registered in `tess_operators`/`tess_eigen` — in the cortex right-click menus.
2. Compute those two eigen files for **Subject01**'s cortex (K=400, Tau=0.5), reusing the
   default subject's HCP-1065 fibers via FreeSurfer-sphere registration.

## Background (verified)

- `tess_eigen.m` variant dispatch (lines 104–124) already maps `'LB-Connectome'` and
  `'Dirac-Connectome'`. `tess_operators.m` has the parallel dispatch.
- `LB-Connectome` → `local_build_connectome_operator` → `tess_connectome`, which uses the
  subject's own Fibers surface if present, else falls back to the **default subject's
  HCP-1065** fibers registered onto the subject cortex (`tess_connectome.m:84–107`).
- `Dirac-Connectome` is a quaternion **lift** of the LB-Connectome base (output = 3×K
  columns); auto-builds/reuses the LB-Connectome operator + eigen base.
- Existing cortex menus in `tree_callbacks.m`: `Compute operator` (1182–1186) and
  `Compute eigenmodes` (1188–1194), each listing Laplace-Beltrami / Connection Laplacian /
  Dirac, guarded by `strcmpi(nodeType,'cortex') && ~bst_get('ReadOnly')`.
- Interactive mode uses default K=400 / Tau=0.5; it only prompts on an exact-spec overwrite.

## Part A — GUI menu (additive, `toolbox/tree/tree_callbacks.m`)

Add two `MenuItem` lines to **each** existing submenu, wired exactly like the neighboring
`Dirac` item (same icon, `bst_call` idiom):

- `Compute operator` (after 1185): `LB-Connectome`, `Dirac-Connectome` → `@tess_operators`.
- `Compute eigenmodes` (after 1193): `LB-Connectome`, `Dirac-Connectome` → `@tess_eigen`.

No backend edits.

## Part B — Compute Subject01 eigen files (live)

1. Start Brainstorm in the MCP MATLAB session; load the protocol holding Subject01.
2. Confirm default subject carries HCP-1065 fibers; locate Subject01's cortex file.
3. Run:
   ```matlab
   tess_eigen(cortexFile, 'LB-Connectome',    'nModes', 400, 'Tau', 0.5);
   tess_eigen(cortexFile, 'Dirac-Connectome', 'nModes', 400, 'Tau', 0.5);
   ```
4. Confirm two `eigen_*` nodes under Subject01's cortex; report nModes + fiber provenance.

## Verification

- **GUI:** reload tree; both items render under a cortex node's two submenus, absent on
  non-cortex surfaces.
- **Data:** `bst_get('EigenFileForSurface', ...)` returns both variants; LB-Connectome
  nModes=400, Dirac-Connectome lifted columns=1200; `Provenance.fibers` =
  `default:HCP-1065 (registered)`.

## Scope guard

No refactoring of the eigen/operator backends. Menu change is two-times-two additive lines.
