# Eigenmodes R2: Per-Hemisphere Computation + Registered-Result Viewer — Design

**Date:** 2026-05-28
**Status:** Approved design pending spec review
**Follows:** `2026-05-28-eigenmode-context-menu-compute-view-design.md` (R1, shipped on `feature/eigenmode-context-menus`)

## Background

R1 added right-click *Compute eigenmodes* / *View eigenmodes* on the cortex. Interactive validation produced five feedback items. Two required investigation:

- **Per-hemisphere (confirmed NO):** `tess_eigenmodes` builds **one** Laplacian on the whole two-hemisphere mesh (`tess_laplacian(Vertices, Faces)`, then one `eigs`, eigenvalues sorted ascending — `tess_eigenmodes.m:153,190`). The two hemispheres are disconnected components, so the operator is block-diagonal and each eigenvector localizes to one hemisphere — but they are **indexed by global eigenvalue**, so left/right modes interleave (the viewer alternates hemispheres). Worse, near-symmetric hemispheres give **near-degenerate eigenvalues**, so `eigs` can return arbitrary left/right *mixtures* — geometrically meaningless modes.
- **Colormap UI broken (root cause):** the R1 viewer renders a hand-built `DataSource.Type='Source'` overlay with no real Results file. `bst_colormaps('SetMaxCustom')` (`bst_colormaps.m:393`) calls `bst_memory('LoadResultsFile', DataSource.FileName)` for a `'Source'` overlay; our FileName is `''`, so it errors before the custom-max dialog opens. The same fake-source pattern caused the R1 render crash. The hand-built overlay fundamentally fights Brainstorm's colormap UI.

`isAbsoluteValues` and `MaxMode` are **global per colormap type**; `'source'` defaults to absolute (→ magnitude). `tess_hemisplit(sSurf)` returns clean left/right vertex sets via the "Structures" atlas. `process_eigenmodes_view` already builds a real modes-as-"time" Results file.

## Goals (the five feedback items)

1. Consolidate the two compute dialogs into **one** (number of modes + mass type).
2. After compute, **auto-open the viewer** instead of a success popup.
3. **Per-hemisphere computation** + **united visualization** (mode k shows left-k and right-k together).
4. Colormap **"Maximum: custom"** works in the viewer.
5. Viewer defaults to a **diverging** colormap (signed +/− lobes).

## Design

### A. Per-hemisphere computation — `tess_eigenmodes` + storage + I/O

`tess_eigenmodes` solves **per connected component** (each disconnected component is a separate manifold — the correct LBO unit), instead of one whole-mesh operator:

1. Validate manifold (unchanged).
2. **Find connected components** of the mesh (vertex adjacency from `Faces`; e.g. `conncomp(graph(...))` or an existing Brainstorm connectivity helper).
3. For **each component**: assemble `L`,`M` on that component's sub-mesh, solve `nModes` eigenmodes, remove the single DC mode, M-orthonormalize.
4. **Assemble** outputs across components:
   - `Vectors [nV × (nModes·nComp)]` — each column one component's eigenmode, zero-padded on other vertices (**same format as today** → downstream basis usage unchanged).
   - `Values [k×1]` — per-mode eigenvalue.
   - **New** `Component [k×1]` — component id per mode; `CompRank [k×1]` — within-component rank (1-based).
   - Columns ordered grouped-by-component then by rank (so `CompRank` pairing is direct). `Values` are per-component-ascending.
5. `nModes` now means **modes per component** (per hemisphere). Cap per component at `nCompVerts − 2`.

Single-component surfaces (one hemisphere, non-cortex meshes) yield `Component` all-ones and behave exactly as before (**backward compatible**).

`out_tess_eigenmodes` / `in_tess_eigenmodes` persist `Component` and `CompRank` (with safe defaults — treat a missing field as a single component — so older files still load).

### B. Registered-result viewer — `view_eigenmodes` rewrite (fixes #3-viz, #4, #5)

The viewer stops hand-building an overlay. Instead it builds a **real Brainstorm Source result** and displays it through the standard surface-data path, so the entire colormap UI works natively.

- **Paired display matrix:** `ImageGridAmp [nV × nModesPerComp]`, where **column k = Σ over components of that component's CompRank==k mode** (disjoint support → left-k on left vertices, right-k on right vertices). Stepping to mode k shows both hemispheres' k-th mode together (#3 united).
- **Display via the standard path** (`Time = 1..nModesPerComp` as the mode axis); the source figure's **←/→ step modes natively**.
- **#4:** custom-max and all colormap menus work because it is a properly-registered result.
- **#5:** set the result's `ColormapType` to a **non-absolute diverging** colormap (e.g. `'stat2'`/`cmap_mandrill`) so signed lobes show by default, **without** globally flipping `'source'`. (Trade-off: diverging types are shared globally like all Brainstorm colormaps; a dedicated `'eigenmode'` colormap type is a possible future refinement, noted not built.)
- **Legend:** `Mode k / K — λ(L)=… , λ(R)=… — wavelength ~ …`, updated as the mode changes.
- **Lifecycle:** the result is created in the subject's intra/Analysis study, displayed, and **auto-removed when the viewer figure closes** (chained `CloseRequestFcn` → delete file + node + refresh tree). *Implementation risk:* the figure/dataset close lifecycle is fiddly; if auto-clean proves unreliable, the documented fallback is a persistent, clearly-named node the user can delete. This is confirmed during interactive validation.
- No-eigenmodes path: friendly error pointing to Compute (unchanged behavior).

### C. Consolidated compute dialog — custom panel (#1)

Replace the two sequential `java_dialog` calls in `ComputeInteractive` with **one modal panel** (built with `gui_river`/`gui_component`, shown via `gui_show_dialog`): a numeric field "Number of eigenmodes per hemisphere" (default 300) and a **mass-matrix dropdown** (barycentric / voronoi / galerkin), with OK/Cancel. Cancel → no-op. Overwrite confirm unchanged.

### D. Auto-open viewer after compute (#2)

`ComputeInteractive` drops the success `msgbox`; after a successful `Compute` + `db_reload_subjects`, it calls `view_eigenmodes(SurfaceFile)` — the rendered modes are the confirmation.

## Data flow

```
Compute: right-click cortex -> ComputeInteractive
  -> custom panel (nModesPerHemi, massType)
  -> Compute(SurfaceFile, ...) -> tess_eigenmodes (per-component) -> out_tess_eigenmodes (+Component/CompRank)
  -> db_reload_subjects -> view_eigenmodes(SurfaceFile)

View: view_eigenmodes(SurfaceFile)
  -> in_tess_eigenmodes (Vectors/Values/Component/CompRank)
  -> BuildPairedGrid: ImageGridAmp[:,k] = sum_c modeOf(component=c, CompRank=k)
  -> create temp Source result (ColormapType=diverging) -> view_surface_data
  -> legend + native arrow stepping
  -> figure close -> auto-remove the temp result
```

## Error handling

- View without eigenmodes → friendly error → Compute.
- Compute dialog / overwrite cancel → no-op.
- Non-manifold (repair off) → existing `Compute`/`tess_manifold` path (interactive prompt).
- Component with fewer vertices than requested modes → cap per component (existing `maxModes` logic, applied per component).
- `in_tess_eigenmodes` on an old file without `Component`/`CompRank` → default to single component (rank = column index).
- Viewer auto-clean failure → leave a named node rather than erroring (fallback).

## Testing strategy

Repo idiom: `dev/tests/*.m` printing `ALL TESTS PASSED`, run via the MATLAB MCP `evaluate_matlab_code` (not `run_matlab_test_file`).

1. **Per-component `tess_eigenmodes`** (`test_eigenmodes_perhemisphere.m`, new): build a synthetic **two-component** manifold mesh (two disjoint patches). Assert: total columns = nModes·2; each column localized to one component (≈0 on the other); `Component`/`CompRank` correct; one DC removed per component; per-component M-orthonormality. Also a **single-component** case → behaves as before.
2. **I/O round-trip** (extend `test_io_eigenmodes_roundtrip.m`): `Component`/`CompRank` persist; an eigenmodes struct missing them loads with single-component defaults.
3. **Paired-grid helper** (pure, in `view_eigenmodes`): given an Eigenmodes struct, assert `ImageGridAmp(:,k)` equals left-rank-k on left vertices + right-rank-k on right vertices, zero elsewhere; K = nModesPerComp.
4. **Downstream regression:** re-run `test_process_eigenmodes_{transform,denoise,dispersion,wavelet,wiener,coeffsfilter}_options`, `test_eigenmodes_{transform,project,filter,filter_gain,noisefloor,dispersion,wavelet,wiener}_pure`, `test_io_eigenmodes_roundtrip`, `test_eigenmodes_manifold_gate`, `test_process_eigenmodes_options` — confirm green.
5. **Interactive (user):** consolidated dropdown dialog; auto-open after compute; united stepping (both hemispheres per mode); custom-max dialog opens/works; diverging default; auto-remove on close.

## Downstream impact & migration

Downstream processes consume `Vectors [nV×K]` / `Values [K]` as a basis and are unaffected *as a basis*. Total basis size grows (nModes per component × components); the transform's auto-mode selection still caps at channel count. Existing on-disk eigenmodes must be **recomputed** to gain per-hemisphere structure + metadata (old files still load via single-component defaults).

## Out of scope (YAGNI)

- A dedicated `'eigenmode'` global colormap type (reuse an existing diverging type).
- Surface types other than cortex for the menu (unchanged from R1).
- Changing the downstream processes' own behavior beyond confirming they still work.
- Per-hemisphere *eigenvalue spectra* analysis features (this round is computation + viewer only).

## Files touched

- `toolbox/anatomy/tess_eigenmodes.m` — per-component solve + `Component`/`CompRank`.
- `toolbox/io/out_tess_eigenmodes.m`, `toolbox/io/in_tess_eigenmodes.m` — persist/load metadata.
- `toolbox/gui/view_eigenmodes.m` — rewrite to registered-result + paired grid + diverging colormap + auto-clean + legend.
- `toolbox/process/functions/process_eigenmodes.m` — `ComputeInteractive`: custom panel dialog + auto-open viewer (drop msgbox).
- `dev/tests/test_eigenmodes_perhemisphere.m` (new); extend `dev/tests/test_io_eigenmodes_roundtrip.m`; new/updated viewer pure-helper test.
