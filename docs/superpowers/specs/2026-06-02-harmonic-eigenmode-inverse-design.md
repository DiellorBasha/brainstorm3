# Harmonic Eigenmode Inverse + Consistent Time Series — Design

**Date:** 2026-06-02
**Author:** Diellor Basha (with Claude)
**Status:** Approved design — ready for implementation plan
**Builds on:** `2026-06-02-eigenmode-timeseries-design.md` (the time series viewer, merged),
`bst_inverse_eigenmodes.m` / `process_eigenmodes_inverse.m` (existing eigenmode inverse),
`bst_eigenmodes_transform.m` (rank-safe unregularized pseudoinverse).

## Goal

Register a **Harmonic** source-imaging method: the eigenmode reconstruction **after noise
whitening but with no source/inverse regularization** (no Tikhonov `λ`, no spatial prior).
It produces a normal results node you display on the cortex, and the eigenmode time series
reads that node's operator so the traces are **exactly consistent** with the displayed map.

This replaces the current launch path, where the time series was opened from a regularized
dSPM/MNE figure and computed its own *unwhitened* lead-field transform — a quantity that did
not correspond to anything shown on the cortex.

## Conceptual model

The eigenmode method is a **filter bank**, not an auto-regularized inverse. The Harmonic node
is the **raw bank**: the unregularized (whitened) eigenmode reconstruction. Spectral filtering
(the `EigenModes` panel band/window) is a *separate*, optional step applied on top for
regularized mapping — out of scope here, and the lever's source-map filtering already exists.

## The math

With constrained lead field `L = Gain` `[nCh × nVert]`, eigenvectors `Φ` `[nVert × K]`, and the
noise whitener `iW` `[nCh × nCh]` (from the noise covariance):

```
M̃ (EigenKernel) = pinv( iW·L·Φ ) · iW         [K × nCh]     eigenmode-space, whitened, unregularized
ImagingKernel    = Φ · M̃                        [nVert × nCh] vertex-space → cortex display
θ(t)             = M̃ · D(t)                     eigenmode time series
u(t)             = ImagingKernel · D(t) = Φ·θ(t) displayed cortex map
```

`Φ·θ(t) = u(t)` holds **by construction**, provided the time series uses the *same* `M̃` (same
whitener, same rank flooring) as the node — hence the node stores `M̃` and the viewer reads it.

`pinv` is the rank-safe SVD pseudoinverse from `bst_eigenmodes_transform` (small singular values
floored — high spatial modes are nearly invisible to the sensors). "Unregularized" means **no
source regularization**; the noise-covariance conditioning that builds `iW` keeps its standard
stabilization (it makes the noise inverse well-posed, it is not a source prior).

## Architecture (4 units)

1. **`harmonic` method (`bst_inverse_eigenmodes.m`).** Add `Method='harmonic'`. After building
   the whitener `iW` (existing code path), compute `Kt = bst_eigenmodes_transform(iW*L, Φ)`
   (= `pinv(iW·L·Φ)`, `[K×nCh]`), then `Kernel = Kt * iW`. `EigenGains = ones(K,1)`,
   `SourcePrior = ones(K,1)`. The `SNR` / `PriorAlpha` options are ignored for this method. The
   returned `Results.ImagingKernel` stays **eigenmode-space** `[K×nCh]` (as for the other
   methods); `Results` already carries `Whitener`, `nModes`, `GoodChannel`, `Eigenvalues`,
   `SurfaceFile`. If `NoiseCovFile` is empty, `iW = I` (unwhitened fallback) — same behavior the
   function already has, surfaced to the user as a warning by the process.

2. **Process option + results node (`process_eigenmodes_inverse.m`).** Add
   `'harmonic'` → "Harmonic (unregularized)" to the method dropdown. When the method is
   `harmonic`, hide/ignore the SNR and prior controls (or leave them inert). For each input,
   save a **kernel-only results node** via the existing `'sources'` output path with:
   - `ImagingKernel = Φ · M̃` `[nVert × nCh]` (vertex-space, displayable);
   - `EigenKernel   = M̃` `[K × nCh]` (new field — the eigenmode-space operator for the viewer);
   - `Function = 'eigenmode_harmonic'`;
   - `GoodChannel`, `SurfaceFile`, `HeadModelFile`, `DataFile`, `nComponents=1`, `Whitener`,
     `Comment = 'Eigenmode HARMONIC (<n|auto> modes)'`, plus the standard `db_template('resultsmat')`
     fields. Registered with `db_add_data` like any results node.

   (When the existing methods `mne/dspm/sloreta` run with `OutputType=sources`, they continue to
   behave exactly as today; only the `harmonic` branch adds `EigenKernel`/`Function`.)

3. **Re-point the viewer (`view_eigenmodes_timeseries.m`).** Entry point becomes
   `view_eigenmodes_timeseries(ResultsFile)` where `ResultsFile` is a Harmonic node. It:
   - loads the node (`in_bst_results`/`bst_memory`) → `M̃` (`EigenKernel`), `GoodChannel`,
     `DataFile`, `SurfaceFile`; errors clearly if `Function ~= 'eigenmode_harmonic'` or
     `EigenKernel` is missing;
   - resolves the recordings dataset `iDS = bst_memory('LoadDataFile', DataFile)` and reads the
     current window `D = GetRecordingsValues(iDS, GoodChannel, 'UserTimeWindow', 0)` (the lazy
     path already built — raw + imported);
   - `Theta = M̃ * D`; caches `M̃` (as `Kernel`), `GoodChannel`, `SurfaceFile`, `DataFile`,
     `Component`/`CompRank` (first `K` from `in_tess_eigenmodes(SurfaceFile)`), window `Theta`/
     `TimeVector`/`WindowTime`, `Band`.
   - Everything downstream — band→L/R trace selection, butterfly/column, `FireModesChanged`
     band tracking, `FireCurrentTimeChanged`/`ReloadFigures` page tracking — is **unchanged**;
     only the source of `M̃` and `Theta` differs (read from node + whitened, instead of computed
     from the head model unwhitened).

   The `SyncWindow` re-read recomputes `Theta = cache.Kernel * D(newWindow)` with the cached `M̃`.

4. **Launch from the DB tree (`toolbox/tree/tree_callbacks.m`).** Add an "Eigenmode time series"
   item to the context menu of a **results node** whose `Function == 'eigenmode_harmonic'`,
   calling `view_eigenmodes_timeseries(ResultsFile)`. **Remove** the 3D-figure-popup launch added
   in `figure_3d.m` (it computed the unwhitened transform from any source figure — superseded).

## Data flow

```
process_eigenmodes_inverse (Method='harmonic', OutputType='sources')
   -> bst_inverse_eigenmodes('harmonic'): M̃ = pinv(iW·L·Φ)·iW
   -> save results node: ImagingKernel = Φ·M̃, EigenKernel = M̃, Function='eigenmode_harmonic'
        |
        |  user displays node on cortex  -> u(t) = Φ·M̃·D(t)
        |  user right-clicks node in tree -> view_eigenmodes_timeseries(ResultsFile)
        v
   read M̃ + GoodChannel + DataFile -> θ(t) = M̃·D(window)   (Φ·θ = u, exact)
   -> band -> L/R traces -> view_timeseries_matrix (butterfly/column, live tracking)
```

## Edge cases & errors

- **No noise covariance** — `iW = I`; the node is the *unwhitened* harmonic reconstruction. The
  process warns ("Harmonic without whitening; compute a noise covariance"). Not a hard error.
- **Node is not Harmonic** (no `EigenKernel`/wrong `Function`) — viewer shows
  `bst_error('Open this from an "Eigenmode HARMONIC" results node.')`.
- **Raw recordings** — handled by the existing lazy current-window read; `GetRecordingsValues`
  loads the current page on demand.
- **Channel set** — `M̃` columns and `GoodChannel` are stored together; `D` is read for exactly
  those channels, so ordering is consistent (`M̃` was built on `Gain(GoodChannel,:)` whitened).
- **Vertex / mode mismatch** — `bst_inverse_eigenmodes` already guards head-model vs eigenmode
  vertex counts and clamps `K ≤ min(nCh, nModes)`.
- **Backward-compat** — pre-existing `mne/dspm/sloreta` eigenmode nodes lack `EigenKernel`; the
  tree menu only appears for `Function='eigenmode_harmonic'`, so they are never opened by the
  viewer.

## Testing

- **`test_inverse_eigenmodes_harmonic_pure.m`** *(new)* — synthetic `Gain`, `Φ`, `iW`: assert
  `harmonic` returns `M̃` with `Φ·M̃` equal to `Φ·pinv(iW·L·Φ)·iW`; assert rank-safety (no blow-up
  when `Gain·Φ` is rank-deficient); assert `SNR`/`PriorAlpha` do not change the result.
- **Reuse** `test_eigenmodes_transform_pure.m` (the rank-safe pinv is already covered).
- **`test_harmonic_timeseries_e2e.m`** *(new, guarded smoke)* — compute a Harmonic node on a
  protocol that has eigenmodes + head model + noise cov + recordings; open
  `view_eigenmodes_timeseries(ResultsFile)`; assert the cached `Theta` equals `M̃ · D(window)`
  read independently, and that `Φ·Theta` (band = all modes) matches the node's
  `ImagingKernel · D(window)` to tolerance. Skips cleanly when no suitable data.

## Out of scope (YAGNI)

- Applying the panel band/window as a regularizing filter to the Harmonic reconstruction (the
  "filter bank" use) — the lever's source-map filtering already exists; this spec only creates
  the raw node + consistent time series.
- Adding `harmonic` to the main "Compute sources" GUI (`panel_inverse_2018`) — rejected during
  brainstorming in favor of extending the existing eigenmode inverse process.
- A free/loose-orientation harmonic inverse — fixed (constrained) orientation only, matching the
  eigenmode basis on the cortex surface.

## File-level summary

| Unit | File | Change |
|------|------|--------|
| `harmonic` method | `toolbox/inverse/bst_inverse_eigenmodes.m` | add `Method='harmonic'` branch |
| Process option + node | `toolbox/process/functions/process_eigenmodes_inverse.m` | add method; save `EigenKernel`/`Function` |
| Viewer reads node | `toolbox/gui/view_eigenmodes_timeseries.m` | entry `(ResultsFile)`; read `M̃`; drop self-transform |
| Tree launch | `toolbox/tree/tree_callbacks.m` | "Eigenmode time series" on Harmonic results node |
| Remove old launch | `toolbox/gui/figure_3d.m` | remove the popup item |
| Pure + e2e tests | `dev/tests/test_inverse_eigenmodes_harmonic_pure.m`, `dev/tests/test_harmonic_timeseries_e2e.m` | **new** |
