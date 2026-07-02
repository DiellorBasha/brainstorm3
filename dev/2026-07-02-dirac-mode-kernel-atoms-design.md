# Dirac mode-kernel atoms — sourcing coefficients from the inverse — design

**Date:** 2026-07-02
**Status:** design (approved sketch + 2 decisions; spec for review)
**Part of:** the Dynamics-portal program. Builds on dimensional atoms (HEAD 371c102b) + the D Dirac
sensor forward + the operator-applicability work.

---

## 1. Motivation

For a **Dirac-dSPM** source the eigen-coefficients are one kernel multiply from the sensors —
`c = ImagingKernelMode · recordings` — which is (a) what the inverse already produces, (b) ~150× cheaper
than projecting the reconstructed field (`nChan≈274` vs `4nV≈40,968`), and (c) expressed in the inverse's
**own** eigenbasis at its **full** mode count. An atom is then just `c_filt = g(λ)·c` — a diagonal
reweighting of the inverse's modes — and reconstruction is needed only for the views.

Today the panel does the opposite for a real source: `GetResultsValues` reconstructs the field, then
`manifold_ft` **re-projects** it onto a fresh **60-mode** `bst_eigen('Axes')` — a round-trip that both
pays for a projection it could get free and truncates the analysis to 60 modes (the interactive scalogram
/ localize / Apply are all capped at 60, while the batch `process_source_frame` runs at 200). This change
sources `c` from the mode kernel instead.

**Proven machinery to reuse:** `view_eigen_timeseries.m` and `view_eigenmode_spectrum.m` already compute
`c(t) = ImagingKernelMode · M(GoodChannel,t)` and load the basis via `in_bst_eigen(DiracEigenFile)` +
`in_bst_operator(OperatorFile)` (Φ `{1×2}[4Vh×K]`, gv, Mass), with the exact pure-imaginary quaternion
embed/stack (`psi(2:4:end)=Jx …`, stacked L-then-R per `ModeHemisphere`). We reuse that pattern.

## 2. Decisions (approved)

- **Unify on the inverse basis.** For a Dirac-dSPM session, BOTH Design (impulse probe) and Apply use the
  inverse's `DiracEigenFile` basis (its modes), so the previewed point-spread has the same resolution as
  the filter actually applied. One basis per session.
- **All Dirac real-source operations** switch to the mode-kernel coefficients: cortex Apply, sensor
  forward, Analyze-window scalogram, Localize-bands (they all consume one shared `c`).
- **Scope-out (unchanged, fall back to today's path):** scalar-magnitude atoms (`|J|`, LB / LB-Connectome
  — nonlinear in `c`, cannot use the linear mode kernel); non-Dirac-dSPM sources (no `ImagingKernelMode`).

## 3. Architecture

### 3.1 Session Dirac basis (the atom's operative axes)

When the linked source is a Dirac-dSPM result (`ImagingKernelMode` non-empty AND `DiracEigenFile` set) and
the operator is `Dirac`, `i_atom_axes(st,'Dirac')` returns the **inverse's** basis instead of a fresh
`bst_eigen('Axes')`:
- `E = in_bst_eigen(R.DiracEigenFile)` → `ax.Phi = E.Phi {1×2}`, `ax.GlobalVertices = E.GlobalVertices`.
- `O = in_bst_operator(E.OperatorFile)` → `ax.Mass = O.Mass {1×2}`, `ax.Operator = O` (carries Registry/
  Frame/Tau for the dimensional-atoms `Fiber`/decode + the sensor-forward Tau).
- `ax.Lambda` per hemi from the stored `R.Eigenvalues` split by `R.ModeHemisphere` (so `g(λ)` uses the
  inverse's eigenvalues, aligned with the mode-kernel rows).
- `ax.SurfaceFile = R.SurfaceFile`, plus `ax.Time`/`nT` as today. Cached per `(DiracEigenFile)`.
- Guard: `numel(Eigenvalues) == ΣK_h` and `size(Phi{h},1)==4·numel(gv{h})`; on mismatch, fall back to the
  canonical `bst_eigen('Axes')` path with a one-line info note (never a silent wrong basis).

Non-Dirac operators, and Dirac with a non-dSPM source, keep the canonical `bst_eigen('Axes')` path.

### 3.2 Mode coefficients (the free projection)

`c = i_mode_coeffs(st, D, iWin)` → `{1×2}` cell of `[K_h × nWin]` (or one stacked `[nMode×nWin]` split by
`ModeHemisphere`):
- `R = in_bst_results(src, 0, 'ImagingKernelMode','GoodChannel','DataFile','ModeHemisphere')`.
- Load the recordings window `Fdat[GoodChannel, iWin]` (same loader `view_eigen_timeseries` uses).
- `cAll = double(R.ImagingKernelMode) * Fdat` → `[nMode × nWin]`; split into `c{h}` by `ModeHemisphere`.
- Cached on `getappdata(0,'DynamicsModeCoeffCache')`, keyed `DiracEigenFile|DataFile|iWin` (mirrors B's
  `DynamicsApplyCache`). Invalidated on source/operator/session change.

### 3.3 Apply (Dirac branch), reconstruction & the four views

`i_atom_apply` Dirac branch, over the cursor window `iWin`:
1. `ax = i_atom_axes(st,'Dirac')` (§3.1, the inverse basis).
2. `c = i_mode_coeffs(st, D, iWin)` (§3.2, free).
3. `c_filt{h} = g(λ_h) .* c{h}` with `g = bst_eigfilter_kernel(kernel, kp)` on `ax.Lambda`.
4. Views from the single `c_filt`:
   - **cortex:** `Uf_h = Φ_h · c_filt{h}` → extract imag 3-vector (rows `2:4:end/3:4:end/4:4:end`) →
     full-surface `V3` + per-vertex magnitude → `SetFilteredField` + quiver `V3` (reuse dimensional-atoms
     decode + Task-6 display).
   - **sensor:** `Dfilt = L_eig · c_filt_stacked`; `L_eig` from `i_dirac_leadfield` built to the **inverse
     basis** (Tau/nModes from `ax`), asserting `CompHM.Eigenvalues == [Eigenvalues]` (reuse the D guard);
     on mismatch skip the sensor view, cortex still previews.
   - **scalogram (Analyze):** `bst_eigenwavelet('Scalogram', ax, gCell, C)` with `C = c` (the mode-kernel
     coefficients) — energy `‖g_m(λ)·c‖²` per band **in coefficient space, no reconstruction**.
   - **localize (JTVAtoms):** on the reconstructed bands `scal.W` (reuse).
5. Static-kernel domain guard for the sensor forward (ts/js bail) — reuse D.

`i_apply_projection` / `OnAnalyzeWindow` / `OnLocalizeBands`: for a Dirac-dSPM session, take `C` from
`i_mode_coeffs` instead of `manifold_ft`-projecting the reconstructed magnitude; else current path.

### 3.4 Design impulse (unified basis)

Because `i_atom_axes(st,'Dirac')` now returns the inverse basis, the existing dimensional-atoms
`bst_eigenfilter('Atom', ax, kernel, kp, seed, seedDir)` automatically realises the impulse in that basis
(delta projection stays cheap — sparse; reconstruction at the inverse's mode count). No separate change;
the impulse preview now matches the Apply resolution.

## 4. Correctness anchor (drives the tests)

The inverse's vertex kernel is the reconstructed mode kernel: `ImagingKernel[3nV×nCh] = Φ_imag · Imaging
KernelMode`. Therefore **reconstructing from the mode-kernel `c` reproduces `GetResultsValues` exactly**
(to numerical precision) when using the DiracEigenFile's full mode set:

```
reconstruct(Φ, ImagingKernelMode · d)   ==   GetResultsValues(src, iWin)   (imag 3-vector)
```

This is the lossless-projection test — it proves the free path drops nothing and that `Φ`, `ModeHemisphere`,
and the mode kernel are mutually aligned. A truncated basis would differ by exactly the dropped modes.

## 5. Components & files

- `toolbox/gui/panel_bst_dynamics.m`:
  - `i_atom_axes` — Dirac-dSPM session-basis branch (§3.1) + guard/fallback.
  - `i_mode_coeffs(st, D, iWin)` — new (§3.2), cached.
  - `i_is_dirac_dspm(D)` — helper: source has `ImagingKernelMode` + `DiracEigenFile`.
  - `i_atom_apply` Dirac branch — consume `c`, filter, reconstruct for the 4 views (§3.3).
  - `i_apply_projection` — Dirac-dSPM branch returns `c` (per-hemi) instead of re-projecting.
  - `i_dirac_leadfield` — build `L_eig` to the inverse basis (Tau/nModes from `ax`), assert vs `Eigenvalues`.
- **Reuse:** `in_bst_eigen`/`in_bst_operator` (basis), the `view_eigen_timeseries` recordings loader,
  `ImagingKernelMode`/`Eigenvalues`/`ModeHemisphere`/`GoodChannel`/`DiracEigenFile` on the result,
  `bst_eigenwavelet('Scalogram'/'JTVAtoms')`, `bst_dirac` (L_eig), the dimensional-atoms decode + quiver
  display, `bst_eigenfilter('Atom'/'Fiber')`.

## 6. Scope & edge cases

- **In:** Dirac-dSPM Apply (cortex+sensor) + scalogram + localize + Design, all in the inverse's basis via
  mode-kernel coefficients.
- **Out:** scalar-magnitude atoms (nonlinear → current path); non-Dirac-dSPM sources (current path); a
  saved filtered-recording/coefficient file (transient, as today); re-deriving `ImagingKernelMode` for
  results that predate it (fall back).
- **Edge:** measure scaling (amplitude/dSPM/sLORETA) is already baked into `ImagingKernelMode` — `g(λ)`
  applies on top, unchanged. GoodChannel selects recording rows. Empty/absent `ImagingKernelMode` →
  fallback. Basis/eigenvalue-count mismatch → fallback + info. Large mode count (200–400): reconstruction
  per view scales linearly (acceptable; projection is now free and cached).

## 7. Testing

- **Headless (matlab -batch):**
  - **Lossless projection (anchor):** `reconstruct(Φ, ImagingKernelMode·d)` imag 3-vector == `GetResults
    Values` field over a window, to `~1e-8` relative (full mode set).
  - `i_mode_coeffs` shape `[K_h×nWin]`, split matches `ModeHemisphere` counts; cache hit returns identical.
  - Filter equivalence: `c_filt = g(λ)·c` then reconstruct == field filtered by the same `g` in the same
    basis (sanity vs a direct `manifold_ift(Φ, g·(Φ'Bψ))` on the reconstructed field), to tolerance.
  - `i_atom_axes` returns the inverse basis (K = inverse's nModes, `Lambda` == `Eigenvalues` split) for a
    Dirac-dSPM source; canonical otherwise; guard falls back on a forced mismatch.
- **Live (controller, MCP):** Dirac-dSPM session → Apply: cortex filtered magnitude + quivers, sensor
  overlay, Analyze scalogram (now spanning the inverse's full `√λ` range, not 60), Localize bands. Confirm
  the scalogram/atom now reach finer scales than the old 60-mode cap. Screenshot.

## 8. Risks / notes

- **Basis alignment is the load-bearing invariant** — `Φ` (from DiracEigenFile), `Eigenvalues`,
  `ModeHemisphere`, and `ImagingKernelMode` rows must share one ordering. The lossless-projection test is
  the guard; add the count/shape asserts in `i_atom_axes` + the D-style eigenvalue assert in
  `i_dirac_leadfield`.
- **Recordings loader** — reuse exactly what `view_eigen_timeseries` uses (raw/epoched handling,
  GoodChannel) rather than re-rolling; a subtle channel/row mismatch would corrupt `c`.
- **Live-session instability** (Apple-silicon `GlobalData` drops) — prefer `matlab -batch` for the anchor
  test; short self-contained live passes.
- Keep the scalar and non-dSPM paths byte-unchanged (guarded branch), so shipped B/C/D behavior is intact.
