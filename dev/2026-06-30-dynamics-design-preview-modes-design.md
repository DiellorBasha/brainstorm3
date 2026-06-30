# Dynamics panel — Design / Preview (Apply) modes — design

**Date:** 2026-06-30
**Status:** design (approved; to be implemented in a fresh session)
**Depends on:** the Dynamics filterbank + Set-operator work (atoms as generators with `G.Operator`, `i_atom_axes`/`i_atom_realise`/`i_atom_preview`, `SetAtomField` overlay) — all committed on `development`.

---

## 1. Motivation

The Dynamics panel currently launches a linked trio (docked panel + a `view_surface_data` figure on the **real** unconstrained source + the recording timeseries) and, on adding an atom, paints the atom's **impulse response** *over* the real source display. That conflates two distinct things. Split them into two explicit modes:

- **Design view (default):** a clean cortex, **no real source data loaded**; the displayed field is the atom's impulse response (the filter applied to a Kronecker delta at the seed, t = 0). This is for *shaping* the atom.
- **Preview view (Apply):** load the real data through the imaging kernel, **filter it through the atom**, and show the filtered source map + the filter's effect in the **sensor/recording** view. This is for *tuning* the atom against data.

Per earlier decisions: the preview operates on the **selected atom** over a **4 s window at the cursor**; the whole-filterbank/frame batch transform stays a separate later step.

## 2. Design view (clean-cortex launch)

`view_dynamics` opens a **clean `figure_3d` cortex** for the design session:
- Use `view_surface` (geometry only), not `view_surface_data` (real source) — no source data painted, no recording-timeseries figure opened.
- Keep the source-result reference on the overlay appdata (`DynamicsOverlay.srcResult`/`srcDS`): it carries the **`ImagingKernel`** and the **`DataFile`** (the recording) that Apply needs, plus the operator hint (`i_launch_operator`).
- The displayed cortex field = the selected atom's **impulse response** via the existing `SetAtomField` path (`bst_eigenfilter('Atom')`). Overlay mode `atom`.
- Time axis: the atom's intrinsic window (4 s @ recording `Fs`, from `i_atom_axes`). No real time cursor over data is required in Design mode.

## 3. Preview view (Apply)

An **Apply** toggle on the east toolbar. When ON (overlay mode `atom-filtered`):

1. **Reconstruct the source over the 4 s window at the cursor:** `F = ImagingKernel × Recording[window]` (window = `[t0, t0+4s]` around the global cursor, clamped to the recording; unconstrained → `F` is `[3·nV × nWin]`).
2. **Filter through the selected atom:** `Ffilt = bst_eigenfilter('Analysis', F, EigenMat, OperatorMat, G.KernelName, G.KernelParams)` on the atom's operator (`G.Operator`). `Analysis` projects `F` onto the modes, scales each by `h(λ)`, and reconstructs — spatial `g(λ)` and, for ts/js kernels, the temporal response over the wavelet's support. (`EigenMat`/`OperatorMat` come from the same eigen/operator files `i_atom_axes` already loads for the variant.)
3. **Paint the filtered source map** on the cortex (magnitude for vector/complex operators), as a function of the window's time (the global cursor scrubs within the window).
4. **Sensor timeseries:** relate `Ffilt` back to the **sensors** via the imaging-kernel relationship and show the **recording/sensor view** — the filter's effect on the data's own axes. *(The exact source→sensor map is the one open implementation detail — see §6.)*
5. Editing the atom's params (or operator/seed) **re-filters live**. Toggling Apply OFF returns to the clean Design impulse view (overlay mode back to `atom`).

## 4. Overlay modes

`view_dynamics`'s `DynamicsOverlay` gains a third `Op`: `none` (native — unused in Design), `atom` (impulse response, Design), `atom-filtered` (real data through the filter, Preview). `SetAtomField` stays the impulse path; add `SetFilteredField(hFig, Ffilt, gv, ...)` for the Preview path (same scatter-paint, different source).

## 5. Operator / source compatibility gate

Because Apply touches the real source, the operator must match the source's structure: read `R.nComponents` (1 = scalar, 3 = unconstrained vector) from the source result. **Disable** incompatible operators in the Set-operator menu — scalar source → only Geometric/Connectomic; vector source → all four (scalar operators act on the magnitude). Don't default to an unavailable/incompatible operator. (This guard is optional in Design mode but **required** in Apply.)

## 6. The source→sensor projection — RESOLVED (2026-06-30): Dirac-eigenbasis kernel

**Decision (user):** the **imaging kernel** route — and specifically the **Dirac dSPM** path. The imaging
kernel already encodes the head model, all the physics, and the priors, and it maps leadfield vectors to
sensor channels. The move is to transform the imaging kernel from `nSources × nChannels` into the **Dirac
eigenbasis** (`nEig × nChannels`) — *exactly what the Dirac dSPM does* — apply the atom's spectral gain
`g(λ)` to the eigenmodes, and forward back to channels:

> `D_filt = L_eig · diag(g(λ)) · K_eig · D`   (eigenbasis forward · atom filter · eigenbasis inverse · sensors)

**Critical constraint:** *only the Dirac dSPM can filter the 3-component leadfield vectors and forward to
the sensors.* Scalar operators (Laplace-Beltrami / LB-Connectome) filter a source **magnitude** — not an
oriented field — so they have **no sensor view**. The filtered-sensor timeseries is therefore a
**Dirac-operator-only** feature.

**Implications for Task 3 (now its own design):**
- The Dynamics session must be launched from the **Dirac dSPM kernel** (`results_DiracEig_KERNEL_*`), not a
  plain unconstrained dSPM — that kernel carries the eigenbasis transform.
- The current **Dirac Apply guard** (Task 2: "scalar-only for now") must be **lifted** for the Dirac path.
- Reuse `bst_dirac` (Transform/Reconstruct) for the eigenbasis leadfield transform (`L_eig`/`K_eig`); see
  the [[dirac-eigenmode-leadfield]] / [[bst-inverse-dirac]] machinery.
- The cortex Preview (Tasks 1+2) is independent and already shipped (scalar magnitude filtering); the sensor
  view is the Dirac-operator extension layered on top.

## 7. Scope

**In:** the clean Design launch; the Apply toggle filtering the *selected* atom over a 4 s window (`bst_eigenfilter('Analysis')`) → filtered source map + sensor timeseries; the operator/source gate; the `atom-filtered` overlay mode.
**Out (later):** whole-filterbank/frame application; full-record batch wavelet transform; multi-atom preview.

## 8. Testing

- Headless: `bst_eigenfilter('Analysis', F, EigenMat, OperatorMat, kernel, kp)` filters a synthetic `F` on each operator (scalar + vector) into the right shape; the window extraction (`Recording[window]` at the cursor) returns the expected sample range; the operator gate disables the right items for `nComponents`.
- Live (controller + user): launch → clean cortex (no source data, no timeseries); add atom → impulse response; **Apply** → the filtered real-source map appears + the sensor timeseries shows the filter's effect; editing params re-filters; Apply off → back to impulse; incompatible operators greyed out.

## 9. Risks / notes

- The source→sensor projection (§6) is the main unknown — resolve first.
- Reconstructing `F = K × Recording[window]` for an unconstrained kernel over 4 s × `Fs` can be large; window-limit keeps it bounded; reuse `bst_memory` source-loading helpers rather than multiplying `K` by hand if a helper exists.
- `bst_eigenfilter('Analysis')` was built for surface source fields — verify its `RowMap` handles the unconstrained (3·nV) layout for each operator; guard like the impulse path.
- Substantial `view_dynamics` change (clean-cortex launch) — keep the impulse/Design path working throughout; build Preview behind the toggle.
