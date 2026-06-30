# Dynamics atom "Set operator" + recording-coupled axes — design

**Date:** 2026-06-30
**Status:** design (approved direction; pending spec review → plan)
**Depends on:** the Dynamics filterbank step-1 (atoms as generators, `BstClusterList`, `i_atom_ensure_axes`/`i_atom_realise`/`i_atom_preview`) — committed on `development`.

---

## 1. Motivation

In the Dynamics panel an atom's eigenbasis is currently hardcoded (`Variant='Laplace-Beltrami'`, `nModes=60`, synthetic `1 s @ 100 Hz`). That's wrong: a Dynamics session is launched from an **inverse kernel**, so the atom should ride that operator's eigenbasis (the **Dirac** eigen file when launched from a Dirac inverse), and its time grid should follow the **recording**, not a synthetic axis. The user also wants to **choose** the operator per atom, the way scouts choose a Function.

## 2. Operator variants

| Menu label | `bst_eigen` Variant | basis type |
|---|---|---|
| Geometric | `Laplace-Beltrami` | scalar |
| Connectomic | `LB-Connectome` | scalar |
| Tangent (connection Laplacian) | `Connection Laplacian` | complex tangent vector |
| *(launch default — not a menu pick)* | `Dirac` (if Dirac-launched) else `Laplace-Beltrami` | 3D embedded vector / scalar |

## 3. Per-atom operator

Each atom carries `G.Operator` (a `bst_eigen` Variant string) on its generator, alongside `KernelName/KernelParams/vertices`. Default at create time = the **launch-derived** variant:
- detect from the source result behind the figure (`getappdata(hFig,'DynamicsOverlay').srcResult`): if it is a **Dirac** inverse result (its comment/method indicates `dirac`), default `'Dirac'`; otherwise `'Laplace-Beltrami'` (LBO, the common case).
- `i_default_atom` sets `G.Operator`; `i_atom_detail` shows it.

"Set operator" changes the **selected atom's** `G.Operator` (mirrors scout "Set function"), then rebuilds + re-previews that atom.

## 4. "Set operator" menu (mirrors scout `CreateMenuFunction`)

A **"Set operator"** submenu under the **Atoms** menu, built dynamically when opened (like `panel_scout`'s `CreateMenuFunction`): `RadioMenuItem`s for **Geometric / Connectomic / Tangent**, with the selected atom's current `G.Operator` checked (a 4th `Dirac` item shown checked-but-also-selectable only when the atom is already on Dirac, so a Dirac atom isn't silently unlabeled). Callback → `OnSetOperator(variant)`.

## 5. Recording-coupled axes (rework `i_atom_ensure_axes`)

Replace the single cached `st.atomAx` with a **per-variant cache** `st.atomAxMap` (a `containers.Map` Variant→ax) + `st.atomBoundsMap`. `i_atom_axes(st, variant)`:
- find-or-build `ax = bst_eigen('Axes', struct('SurfaceFile',surf, 'Variant',variant, 'nModes',60, 'TimeWindow',[0 (nF-1)/Fs], 'SampleRate',Fs))`, where **`Fs`** is the recording's sample rate (from `srcResult`'s Time vector via `bst_memory('GetTimeVector', srcDS, srcResult)`; fall back to 100 Hz) and **`nF = round(4*Fs)`** — a **4-second** window (enough cycles for low-frequency oscillations, for visualization; a batch wavelet transform over the full record comes later).
- cache by variant; reuse on repeat.

`i_atom_realise`/`i_atom_preview` use the selected atom's `G.Operator` → `i_atom_axes(st, G.Operator)`.

## 6. Preview for vector / complex operators

`Laplace-Beltrami`/`LB-Connectome` are scalar → painted as today. `Connection Laplacian` (complex tangent) and `Dirac` (3D vector) realise a **non-scalar** field; the source overlay is scalar, so the preview paints the field **magnitude** `|W|` (one-signed → sequential colormap). `i_atom_realise` returns a scalar (magnitude for vector/complex variants), so `SetAtomField` is unchanged.

## 7. Scope

**In:** the per-atom `G.Operator` (+ launch-derived default), the "Set operator" Atoms-submenu (Geometric/Connectomic/Tangent), the recording-coupled per-variant axes (4 s @ recording Fs), and magnitude preview for vector/complex operators.
**Out / later:** the batch wavelet transform over the full recording; sensor-timeseries projection; the Dirac item being a first-class menu pick.

## 8. Risks / notes

- The atom realiser (`bst_eigenfilter('Atom')`) was built on scalar eigenbases. Realising on `Dirac`/`Connection Laplacian` (vector/complex `Phi`) must be **verified during implementation**; if a variant doesn't realise into a paintable field, guard it (warn + skip the preview) rather than crash. This is the main implementation risk.
- A variant's eigen file may not exist for the surface; `bst_eigen('Axes')` find-or-creates it (can be slow on first use — wrap in `bst_progress`). The Dirac file exists from the launch.
- 4 s at a high `Fs` (e.g. 1200 Hz → 4800 frames) is heavier than the old 100-frame axis; acceptable for one preview, watch responsiveness.

## 9. Testing

- Headless: `i_default_atom` sets `G.Operator`; the launch-derived default picks `Dirac` vs `Laplace-Beltrami` from a result comment; `i_atom_axes` caches per variant and uses the given `Fs`/4 s window; `OnSetOperator` updates the selected atom's `Operator`.
- Live (controller + user): in a Dirac-launched session, a created atom defaults to `Dirac`; Atoms → Set operator → Geometric/Connectomic/Tangent switches the basis and re-previews (magnitude for Tangent); the readout shows the operator.
