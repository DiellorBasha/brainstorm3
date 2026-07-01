# Dynamics panel — Preview completion (sub-project D) — design

**Date:** 2026-07-01
**Status:** design (approved; to be implemented in a fresh session)
**Depends on:** A (cleanup), B (frame core), C (analyze) — all shipped + pushed to `origin/development`.
**Part of:** the Dynamics-portal program (A→B→C→D). See [[dynamics-portal-optimization]].
**Supersedes:** the resolved Task-3 sketch in `dev/2026-06-30-dynamics-design-preview-modes-design.md` §5-6.

---

## 1. Motivation

The Apply/Preview loop is incomplete for the **Dirac** operator: `i_atom_apply` currently guards Dirac with
"scalar-only for now". D completes it. Lifting the guard gives the filtered eigenmode coefficients
`c_filt = diag(g(λ)) · K_eig · D`; from those, both a **cortex** filtered-source map and — the headline —
a **filtered-sensor** view, obtained by forwarding the filtered coefficients back through the head model
to the channels (`D_filt = L_eig · c_filt`). This is the "does the atom's spatial-scale filter change the
data the way I expect?" validation on the sensors themselves. Plus the operator/source **compatibility
gate** so the Set-operator menu reflects what the linked source supports.

**Resolved math (from the preview-modes doc §6):** `D_filt = L_eig · diag(g(λ)) · K_eig · D` — eigenbasis
inverse (`K_eig`) · atom filter (`g(λ)`) · eigenbasis forward (`L_eig`) · sensors (`D`). Only the Dirac
dSPM path can filter the 3-component leadfield vectors and forward to sensors; scalar operators filter a
magnitude and have **no** sensor view.

## 2. Decisions (from brainstorming)

- **Filtered-sensor display = overlay filtered vs raw on the recording time series** (a custom overlay on
  `figure_timeseries`, analogous to the cortex `DynamicsOverlay`), for direct raw-vs-filtered comparison.
- **Compatibility gate = grey out incompatible operators by `nComponents` + keep the Apply guard** as a
  backstop.
- **Dirac Apply = both cortex map + sensor overlay** from the shared `c_filt` (one compute, two views).

## 3. The shared filtered-coefficient pipeline (`i_atom_apply` Dirac branch)

When the atom's operator is `Dirac`, a Dirac-dSPM source is linked, and Apply is ON — over the 4 s cursor
window `iWin`:
1. `J = double(bst_memory('GetResultsValues', D.srcDS, D.srcResult, [], iWin, 0))` → Dirac source
   **3-vector** field `[3nV × nWin]` (keep the vector; do NOT reduce to magnitude).
2. `c = manifold_ft(Phi_d, B_d, J_embedded)` → project onto the Dirac quaternion eigenbasis → `[nEig × nWin]`
   coefficients. (Reuse `bst_eigenfilter`'s Dirac `RowMap`, which already maps the `3nV`→`4nEig` quaternion
   layout, rather than hand-rolling the embedding.)
3. `c_filt = g(λ) .* c` where `g = bst_eigfilter_kernel(kernel, kp)` on the Dirac `Lambda`.

From this single `c_filt`:
- **Cortex:** `J_filt = manifold_ift(Phi_d, c_filt)` → `i_paintable_scalar` per-vertex magnitude →
  `view_dynamics('SetFilteredField', …)` (the existing paint path).
- **Sensor:** `D_filt = L_eig · c_filt` `[nCh × nWin]` → `view_dynamics('SetFilteredSensors', …)` overlay.

**Reuse note:** steps 2-3 + the cortex reconstruction are exactly what `bst_eigenfilter('Analysis', J,
EigenMat_d, OperatorMat_d, kernel, kp)` computes for the Dirac variant. To get BOTH the reconstructed
cortex field AND the intermediate `c_filt` (needed for the sensor forward), either (a) call `Analysis`
for the cortex field and recompute `c_filt` once for the sensor, or (b) add a thin `bst_eigenfilter`
path that returns `c_filt` alongside the field. The plan picks one; the math is identical either way.

## 4. `L_eig` sourcing & cache

`L_eig` (Dirac-eigenbasis leadfield, `[nCh × nEig]`) = the Gain of `bst_dirac(HeadModel, 'nModes', K,
'Tau', tau)` (TRANSFORM verb — composed Dirac eigenmode head model) on the session's head model. Its
`Tau`/`nModes` **must match** the atom's Dirac eigenbasis from `i_atom_axes(st,'Dirac')`, so the eigenmode
ordering aligns with `c_filt` — assert/guard this. Built once per `(HeadModelFile, nModes, tau)` and
cached on `getappdata(0,'DynamicsDiracFwd')` (mirrors B's `DynamicsApplyCache`); invalidated on
operator/session change. The head model file comes from the source result's study
(`bst_get('HeadModelForStudy', …)` or the result's `HeadModelFile`).

## 5. Filtered-sensor overlay on `figure_timeseries`

- `view_dynamics('SetFilteredSensors', hRec, Dfilt, iChan, iWin)` stashes `Dfilt` + the window on the
  recording figure's appdata (`FilteredSensorsOverlay`) and triggers a redraw.
- A minimal hook in `figure_timeseries` draws the overlay: for the butterfly/time-series axes, add a
  colored second trace set (`Dfilt` mapped to the same channels/time as the raw), distinct color, so raw
  vs filtered read together. The global time cursor scrubs within the window (map cursor→window column).
- Clearing: `view_dynamics('ClearFilteredSensors', hRec)` removes the overlay traces; called when Apply
  turns off, the operator leaves Dirac, or the session closes.
- Find the recording figure via the existing pattern (`bst_figures('GetFiguresByType','DataTimeSeries')`
  matched to the session's `DataFile`); if none is open, open it (`view_timeseries`) so the overlay has a
  host — the Design launch keeps a clean cortex, but Dirac Preview needs the recording figure.

## 6. Operator/source compatibility gate

- On `SetTarget` (and when the linked source changes), read `R.nComponents` (`in_bst_results(srcFile, 0,
  'nComponents')`; 1 = scalar, 3 = unconstrained vector).
- Grey out (`jOpItems(k).setEnabled(false)`) the Set-operator radio items that don't fit:
  - scalar source (`nComponents==1`) → disable **Dirac** and **Tangent (Connection Laplacian)**;
  - vector source (`nComponents==3`) → all four enabled.
- Never auto-select a disabled operator (if the launch default would be disabled, fall back to the first
  enabled one). Keep the Apply-time "scalar-only for now" message as a backstop.
- `i_gate_operators(st)` centralizes this (reads `nComponents`, sets `jOpItems` enabled state).

## 7. Components & files

- **`toolbox/gui/panel_bst_dynamics.m`:**
  - `i_atom_apply`: lift the Dirac guard; add the Dirac branch (steps §3) → cortex `SetFilteredField` +
    sensor `SetFilteredSensors`.
  - `i_dirac_forward(st, D, iWin)` → `[Dfilt, cortexField]` or `c_filt` (computes/caches `L_eig`, projects,
    filters, forwards).
  - `i_gate_operators(st)` → the compat gate; call from `SetTarget`.
- **`toolbox/gui/view_dynamics.m`:** `SetFilteredSensors` / `ClearFilteredSensors` hooks (stash + trigger
  `figure_timeseries` redraw; find/open the recording figure).
- **`toolbox/gui/figure_timeseries.m`:** a small, guarded overlay-draw hook (draw `Dfilt` traces / clear),
  kept minimal and behind an appdata check so non-Dynamics figures are unaffected.
- **Reuse:** `bst_dirac` (TRANSFORM → `L_eig`), `bst_eigenfilter` Dirac `RowMap`/`Analysis`,
  `i_atom_axes(st,'Dirac')`, `i_cursor_window`, `i_paintable_scalar`, `SetFilteredField`, `bst_memory`.

## 8. Scope & edge cases

- **In:** Dirac Apply → cortex filtered magnitude + filtered-sensor overlay; the compatibility gate.
- **Out:** whole-frame/multi-atom sensor forward (single selected atom only, per the Apply model); sensor
  views for scalar operators (they filter magnitude — no oriented sensor view; cortex-only, as shipped in
  B/C); a saved filtered-recording file (the overlay is transient).
- Requires a **Dirac-dSPM source** (`results_DiracEig_KERNEL_*`) — carries `K_eig`. If the linked source
  isn't Dirac-dSPM, the sensor view is unavailable (info message); the cortex Dirac filter still previews.
- `L_eig` and `c_filt` eigenmode ordering must match (same `Tau`/`nModes`) — assert; on mismatch, bail the
  sensor view with a message (cortex still works).
- Dynamic (ts/js) kernels on the Dirac operator: the sensor forward uses the static `g(λ)`; a dynamic
  kernel's temporal response is out of scope for the sensor overlay (guard to static, or apply `g(λ,·)`'s
  spatial part — the plan picks; default: static gain for the sensor forward).

## 9. Testing

- **Headless (controller-run):** `bst_dirac(HeadModel)` returns `L_eig [nCh × nEig]`; **forward
  consistency** — for a random per-vertex 3-vector field `J`, `L_eig · (Phi_d' B_d J_embedded)` reproduces
  `Leadfield · J` on the represented modes (to a tolerance set by the mode truncation); `i_gate_operators`
  disables the right items for `nComponents ∈ {1,3}`.
- **Live (controller + user, MCP):** Dirac session (sub-MTL0002 `results_DiracEig_KERNEL_*`) → Apply →
  the cortex shows the filtered Dirac magnitude AND the recording time series shows the filtered-sensor
  overlay (raw vs filtered); editing a kernel param re-filters both; Apply OFF clears both; loading a
  scalar source greys out Dirac/Tangent. Screenshots.

## 10. Risks / notes

- **§5 (the `figure_timeseries` overlay) is the riskiest** — everything else reuses existing Dirac/
  eigenfilter machinery. Keep the hook minimal and behind an appdata guard so ordinary recording figures
  are untouched; verify the overlay tracks the time cursor and the Smooth/montage changes don't orphan it
  (cf. the [[overlay-on-displayed-surface]] gotcha for the cortex).
- `L_eig` build (`bst_dirac`) can be heavy — cache per session; reuse the atom's already-built Dirac
  eigenbasis (`i_atom_axes`) rather than recomputing the eigendecomposition.
- Keep the scalar Apply path (B/C) working throughout; the Dirac branch is additive behind the operator +
  Dirac-dSPM guards.
