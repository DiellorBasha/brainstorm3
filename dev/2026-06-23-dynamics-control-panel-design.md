# Dynamics computational control panel — design

**Goal:** Turn the Dynamics panel into an interactive *atom-coordinate selector*: a computational
control panel on top of the Atoms section that filters a field along the conjugate axes
(time↔frequency, region↔scale), driving a linked set of figures so each manipulation is visible.
The filtered field's features later become atoms.

## The atom, crystallized

An atom is one point pinned in four coordinates, set by two conjugate "windows" each:

| Axis | Window the control sets | Operator (how the field is shaped) |
|------|-------------------------|------------------------------------|
| Time | time window | (later) burst/window detector |
| Frequency | **band** δ/θ/α/β/γ | temporal bandpass (Filter-tab DSP) |
| Space | region | differential — `bst_operators` (grad/div/curl) |
| Scale | spatial **eigen-scale** | eigenfilter smoothing — `bst_eigenfilter` |

## Linked figures (the launch)

`view_dynamics('FromResult', ResultsFile)` (and the direct `view_dynamics(file)` path when the
table carries provenance) opens, time-linked by the global cursor:

- the **Dynamics/Atoms panel** (docked tools tab) — controls + atom table,
- **`view_timeseries`** on `T.DataFile` — where temporal/frequency filtering is seen,
- **`view_helmholtz`** 3D on the source result — where spatial filtering is seen; the atom
  **markers are drawn on this cortex** (no separate static surface).

Atom selection already drives `panel_time('SetCurrentTime')`, so both figures follow the cursor.

## Control panel (above the Atoms section)

```
┌ Frequency ───────────────────────────┐   Increment 1 (new)
│  [δ] [θ] [α] [β] [γ]                   │   mutually-exclusive; re-click active = off (raw)
└───────────────────────────────────────┘
┌ Space ────────────────────────────────┐   Increment 2 (fold view_helmholtz controls in)
│  Smooth [eigen cutoff]  Op |J|·Φ·Ψ·∇   │
└───────────────────────────────────────┘
┌ Atoms ────────────────────────────────┐   built
│  stack tree  │  per-window atom list    │
└───────────────────────────────────────┘
```

### Frequency section (Increment 1)
Five mutually-exclusive band toggles (δ 2–4, θ 4–8, α 8–13, β 13–30, γ 30–60 Hz). Selecting a band
calls `panel_filter('SetFilters', LPon=1, LP=hi, HPon=1, HP=lo, sinOff, [], mirror=0, FullSources=1)`,
which bandpasses the displayed recording **and** the source results. Verified load-bearing fact:
`view_helmholtz` reads its source via `bst_memory('GetResultsValues')` (view_helmholtz.m:105), which
returns the display-filtered source when `FullSourcesEnabled` is on — so one band button updates both
the time series and the Helmholtz 3D. Re-clicking the active band (or none) disables the filter.
The selected band is the atom's **frequency coordinate**.

### Space section (Increment 2 — DONE)
Bordered "Space" section between Frequency and Atoms, driving the linked Helmholtz figure via
`view_helmholtz('SetSmoothing'/'SetComponent', hFig, …)`:
- **Smooth** = eigenfilter low-pass (the **scale** axis): a checkbox + kernel combobox + a parameter
  slider built with the SAME `panel_eigenfilter_design` machinery the Helmholtz panel uses (the 'heat'
  kernel's `t` slider = *Scale*). Built at `SetTarget` from `St.Lambda` read off the Helmholtz figure.
- **Op** = differential operator: Φ (Potential / divergence) and Ψ (Stream / curl) toggles, mutually
  exclusive, all-off = Total |J| — `view_helmholtz('SetComponent', …, 'Irrot'/'Solen'/'Total')`.
- The standalone Helmholtz panel is hidden on launch (`gui_hide('Helmholtz')`); the figure stays live
  and the Dynamics panel's Space section is its sole controller. The chosen operator/scale are stashed
  on the atom target (`st.curOp`, `st.curScale`) as the atom's space/scale coordinates.
- Extending the operator palette to gradient/Laplacian (full `bst_operators`) is a later refinement —
  div/curl (Φ/Ψ) are the core differential quantities for a source vector field.

## Out of scope (later)
- Atom storage: writing the manipulated field's extrema at the chosen (band, scale, operator) as atoms.
- A shared ephys-band constant (δ/θ/α/β/γ duplicated across process_source_atoms / _helmholtz_events /
  _evt_refphase and now the panel) — candidate for a `bst_ephys_bands` helper.

## Increment 1 deliverables
1. `view_dynamics` linked-trio launch (Helmholtz + timeseries + panel, markers on the Helmholtz cortex),
   with graceful fallback to the static surface when the table lacks DataFile/ResultsFile provenance.
2. A bordered **Frequency** section in `panel_bst_dynamics` with the 5 band toggles driving
   `panel_filter('SetFilters', …, FullSources=1)`.
3. Regression: existing suite stays green; new check that a band toggle sets `GlobalData.VisualizationFilters`.
