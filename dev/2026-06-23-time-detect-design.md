# Time-window detection in the Dynamics panel — design

**Goal:** Make `process_evt_refphase` (the time-window detector) the panel's **entry point**. Selecting a
frequency band and hitting **Detect** runs refphase for that band and writes the *temporal skeleton* of
the atom table — the band-window stack + the per-cycle phase markers. This is the FIRST atom creation;
all subsequent analysis (Space/Record) hangs off it. **Decision: temporal only** — no source localization
in this click (Record adds space later).

## Time ↔ Frequency, one band selection
The Frequency section's δ/θ/α/β/γ already names the band (and drives the live band-pass). A **Detect**
button beside it uses the SAME band to run the window detector — one frequency choice for both the
display band-pass (Frequency) and the window detection (Time).

## Detect action (`panel_bst_dynamics` `OnDetect`)
1. Require a selected band (`st.curBand`/`curBandName`); else warn.
2. Load the recording sensors `F` + `TimeVector` for `T.DataFile` (MEG channels), data or raw.
3. `OPTIONS = process_evt_refphase('Compute')`; `OPTIONS.freqRange = st.curBand`.
   `[evt, markers] = process_evt_refphase('Compute', F, TimeVector, OPTIONS, validMask)`.
4. Build the **temporal** atom groups (replacing any prior groups for this band):
   - **band-window** extended group `"<band> (lo-hi Hz)"`: `times = evt` [2×N], `band/bandName`, color.
   - **4 phase-marker** simple children `"<band>_peak/_trough/_rising/_falling"`: `parent`=window,
     `phase`, `times = markers.<phase>` [1×M], `band/bandName`. **No** `vertices/pos/hemi` — time only.
5. `i_apply` (refresh tree; no cortex markers yet — these are temporal) + auto-save.

## Tree / right-list
Unchanged shape: the band-window stack expands to its window leaves; selecting a window lists the phase
markers within it. `i_window_atoms` is extended to handle **no-vertex** markers — iterate by occurrence
time, show `time · phase` (vertex column omitted when absent); selecting one still jumps the time cursor
(no cortex highlight until the atom is localized).

## Entry flow (Detect becomes the first step)
`AtomsFromResult` no longer auto-runs `process_source_atoms`. It reuses an existing `dynamics_*` table,
else creates an **empty** one carrying provenance (`DataFile`, `SurfaceFile`, `ResultsFile`). `view_dynamics`
opens the linked trio + panel on an empty table (empty tree, no markers). The user then: **pick band →
Detect** (Time) → step to a marker → shape the field (Space/Scale) → **Record** (space). Both creation
paths feed one table; `process_source_atoms` stays as the standalone batch process.

## Out of scope (later)
- Source localization on Detect (kept separate, in Record).
- Bad-segment mask / time-window restriction in the panel Detect (use refphase defaults for now).
- Phase tagging of recorded spatial atoms from these markers.

## Deliverables
1. `panel_bst_dynamics`: a **Detect** button in the Frequency section + `OnDetect` (load recording → refphase
   → band-window stack + 4 phase-marker children, temporal; replace-by-band; refresh + auto-save).
2. `i_window_atoms`: handle no-vertex temporal markers (`time · phase`).
3. `view_dynamics` + `AtomsFromResult`: open/empty-table path so Detect is the first creation step.
4. Regression: suite green; new check that Detect on the test recording creates the `<band> (lo-hi Hz)`
   stack + the 4 phase-marker children with times and no vertices.
