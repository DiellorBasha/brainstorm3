# Navigate / Detect / Save contract — design

**Goal:** Make the atom panel's persistence boundary explicit: **Navigate** writes only the transient
cursor, **Detect** stages a time skeleton as Brainstorm **events** (rendered + navigable, not yet atoms)
parameterized by a live frequency band, and **only Save writes the atom table** — committing atoms whose
**time and frequency are recorded as numeric tensor indices**.

**Status:** design, ready to plan. Phase 4 of the atom-tensor architecture
(`docs/superpowers/specs/2026-06-24-atom-tensor-architecture-analysis.md`, §4 / §6 step 4). Builds on
Phase 1 (`bst_atom`), Phase 2 (`bst_geodesic_tool`), Phase 3 (the navigator panel).

**Author:** Diellor Basha, 2026-06-24

---

## 1. Motivation

Today all three panel actions persist immediately: `OnDetect` writes the band-window + phase children
into `st.T` and auto-saves; `OnRecord`/`OnCaptureRegion` likewise commit + auto-save. There is no
separation between "looking at a detection" and "committing atoms." The atom-tensor model wants a clean
**Navigate → Detect → Save** contract (analysis §4): Navigate and Detect never touch the saved table; only
Save does.

The realization piggybacks on existing Brainstorm infrastructure, because for the **time** and
**frequency** axes an atom's coordinates are exactly what Brainstorm already renders:
- **Time** ⇄ Brainstorm **Events** (auto-rendered on the time-series view; navigable via the global cursor,
  which also drives figure_3D).
- **Frequency** ⇄ Brainstorm's **filter panel** band-pass (dynamically re-filters the time-series view) —
  already driven by the navigator's Frequency block.

## 2. Decisions (locked in brainstorming)

1. **Three buffers, one writer set.** `st.nav` (transient cursor; Navigate) · **detection events** (the
   Detect staging — a real Brainstorm event group on the recording, *not yet atoms*) · `st.T` (the saved
   atom table; the only persisted buffer). Writers that commit to `st.T`: **Record** (extrema), **Capture**
   (region), **Save cursor** (one atom from `st.nav`), **Save detection** (events → atoms). Detect writes
   only events; Navigate only `st.nav`.
2. **Detect stages as Brainstorm events** (not a custom overlay). The refphase result becomes a recording
   event group → auto-rendered on the time series, listed in the Record-panel events GUI, mirrored in the
   Atoms tree, navigable for free.
3. **Both ephemeral axes piggyback symmetrically.** Time = events; frequency = the filter-panel band-pass.
   The detection state is a `(time, frequency)` pair; the user can try different frequency windows before
   saving (change the Frequency block → re-Detect → re-filter + regenerate events).
4. **Frequency is recorded on the atom as a numeric tensor index** `(center, extent)` via `bst_atom`, not
   merely a `bandName` label. Save stamps it from the navigator's current band.
5. **Save cursor as one atom** — commit the literal 4-D cursor (`st.nav`) + the measured descriptor as one
   atom (distinct from Record's extrema search).
6. **Load atom into navigator** — an explicit action populates the four blocks + `st.nav` from a saved atom
   (round-trips Navigate ↔ saved). Not auto-on-select (won't clobber the cursor).
7. **Source detection stays manual this phase.** Building an atom's source = navigate to a detection event
   → Region tool / Capture / Record. The auto source-feature detector (Φ/Ψ extrema + saddles via
   `bst_operators`/`bst_eigen`) is a **separate, dedicated future phase** (plan + develop + validate).

## 3. Architecture

### 3.1 Detect → detection events (no st.T write)

`OnDetect` (rewritten):
1. Require a navigator frequency band (`st.nav` freq / `st.curBand`), which is already applied as the live
   band-pass via the Frequency block (filter panel).
2. Load MEG + run `process_evt_refphase('Compute', F, TimeVector, OPTIONS)` for that band (as today).
3. **Write the result as a Brainstorm event group on the recording** (reusing the events API, e.g. the same
   event structure `process_evt_refphase('Run')` produces): an extended event `<band> (lo-hi Hz)` +
   simple events `<band>_peak/_trough/_rising/_falling` with the marker times. Refresh the time-series so
   the events render. Re-running for a band replaces that detection group.
4. **No `AddGroup(st.T,…)`, no auto-save.** The events are the staging.

The Atoms tree mirrors the detection event group under a distinct `Detection (events)` node (read from the
recording's events), so the user sees it in the panel too; selecting a marker jumps the time cursor.

### 3.2 Save detection → atoms (the writer)

`OnSaveDetection` (new button in a Save area):
1. Read the detection event group from the recording (the `<band> (lo-hi Hz)` extended event + the 4 phase
   simple events).
2. Build atom groups in `st.T` exactly as today's `OnDetect` did — the band-window extended group + 4
   phase children (times from the events) — but **stamp each group's frequency Localization numerically**
   from **the band the detection was run at** (recorded with the detection event group — equal to
   `st.nav` freq at detect time, but read from the event group so it can't drift if the user changed the
   navigator band without re-detecting): `G = bst_atom('Set', G, 'freq', 1, freqLoc)` where
   `freqLoc.center/extent` derive from that band (`band=[center−extent, center+extent]`), plus `bandName`.
3. `i_remove_band(st.T, bandName)` first (replace any prior saved windows for that band); `i_apply`;
   auto-save `st.T`.

`OnClearDetection` removes the detection event group (discard without committing).

### 3.3 Save cursor → one atom (the writer)

`OnSaveCursor` (new): snapshot the current cursor as one atom.
- Read `st.nav` (the 4-D coords already written by the navigator blocks via `bst_atom`).
- Measure the descriptor at the cursor: if source is localized and a Helmholtz field is present, sample the
  active operator's scalar (Φ/Ψ/|J|) at the seed → `strength` (+ `charge` sign); else descriptors empty.
- Find-or-create the `(bandName, Function)` group (as `OnRecord` does), append **one** occurrence carrying
  the cursor's time / freq / source(seed,region,radius) / scale coords (copied from `st.nav` via the group
  fields) + the descriptor. `i_apply`; auto-save.

### 3.4 Load atom → navigator (round-trip)

`OnLoadAtom` (new; explicit — a "Load into navigator" Atoms-menu item / double-click on a saved
occurrence): for the selected saved atom occurrence `(g, o)`, for each axis
`loc = bst_atom('Get', st.T.Groups(g), axis, o)`; write `loc` into `st.nav` (`bst_atom('Set', st.nav, …)`)
and populate the four block fields from `loc` (center/window text; the freq band combobox set to the
matching preset or 'custom'); then drive the viewers (`i_drive` per axis) so the navigator reflects the
atom. Source loads the seed/radius into the Source block (no Region-tool activation needed).

## 4. Data flow

```
NAVIGATE  block edits        -> st.nav (bst_atom)         -> viewers (live); no persist
DETECT    OnDetect           -> refphase at st.nav freq   -> Brainstorm EVENTS on the recording
                                                            (time series + tree + nav); freq = live filter
          (try other bands: change Frequency block -> re-Detect -> re-filter + new events)
SAVE      OnSaveDetection    -> events -> atom groups in st.T (freq stamped numeric via bst_atom) + disk
          OnSaveCursor       -> one atom from st.nav (+ descriptor) in st.T + disk
          OnRecord/OnCapture -> commit extrema / region in st.T + disk   (existing writers)
LOAD      OnLoadAtom         -> saved atom -> st.nav + blocks + viewers  (round-trip)
```

Only the SAVE/Record/Capture paths write `st.T`. Detect writes events; Navigate writes `st.nav`.

## 5. UI changes

- The Actions row gains a **Save** grouping: `Save detection`, `Save cursor`, `Clear preview` alongside the
  existing `Detect`, `Record`, `Capture`. (`Detect` now stages events; the three Save-family buttons are the
  explicit writers.)
- An Atoms-menu **Load into navigator** item (and/or double-click on an occurrence).
- The Atoms tree shows a distinct **Detection (events)** node mirroring the recording's detection event
  group (unsaved), separate from the saved atom stacks.

## 6. Testing

- **`dev/test_detect_save.m`** (new, headless under `brainstorm nogui`):
  - `OnDetect` at the alpha band writes a detection **event group** on the recording (assert the events
    exist via the recording's event list) and adds **nothing** to `st.T`.
  - `OnSaveDetection` promotes them: `st.T` gains the `alpha (8-13 Hz)` extended group + 4 phase children
    (times match the events) and each carries a **numeric freq Localization** (`bst_atom('Get',G,'freq')`
    center≈10.5 extent≈2.5); auto-saved.
  - `OnClearDetection` removes the detection event group; `st.T` unchanged.
  - `OnSaveCursor`: set `st.nav` via the blocks, call it → exactly one new occurrence in `st.T` with the
    cursor's coords.
  - `OnLoadAtom`: select a saved atom → the four block fields + `st.nav` round-trip its coordinates.
- **`dev/test_dynamics_atoms.m`** — **T5 splits**: today it asserts `OnDetect` populates `st.T`; under the
  new contract `OnDetect` populates events and `OnSaveDetection` populates `st.T`. Update T5 to detect →
  assert events + empty `st.T` band, then save → assert the band-window + 4 children (and the numeric freq).
  T1–T4, T6–T8 stay green (Record/Capture/region/filter unchanged).
- **`dev/test_nav_panel.m`** — unchanged (T1–T6 still pass; the navigator blocks are untouched).

## 7. Out of scope (later phases)

- **Source-feature detector** — Φ/Ψ extrema + saddles on the cortex (the source analog of refphase),
  using `bst_operators`/`bst_eigen` to smooth/filter on the scale axis. A dedicated future phase: plan,
  develop, validate. This phase adds source manually.
- **Full Scale activation** (eigenvalue-spectrum navigation, eigen-band, eigen-salience) — Phase 5.
- **Atlas presets beyond Frequency** (time-event presets, anatomical source atlases) — future (§2.5).

## 8. Risks / notes

- **Events API**: writing/reading the detection event group on the loaded recording (raw link or imported)
  + refreshing the time-series. Reuse the established event mechanism (`process_evt_refphase('Run')` already
  writes events; the panel path mirrors it). Confirm the API for adding events to the in-memory recording
  and triggering the time-series redraw during implementation.
- **Detect↔Save band consistency**: Save stamps the freq from `st.nav` at save time. If the user changes the
  band between Detect and Save without re-detecting, the events (old band times) and the stamped freq could
  disagree. Mitigation: `OnSaveDetection` reads the band from the detection event group's name/stored band,
  and warns if it differs from `st.nav` freq — or simply re-derives freq from the event group. Decide in the
  plan; default: stamp from the detection event group's band (the band it was detected at), not `st.nav`.
- **Docked-panel tests require `brainstorm nogui`** (GuiLevel 0), not `server` (-1).
- **Reorganization of existing logic**: `OnDetect`'s group-building moves into `OnSaveDetection`; `OnDetect`
  now writes events. `OnRecord`/`OnCaptureRegion` are unchanged. The risk is wiring + the events API, not
  new numerics.
