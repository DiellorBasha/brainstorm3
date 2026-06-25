# Atom navigator panel — design

**Goal:** Replace `panel_bst_dynamics`'s feature-organized control sections (Frequency / Space / Record)
with **four uniform `(center, extent)` axis blocks** (Time / Frequency / Source / Scale) + a Measurement
selector, all driven through the `bst_atom` accessor into a transient cursor atom — the live **Navigate**
state of the atom-tensor model.

**Status:** design, ready to plan. Phase 3 of the atom-tensor architecture
(`docs/superpowers/specs/2026-06-24-atom-tensor-architecture-analysis.md`, §6 step 3). Builds on Phase 1
(`bst_atom` accessor, merged) and Phase 2 (`bst_geodesic_tool`, merged).

**Author:** Diellor Basha, 2026-06-24

---

## 1. Motivation

Today `panel_bst_dynamics`'s control area is organized by feature, not by the four axes: a **Frequency**
section (band toggles + Detect), a **Space** section that mixes the scale axis (eigenfilter smoothing)
with the operator (Φ/Ψ measurement), and a **Record** section. Source has no block (it was done via the
Scout tool), Time has no block (just the global cursor). The atom-tensor model says every axis localizes
the same way — `(center, extent)` — so the panel should present **four identical blocks**, one per axis,
with the axis / atlas / measurement layers kept distinct (analysis §2.7).

## 2. Decisions (locked in brainstorming)

1. **Replace** the Frequency / Space / Record-navigation sections with four uniform blocks + a Measurement
   row. The Atoms table and the Detect / Record / Capture **actions** are kept (the Region-tool toggle
   moves into the Source block).
2. **Symmetric layout** (the rule): every block is `center [ ] · window [± ]` on the LEFT (uniform,
   including Source) and **the selector/preset in the SAME right-hand slot** — so the user learns "the
   selector is always on the right." Source uses center+window like the rest (center = seed vertex,
   window = radius), with the Region tool in the right slot — the position the band combobox holds for
   Frequency.
3. **Basic Scale block now**; full eigenvalue-spectrum navigation + eigen-salience detection is Phase 5.
4. **Navigate-only**: the blocks drive the linked viewers live; nothing persists. Detect / Record /
   Capture remain the only writers, unchanged this phase. The Navigate / Detect / Save contract is Phase 4.
5. **Axis vs Atlas vs Measurement kept distinct** (§2.7): the numeric center/window IS the pure axis; the
   right-slot selector (band combobox = frequency-atlas preset; Region tool = source picker) is the
   selection affordance; the operator (Φ/Ψ/|J|) is a Measurement descriptor, not an axis.

## 3. Layout

```
┌ TIME ──────────────────────────────────────────┐
│ center [ 0.911 s ]  window [± 0.05 s ]   [  —  ]│   right slot: time-events preset (future)
├ FREQUENCY ──────────────────────────────────────┤
│ center [ 10.5 Hz ]  window [± 2.5  ]   [alpha ▾]│   right: band combobox (frequency-atlas preset)
├ SOURCE ─────────────────────────────────────────┤
│ center [ v12043  ]  window [ 6 mm  ]   [◉Region]│   right: Region tool (picks center, grows window)
├ SCALE ──────────────────────────────────────────┤
│ center [ λ …     ]  window [± λ-bw ]   [  —  ]  │   right: eigen preset (future); center reserved (Ph 5)
└────────────────────────────────────────────────┘
  Measurement:  [Φ] [Ψ] [|J|]      ← descriptor (view_helmholtz SetComponent), NOT an axis
  Peaks [3]   [ Detect ]  [ Record ]  [ Capture ]   ← actions, unchanged this phase
```

- **Left two columns are always numeric** `center`/`window` (typeable + driven).
- **The right slot is always "the selector"** for that axis: empty for Time/Scale (until their presets
  exist), the band combobox for Frequency, the Region-tool toggle for Source.

## 4. Architecture

### 4.1 The block-builder

One builder produces the symmetric row; all four axes call it identically — only the right selector and
the engine driver differ.

```matlab
% Build a uniform axis block: [center field] [window field] [right selector slot].
%   parent       the jCtrl panel
%   axis         'time'|'freq'|'source'|'scale'
%   labels       struct('center', 'Hz'|'s'|'vertex'|'λ', 'window', '± Hz'|...) unit hints
%   rightSel     the axis-specific selector component (combobox / toggle), or [] for none
% Returns h = struct('center', jCenterField, 'window', jWindowField).
function h = i_axis_block(parent, axis, labels, rightSel)
```

Editing a field or operating the right selector → `OnAxisChange(axis)`.

### 4.2 The transient cursor atom + OnAxisChange

The panel holds `st.nav` = a transient single-occurrence group (`bst_dynamics('NewGroup', 'cursor')`),
the atom-under-construction. Every block writes its axis through the Phase-1 accessor:

```matlab
function OnAxisChange(axis)
    [ctrl, st] = i_cs();
    loc = i_read_block(ctrl, axis);                  % (center, extent) from the block's fields
    st.nav = bst_atom('Set', st.nav, axis, 1, loc);  % uniform write, all 4 axes
    setappdata(0, 'DynamicsTarget', st);
    i_drive(axis, loc);                              % thin per-axis engine driver
end
```

### 4.3 Per-axis drivers (reuse existing engines — no new DSP)

| Axis | `i_drive(axis, loc)` | reuses |
|------|----------------------|--------|
| time | `panel_time('SetCurrentTime', loc.center)` | global cursor |
| freq | `panel_filter('SetFilters', 1, loc.center+loc.extent, 1, loc.center-loc.extent, 0, [], 0, 1)` | today's `OnBand` |
| source | `bst_geodesic_tool` (Toggle via the right-slot Region toggle; GetState seed/radius → the center/window fields) | Phase 2 |
| scale | `loc.extent > 0` → `view_helmholtz('SetSmoothing', st.hFig, 1, 'heat', params(loc.extent))`; else off | today's `OnSpaceSmooth` |

- **Frequency right slot** = a band combobox (δ/θ/α/β/γ/custom = `i_bands()`); selecting a band fills the
  freq center/window fields and calls `OnAxisChange('freq')`. Custom = leave the numeric fields editable.
- **Source right slot** = the Region-tool toggle → `bst_geodesic_tool('Toggle', state)`. After a pick or
  scroll, `bst_geodesic_tool('GetState')` populates the Source center (seed vertex id) + window (radius mm)
  fields (a small sync called from the panel, e.g. when the field is focused or via the existing draw path).
- **Scale**: the eigenvalue **center** field is present for symmetry but **reserved** — its band-pass
  effect and the eigen-spectrum navigation are Phase 5. The **window** drives a low-pass heat eigenfilter
  via `view_helmholtz('SetSmoothing')`. The kernel/slider UI from the old Space section is dropped (a
  Phase-5 concern); the kernel defaults to the low-pass heat kernel.
- **Measurement** (Φ/Ψ/|J|): `view_helmholtz('SetComponent', st.hFig, 'Irrot'|'Solen'|'Total')` — today's
  `OnSpaceComp`, relabeled and moved out of the Space section. Mutually exclusive toggles; all-off = Total.

### 4.4 Controls struct + state

New handles replace the old `jBands` / `jSpace*`:
`jTimeC`, `jTimeW`, `jFreqC`, `jFreqW`, `jFreqBand`, `jSrcC`, `jSrcW`, `jRegionTool`, `jScaleC`, `jScaleW`,
`jMeasPot`, `jMeasStr`. Kept: `jTree`, `jListOccur`, `jMenuFile`, `jMenuAtoms`, `jPhaseItems`, `jPeaks`.
`st.nav` is added to the `DynamicsTarget` appdata struct. Removed callbacks: `OnBand` (folded into the freq
block), `OnSpaceSmooth`/`OnSpaceKernel` (folded into the scale block, kernel fixed to heat). `OnSpaceComp`
becomes `OnMeasurement`. The Region-tool toggle callback (`ctrl_region_state` → `bst_geodesic_tool`) moves
into the Source block.

## 5. Data flow (Navigate)

```
user edits a block field / operates its right selector
   → OnAxisChange(axis)
      → loc = read (center, extent) from the block
      → st.nav = bst_atom('Set', st.nav, axis, 1, loc)
      → i_drive(axis, loc)   → panel_time / panel_filter / bst_geodesic_tool / view_helmholtz
   → the linked time-series + Helmholtz 3-D viewers update live
(nothing persists; Detect / Record / Capture remain the only writers)
```

## 6. Testing

- **`dev/test_nav_panel.m`** (new, headless under `brainstorm nogui` — GuiLevel 0, NOT server mode):
  - `OnAxisChange('freq')` with center=10, ext=2 → `st.nav` freq Localization round-trips
    (`bst_atom('Get', st.nav, 'freq')` center≈10 ext≈2) and the band-pass display filter is set.
  - The freq band combobox preset (select 'alpha') fills the freq center/window fields to 10.5 / 2.5.
  - `OnAxisChange('source')` after `bst_geodesic_tool('Seed', Surf, vi)` syncs the source center field to
    `vi` and window to the radius; `st.nav` source Localization reflects them.
  - `OnMeasurement('Irrot')` calls `view_helmholtz('SetComponent', …, 'Irrot')` (verify via the figure's
    HelmholtzState component, as the existing Space test did).
  - Panel constructs without error; the controls struct has the new handles and not the old `jBands`.
- **`dev/test_dynamics_atoms.m`** — T1–T8 must stay green: the Atoms table, Detect, Record, Capture, and
  the Show-phases filter are untouched. (Note: tests that drive the panel reference control handles; any
  T1–T8 step that used `jBands`/`jSpace*` is updated to the new block handles — e.g. T4 selected the alpha
  band via `jBands(3)`; it now selects via the freq band combobox.)

## 7. Out of scope (later phases)

- **Navigate / Detect / Save contract** (Detect/Capture become non-persisting guidance; explicit Save) —
  Phase 4.
- **Full Scale activation** — eigenvalue-spectrum navigation, eigen-band band-pass, eigen-salience
  detection — Phase 5. (This phase: basic low-pass via the window; center reserved.)
- **Atlas presets beyond Frequency** — time events, source anatomical atlases (Desikan-Killiany centroids),
  eigen presets — future (§2.5).
- **Loading a selected list atom back into the navigator blocks** — a Phase-4 nicety; Navigate is free
  exploration this phase.

## 8. Risks / notes

- **T1–T8 handle churn**: the existing atom suite drives the panel via `jBands`/`jSpace*` handles in a few
  steps (e.g. T4 `jBands(3).doClick()`, the Space `jSpaceStr`/`jSpacePot` toggles). These steps must be
  retargeted to the new block handles (freq band combobox; the `jMeasStr`/`jMeasPot` measurement toggles)
  in the same change, or T4/T6 break. Plan the panel rebuild + the test retarget as coupled.
- **Docked-panel tests require `brainstorm nogui`** (GuiLevel 0), not `brainstorm server` (-1) — a known
  env gotcha; `bst_get('PanelControls','Dynamics')` is `[]` under server mode.
- **Source field sync**: the seed/radius fields reflect `bst_geodesic_tool('GetState')`. The sync point
  (when the fields update after an interactive pick/scroll) should be lightweight — update on the Region
  tool's draw or when the Source block gains focus; do not poll.
- **Reorganization, not new numerics**: every driver calls an existing engine verbatim; the risk is wiring,
  not algorithms. The `bst_atom` round-trip is the correctness anchor for each block.
