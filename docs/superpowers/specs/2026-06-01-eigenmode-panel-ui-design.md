# EigenModes Panel — UI Refinement Design

**Date:** 2026-06-01
**Author:** Diellor Basha (with Claude)
**Status:** Approved design — ready for implementation plan
**Builds on:** the eigenmode lever (`panel_eigenmodes.m`, merged)

## Goal

Make the `EigenModes` panel's controls intuitive. Today it exposes the lever's
internals — two unlabeled band-edge sliders (`lo`/`hi`) and a cryptic `lo=1 c=15 hi=30`
label, with `Single` mixed in as a 4th "Window" radio. Replace this with a
**center mode + width + shape** control surface. Only the panel's control surface
changes; the lever engine (`BuildWeights`, `ApplyToColumn`, paired-rank, the viewer)
is untouched.

## Mental model

A selection = **a center mode + a width + (when width > 0) a weighting shape**:

- **Center mode** — which spatial frequency (paired rank) you're focused on; the slider
  thumb. Arrow stepping moves it.
- **Width (± modes)** — how many modes around the center. `0` = a single mode.
- **Shape** — how modes in the band are weighted: `Box` (flat), `Taper` (cosine edges,
  anti-ringing), `Gauss` (bell). Meaningless at width 0, so the shape row **greys out**
  when width = 0. "Single" is no longer a shape — it is simply width 0.

## Layout

```
+- Spatial scale (eigenmodes) ----------+
| [x] Active                            |   (shown only in source-map context)
|                                       |
| Center mode                           |
| 1 ----[-----o-----]------------ 600   |   one JSlider: thumb = center;
|       30   42    55                    |   axis tick-labels (1..K) + [ ] window
|                                       |    markers at center-width / center+width
| Width (+/- modes):  [ 13 ]            |   text field; 0 = single mode
|                                       |
| Shape:  (o) Box  ( ) Taper  ( ) Gauss |   greyed when Width = 0
|                                       |
| Keeping 27 modes      lambda 9 - 61   |   readout
+---------------------------------------+
```

- **Center `JSlider`** over `[1, K_paired]`. `setPaintLabels(true)` with a **label table**
  carrying (a) a few numeric **axis labels** (`1`, ¼K, ½K, ¾K, `K`) so position in the
  spectrum is visible, and (b) `[` / `]` **window markers** at `center−width` /
  `center+width`. The label table is rebuilt whenever center or width changes; at width 0
  no markers are drawn (just the thumb).
- **Width** text field (`JTextField`), integer ≥ 0; clamped to `[0, K−1]`. `0` = single mode.
- **Shape** radios `Box` / `Taper` / `Gauss` (the lever's `box`/`tapered`/`gain`). Disabled
  when width = 0.
- **Active** checkbox and **readout** label unchanged in role.

## Control → lever-state mapping (engine unchanged)

The panel translates Center + Width + Shape into existing lever verbs:

- **Center slider released** → `SetBand(center−width, center+width)` (this sets
  `iCurrentMode = center` and slides the band; `width` from the field).
- **Width field committed** → `SetBand(center−width, center+width)` around the current center.
- **Shape radio** → `SetWindowShape('box'|'tapered'|'gain')`.
- **Width = 0** → `SetWindowShape('single')` and grey the shape radios.

A pure helper `BandFromCenterWidth(center, width, K)` returns the clamped `[lo, hi]` band
and is unit-tested.

While dragging the center slider, a lightweight `StateChanged` handler updates the markers
and readout (preview) without committing; `MouseReleased` commits via `SetBand` (so the
expensive lever recompute / synthesis fires once on release) — the `panel_freq`
quick-preview pattern.

## RefreshControls

Reflects lever state back into the controls: thumb ← `iCurrentMode`; width field ←
half-width derived from the band (`round((hi−lo)/2)`); rebuild the axis + `[ ]` marker
label table; grey the shape radios when width = 0; set the readout
(`Keeping N modes   lambda a - b`, or `Mode k   lambda x` when width 0).

## Control struct changes

Replace `jSliderLo`, `jSliderHi`, `jLabelBand`, `jRadioSingle` with `jSliderCenter`,
`jTextWidth`. Keep `jCheckActive`, `jLabelReadout`, `jRadioBox`, `jRadioTaper`,
`jRadioGain` (relabel the Gain radio text to "Gauss"; the lever shape string stays
`'gain'`). `UpdatePanel`/`SetSelectEnabled` update their referenced control names; the
selection-control set becomes `{jSliderCenter, jTextWidth, jRadioBox, jRadioTaper,
jRadioGain}` (the Active toggle stays gated separately by context).

## Error handling

- **Non-numeric / out-of-range width** → parse failure or out-of-range falls back to the
  last valid width (clamp to `[0, K−1]`); never crash.
- **Width ≥ K** → clamped; band becomes `[1, K]`.
- **No eligible view** (lifecycle, unchanged) → controls disabled, readout "no eigenmode view".

## Testing

- **Pure (headless):** `BandFromCenterWidth(center, width, K)` — clamping at both ends,
  `width = 0 → [c, c]`, symmetric interior, `width ≥ K → [1, K]`.
- **Panel smoke (MCP):** `CreatePanel` builds; the control struct exposes
  `jSliderCenter`, `jTextWidth`, `jRadioBox/Taper/Gauss`, `jCheckActive`, `jLabelReadout`.
- **Regression:** the existing lever/viewer suite still passes (engine untouched).
- **e2e (live, best-effort):** open a source map, set center + width via the panel, confirm
  the displayed band matches; width 0 greys the shape row.

## Out of scope

- A custom Java two-handle range slider (decided against — keep it pure-MATLAB).
- Any change to the lever engine, paired-rank logic, or the viewer.
