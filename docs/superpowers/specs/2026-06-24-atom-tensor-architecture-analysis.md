# Atom tensor architecture — analysis & target model

**Status:** architecture analysis / design only. No implementation footprint is committed here; the
final migration strategy is decided at planning time. This document defines the target conceptual
model for the spatiotemporal **atom** and the `panel_bst_dynamics` navigator, and recommends a phased
path from the current code.

**Author:** Diellor Basha, 2026-06-24

---

## 1. Motivation

Today an atom is, in practice, a **labelled Event that is also localized on the cortex** — the
`atomgroup` struct (`db_template`) carries Events-style `times`, plus per-occurrence `vertices`/`pos`/
`region`, plus group-level `band`/`scale`. The four axes are represented **heterogeneously** (time as
onset/offset, frequency as `[fLo fHi]`, source as seed+region vertices, scale as `[k1 k2]`), and the
panel's controls are organized by feature (Frequency band toggles, Space smoothing/operator, Record)
rather than by a single principle.

The insight that unifies them: **every axis localizes the same way — a center plus a window (extent),
or in the simplest case just a center (a point).** Once that is explicit, the atom becomes a single
index into a 4-D `(time, frequency, source, scale)` box, and the panel becomes a uniform navigator
over that box. This is the target model.

## 2. The core model

### 2.1 Localization — the one primitive

```
Localization = (center c, extent w, weighting g)   a kernel / region on its axis
  point        when w = 0                  a single t / f / vertex / eigenvalue
  window       when w > 0                  interval / band / geodesic disk / eigen-band
  unlocalized  when undefined              the axis is not yet pinned
  weighting g  hard (default) | soft       sharp cutoff vs a smooth (wavelet) decay — see 2.6
```

The **three-state** per axis (unlocalized / point / window) is essential: it is what lets a detector
produce an atom localized in **time + frequency only** (source and scale still unlocalized), and a
later step pin **source** (and **scale**). An atom's effective dimensionality is just how many of its
axes are localized.

### 2.2 The four axes

| Axis | `center` | `extent` (window) | Kernel it drives | Engine |
|------|----------|-------------------|------------------|--------|
| **Time** | central time | half-duration | boxcar → wavelet envelope | global time cursor / `panel_filter` |
| **Frequency** | center frequency | bandwidth/2 | band-pass | `panel_filter` (display filter) |
| **Source** | seed vertex (+ 3-D pos) | geodesic radius | heat-distance disk | `tess_scout_area` |
| **Scale** | center eigenvalue | eigenvalue bandwidth | eigenfilter | `panel_eigenfilter_design` |

Each kernel above is the **hard** form; every axis also has a **soft / wavelet** form (§2.6) using the
same `center`+`extent`, differing only by the weighting. Time↔Frequency is the **wavelet / Heisenberg
box**; Source↔Scale is its **graph-Fourier analog**. The four axes are therefore two
`(center,extent) × (center,extent)` conjugate planes — which is what makes joint **wavelet markers**
(a time+frequency box) and **spatial-scale markers** (a source+scale box) fall out of the same grammar
later, with no new abstraction.

### 2.3 Three concerns, kept separate

The current `atomgroup` mixes these; the target model separates them cleanly.

1. **Coordinates** — the four `Localization`s. *The only thing the navigator moves.*
2. **Descriptors** — measured *at* the atom: `strength`, `charge`, `chirality`, `persistence`, and the
   `Function`/**operator** (Φ potential / Ψ stream / |J| magnitude; later grad / Laplacian) that
   produced the measured scalar field. **The operator is a measurement selector, not an axis.**
3. **Provenance** — `DataFile` (recording), `ResultsFile` (inverse/kernel), `SurfaceFile` (cortex).

### 2.4 Atom vs group

- An **atom** is one occurrence — one index into the 4-D box (each axis unlocalized / point / window),
  with its descriptors and provenance.
- A **group** is a labelled collection of atoms that **share some coordinates** (e.g. `alpha_peak` =
  fixed frequency band + fixed time-phase, varying time and — once recorded — source). The existing
  per-occurrence arrays already *are* a set of atoms sharing the group's coordinates; the model just
  names that structure.

### 2.5 Atlas = mask = labelling

`center` and `extent` are **numeric**; an **atlas** is an *arbitrary labelling* of an axis — a partition
into named regions. **An atlas and a mask are the same object**: a labelled extent on an axis, with no
measurement attached. A scout (a cortical region), a frequency band, and an alpha time-window are all the
*same kind of thing* — a labelled region of their axis.

An atlas has two sources:

- **Presets** (supplied externally): frequency bands (`(10 Hz, ±4)` → "alpha"), cortical anatomical
  atlases (Desikan-Killiany; a region centroid is a preset `center`), time events (later); an eigenvalue
  band → a spatial scale ("column / gyrus / lobe").
- **Data-driven** (produced by detection, §4.2): alpha windows = a *time mask* from GFP-power thresholding;
  a cortical region = a *source mask* from `|J|`-norm thresholding; a band = a *frequency mask* from
  spectral peak-picking. The thresholding produces a mask exactly as GFP power produces a time mask.

An atlas **guides navigation** (snap `center` to a region or its centroid) and labels stored atoms; it
carries no field measurements. The full atlas build-out (anatomical source atlases, `|J|` masks, time
events) is **future**; the model reserves the optional-label slot and this section fixes the concept so
later phases stay aligned.

### 2.6 Weighted windows — wavelets as soft localizations

The kernels in §2.2 are **hard** windows: sharp cutoffs (a boxcar interval, a geodesic disk truncated
at an isoline, a rectangular band, an eigenfilter passband). A **wavelet** is the *same* localization
with a **weighting function** `g` that decays smoothly from the center instead of cutting off — so the
primitive is `(center, extent, weighting)` where the `extent` sets the kernel scale and `g ∈ {hard
(default), soft}` sets sharp-vs-smooth. **It is the same engine as center+extent, plus a weighting
function.**

This is already standard on the time–frequency plane: a Morlet / Gabor wavelet is a center-frequency +
bandwidth with a Gaussian envelope, smoothly localized in **both** time and frequency at once (this is
exactly what Brainstorm's time-frequency wavelet analysis does). The same construction extends to the
source–scale plane:

- **Source** — the geodesic **heat-distance** field is *already* a smooth scalar. Truncating it at an
  isoline gives the hard disk; **not** truncating it gives a heat-kernel-weighted region — a spatial
  wavelet centered at the seed (the geodesic-area engine with the weighting kept instead of thresholded).
- **Scale** — `bst_eigen` / `bst_eigenwavelet` define **spectral wavelets** on the eigenvalue axis: a
  smooth eigenfilter response (itersine tight frame) rather than a hard eigen-band.

So every axis has a hard form and a soft (wavelet) form from one `center`+`extent`, differing only by
`g`. **Joint wavelets** then follow with no new abstraction: a single smooth kernel localized
*simultaneously* on all four axes (time × frequency × source × scale) — the tensor product of the
per-axis weightings. A joint-wavelet atom is the soft-`g` generalization of the hard-window atom: same
4-D index, same engines, weighted instead of truncated. **Future development**, but the primitive is
designed for it now (the `weighting` slot).

### 2.7 Axis / Atlas / Measurement — three orthogonal layers

Each axis carries three *separate* concerns, and conflating them muddles the panel. They must stay
orthogonal:

| Layer | What it is | time / frequency / source / scale |
|-------|-----------|-----------------------------------|
| **Axis** (navigation) | pure numeric `(center, extent, weighting)` — free continuous navigation | the cursor / a frequency / a vertex+radius / an eigenvalue |
| **Atlas** = mask (§2.5) | an arbitrary labelling of the axis into named regions; **no measurement** | alpha windows / α-β-γ bands / Desikan-Killiany regions / eigen-bands |
| **Measurement** | extrema & zero-crossings of a *reference field* along the axis (differential analysis) | refphase peak/trough/rising/falling / spectral peaks / Φ-Ψ extrema / eigen-peaks |

Two consequences that the navigator and detectors must respect:

- **Detection produces a data-driven atlas, and (optionally) the measurements within it.** `process_evt_refphase`
  is the canonical example: the alpha *window* is the **atlas** (a time mask from GFP power); the phase
  *markers* are the **measurements** (extrema/zero-crossings of the reference field) inside that mask. The
  cortical analog: threshold the `|J|`-norm → a **mask** (atlas, no measurement), then take Φ/Ψ extrema →
  **measurements**. Frequency: pick the band where spectral peaks live → a mask; the peaks are the
  measurements.
- **Each axis owns an analysis/visualization toolbox** for its measurement step — power-spectral density on
  frequency, `process_evt_refphase` on time, the Helmholtz/vortex detectors on source, the eigen-spectrum
  tools on scale. The navigator *invokes* these; it does not reimplement them.

In the panel (§3), the pure axis is just the numeric `center`/`extent`; **atlas presets are a guidance
overlay** that snaps the center to a labelled region (the frequency band combobox is a *frequency-atlas
preset selector*, not the axis); **detection/measurement are a separate action** that writes a data-driven
atlas + measurements.

### 2.8 The geodesic Area tool belongs in the dynamics suite

The geodesic Area tool (heat-distance disk: a seed `center` + a radius `extent`) was added as a custom
extension to `panel_scout`, but it is **not** standard Brainstorm scout behaviour. It belongs in the
dynamics suite, because there a seed + radius is *natively* the **source-axis point+extent primitive** of
the tensor — not a bolt-on to scouts. Moving it (the `tess_scout_area` engine and the tool's
UI/figure-interaction) into the dynamics suite is the prerequisite refactor for the source-axis navigator
block. (Standard scout atlases — Desikan-Killiany and friends — remain in Scouts; they are the *source
atlas presets* of §2.5, a different concern.)

## 3. The panel as a 4-axis navigator

With one `Localization` primitive, `panel_bst_dynamics` collapses into **four identical control
blocks** — each a *center selector* + a *window control* + the kernel it drives:

```
┌ TIME ─────── center [cursor / detected marker] · window [± duration] → boxcar / wavelet
┌ FREQUENCY ── center [band / frequency]          · window [± bandwidth] → band-pass
┌ SOURCE ───── center [seed vertex]               · window [geodesic R]  → heat-distance disk
┌ SCALE ────── center [eigenvalue]                · window [± eigen-bw]  → eigenfilter
   Operator (descriptor): [Φ] [Ψ] [|J|]   ← measurement, shared by source/scale; NOT an axis
```

**The panel state *is* the atom-under-construction** — a live cursor in the 4-D box. This is a
*consolidation* of today's machinery, not new engine code: the Frequency band toggles, the Space
smoothing + operator, and the global time cursor already drive these engines; the redesign reorganizes
them under one grammar and adds the missing blocks (an explicit Time window; the Scale center/extent).

## 4. Operational model — three separated states

The single most important contract change: **Navigate and Detect never touch storage; only Save
does.** Today the panel conflates them (its "Detect windows" writes the table immediately; "Record"
saves).

```
 NAVIGATE ───────────→ DETECT ───────────→ SAVE
 live, no persist      guidance, no         commit
                       persist
```

### 4.1 Navigate (default)

Move the `center` on any axis; the **time-series viewer and 3-D figure update live**. Pure exploration
of the 4-D box — nothing is written. The window on each axis sets the active kernel (band-pass extent,
geodesic radius, eigenfilter band, time extent).

### 4.2 Detect (intermediate guidance)

Detection **produces a data-driven atlas (a mask) on an axis, and optionally the measurements within it**
(§2.7) — giving **immediate visual feedback but persisting nothing**. It steers the user toward relevant
activity. Per axis, the *mask* and the *measurement* are distinct outputs:

- **Frequency** → mask: the band where spectral peaks live (or a preset band); measurement: spectral peaks
  (PSD). Bands are also available as presets (§2.5).
- **Time** → mask: band-power periods (GFP power threshold = the alpha-window atlas); measurement:
  `process_evt_refphase` phase markers (rising / peak / falling / trough) within the mask. (Markers come
  from the monotonic GFP phase + one-per-period polarity seed — see the refphase fix; numeric phase value.)
- **Source** → mask: threshold the `|J|`-norm (or an eigenmode-smoothed field) → a cortical region (the
  same kind of object as a scout); measurement: Φ/Ψ extrema & prominences inside the mask
  (`bst_dynamics('Extrema')`, vortex/prominence detectors).
- **Scale** → mask: an eigen-band (future); measurement: eigenvalue-spectrum peaks (future).

"Snap `center` to the mask region / its measured extrema"; navigation steps among the candidates.
Detection parameters are exactly the axis windows + operator already set during navigation. **A mask is an
atlas (§2.5); the measurement is the differential-analysis step** — keep them separate.

### 4.3 Save (explicit commit)

When satisfied, persist — either:

- **this index** as one atom (snapshot the current 4-D cursor + measured descriptors), or
- **batch-run the same detection parameters over the whole recording** to harvest a group of atoms.

Only here does the `dynamics_*` table change. The batch path reuses the interactive parameters, so
"what you tuned is what you harvest."

## 5. Relation to the current code (gap analysis)

| Target concept | Today | Gap |
|----------------|-------|-----|
| Localization `(center,extent)` per axis | heterogeneous fields (`times` onset/offset, `band`, `region`, `scale`) | no uniform accessor; storage materializes windows but presents them differently per axis |
| 3-state axis (unlocalized/point/window) | implicit (empty vs filled fields; NaN seeds for unlocalized occurrences) | not named or enforced; works ad-hoc |
| Coordinates / descriptors / provenance split | all flat in `atomgroup` | conceptually mixed; descriptors and coords interleaved |
| 4-axis navigator panel | Frequency + Space + Record sections | Time-window block and Scale center/extent block missing; sections not uniform |
| Navigate / Detect / Save separation | Detect **writes**; Record **saves** | Detect must become non-persisting guidance; Save must be the only writer |
| Operator = descriptor | Space "Op" toggles drive the figure | already correct; keep, relabel as measurement |
| Atlas labels | `bandName` / `scaleName` strings | partial (frequency only); generalize later |

**No existing engine is missing** — `panel_filter`, `tess_scout_area`, `panel_eigenfilter_design`,
`process_evt_refphase`, `bst_dynamics('Extrema')`, the vortex detectors all exist and already back the
four axes. The work is *reorganization under one grammar* + the **Navigate/Detect/Save** contract,
not new numerics.

## 6. Recommended phased path (decided at planning time)

A migration that preserves the just-built region/phase/filter/detect work and the Events parity:

1. **Localization accessor (non-breaking).** ✅ DONE (`toolbox/dynamics/bst_atom.m`, merged 78f72a0d).
   Uniform `(center, extent, weighting)` Get/Set over the existing `atomgroup` fields (time ⇄ onset/offset,
   frequency ⇄ `[fLo fHi]`, source ⇄ seed + new `radius`, scale ⇄ `[k1 k2]`). No storage change.
2. **Move the geodesic Area tool into the dynamics suite (§2.8).** Relocate the `tess_scout_area` engine
   and the Area-tool UI/figure-interaction out of `panel_scout` into the dynamics suite, so seed + radius
   become the native source-axis point+extent primitive. Prerequisite for the source-axis navigator block.
3. **Navigator panel.** Reorganize `panel_bst_dynamics` into four uniform `(center,extent)` blocks
   driven through the accessor; add the missing Time-window and Scale blocks. Relabel the operator as a
   measurement selector. **Keep Axis / Atlas / Measurement separate (§2.7):** the axis blocks are pure
   numeric navigation; atlas presets (frequency bands) are a guidance overlay; detect/measure stay actions.
4. **Navigate / Detect / Save separation.** Make Detect a non-persisting guidance overlay that produces a
   data-driven atlas (mask) + measurements (§2.7); route all writes through an explicit Save (single index
   or batch).
5. **Scale axis activation.** Wire the eigenvalue center/extent to the eigenfilter + an eigen-salience
   detector (currently `scale`/`scaleName` are reserved/unused).
6. **(Future) Atlas system** — presets (anatomical source atlases, ephys band presets) + data-driven masks
   as first-class labellings (§2.5); **(future) joint wavelet (time+freq) and spatial-scale (source+scale)
   box markers** (§2.6) from the conjugate-plane grammar.

Each phase is independently testable and reversible; phase 1 (the accessor) is the keystone that makes
every later phase uniform. Phase 2 (geodesic-tool move) unblocks the source-axis navigator block.

## 7. Open questions for the planning stage

- Canonical time representation: keep onset/offset (Events parity) as storage with `(center,extent)`
  as the accessor view (recommended), or migrate storage to center/half-width.
- Navigator replace decision: **resolved** — the four blocks replace the current Frequency/Space/Record
  sections; Detect/Record/Capture stay as actions until the Phase-4 contract.
- Geodesic-tool move: where the `tess_scout_area` engine lands (a dynamics-suite home) and how the
  figure-interaction (currently `panel_scout`-routed click/scroll) is re-homed — decided in that phase's spec.

## 8. Out of scope (here)

- Any code changes — this is analysis only.
- The atlas/labelling layer (future, §2.5).
- Soft/weighted-window **wavelets** and **joint wavelets** (future, §2.6) — the model and the
  `weighting` slot are designed for them, but implementation is deferred.
- Scale-axis detector science (eigen-salience) beyond noting where it plugs in.
