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

### 2.5 Atlas — labelling layer (future)

`center` and `extent` are **numeric**; each axis localization may carry an **optional human-readable
label**. A tensor **atlas** maps regions of any axis to names:

- frequency `(10 Hz, ±4)` → "alpha";
- an eigenvalue band → a spatial scale: "column / gyrus / lobe";
- a geodesic region → an anatomical label.

The atlas is a *labelling layer over the numeric tensor* — deferred, but the model reserves the
optional-label slot so it can be added without a schema change. **Future development.**

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

Each axis runs a detector that **highlights salient centers** for the current window/params, giving
**immediate visual feedback but persisting nothing**. It steers the user toward relevant activity:

- **Frequency** → known ephys bands (alpha / beta / gamma / …) — already guiding.
- **Time** → `process_evt_refphase` phase markers (rising / peak / falling / trough) inside band-power
  periods. (Phase markers now come from the monotonic GFP phase + one-per-period polarity seed —
  see the refphase fix; they carry a numeric phase value.)
- **Source** → eigenmode-smoothed Φ/Ψ peaks & prominences: the **scale** axis localizes (eigenfilter
  smoothing), then the **operator's** extrema are the candidate seeds (`bst_dynamics('Extrema')`,
  vortex/prominence detectors).
- **Scale** → (future) eigenvalue salience / dispersion features.

Detection = "snap `center` to the salient features on this axis"; navigation steps among the
candidates. Detection parameters are exactly the axis windows + operator already set during navigation.

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

1. **Localization accessor (non-breaking).** Add a uniform `(center, extent)` Get/Set layer over the
   existing `atomgroup` fields (time ⇄ onset/offset, frequency ⇄ `[fLo fHi]`, source ⇄ seed+radius,
   scale ⇄ `[k1 k2]`). No storage change; both representations stay in sync. This delivers the
   conceptual unification and lets the panel treat all axes uniformly.
2. **Navigator panel.** Reorganize `panel_bst_dynamics` into four uniform `(center,extent)` blocks
   driven through the accessor; add the missing Time-window and Scale blocks. Relabel the operator as
   a measurement selector.
3. **Navigate / Detect / Save separation.** Make Detect a non-persisting guidance overlay (per-axis
   salient-center highlighting); route all writes through an explicit Save (single index or batch).
4. **Scale axis activation.** Wire the eigenvalue center/extent to the eigenfilter + an eigen-salience
   detector (currently `scale`/`scaleName` are reserved/unused).
5. **(Future) Atlas labelling layer** over the numeric tensor; **(future) joint wavelet (time+freq)
   and spatial-scale (source+scale) box markers** from the conjugate-plane grammar.

Each phase is independently testable and reversible; phase 1 is the keystone (the accessor is what
makes every later phase uniform).

## 7. Open questions for the planning stage

- Canonical time representation: keep onset/offset (Events parity) as storage with `(center,extent)`
  as the accessor view (recommended), or migrate storage to center/half-width.
- Source `extent`: store the geodesic **radius** as the window parameter (materialized region cached),
  so the source window is a single dialable scalar like the other axes.
- Whether the navigator replaces or coexists with the current sections during the transition.

## 8. Out of scope (here)

- Any code changes — this is analysis only.
- The atlas/labelling layer (future, §2.5).
- Soft/weighted-window **wavelets** and **joint wavelets** (future, §2.6) — the model and the
  `weighting` slot are designed for them, but implementation is deferred.
- Scale-axis detector science (eigen-salience) beyond noting where it plugs in.
