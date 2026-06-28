# Geometry-based epicenter detection of cortical amyloid/tau — design

**Date:** 2026-06-28  **Branch:** `experimental/pet` (analysis layer; may fork `experimental/pet-analysis`)
**Status:** Approved

## Goal / scientific question
Detect, per subject, the cortical **concentration foci (epicenters)** of amyloid (18F-NAV4694) and
tau (18F-flortaucipir) deposition from the validated surface SUVR maps, using differential-spectral
geometry — a **dominant epicenter + persistence-ranked secondary foci**. Scientific anchor: do the
geometric epicenters land in the **known seeding/vulnerable regions** (amyloid: precuneus / posterior
& isthmus cingulate / medial-orbitofrontal; tau: entorhinal / inferior-temporal / fusiform)?

This is an **analysis layer on the pipeline outputs**, scientifically distinct from the preprocessing.

## Scope (converged with the user)
- **Epicenter = current concentration focus** — where the SUVR gradient field *converges*
  (Hodge / divergence). Model-free, per-subject, static.
- **Hierarchical** — dominant epicenter + persistence-ranked secondary foci (captures multifocality).
- **Not** an inferred seeding origin (no spreading model); **not** longitudinal.

## Method — variant B: covariant gradient → Morse-Smale → persistence
Input: a subject's surface SUVR map (vertex 0-form), on the cortex surface.
1. **Regularize** — LBO heat smooth via `tess_operators` `(M + t·K)⁻¹` solve: suppress voxel noise
   *tangentially* without blurring across the boundary (the smoother built for `mri_bbregister`).
2. **Covariant gradient** — `bst_gradient` / `tess_operators` → per-face vector field `∇(SUVR)`
   (the "deposition flow"); it converges at maxima.
3. **Critical points (foci)** — the surface Laplacian `Δ(SUVR)` (= divergence of the gradient,
   `tess_operators`) localizes maxima; foci = vertices where the discrete gradient vanishes with a
   concave (negative-Laplacian) neighborhood. The Hodge view: foci are where the flow converges.
4. **Basins** — discrete Morse-Smale / watershed: assign each vertex to the maximum its gradient
   ascent reaches → partition the cortex into catchment basins, one per focus.
5. **Persistence ranking** — topological persistence on the super-level-set filtration of the SUVR:
   each focus gets a persistence (birth→death over the threshold sweep); rank dominant→secondary and
   drop low-persistence (noise) foci. Adapts the persistence concept from `bst_vortex_persistence`
   to scalar critical points.

**Output per subject:** ranked foci `[vertex, peakSUVR, persistence, basinArea]` → dominant
epicenter + secondaries (+ a basin label map on the surface).

## Group analysis & validation
1. **Location validation (scientific anchor).** Map each subject's dominant (and secondary) foci to
   the Desikan region; test concentration in the expected regions:
   - amyloid: precuneus, posterior/isthmus cingulate, medial-orbitofrontal (Centiloid/DMN);
   - tau: entorhinal, inferior-temporal, fusiform, amygdala (Braak I–IV).
   Metrics: fraction of subjects whose dominant focus is in the expected region; a cohort
   **focus-density map** (where foci concentrate across subjects).
2. **Robustness vs naïve SUVR-argmax.** Compare the geometric foci to the raw vertex-argmax:
   reproducibility (under added noise / sharp-vs-smoothed map) and anatomical plausibility.
3. **Topology (multifocality).** Number of persistent foci and the 2nd/1st persistence ratio —
   expect amyloid multifocal/diffuse vs tau focal/temporal.
4. **Persistence as a positivity/severity marker (by-product).** Dominant-focus persistence/strength
   across all 66 vs global cortical SUVR / amyloid+ status — does focus persistence index positivity?

Emphasis: **locations** on the ~7 amyloid+ and any tau-positive; **persistence-marker** on all 66;
**topology** amyloid vs tau.

## Data scope & honest constraints
- All 66 subjects, both tracers. **Foci are detected in each subject's own cortex space** (on the
  per-subject surface SUVR), and **region-assigned via that subject's Desikan atlas**. For the cohort
  **focus-density map**, the per-subject foci are mapped to the default-anatomy template via the
  registration spheres (`tess_interp_tess2tess`).
- **Cross-sectional** → current concentration foci, not seeding origin; no dynamics.
- **Small positive n (~7 amyloid+)** → method development + descriptive validation, not clinical claims.
- The epicenter of a diffuse field is inherently fuzzy → **persistence + basins make it well-defined**
  and the dominant/secondary ranking honest about multifocality.

## Deliverables
- **`pet_epicenter.m`** (analysis function): `[foci, basinLabel, info] = pet_epicenter(SurfaceFile,
  suvrMap, Opts)` — heat-smooth → covariant gradient → critical points → Morse-Smale basins →
  persistence ranking. Reuses `tess_operators` (LBO + Laplacian), `bst_gradient`, the heat solve.
  Opts: `.HeatT`, `.MinPersist`, `.nFociMax`.
- **Cohort analysis script** (`dev/benchmarks/`): run on all 66; the group validation + figures
  (per-subject foci on the inflated surface; cohort focus-density map; known-region validation table;
  amyloid-vs-tau topology).

## Alternatives considered
- **Variant A** (Laplacian peaks + watershed, no explicit gradient field): simpler, less principled — rejected.
- **Variant C** (scale-space persistence over heat scales): most robust, more work — deferred as a refinement.
- **Inferred seeding origin** (diffusion-model deconvolution): rejected (cross-sectional; user chose
  current concentration focus).
- **Single epicenter only**: rejected (hierarchical chosen to capture multifocality).

## Out of scope (future)
Longitudinal spreading dynamics; connectome-based (white-matter) propagation; the inferred-origin
diffusion model; the geometric-eigenmode reconstruction question (a separate study).
