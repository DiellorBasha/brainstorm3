# Frame-Based Source Imaging & Dynamics — Research Roadmap

> **Status:** exploratory roadmap (2026-07-02). Captures the conceptual framework and a phased set of
> experiments for using the atom/eigenfilter **frames inside the inverse solution** (and for the dynamics
> that motivate them). Grounded in existing modules so each phase is a concrete next step, not a rewrite.
> This is a *research* document — it lives on `development`, is not a commitment to ship any given piece.

---

## 0. North star

Replace the **arbitrary L2 smoothness** of wMNE/dSPM with a **frame prior** that is (i) physiologically
explicit (activity = a few coherent multi-scale, possibly *dynamical* structures), (ii) numerically
well-conditioned by the frame bounds rather than an ad-hoc `λ`, and (iii) **observability-aware** — aligned
with what the leadfield can actually see. The same atom frames then double as a **dynamical-regime
classifier** (diffusion vs wave vs vortex) and a **localization** engine, jointly, from sensors.

**Honest ceiling (design invariant).** None of this adds *observability*: the leadfield null space is
unrecoverable regardless of prior. This is already established in `bst_inverse_dirac` (three-axis
noise/observability/geometry picture; "observability ceiling"; Matérn gave only marginal gains). A frame
prior improves **conditioning, plausibility, interpretability, and legible failure** — not resolution
beyond the physics. Every phase below is judged against that boundary, not against ground truth it cannot see.

---

## 1. Conceptual foundation (the vocabulary the roadmap uses)

**Localization axes.**

| | scale (λ) | time | cortex |
|---|---|---|---|
| **Eigenfilter** `g(λ)` | ✅ | ✗ | ✗ (global spectral operator) |
| **Atom** (seeded, windowed) | ✅ | ✅ | ✅ (a cortical–spectral–temporal wave packet) |

An eigenfilter is the spectral-domain filter; an **atom is its wave packet** — the cortical analogue of a
Morlet wavelet, with `√λ` playing frequency's role. Dispersive kernels (`wave`, `dampedwave`,
`kleingordon`, `travwave`) tie `√λ ↔ ω`; global as eigenfilters, they become **traveling packets** as
seeded/windowed atoms.

**The enabling primitive (already built).** Composing an atom `g_m(λ)` with the imaging kernel gives a
tiny, data-independent **sensor→mode operator**

```
A_m = diag(g_m(λ)) · Φᵀ M P K_imaging          [K x nCh]      (i_vector_modekernel is the m-agnostic A = Φᵀ M P K)
```

- **Analysis** (linear reparametrization): `c_m = A_m · data` → per-scale/per-atom coefficients; energy
  `‖c_m‖²` (Parseval / per-hemi Gram) is a cheap reduction; per-vertex localization `Φ·c_m` is reconstructed
  on demand, one frame at a time. Nothing `[nVertex × nTime × nScale]` is ever stored (cf. how
  `bst_timefreq` keeps source TF in kernel space). *This is the cheap, exact regime and is essentially done.*
- **Synthesis** (sparse dictionary inversion): atoms as a **source model**, forward-projected
  `d_{s,m} = L·ψ_{s,m}`, sparse-fit to sensors. *This is the research program below.* Note the MP correlation
  step `⟨residual, d_{s,m}⟩` **is** an `A_m`-type analysis operator — the composition we built is the
  correlation engine for the sparse solve.

**Two things "nonlinear" bundles:** *non-stationarity* (moving pattern — a linear space-scale-time frame
handles it via time-localized atoms) vs *nonlinearity proper* (self-advection, mode coupling — needs an
operator-theoretic layer, Phase 5).

---

## 2. Existing infrastructure this builds on

- `bst_inverse_dirac` — whitened MNE in the Dirac mode basis; three-axis framework; the observability
  characterization is the yardstick for every inverse experiment here.
- `bst_eigen` / `tess_eigen` — operator eigenbases: Laplace–Beltrami, LB-Connectome, Connection Laplacian
  (phase-bearing), Dirac, Dirac-Connectome, Hodge-Face.
- `bst_eigenfilter` (`Design/Evaluate/Analysis/Synthesis/Bounds/RowMap/Atom/Fiber`) and `bst_eigenwavelet`
  (`Scalogram`, itersine tight frame) — the frame machinery + tight-frame bounds.
- `i_vector_modekernel` / `i_vector_coeffs` (panel_bst_dynamics) — the `A = Φᵀ M P K` composition, cached
  (2026-07-02).
- `process_helmholtz` + `bst_divergence`/`bst_curl` + Hodge — gradient (source–sink) vs rotational (vortex)
  decomposition; the residual-diagnosis tool.
- `bst_vortex_persistence` + `bst_vortex_track` — core detection + cross-frame trajectories.
- Atom kernel banks — `diffusion/heat/mexhat/matern/log/tikhonov/...` (static) and
  `wave/dampedwave/kleingordon/travwave/gabor/resonator/stmatern` (dynamic, dispersive).
- Validation assets — `dev/benchmarks` synthetic-on-real-cortex harness; posterior-alpha test segment
  (S01_AEF_01_notch, 20–25 s, 7–13 Hz, ~10.55 Hz right parieto-occipital burst ~22.6 s); Dirac
  source-vortex test block; real M100/MMN for inverse validation.

---

## 3. Phased program

### Phase 1 — Observability-aware, per-scale *linear* frame inverse
**Goal:** the cheapest principled improvement over global Tikhonov — regularize **per scale by what the
leadfield sees**, not by one global `λ`.

**Approach.** The leadfield ill-posedness is spectral (fast-decaying singular values, smooth right-singular
vectors) and the frame is spectral — align them. For each frame scale `m`, score its observability against
the head-model leadfield SVD (energy of `d_{s,m}=L·ψ_{s,m}` retained in the well-conditioned SVD subspace),
then damp noise-dominated scales and keep well-observed ones. Effectively a **diagonal, scale-resolved
regularizer** in the eigenfilter frame.

**Steps.**
1. Per-surface: leadfield SVD (or reuse the observability axis already in `bst_inverse_dirac`).
2. `observability(m)` = fraction of `d_{·,m}` energy in the top-`r` SVD subspace (sweep `r`).
3. Build a per-scale gain `w(m)` (e.g. Wiener-like `obs/(obs+noise)`), apply in the eigenfilter frame.
4. Reconstruct; compare resolution kernels / point-spread vs wMNE on the synthetic harness.

**Deliverable:** a `bst_eigenfilter`-based regularizer + an inverse variant option.
**Validation:** point-spread width, localization error, depth bias on `dev/benchmarks`; must not exceed the
`bst_inverse_dirac` observability ceiling (sanity: it shouldn't).
**Risk:** low — it's a reparametrized linear inverse; worst case it equals wMNE.

---

### Phase 2 — Structured-sparse dictionary inverse (MxNE-style)
**Goal:** focal, physiologically-plausible estimates via **sparsity in the atom frame** instead of L2 spread.

**Approach.** `min ‖b − L J‖²_Σ + μ·Ω(α)`, `J = Σ α_{s,m} ψ_{s,m}`, `Ω` = L1 or group-L21 (group atoms by
seed or by band). Extended multi-scale atoms cure the spikiness/instability of point-source L1 (MCE). MP/L1
solve; the `A_m` operators are the correlation step.

**Steps.**
1. Assemble the source dictionary `Ψ` (seeds × scales) from `bst_eigenfilter('Atom', …)`; forward-project
   `D = L·Ψ` (reuse the face/vertex leadfields).
2. Baseline solver: irMxNE (reweighted L21) or greedy OMP with the `A_m` correlation engine.
3. Coherence audit: `μ_coh(D)` and per-seed coherence maps — *quantify* the "adjacent seeds are
   indistinguishable" limit before trusting any focal claim.
4. Reconstruct focal test sources; compare vs Phase 1 and wMNE.

**Deliverable:** a dictionary-sparse inverse prototype + a coherence-audit report per head model.
**Validation:** focality + stability across data splits on synthetic focal/patch sources; MMN/M100 sanity.
**Risk:** medium — solver stability, `μ` selection (→ Phase 3 fixes this), coherence-driven instability
(quantified in step 3, not assumed away).

---

### Phase 3 — Hierarchical-Bayes / ARD frame inverse (learn the prior)
**Goal:** remove the *arbitrary* from the prior entirely — let the **data** choose which frame atoms are
active via evidence maximization. This is the principled apex of "no arbitrary smoothing."

**Approach.** Sparse Bayesian learning (ARD / Champagne lineage): each atom (or atom group) is a
**covariance component** with its own hyperparameter; maximize the model evidence to infer relevances. No
hand-set `λ`; irrelevant scales/locations are pruned automatically.

**Steps.**
1. Cast the atom dictionary as covariance components `Σ_source = Σ_γ γ_{s,m} ψ_{s,m} ψ_{s,m}ᵀ`.
2. Evidence maximization (fixed-point / EM) for `{γ}`; warm-start from Phase 2's active set.
3. Compare evidence across **restricted** dictionaries (bands only, seeds only, full) — the evidence *is*
   the model comparison.

**Deliverable:** an ARD inverse over the atom frame; per-atom relevance maps + model evidence.
**Validation:** evidence-vs-truth on synthetic; sharper alpha reconstructions than wMNE *within*
observability; graceful (low-evidence, dense) behavior under mismatch.
**Risk:** medium-high — ARD nonconvexity/local optima, compute cost; mitigate with Phase-2 warm starts and
a small seed grid.

---

### Phase 4 — Dynamics-typed inversion + model selection by coding cost
**Goal:** fuse localization with **dynamical-regime identification** — the atom bank the data codes most
sparsely is evidence for that regime (diffusion vs wave vs rotational/vortex), localized in space/scale/time.

**Approach.** Make the **dynamic atoms** (dispersive `wave/dampedwave/travwave`, diffusive `heat`,
rotational/curl from the Helmholtz Ψ side, complex from the Connection-Laplacian) the covariance components
of the Phase-3 inverse. Run per-family; compare **coding cost** (residual at fixed sparsity, or evidence).
The **Helmholtz residual** after coding is the discovery channel: net-curl residual ⇒ add rotational atoms.

**Steps.**
1. Curate a **union dictionary** across families (don't commit to one bank).
2. Per-family evidence/coding-cost on the alpha-vortex segment; rank.
3. Residual Helmholtz split → identify unmodeled (gradient vs rotational) structure → augment.

**Deliverable:** a dynamics-typed inverse that returns *where + at what scale + which dynamical regime*.
**Validation:** on the alpha burst and Dirac source-vortex blocks — does the rotational family win where you
observe vortices? Does a stationary focal source correctly select diffusion, not wave?
**Risk:** medium — dictionary coherence *across* families; interpretability of mixed selections.

---

### Phase 5 — Nonlinear / advecting dynamics (Koopman + phase singularities)
**Goal:** capture advecting vortices during oscillatory bursts (your alpha observation) — the genuinely
non-stationary / weakly-nonlinear regime a fixed linear dictionary only *describes*.

**Approach — three layers, escalating.**
1. **Complex phase-bearing basis.** Code alpha vortices in the Connection-Laplacian / Dirac complex basis,
   where a vortex is a **single phase-singularity atom** (winding phase), not a sum of real counter-rotating
   pieces. Helmholtz isolates the rotational part.
2. **Kinematics = a sparse path.** A moving structure is a time-ordered *sequence* of vortex atoms —
   exactly `bst_vortex_persistence` + `bst_vortex_track`. Gives trajectory/velocity/persistence (description).
3. **Dynamics = Koopman/DMD on atom coefficients.** Treat frame coefficients `c(t)` as **observables**, fit
   `c(t+1) ≈ M c(t)`. An advecting oscillating vortex = a **Koopman mode** with a complex eigenvalue
   (frequency = alpha) and a spatial phase gradient (advection velocity). Optionally in a **co-moving
   (Lagrangian) frame** (demodulate carrier, track core → quasi-stationary; cf. the advection montage).

**Steps.**
1. Complex-basis vortex atoms + phase-singularity detector on the analytic signal.
2. DMD on `c(t)` over bursts; extract traveling/rotating Koopman modes; validate frequency = alpha.
3. Co-moving frame demodulation; check the structure becomes stationary (few atoms) in the moving frame.

**Deliverable:** a Koopman-mode / phase-singularity analysis over the atom coefficients; co-moving montage.
**Validation:** alpha segment (10.55 Hz, right parieto-occipital burst ~22.6 s) — recover an advecting
rotating mode with alpha frequency and a coherent velocity.
**Risk:** high — DMD noise-sensitivity (works on clean coherent bursts, fails legibly otherwise); requires
sufficient window length vs the ~10 Hz carrier.

---

## 4. Cross-cutting

**Validation ladder (applies to every phase).**
1. Synthetic-on-real-cortex (`dev/benchmarks`): known sources, measure localization error / point-spread /
   depth bias / stability.
2. Ceiling check: compare against `bst_inverse_dirac` observability — flag any claim that would require
   information outside the observable subspace.
3. Real data: M100/MMN (focal, well-characterized) → then the alpha-vortex segment (the hard target).
4. Failure legibility: under deliberate dictionary mismatch, confirm the method reports it (dense code / low
   evidence / high residual) rather than silently blurring.

**Dependencies.** Phase 1 → 2 → 3 are a stack (each reuses the prior's operators/active sets). Phase 4
sits on Phase 3. Phase 5 is largely independent of 1–4 (it analyzes coefficients) but is *most* useful once
the complex/rotational atoms of Phase 4 exist.

**Non-goals / guardrails.**
- Not chasing super-resolution beyond observability — the frame improves the *estimate within* the
  observable+prior span, never the null space.
- Not committing to a single kernel family — mismatch is handled by *union dictionaries + evidence*, not by
  guessing the dynamics.
- Keep the enabling primitive honest: everything routes through `A_m = diag(g)·Φᵀ M P K` and on-demand
  per-frame reconstruction; never materialize `[nVertex × nTime × nScale]`.

**First concrete experiment (lowest cost, highest information).** Take the existing Connection-Laplacian /
Dirac atom frame, compute per-scale observability against the head-model leadfield SVD (Phase 1), and use
those as ARD covariance components (Phase 3) in the `bst_inverse_dirac` solve on the alpha-vortex segment.
One run tests whether an evidence-learned frame prior sharpens the reconstruction over wMNE **without**
exceeding the observability already characterized — and tells us whether to invest in Phases 2/4/5.

---

## 5. Math appendix (one-line reference)

- Mode kernel: `A = Φᵀ M P K_imaging` `[K×nCh]`; atom: `A_m = diag(g_m(λ)) A`; coeffs `c_m = A_m b`.
- Energy (Parseval, orthonormal Φ): `E_m(t) = ‖c_m(:,t)‖²`; per-hemi via Gram `G_h = Φ_hᵀ M_h Φ_h`.
- Sparse synthesis inverse: `min_α ‖b − L Ψ α‖²_Σ + μ Ω(α)`, `J = Ψ α`, `Ψ`=atoms, `Ω`=L1/L21.
- ARD: `Σ_src = Σ γ_{s,m} ψψᵀ`; maximize evidence `log|Σ_b| + bᵀ Σ_b⁻¹ b`, `Σ_b = Σ_noise + L Σ_src Lᵀ`.
- Koopman/DMD: `c(t+1) ≈ M c(t)`; advecting oscillation ⇒ complex eigenpair (freq = Im, advection =
  spatial phase gradient of the mode).
- Observability gate (Phase 1): `w(m) = obs(m)/(obs(m)+noise)`, `obs(m)` = `‖U_rᵀ d_{·,m}‖²/‖d_{·,m}‖²`.
