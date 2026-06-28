# Longitudinal Aβ/tau spread — generative model + reaction-diffusion coupling inversion (design)

**Date:** 2026-06-28  **Branch:** `experimental/pet` (analysis layer)
**Status:** Approved

## Goal / scientific question
Recover the **Aβ→tau spatiotemporal coupling** (the amyloid-drives-tau cascade) from a *longitudinal*
surface-PET series, by inverting a coupled reaction-diffusion model on the cortical manifold.
PREVENT-AD PET is **cross-sectional**, so the study is built on a **synthetic generative model**
(ground truth) and cross-checked against the real cohort as a **pseudo-longitudinal** trajectory.

This is the longitudinal extension of the epicenter work, exercising the synergy between Brainstorm's
spatiotemporal environment (a longitudinal PET series IS a `[nVertex × nTime]` source-map) and the
differential/spectral geometry toolkit (LBO diffusion via `tess_operators`).

## Scope (converged with the user)
- **Minimal compelling study**: the generative model + the coupling estimator + validation.
- Coupling estimator = **reaction-diffusion model inversion** (recover global κ and the diffusion/
  growth parameters), implemented as a **linear-in-parameters regression**.
- Headline validation = the **κ-recovery curve** (sweep κ, confirm κ_est tracks κ_true).
- Substrate = **surface LBO** (white-matter connectome deferred to future work).
- **Not** in scope (deferred): spreading-front tracker, Hodge source/flux/sink decomposition,
  eigenmode-dynamics, full Brainstorm source-series/atoms integration.

## Component A — `pet_spread_simulate` (generative model)
Productionized from the prototype `dev/benchmarks/proto_pet_spread.m`. Two coupled Fisher-KPP
reaction-diffusion fields on the LH cortical manifold; LBO diffusion (`tess_operators`
Mass `M` / stiffness `K`), implicit-Euler diffusion + explicit logistic reaction:
```
∂a/∂t = r_a · a·(1−a)            + D_a · Δa          (Aβ,  seed = precuneus)
∂τ/∂t = r_τ · (1 + κ·a)·τ·(1−τ)  + D_τ · Δτ          (tau, seed = entorhinal; growth gated by Aβ)
```
- Seeds from the validated epicenters (Desikan `precuneus L`, `entorhinal L`).
- `[a, tau, info] = pet_spread_simulate(SurfaceFile, Opts)` with `Opts` = `{nT, dt, Da, Dt, ra, rt,
  kappa, seedAmp, seedA, seedT}`; returns `a`,`tau` as `[nLH × nT]`.
- Prototype-confirmed: spread + Aβ-precedes-tau cascade + a strong, tunable coupling
  (t=24 tau coverage 61% at κ=3 vs 9% at κ=0).

## Component B — `pet_spread_invert` (coupling estimator)
The PDE is **linear in its parameters**, so inversion is a regression. For tau:
```
dτ/dt = r_τ·[τ(1−τ)] + (κ·r_τ)·[a·τ(1−τ)] + D_τ·[Δτ]
        = β1·X1       + β2·X2              + β3·X3
```
- `dτ/dt` from finite differences (central where possible); `Δτ = X3` from the LBO. Use the **same
  Laplacian and sign convention as the simulator** — `Δ = −M⁻¹K` (heat-diffusion sign) — so the
  recovered `D_τ = β3 > 0`. A sign mismatch flips `D`; the validation must check `D_est > 0`.
- Stack all `(vertex, time)` samples, solve the linear least-squares for `β = [β1 β2 β3]`.
- Recover: `r_τ = β1`, `D_τ = β3`, and **`κ = β2/β1`**. Aβ fit (no coupling term) recovers `r_a, D_a`.
- `[est, info] = pet_spread_invert(SurfaceFile, a, tau, dt)` → `est = {ra, Da, rt, Dt, kappa}`,
  `info` = design-matrix condition number + residuals.
- **Identifiability handling:** report `cond(X)`; apply light ridge regularization if ill-conditioned;
  heat-smooth the fields before differencing to control finite-difference noise.

## Component C — validation (`dev/benchmarks/validate_pet_spread_coupling.m`)
1. **Parameter recovery** — simulate with known `{ra,Da,rt,Dt,κ}` (κ=3) → invert → recover each within
   tolerance (e.g. κ within ±20%). Seeds recoverable from the early-time argmax.
2. **κ-recovery curve (headline)** — sweep κ ∈ {0, 1, 3, 5} → plot κ_est vs κ_true; expect a monotone,
   near-identity relationship. The figure that proves the interaction is recoverable.
3. **Robustness** — (a) additive noise on the maps; (b) **subsample to few timepoints** (3–5, the real
   PET regime) → κ recovery degrades gracefully; report the degradation.
4. **Pseudo-longitudinal cross-check** — order the real 66-subject cohort (surface SUVR, amyloid and
   tau) by global SUVR severity to form a pseudo-time series; fit the model; sanity-check the recovered
   dynamics are physiologically sensible. Honestly labeled pseudo-time, not true longitudinal.

## Data scope & honest constraints
- Synthetic data on the real cortex (LH) for development/validation; real cohort only as the
  pseudo-longitudinal cross-check.
- **Surface geodesic diffusion ≠ white-matter connectome** — tau spreads trans-synaptically; surface
  diffusion is a first-order approximation. Connectome substrate (and the connection Laplacian for
  directed spread) are future work.
- **Finite-difference derivatives are noise-sensitive** — heat-smooth + regularize; report identifiability.
- **Real longitudinal PET is sparse in time** (2–4 timepoints over years) — the sparse-sampling test is
  required, not optional; it sets expectations for real-data application.
- **Synthetic circularity** — the inversion could recover what was built in; mitigated by the κ-sweep
  (a *curve*, not one point) and the independent pseudo-longitudinal cross-check.

## Deliverables
- `toolbox/anatomy/pet_spread_simulate.m` — generative model (coupled Fisher-KPP on the LBO).
- `toolbox/anatomy/pet_spread_invert.m` — reaction-diffusion coupling inversion (linear regression).
- `dev/benchmarks/validate_pet_spread_coupling.m` — recovery, κ-sweep, robustness, pseudo-longitudinal,
  with figures (κ-recovery curve; spread dynamics; pseudo-longitudinal cross-check).

## Alternatives considered
- **Model-free lagged regression** (Granger-like coupling map): spatially-resolved, model-agnostic —
  not chosen (user chose model inversion, which recovers the generative κ directly).
- **Forward-simulation fit** (optimize κ by repeated forward PDE sims): robust but slow and
  optimizer-dependent — rejected in favor of the linear-in-parameters regression (faster, well-posed).
- **Broader toolkit / Brainstorm integration**: deferred (scope = minimal study).

## Out of scope (future)
Connectome/connection-Laplacian directed spread; spreading-front tracker; Hodge source/flux/sink of
∂/∂t; eigenmode-dynamics; Brainstorm source-series + Dynamics-atoms integration; real longitudinal
data (ADNI/A4 or PREVENT-AD follow-up).
