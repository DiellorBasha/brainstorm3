# Eigen / JTV-atom architecture — design

**Date:** 2026-06-30
**Status:** design (approved direction; pending spec review → plan)
**Topic:** organise the `eigen` module and retire `dynamics/bst_atom` so single **atoms** (the atom designer's output) and **joint time-vertex (JTV)** wavelets share one domain-aware engine, named by a consistent vocabulary.

---

## 1. Motivation

Three partially-overlapping implementations of "propagate a kernel on an eigenbasis" exist:

- `bst_eigenwavelet` — static spectral wavelet **frames** (Design/Evaluate/Bounds/Analysis/Synthesis) *plus* a bolted-on JTV path (`JTVAnalysis`/`JTVAtoms`), an `Atom` verb, and the Dirac/quaternion vector specialisation (Steer/ToVec/ToQuat).
- `bst_eigenfilter` — a simpler **static** filter engine (Design/Evaluate/Analysis/Synthesis), *no* `Atom`, *no* time/frequency axis.
- `dynamics/bst_atom` — a marker **localiser** (Axes/Get/Set on time/freq/source/scale) fused with a *ts-only* single-atom realiser (`Evaluate`) and a level-set bridge (`Levelset`).

The atom designer (`view_atom_designer`) calls `bst_atom('Evaluate')`, which is ts-only: it cannot realise a `js` (frequency-eigenspectrum) kernel, and the terminology ("jtv" / "dynamic" / "atom") is used loosely. We want a clean, vocabulary-aligned architecture in which the designer produces **atoms** through a single **domain-aware** realiser, and the dynamics layer is purely markers + GUI.

## 2. Vocabulary (the naming contract)

These terms are used precisely throughout, and each names exactly one tier:

| term | meaning |
|---|---|
| **filter** / **kernel** | one spectral shape `g(λ)` (static) or `g(λ,t)`/`g(λ,ω)` (JTV). The `bst_eigfilter_*` library. |
| **atom** | one **localised** filter — `g(L)·δ_seed` — a single vector/field. *This is what the atom designer makes.* |
| **filterbank** | a *set* of filters (several shapes/scales). |
| **frame** | all atoms from a filterbank, with bounds `A,B` → stable, invertible Analysis↔Synthesis. |
| **wavelet** | a *band-pass, doubly-localised, dilation-generated* filter/atom (admissibility `g(0)=0`). A **special kind** of filter, marked in the registry by `bandpass=true`. Wavelet ⊂ filter. |
| **static vs JTV** | a *domain* qualifier, orthogonal to the tiers: static `g(λ)`, JTV `g(λ,t)` (ts) / `g(λ,ω)` (js). |
| **marker** | the stored dynamics reference (Scout + Event), *derived from* an atom via level-set. Lives in `bst_dynamics`. |

Consequences that drive the design:
- "atom" is a **verb/output**, never a module name → `bst_atom` is misnamed and is retired.
- wavelet ⊂ filter → `bst_eigenwavelet` is a **specialisation** of `bst_eigenfilter`, not a parallel engine.
- domain (static/JTV) is handled **inside** the engine via an axes struct, not by separate modules or `JTV*` verbs.

## 3. Target architecture

```
eigen/
  bst_eigen.m            ORCHESTRATOR + file I/O + AXES source        [+ 'Axes' verb]
  bst_eigenfilter.m      GENERAL engine (filter/atom/filterbank/frame), domain-aware
                           target verbs: Design · Evaluate · Atom · Analysis · Synthesis · Bounds
                           (today: Design/Evaluate/Analysis/Synthesis static-only; Phase 1 adds Atom;
                            Bounds + frame consolidation are a later phase — see §6)
  bst_eigenwavelet.m     WAVELET specialisation on bst_eigenfilter
                           band-pass dilation families, mm/Hz scale calibration,
                           tight-frame design (itersine), Dirac/quaternion vector wavelets
  bst_eigenspectrum.m    (unchanged)
  eigfilter/             FILTER LIBRARY (shared substrate)            [keep]
                           bst_eigfilter_kernel (registry: domain, bandpass, separable)
                           bst_eigfilter_design_* · bst_eigfilter_evaluate
                           bst_eigfilter_jtv_evaluate (ts↔js bridge)
math/
  manifold_ft / ift      cortex (graph) Fourier transform            [keep]
  manifold_jft / ijft    joint time-vertex Fourier transform         [keep]
dynamics/
  bst_dynamics.m         MARKERS: atomgroup data-model + Levelset (atom→Scout/Event)
                           + the marker-localisation accessor (absorbed from bst_atom)
  panel_bst_dynamics.m   GUI; eventually absorbs the atom designer
  (bst_atom.m            RETIRED)
gui/
  view_atom_designer.m   atom designer → calls bst_eigenfilter('Atom'); later folds into panel
```

**`bst_eigen` is the axes source** (no new `bst_eigen_axes` file). `bst_eigen` already owns eigen-file resolution and loading (`GetEigenBasis` → `in_bst_eigen`/`in_bst_operator`); it is the natural place to assemble the canonical axes struct, exposed as a verb:

```
ax = bst_eigen('Axes', SurfaceFile, OPTIONS)
```

### 3.1 The canonical axes struct (the keystone)

Built once by `bst_eigen('Axes', …)` and consumed by every realiser/analyser. Replaces the three hand-built `ax` structs scattered across the designer, `bst_dynamics`, and `JTVAnalysis` (none of which carry `ω` today).

```
ax = struct(
  % --- spatial (from GetEigenBasis) ---
  'Phi',{...}, 'Lambda',{...}, 'Mass',{...}, 'GlobalVertices',{...},  % per-hemisphere cells
  'Variant', '...', 'nModes', K,
  % --- temporal ---
  'tlag',  [1 x nT],   'nT', nT, 'dt', dt,        % time-lag axis (s)
  % --- spectral (temporal frequency) ---
  'omega', [1 x NFFT], 'NFFT', NFFT, 'df', df     % temporal-frequency axis (Hz)
)
```

Consistency invariant: `omega = fft frequencies(nT, 1/dt)` so ts↔js round-trips are exact (Parseval via the `1/sqrt(NFFT)` normalisation already in `manifold_jft`).

`OPTIONS` extends the existing eigen OPTIONS (`Variant`, `nModes`, `Tau`, `EigenFile`) with temporal fields (`TimeWindow`, `SampleRate`); spatial-only callers omit them and get a spatial-only axes (no `tlag`/`omega`).

### 3.2 The domain-aware `Atom` realiser

`bst_eigenfilter('Atom', …)` is the single realiser the designer (and `bst_dynamics`) call. It dispatches on the kernel's `domain`:

```
[W, gv] = bst_eigenfilter('Atom', ax, KernelName, KernelParams, seedVert)

c0 = manifold_ft(Phi, Mass, δ_seed);                       % seed in the eigenbasis [K x 1]
g  = bst_eigfilter_kernel(KernelName, KernelParams);
switch bst_eigfilter_kernel('info',KernelName).domain
  case 'static' : W = repmat(manifold_ift(Phi, g(Lam).*c0), 1, nT);          % const in time
  case 'ts'     : W = manifold_ift(Phi, g(Lam, ax.tlag) .* c0);              % today's designer path
  case 'js'     : Gjs = bst_eigfilter_jtv_evaluate(g,'js',Lam,ax.omega,NFFT);% NEW
                  W   = real(manifold_ijft(Phi, (c0*ones(1,NFFT)) .* Gjs, nT));
end
W :: [nV x nT]
```

This reuses the same **domain dispatch** `bst_eigenwavelet.JTVAnalysis` performs, lifted into the single-atom path so the designer inherits full ts/js support. Note the direction differs: `JTVAnalysis` *analyses* a signal (inner product → `conj(Gjs)`), whereas `Atom` *realises* the kernel's own response to a seed (synthesis → `Gjs`, no conjugate). (`static` is the absent-`domain` case; `bandpass` is irrelevant to realisation — it only marks whether the filter is a *wavelet*.)

### 3.3 `bst_atom` retirement map

```
bst_atom('Evaluate', G, occ, ax, kernel, kp)  →  bst_eigenfilter('Atom', ax, kernel, kp, seed)
bst_atom('Levelset', W, gv, thr)              →  bst_dynamics('Levelset', …)   (atom → Scout/Event)
bst_atom('Axes'/'Get'/'Set' localisation)     →  bst_dynamics  (marker-localisation accessor)
```

The marker **data-model** (`atomgroup` template: `times`/`band`/`vertices`/`scale`) is unchanged and stays in `bst_dynamics`; only the *accessor API* moves out of `bst_atom`. An atom's localisation becomes a **derived** quantity (compute via `Levelset`), not a hand-set property — Brainstorm's native time-stepping / scouts / frequency tools consume the resulting markers.

## 4. Decided design questions

- **Two engines, not one.** `bst_eigenfilter` = general (any kernel, incl. low-pass/propagator); `bst_eigenwavelet` = band-pass wavelet specialisation that *uses* it. Mirrors GSPBox (`gsp_filter_*` vs the wavelet designers) and keeps `bst_eigenwavelet` honestly about wavelets.
- **`Design` (dynamic) = one canonical atom now.** A *designed JTV frame* (tight filterbank over `(λ,ω)` with bounds) is deferred — see §6.

## 5. Scope of THIS spec (Phase 1 — single implementation plan)

The plannable, no-regression chunk:

1. **`bst_eigen('Axes', SurfaceFile, OPTIONS)`** — assemble the canonical axes struct (spatial from `GetEigenBasis`, temporal `tlag`, spectral `omega`). Spatial-only when temporal OPTIONS absent.
2. **`bst_eigenfilter('Atom', …)`** — add the domain-aware single-atom realiser (static/ts/js) per §3.2; reuse `bst_eigfilter_jtv_evaluate` + `manifold_(j)ift`.
3. **Retire `bst_atom`** per §3.3: `Evaluate`→`bst_eigenfilter('Atom')`; `Levelset` + localiser → `bst_dynamics`; delete `dynamics/bst_atom.m`.
4. **Point `view_atom_designer` at `bst_eigenfilter('Atom')`** via `bst_eigen('Axes')`; remove its inline `i_eval_atom`/`ax` assembly. Behaviour preserved for static + ts; **js kernels now realisable** in the designer.

**Success criteria:** the designer renders identical output for existing kernels (diffusion/wave/heat/mexhat); a `js` kernel (when one exists) realises without error; `bst_atom.m` is gone with no dangling callers; no change to the marker DB schema.

## 6. Out of scope (documented, not specced here)

- **Later phase — engine consolidation:** move the general JTV / `Bounds` / frame `Synthesis` machinery *down* from `bst_eigenwavelet` into `bst_eigenfilter`, leaving `bst_eigenwavelet` as the band-pass + Dirac specialisation. (Today JTV lives in `bst_eigenwavelet`; Phase 1 only *reads* `jtv_evaluate`, it does not relocate code.)
- **Later phase — GUI fold-in:** absorb `view_atom_designer` into `panel_bst_dynamics` as its Design front-end.
- **Research frontier — JTV tight frame:** `bst_eigenfilter('Bounds'/'Synthesis')` for *joint time-vertex* wavelets (designing a tight filterbank over `(λ,ω)`, not just filtering). A genuine band-pass-in-both-planes JTV wavelet family does not yet exist (`diffusion`/`wave` are low/all-pass in `λ`). Explicitly deferred; must not block Phase 1.

## 7. Risks / notes

- `bst_eigen`'s entry is `bst_eigen(Data, OPTIONS)`, not macro-dispatch; the `'Axes'` verb needs a leading-string branch (small refactor of the entry, not a new file).
- `tess_callbacks`/CRLF and the throwaway-results-file designer pattern are unaffected.
- Keep the shared substrate single-sourced: both engines and the designer pull kernels from `bst_eigfilter_kernel` and transforms from `manifold_(j)ft/(j)ift` — no fourth copy of the propagation math.
