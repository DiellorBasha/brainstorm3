# Joint time-vertex (time-cortex) dynamic wavelets — design

**Goal.** Extend the eigenfilter/eigenwavelet system with JOINT TIME-VERTEX dynamic wavelets à la GSPBox:
kernels `g(λ, ω)` over BOTH the cortical eigen-spectrum (`λ`, from `tess_eigen`; LBO or LB-Connectome)
and the recording's temporal frequency (`ω`, Hz, from the data `Time`/`Fs`). This turns the static
wavelets (heat/mexhat/itersine) into a *registry* that also holds time-varying wavelets — damped waves
with a physical speed via a dispersion relation, reaction-diffusion, Klein-Gordon — so neurophysiological
signals can be analysed as they actually **propagate over geometry+connectome with physical dynamics**.
The canonical `(cortex, time)` axis source-of-truth is the **atoms** system (already half-built for it).

## The mapping — GSPBox JTV ↔ Brainstorm (near-exact; this is an expansion, not a rewrite)
| GSPBox | Brainstorm equivalent |
|---|---|
| graph eigenvalues `λ` (`gsp_compute_fourier_basis`) | cortex `Lambda` (`tess_eigen`; LBO / **LB-Connectome**) — mm-calibrated `2π/√λ` |
| `G.jtv` : `T`, `fs`, `omega`, `NFFT` | recording `.Time` + `Fs = 1/dt` → `omega` (FFT grid, Hz) |
| `gsp_jft = U_G' X conj(U_T)` | `manifold_ft` (`Φ'·B·F`) along cortex × `fft` along time |
| filtertype `'ts'` `g(λ,t)` / `'js'` `g(λ,ω)` | **eigen-time** / **eigen-frequency** (the user's distinction) |
| separable `g₁(λ)g₂(ω)` / non-separable | same |
| `gsp_jtv_filter_analysis` (joint-domain multiply) | new `bst_eigenwavelet('JTVAnalysis', …)` |
| `gsp_jtv_design_*` (damped wave, …) | new `bst_eigfilter_design_*` dynamic kernels |

Our `bst_eigenfilter`/`bst_eigenwavelet` already expose `Design/Evaluate/Bounds/Analysis/Synthesis`
verbs (GSPBox-style), so the joint (time) axis slots straight in.

## A. Canonical JTV axis — the atoms expansion (source of truth)
`db_template('atomgroup')` ALREADY references the entire joint structure:
- `.SurfaceFile` → the cortex (**space axis**)
- `.DataFile` / `.ResultsFile` → the recording (**time axis**: `.Time`, `Fs = 1/(Time(2)-Time(1))`)
- `.scale [k1 k2]` + `.scaleName` → **eigen-spectrum** sub-band ("reserved for bst_eigen")
- `.band [fLo fHi]` + `.bandName` → **temporal-frequency** sub-band

So the atom is already the canonical `(cortex, time, eigen-scale, temporal-band)` anchor. The expansion is
a resolver turning an atom (or a results file + an operator variant) into a concrete JTV struct — the
direct analog of `gsp_jtv_graph(G, T, fs)`, but with the `(cortex, time)` sources taken from the atom:
```
JTV = bst_jtv('Resolve', atomOrResultsFile, OperatorVariant)
   .EigenMat  % tess_eigen(.SurfaceFile, OperatorVariant) -> Phi, Lambda   (CORTEX axis)
   .Time .Fs  % from .DataFile/.ResultsFile .Time                          (TIME axis)
   .omega     % FFT grid (Hz) from numel(Time), Fs                         (TEMPORAL-FREQ axis)
   .lambda    % Lambda, with mm-calibration 2*pi/sqrt(lambda)              (EIGEN-FREQ axis)
   .Atom      % optional; .scale/.band define the joint FOOTPRINT in (lambda,omega)
```
A "dynamic atom" is then literally a wavelet **localized in both cortex-scale and temporal-band** — the
`.scale`×`.band` rectangle of the joint domain.

## B. The dynamic kernel registry (Domain × Separable taxonomy)
Extend `toolbox/eigen/eigfilter/` with a REGISTRY; every kernel is tagged on the two axes the user named:
- **Domain**: `eigen-frequency` (js, `g(λ,ω)`) | `eigen-time` (ts, `g(λ,t)`)
- **Separable**: `true` (`g₁(λ)·g₂(·)`) | `false` (dispersion-coupled)
```
bst_eigfilter_registry entry := { Name, Domain, Separable, Params, Kernel(lambda,·,params,JTV),
                                  Dispersion, Units, Refs }
```
**Static** (current) kernels → `Domain='eigen-frequency'`, ω-independent (the time axis is trivial):
- `heat` `exp(-tλ)`, `mexhat` `(tλ)exp(-tλ)`, `itersine` (tight half-cosine).

**Dynamic** kernels (ported from GSPBox `gsp_jtv_design_*`, eigenbasis-adapted; `λmax`=spectrum max):
| Name | kernel | Domain | Separable | dynamic |
|---|---|---|---|---|
| `diffusion` (physical heat) | `exp(-τ|t|λ/(λmax·Fs))` | ts | yes | the tau→seconds bridge |
| `wave` | `cos(t·acos(1-α²λ/(2λmax)))` | ts | no | propagation, speed α |
| `damped-wave` | `exp(-β|t|)·cos(α·t·acos(1-λ))` | ts | no | speed α + damping β (dispersion) |
| `klein-gordon` | `cos(t·acos(1-α²(λ+μ)/(2λmax)))` | ts | no | mass μ |
| `meyer-jtv` | `g₁(λ)·g₂(|ω|)` | js | yes | tight joint frame |
| `reaction-diffusion` | linearized Turing dispersion (open q.) | ts | no | growth+diffusion |

## C. Physical-time heat (the tau→seconds bridge)
Today `tau` is an abstract spectral diffusion scale (`exp(-tλ)`). With the JTV `Fs`, the `diffusion`
kernel `exp(-τ|t|λ/(λmax·Fs))` makes `τ` a PHYSICAL diffusivity acting over real time `|t|` (seconds).
Combined with the mm-calibration (`λ→2π/√λ`), the dynamic kernels acquire **physical units**: spatial
scale (mm), temporal scale (s), and — for `wave`/`damped-wave` — a **wave speed in m/s** (the α
parameter, since temporal-freq `∝ α√λ` and spatial-freq `∝ √λ` ⇒ speed `∝ α·Fs`). A cortical
traveling-wave speed (≈0.1–10 m/s) maps directly to α. This is the headline upgrade: physically
calibrated traveling-wave wavelets on the cortex+connectome.

## D. JTV Analysis/Synthesis verbs (mirror `gsp_jtv_filter_analysis`)
```
W = bst_eigenwavelet('JTVAnalysis', F, JTV, frame)   % F = [nV x nT] time-vertex signal
   per hemisphere h:
     C    = Phi{h}'·(B{h}·F(gv,:))     % graph FT                          [K x nT]
     Chat = fft(C, NFFT, 2)            % temporal FFT  -> joint spectrum    [K x NFFT]
     for each kernel g_m:
        G  = Evaluate(g_m, Lambda{h}, omega)   % ts->js via fft-along-time if Domain='ts'  [K x NFFT]
        W_m(gv,:) = Phi{h}·real(ifft(Chat .* conj(G), NFFT, 2))             % inverse joint FT [|gv| x nT]
```
`Synthesis` is the dual (sum over m, no conjugate). `'ts'`↔`'js'` conversion lives in `Evaluate`
(`fft`/`ifft` along the time axis), exactly as `gsp_jtv_filter_evaluate` does. For the **whole-brain**
LB-Connectome eigenbasis the hemisphere loop collapses to one block (`GlobalVertices{2}=[]`).

## Reuse map (minimal new code)
| Need | Exists | New |
|---|---|---|
| eigenbasis Phi,Lambda | `tess_eigen` (incl. whole-brain LB-Connectome) | — |
| recording Time,Fs | `in_bst_results`/`.Time`, `bst_memory('GetTimeVector')` | — |
| (cortex,time) anchor | atom `.SurfaceFile`/`.DataFile`/`.scale`/`.band` | `bst_jtv('Resolve')` |
| verb orchestrator | `bst_eigenwavelet`/`bst_eigenfilter` (Design/Evaluate/Analysis/Synthesis) | JTV verbs |
| static kernels | heat/mexhat/itersine | Domain/Separable tags |
| dynamic kernels | — | diffusion/wave/damped-wave/klein-gordon/meyer-jtv |
| joint transform | `manifold_ft` (`Φ'B`) | + temporal `fft` wrapper |

## Phased plan
1. **JTV axis resolver `bst_jtv('Resolve')`** — atom/results + operator variant → `{EigenMat, Time, Fs,
   omega, lambda}`. The canonical source of truth. (Smallest, unblocks everything.)
2. **Kernel registry + physical-time heat** — the registry schema (Domain/Separable) +
   `bst_eigfilter_design_diffusion` (tau→seconds); validate the static heat is the `ω`-independent limit.
3. **JTV Analysis/Synthesis** — the joint apply path in `bst_eigenwavelet`/`bst_eigenfilter`; ts/js
   handling; reconstruction/Parseval (frame-bounds) check.
4. **Dynamic wavelet library** — `wave`, `damped-wave`, `klein-gordon`, `meyer-jtv`; validate on the
   posterior-alpha traveling-wave segment (does `damped-wave` localize the propagation + recover its
   speed in m/s?), on the LB-Connectome eigenbasis (geometry+connectome).
5. **Reaction-diffusion wavelet** — linearized/Turing dispersion form; tie to the PET reaction-diffusion
   work.

## Open questions
- **Reaction-diffusion is NONLINEAR** — no direct GSPBox analog. Options: (a) a *linearized* (Turing)
  dispersion kernel for the registry; (b) a forward-sim "dynamic atom" (a propagator, not a spectral
  filter). Recommend (a) in the registry + (b) as a separate dynamics tool.
- **Whole-brain**: using the LB-Connectome eigenbasis as the cortex axis makes the dispersion act over
  geometry+connectome — the unification (signals propagate locally *and* via fibers, with a speed).
- **Atom footprint**: should `bst_jtv` enforce the atom's `.scale`×`.band` as a hard mask, or expose it
  as a default analysis window? (Recommend: default window, overridable.)
- **DCT vs DFT** temporal transform (GSPBox supports both via `G.jtv.transform`) — start with DFT.
