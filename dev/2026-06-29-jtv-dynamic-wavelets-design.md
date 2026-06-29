# Dynamic wavelet atoms — joint time-cortex wavelets (GSPBox JTV, Brainstorm-native)

**Goal.** Turn the **atom** into a joint time-cortex **wavelet** filter. The atom already localizes on four
axes (`time`, `freq`, `source`, `scale`) with a uniform `center`/`extent`/`weighting` struct, where
`weighting='hard'` is an indicator mask (today) and `weighting='soft'` is documented as *"wavelet decay;
reserved, future"*. We fill that reserved `'soft'` case with **dynamic eigenfilters** so an atom becomes a
smooth, time-localized, cortex-localized wavelet — including **time-varying** wavelets (damped waves with a
physical speed via a dispersion relation) over the geometry **or** the connectome. No new axis machinery:
the **results file already IS the joint cortex-time axis** (`SurfaceFile`+`Time`+`DataFile`+`Fs`).

## Key realization (drives the whole design)
1. **The results/source architecture already binds cortex + time.** `db_template('resultsmat')` carries
   `SurfaceFile` (cortex), `Time` (`Fs = 1/(Time(2)-Time(1))`), and `DataFile`. That is the joint
   time-vertex axis, complete with file I/O. **No `bst_jtv` is needed** — `bst_dynamics` adopts this same
   architecture so atoms inherit the axes.
2. **The atom IS the joint filter.** `bst_atom` already abstracts all four axes as `(center, extent,
   weighting)`. Hard = indicator (current); **soft = wavelet** (this work). The atom's output — a
   time-localized, cortex-localized field with scale + frequency definition — is exactly what a dynamic
   wavelet produces.
3. **A filterbank = a bank of atoms** (the GSPBox bookkeeping, in our object).

## Mapping — GSPBox JTV ↔ Brainstorm (an expansion, not a rewrite)
| GSPBox | Brainstorm |
|---|---|
| graph eigenvalues `λ` | cortex `Lambda` (`tess_eigen`; LBO / **LB-Connectome**) — mm-calibrated `2π/√λ` |
| `G.jtv`: `T`,`fs`,`omega` | **results** `Time` + `Fs` → `omega` (Hz). No new struct. |
| `gsp_jft = U_G' X conj(U_T)` | `Φ'·B·F` (cortex) × `fft` (time) |
| `'ts'` `g(λ,t)` / `'js'` `g(λ,ω)` | **eigen-time** / **eigen-frequency** kernel domain |
| separable `g₁(λ)g₂(ω)` / non-separable | same |
| filterbank `g{m}` | **bank of atoms** |
| `gsp_jtv_design_*` | `toolbox/eigen/eigfilter/bst_eigfilter_design_*` (dynamic) |

## A. Axes — `bst_dynamics` adopts the results architecture (no `bst_jtv`)
`bst_dynamics` binds, like a results file, to:
- `SurfaceFile` → eigenbasis via `tess_eigen(variant)` → `Lambda` (CORTEX / eigen-freq axis, mm),
- `DataFile`/`ResultsFile` → `Time`,`Fs` → `omega` (TIME / temporal-freq axis, Hz).
These axes are *resolved on demand* from the bound files (exactly how a results file reconstructs its time
vector and references its cortex). The atom's `.scale [k1 k2]` and `.band [fLo fHi]` are then *windows on
these axes* — the atom's footprint in the joint `(λ, ω)` domain.

## B. The atom = the joint filter; `'soft'` weighting = wavelet
`bst_atom`'s `(center, extent, weighting)` per axis becomes the wavelet spec:
- **source axis**: `center`=seed vertex, `extent`=geodesic radius; `hard`=disk indicator → `soft`=spatial
  eigenfilter kernel on the cortex eigenbasis (heat/mexhat/damped-wave) centred at the seed.
- **time axis**: `center`/`extent`=window; `hard`=box → `soft`=temporal wavelet envelope.
- **scale axis**: `[k1 k2]` selects the eigen-band the spatial kernel lives in.
- **freq axis**: `[fLo fHi]` selects the temporal-frequency band.
**Separable** atom = soft-source ⊗ soft-time (independent). **Non-separable** atom = the source and time
axes COUPLE through a dispersion (damped-wave: the cortex scale `λ` sets the temporal oscillation) — this
is the genuinely *dynamic* wavelet. `bst_atom` gains a `Kernel`/`Evaluate` path that turns a soft
localization into the actual weighted field via the eigenfilter library.

## C. Dynamic eigenfilter library (`toolbox/eigen/eigfilter/`, Domain × Separable)
The static library already exists (`heat`, `mexhat`, `log`, `tikhonov`, `diffgauss`, `inverse_heat`,
`ideal`, `matern`, `power`, `flat`, + `kernel`/`compose`/`evaluate`). Add **dynamic** designs, each tagged
`Domain` (`eigen-time` ts `g(λ,t)` | `eigen-frequency` js `g(λ,ω)`) and `Separable`:
| design | kernel | Domain | Sep. | dynamic |
|---|---|---|---|---|
| `bst_eigfilter_design_diffusion` | `exp(-τ|t|λ/(λmax·Fs))` | ts | yes | physical heat (tau→seconds) |
| `bst_eigfilter_design_wave` | `cos(t·acos(1-α²λ/(2λmax)))` | ts | no | propagation, speed α |
| `bst_eigfilter_design_dampedwave` | `exp(-β|t|)·cos(α·t·acos(1-λ))` | ts | no | speed α + damping β |
| `bst_eigfilter_design_kleingordon` | `cos(t·acos(1-α²(λ+μ)/(2λmax)))` | ts | no | mass μ |
| `bst_eigfilter_design_meyerjtv` | `g₁(λ)·g₂(|ω|)` | js | yes | tight joint frame |
`bst_eigfilter_evaluate` gains the ts↔js conversion (`fft`/`ifft` along time, as `gsp_jtv_filter_evaluate`
does) so a kernel stored in either domain applies in the joint spectral domain.

## D. Orchestration — `bst_eigenfilter` / `bst_eigenwavelet`
The verb orchestrators gain JTV-aware Analysis/Synthesis that read the axes from the bound results/dynamics
structure and the kernel from the library:
```
W = bst_eigenwavelet('Analysis', F, EigenMat, OperatorMat, frame, JTVaxes)
   C    = Φ'·(B·F)                 % cortex FT          [K x nT]
   Chat = fft(C, NFFT, 2)          % time FT            [K x NFFT]
   G    = bst_eigfilter_evaluate(kernel, Lambda, omega) % ts→js if needed
   W    = Φ·real(ifft(Chat .* conj(G), NFFT, 2))        % inverse joint FT
```
The output is written back as a **bank of atoms** (one atom per band/scale, soft-weighted). For the
**whole-brain LB-Connectome** eigenbasis the per-hemisphere loop collapses to one block, so dispersion acts
over geometry **and** connectome.

## Physical units (the payoff)
`λ→mm` (`2π/√λ`), `ω→Hz` (`Fs`), and the wave parameter `α→m/s` (temporal-freq `∝ α√λ`, spatial-freq `∝
√λ` ⇒ speed `∝ α·Fs`). A cortical traveling-wave speed (~0.1–10 m/s) maps to an atom. Atoms thus carry
real spatial (mm), temporal (s), and speed (m/s) units.

## The capability this unlocks
A first-class, interactive field algebra on the cortex: **drop an atom on a vertex and propagate a damped
wave through the connectome**; analyze/synthesize/steer fields defined on geometry or connectome; select a
`(scale, band)` window and read the localized wavelet response — all as banks of atoms, visualizable and
composable.

## Reuse map
| Need | Exists | New |
|---|---|---|
| joint cortex-time axis + I/O | **results architecture** (`SurfaceFile`/`Time`/`DataFile`/`Fs`) | `bst_dynamics` adopts it |
| 4-axis localization + hard/soft hook | **`bst_atom`** (`weighting='soft'` reserved) | fill `'soft'` (wavelet) |
| eigenbasis | `tess_eigen` (LBO / whole-brain LB-Connectome) | — |
| static eigenfilters | full library in `eigfilter/` | + dynamic ts/js designs |
| verb orchestrators | `bst_eigenfilter`/`bst_eigenwavelet` | JTV Analysis/Synthesis + ts/js evaluate |
| filterbank bookkeeping | atom groups (`bst_dynamics`) | bank-of-atoms output |

## Plan
1. **`bst_dynamics` ← results architecture** — bind `SurfaceFile`/`DataFile`/operator-variant; resolve the
   `(Lambda, Time, Fs, omega)` axes on demand (mirrors how results reconstruct time + reference cortex).
2. **Dynamic eigenfilter library** — `diffusion` (physical-time heat, the bridge) + `wave`/`dampedwave`/
   `kleingordon`/`meyerjtv`; `bst_eigfilter_evaluate` ts↔js conversion. Tag Domain × Separable.
3. **`bst_atom` soft weighting** — `Kernel`/`Evaluate` path: a soft localization → the weighted field via
   the library (separable source⊗time; non-separable dispersion).
4. **Orchestrator JTV Analysis/Synthesis** — joint Fourier (`Φ'B` × `fft`), output as a bank of atoms;
   frame-bounds/Parseval check.
5. **Validate** — drop an atom on a posterior-occipital vertex, propagate `dampedwave` over the
   **LB-Connectome** eigenbasis, recover the alpha traveling-wave speed (m/s); confirm the static heat is
   the `ω`-independent limit.

*(Reaction-diffusion: OUT of scope for now — nonlinear, no GSPBox analog; add later as a forward-sim
dynamic atom.)*

## Open question
- `bst_atom` soft kernel storage: store the *design* (name+params) and evaluate lazily, vs cache the
  evaluated field. Recommend lazy (store design; evaluate against the resolved axes) — keeps atoms light
  and re-resolvable across surfaces/recordings.
