# Joint-spectral (js) eigfilter kernels — design

**Date:** 2026-06-30
**Status:** design (approved direction; pending spec review → plan)
**Depends on:** the Phase-1 JTV-atom architecture (`bst_eigenfilter('Atom')` js branch, `bst_eigen('Axes')`), shipped to origin/development `f30bdffe`.

---

## 1. Goal

Add the first **joint-spectral (`domain='js'`) kernels** to the eigfilter library — filters `g(λ,ω)` over (cortex eigenvalue × temporal frequency) — and wire them into the atom designer. This populates the js branch (today plumbed but unexercised by any real kernel) with four electrophysiology-motivated atoms, and adds the temporal-frequency (Hz) controls the designer was missing.

The four (approved): **gabor**, **travwave**, **resonator**, **stmatern**. No engine changes — they ride the existing `bst_eigenfilter('Atom')` js branch.

## 2. The realness convention (shared by all js kernels)

The `Atom` js branch evaluates `G = g(Λ, ω)` on the FFT grid `ω = (0:N-1)·Fs/N` (0..Fs, Hz, from `ax.omega`) and realises `W = manifold_ift(Phi, real(ifft(G,N,2))(:,1:nT) .* c0)`. For the atom to be real, `g` must be **conjugate-symmetric in the FFT sense**. Each js design file enforces this by folding the grid to **signed frequencies** and defining `g` as an even / Hermitian function of them:

```
Fs = numel(w) * (w(2)-w(1));            % infer Fs from the grid
ws = w - Fs .* (w >= Fs/2);            % 0..Fs  ->  signed [-Fs/2, Fs/2)
```

Even kernels (`gabor` mirror, `travwave |ws|`, `stmatern ws²`) and Hermitian kernels (`resonator` H(-ω)=conj H(ω)) then give `max|imag(atom)| ~ 1e-16`. Spatial frequency uses `k = √λ` (rad/m) with the existing `2π/√λ → mm` calibration.

Each design file follows the registry convention: return the `meta` struct on `'meta'`, else parse params (with `lmax` default) and return `out = @(l,w) i_eval(l, w, <params>)`, where `i_eval` is a local function doing the fold + formula (anonymous handles may call file-local functions).

## 3. The four kernels

`l` = λ column `[K×1]`, `w` = ω row `[1×N]` (Hz), `ws` = signed fold of `w`, `g → [K×N]`.

### 3.1 gabor — joint scale×frequency packet (bandpass wavelet)
```
g(λ,ω) = exp(-(√λ - k0)² / (2σk²)) .* ( exp(-(ws-f0)²/(2σf²)) + exp(-(ws+f0)²/(2σf²)) )
```
- `k0 = 2π/(Scale_mm/1000)` (rad/m); `σk = k0/3` (fixed relative spatial bandwidth).
- `f0` (Hz) temporal centre; `σf` (Hz) temporal bandwidth. The `±f0` pair makes it Hermitian → a real cosine packet at `f0`, spatially localised at scale `Scale`.
- meta: `bandpass=true`, `separable=true`, `priorAdmissible=false`. params: `f0, scaleMM, bwHz`.

### 3.2 travwave — dispersive traveling wave (non-separable, the showcase)
```
g(λ,ω) = exp( -(|ws| - c·√λ/(2π))² / (2σ²) )
```
- Ridge centre frequency per mode `fr(λ) = c√λ/(2π)` (Hz), `c` = phase speed (m/s); `σ` (Hz) ridge width. Even in `ws` → real. Couples temporal frequency to spatial scale.
- meta: `bandpass=true`, `separable=false`, `priorAdmissible=false`. params: `c, widthHz`.

### 3.3 resonator — damped harmonic oscillator (Lorentzian, bandpass)
```
g(λ,ω) = f0² / (f0² - ws² + 1i·ws·f0/Q)
```
- Hermitian (`g(-ws)=conj g(ws)`) → real atom = a decaying oscillation at `f0` (Hz) with quality `Q` (envelope time-constant ~ `Q/(π f0)`). λ-independent here (same resonance everywhere); a dispersive variant (`f0=c√λ`) is a future option.
- meta: `bandpass=true`, `separable=true`, `priorAdmissible=false`. params: `f0, Q`.

### 3.4 stmatern — spatiotemporal 1/f (Whittle–Matérn, low-pass / aperiodic)
```
g(λ,ω) = ( κ² + λ + (ws/v)² )^(-ν)
```
- Spectral density of a spatiotemporal Matérn field: the aperiodic 1/f background / a principled spatiotemporal prior. `κ = 2π/(Corr_mm/1000)` (rad/m) correlation length; `ν` smoothness; `v` (m/s) the space–time coupling speed — **fixed to a default** (`v = 1.0 m/s`), only `Corr` + `ν` exposed (per review). Even in `ws`, positive → real.
- meta: `bandpass=false`, `separable=false`, `priorAdmissible=true`. params: `corrMM, nu` (`v` internal default).

## 4. Designer control wiring (per-kernel, minimal)

Generalise today's morphing (which relabels only slider 1) to **all three sliders** via one per-kernel spec covering existing + new kernels:

```
cfg = i_kernel_sliders(kernel)   % 1x3 struct array: {label, unit, min, max, val, enable}
```

`SyncControls` applies `cfg` to the three sliders (relabel + re-range + enable/disable); `i_phys2kernel` maps the three slider values → kernel params per kernel. This subsumes the current `i_spatial_mode` (slider 1) and the fixed Speed/Decay handling. Existing kernels keep their present mappings, expressed in the same table. New js mappings:

| kernel | slider 1 | slider 2 | slider 3 |
|---|---|---|---|
| gabor | Scale (mm) | Freq (Hz, 0–50) | Bandwidth (Hz) |
| travwave | Speed (m/s) | Ridge width (Hz) | — |
| resonator | Freq (Hz, 0–50) | Q | — |
| stmatern | Corr (mm) | Smoothness ν | — |

Frequency sliders span **0–50 Hz** (Nyquist of the 1 s @ 100 Hz axis; covers δ/θ/α/β/low-γ). js kernels declare a `domain`, so they appear in the dropdown's **dynamic** group automatically. Spatial ranges (`Scale`/`Corr`) derive from the cortex spectrum (`scaleMin/MaxMM`).

## 5. Testing (TDD) — finally exercises the js branch end-to-end

Per kernel, on a synthetic `ax`, realise via `bst_eigenfilter('Atom')` and assert:
1. shape `[nV × nT]`;
2. **real**: `max|imag| < 1e-9` (the realness convention holds);
3. the defining property:
   - **gabor** — temporal power spectrum of the seed-row peaks at `f0` (± one bin);
   - **travwave** — peak temporal frequency increases with `√λ` (compare two `λ` bands);
   - **resonator** — seed-row oscillates at `f0`, and a higher `Q` yields a **longer-lasting envelope** than a lower `Q` (robust monotonic comparison of two atoms, not a brittle exact decay match);
   - **stmatern** — atom is positive (low-pass), temporal PSD is monotone-decreasing (1/f slope).
4. designer wiring — `i_kernel_sliders('gabor')` returns 3 labelled sliders; `i_phys2kernel` (driven by the slider state) yields the expected param struct; each new kernel is in `bst_eigfilter_kernel('list')` and in the designer's dynamic group.

Plus an end-to-end smoke on the real cortex (`sub-MTL0002`, Laplace-Beltrami): each js kernel realises a finite, real atom.

## 6. Files

- Create: `toolbox/eigen/eigfilter/bst_eigfilter_design_{gabor,travwave,resonator,stmatern}.m` (auto-registered by the `dir` discovery).
- Modify: `toolbox/gui/view_atom_designer.m` (per-kernel `i_kernel_sliders` spec + `i_phys2kernel` cases; generalise slider morphing to all three; existing kernels folded into the table unchanged).
- Tests: `dev/tests/test_js_kernels.m`, `dev/tests/test_designer_controls.m`.

## 7. Out of scope

- The **meta-driven** generic slider architecture (kernels declare params, panel renders them) — deferred; this phase keeps the per-kernel hand-wired mapping.
- Dispersive resonator (`f0 = c√λ`), chirplet, transient/edge kernels — future kernels; the design/test pattern here generalises to them.
- A faster time axis for high-γ (`f0 > 50 Hz`) — the 1 s @ 100 Hz axis is fixed this phase.
