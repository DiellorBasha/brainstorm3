# bst_eigfilter: a spectral-filter library for LBO eigenmodes

**Date:** 2026-06-02
**Author:** Diellor Basha
**Status:** Design approved (v1 scope), pending implementation plan

## Motivation

LBO cortical eigenmodes diagonalize any function of the Laplacian, so every
isotropic spatial operator is a per-mode multiplier `g(λ_k)` on the eigenvalues.
That single object — a spectral filter `g(λ)` — currently appears in **two
separate, non-shared switch statements**:

- `toolbox/math/bst_eigenmode_prior.m` — `{flat, power, log}`, hardcoded, used as
  the source-covariance prior `R` in the eigenmode inverse.
- `toolbox/math/bst_eigenmodes_filter_gain.m` — `{lowpass, highpass, bandpass,
  heat, inverse_heat, tikhonov, custom}`, a centralized gain switch consumed by
  `bst_eigenmodes_filter` (spatial smoothing) and `process_eigenmodes_coeffsfilter`
  (coefficient filtering).

The heat kernel, power law, etc. are therefore at risk of being defined twice.
This work introduces a single **named-filter library**, `bst_eigfilter`, that
both the prior step and the analysis/coefficient filters draw from. Each filter
is a standalone analytic factory returning a function handle `g = @(λ) …`,
following the Graph Signal Processing Toolbox (GSPBox) `filters/` module layout
(design / evaluate / compose, with a reserved joint time-vertex family).

The name pairs deliberately with the existing `bst_freqfilter` ("apply a
frequency filter to a signal"): `bst_freqfilter` filters along the temporal-
frequency axis; `bst_eigfilter` filters along the eigenvalue (spatial-frequency)
axis.

## Scope (v1)

In scope:
- The `bst_eigfilter` design library under `toolbox/math/eigfilter/`: one factory
  file per kernel, an auto-discovering dispatcher, an evaluator, and a composer.
- Port the existing kernels (`flat, power, log, heat, inverse_heat, tikhonov,
  ideal` masks) onto the library and add a principled core (`matern, mexhat, dog`).
- Unify `bst_eigenmode_prior` and `bst_eigenmodes_filter_gain` so they both
  delegate to the library (single source of truth).
- Filterbank support: a vector-valued scale parameter yields a cell array of
  handles (design level only).
- Pure tests, including numerical parity guards proving the refactor changed no
  behavior.

Out of scope (documented future work):
- New GUI options exposing the new kernels (no changes to process radios yet).
- A filter-design GUI (live `g(λ)` plot + on-surface vertex delta response).
- The joint time-vertex family `bst_eigfilter_jtv_design_*` (damped-wave / neural
  field, diffusion, wave, Klein-Gordon).
- Neural-field spectrum kernel.

## Design decisions (resolved in brainstorming)

1. **Layout:** one file per kernel, GSPBox-style, under `toolbox/math/eigfilter/`.
2. **Naming:** `bst_eigfilter_*` prefix with GSP verbs (`_design_<name>` factories,
   `_evaluate`, `_compose`); pairs with `bst_freqfilter`. A separate
   `bst_eigfilter_kernel` registry provides string→factory dispatch + `list`/`info`
   (GSPBox has no such dispatcher; factories there are called directly by name).
3. **Unify both** existing `g(λ)` switches onto the library.
4. **Kernels are prior-agnostic.** The prior owns admissibility, not the kernel.
5. **Filterbanks supported** (vector scale → cell of handles); the prior consumes
   single-filter kernels only.
6. **Reuse existing apply-operators.** Because we have explicit eigenvectors Φ,
   filters are applied exactly as `Φ·diag(g(λ))·Φ'·M` (already implemented in
   `bst_eigenmodes_filter`/`_project`/reconstruct). No Chebyshev approximation /
   frame machinery (unlike GSPBox, which avoids computing eigenvectors).

## Architecture

```
toolbox/math/eigfilter/
  Contents.m                      % module index (mirrors GSPBox Contents.m)
  bst_eigfilter_kernel.m          % registry: name string -> factory + list/info (auto-discovers)
  bst_eigfilter_design_flat.m     % g = @(l) ones(size(l))
  bst_eigfilter_design_power.m    % g = @(l) l.^(-alpha)
  bst_eigfilter_design_log.m      % g = @(l) -log(l)      (raw; prior owns scaling)
  bst_eigfilter_design_heat.m     % g = @(l) exp(-t*l)
  bst_eigfilter_design_inverse_heat.m  % g = @(l) min(exp(+t*l), MaxGain)
  bst_eigfilter_design_tikhonov.m % g = @(l) 1./(1+beta*l)
  bst_eigfilter_design_ideal.m    % g = @(l) indicator on a lambda interval
  bst_eigfilter_design_matern.m   % g = @(l) (kappa^2 + l).^(-nu)
  bst_eigfilter_design_mexhat.m   % g = @(l) (t*l).*exp(-t*l)        (band-pass)
  bst_eigfilter_design_dog.m      % g = @(l) exp(-t1*l) - exp(-t2*l) (band-pass)
  bst_eigfilter_evaluate.m        % h = evaluate(g, lambdas)  (handles cell banks)
  bst_eigfilter_compose.m         % g = compose(g1, g2, ...)  (product = serial)
```

### Factory contract (`bst_eigfilter_design_<name>`)

Factories are GSP-faithful: each is a standalone function, called **directly** by
name when the caller knows the kernel (`g = bst_eigfilter_design_heat(params)`),
or via the registry when the kernel arrives as a string (see below).

```matlab
function out = bst_eigfilter_design_<name>(params)
%   out = bst_eigfilter_design_<name>()           -> handle with defaults
%   out = bst_eigfilter_design_<name>(params)     -> handle (or cell of handles)
%   out = bst_eigfilter_design_<name>('meta')     -> metadata struct
```
- `params` is a struct (or name/value list) of the kernel's hyperparameters.
- Returns a function handle `g = @(lambda) ...` evaluating the analytic shape.
- If a scale parameter is a **vector**, returns a **cell array** of handles
  (a filterbank), exactly as `gsp_design_heat`/`gsp_design_mexican_hat` do.
- `'meta'` returns `struct('name', 'heat', 'display', 'Heat / diffusion',
  'params', struct('t', struct('default', 0.01, 'range', [0 Inf])), 'bandpass',
  false, 'priorAdmissible', true)`. The `priorAdmissible` field is advisory
  metadata for UIs; the prior still validates numerically (see below).
- Optional `params.lmax`: if present, the kernel normalizes the spectrum (e.g.
  `exp(-t*l/lmax)`) so scale params are dimensionless and portable across
  surfaces, matching GSPBox. Absent → raw eigenvalues.

### Registry (`bst_eigfilter_kernel`) — separate from the factories

A thin, distinctly-named helper for the Brainstorm-specific need that GSPBox does
not have: turning a kernel **name string** (from a process option or the future
GUI) into a factory call, plus programmatic discovery. The factories themselves
stay pure GSP-style; this just routes to them.

```matlab
%   g     = bst_eigfilter_kernel(name, params)    -> handle via the named factory
%   names = bst_eigfilter_kernel('list')          -> available kernel names
%   meta  = bst_eigfilter_kernel('info', name)    -> metadata for one kernel
```
- `bst_eigfilter_kernel(name, params)` dispatches via
  `feval(['bst_eigfilter_design_' name], params)`. Unknown name → clear error.
- `'list'` scans `bst_eigfilter_design_*.m` in the module folder and strips the
  prefix → kernel names. Adding a kernel needs no registry edit (auto-discovery).
- `'info'`/`name` returns that factory's `'meta'` struct.

(GSPBox has no equivalent — there `gsp_design_heat` etc. are called directly by
name; there is no `gsp_design` dispatcher. The registry is an explicit Brainstorm
add-on for string-keyed callers, named separately so it does not visually shadow
the factories.)

### Evaluate and compose

- `bst_eigfilter_evaluate(g, lambdas)` → `[K x 1]` for a single handle, `[K x Nf]`
  for a cell-array bank. Thin (`g(lambdas(:))`), but centralizes the cell-vs-handle
  handling so consumers do not each reimplement it.
- `bst_eigfilter_compose(g1, g2, ...)` → `@(l) g1(l).*g2(l).*...` (product =
  serial filtering; mirrors `gsp_multiply_filters`). Single handles only.

## Consumers after unification

### `bst_eigenmode_prior.m` (owns admissibility)

Keeps its public signature `R = bst_eigenmode_prior(lambdas, K, priorType, alpha)`.
Internally:
1. Handles its own conventions: DC-mode swap, the `log` millimetre rescaling
   (`λ_mm = λ·1e-6`), and final `max(R)=1` normalization — these stay in the
   prior, not the kernel.
2. Builds the kernel from the (string) prior type: `g = bst_eigfilter_kernel(priorType, params)`.
3. Evaluates: `R = bst_eigfilter_evaluate(g, lambdas_scaled)`.
4. **Admissibility check:** a prior is a covariance, so `R` must be finite and
   `>= 0`. If a caller selects a kernel whose evaluated gain goes negative (e.g.
   `dog`), error clearly: `"Filter '<name>' is not admissible as a prior
   (negative spectral density). Use a non-negative kernel."` Banks (cell arrays)
   are rejected for the prior. This keeps kernels prior-agnostic while the prior
   enforces what it needs.

`power`/`log`/`flat` map to the same-named factories. The numeric result for the
three current prior types must be identical to today (parity test).

### `bst_eigenmodes_filter_gain.m` (delegates)

Keeps its signature `h = bst_eigenmodes_filter_gain(lambdas, FilterType, varargin)`
and its option names, so `bst_eigenmodes_filter` and `process_eigenmodes_coeffsfilter`
are untouched. Internally:
- Analytic λ-kernels (`heat, inverse_heat, tikhonov`) → delegate to
  `bst_eigfilter_kernel` + `bst_eigfilter_evaluate`. These are exact λ-functions,
  so the delegation is numerically identical.
- Index masks (`lowpass, highpass, bandpass` keyed on `CutoffMode`/`ModeRange`)
  stay **index-based, computed in place** — they select by mode *index*, not by a
  λ threshold, so routing them through the λ-based `ideal` kernel could drift at
  ties/off-by-one. Keeping them in place guarantees exact parity.
- `custom` (caller-supplied handle) → unchanged passthrough.
The `ideal` library kernel still exists for λ-threshold use elsewhere/future, but
`filter_gain` does not use it for its index masks. The existing
`test_eigenmodes_filter_gain_pure.m` must still pass unchanged.

## Kernel set (v1)

| Name | `g(λ)` | Role | Prior-admissible |
|------|--------|------|------------------|
| `flat` | `1` | baseline | yes |
| `power` | `λ^{-α}` | smoothness / 1-f | yes |
| `log` | `-log(λ)` (prior rescales) | GBF 2026 prior | yes |
| `heat` | `e^{-tλ}` | low-pass / diffusion | yes |
| `inverse_heat` | `min(e^{+tλ}, cap)` | sharpening | yes (>0) |
| `tikhonov` | `1/(1+βλ)` | low-pass | yes |
| `ideal` | `1[λ∈interval]` | brick-wall low/high/band | yes (0/1) |
| `matern` | `(κ²+λ)^{-ν}` | SPDE/GP prior family | yes |
| `mexhat` | `(tλ)e^{-tλ}` | band-pass (vector t → bank) | no (zero at λ=0) |
| `dog` | `e^{-t₁λ}-e^{-t₂λ}` | band-pass | no |

`mexhat`/`dog` are admissible as *analysis* filters but not as priors (zero or
negative density), which the prior's check enforces.

## Testing

- **`dev/tests/test_eigfilter_pure.m`** (new): for each factory — correct handle
  output, analytic values at sample λ, limits (`heat t→0 ≈ 1`; `power`/`log`
  monotone decreasing; `matern` monotone; `mexhat`/`dog` zero at λ=0, peak at
  intermediate λ, →0 at large λ; `ideal` exact 0/1; `inverse_heat` clamps at
  `cap`); vector-scale returns a cell of the right length; dispatcher `list`
  includes every factory; `info`/`meta` returns a valid spec; unknown name
  errors; `compose` equals the pointwise product; `evaluate` handles both a
  handle and a bank.
- **Parity guards (critical):** `test_eigenmode_prior_pure.m` and
  `test_eigenmodes_filter_gain_pure.m` must pass **unchanged** after the refactor,
  proving `flat/power/log` and the seven filter-gain types are numerically
  identical to pre-refactor.

## Files

**New:**
- `toolbox/math/eigfilter/Contents.m`
- `toolbox/math/eigfilter/bst_eigfilter_kernel.m`  (registry: string -> factory + list/info)
- `toolbox/math/eigfilter/bst_eigfilter_design_{flat,power,log,heat,inverse_heat,tikhonov,ideal,matern,mexhat,dog}.m`
- `toolbox/math/eigfilter/bst_eigfilter_evaluate.m`
- `toolbox/math/eigfilter/bst_eigfilter_compose.m`
- `dev/tests/test_eigfilter_pure.m`

**Modified:**
- `toolbox/math/bst_eigenmode_prior.m` (delegate; keep scaling/DC/normalize/admissibility)
- `toolbox/math/bst_eigenmodes_filter_gain.m` (delegate; index→λ conversion for masks)

**Verify in the plan:** that Brainstorm's startup path setup adds the new
`toolbox/math/eigfilter/` subdirectory (it normally adds the toolbox tree
recursively; confirm, else add explicitly).

## Future work (documented, not built)

- **Filter-design GUI:** choose a kernel, tune hyperparameters with a live plot of
  `g(λ)` on the eigenvalue axis, and render the **vertex delta response** — apply
  the filter to an impulse at a chosen cortical vertex and display the resulting
  point-spread on the surface — so the spectral and spatial effects are visible
  together. The factory `'meta'` specs + handle API are designed to feed this GUI.
- **Joint time-vertex family** `bst_eigfilter_jtv_design_*`: damped-wave (neural
  field), diffusion, wave, Klein-Gordon — mirroring GSPBox `gsp_jtv_design_*`.
- **Neural-field spectrum** kernel: a physiologically-derived `g(λ)` from
  corticothalamic field theory, as a principled prior.
- **Library expansion** of analysis-side wavelet/scaling banks once the pattern is
  established and wired.
