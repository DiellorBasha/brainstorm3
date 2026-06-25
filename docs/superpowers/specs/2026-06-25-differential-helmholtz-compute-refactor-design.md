# Differential module refactor — Helmholtz as a stateless `Compute` algorithm

- **Date:** 2026-06-25
- **Status:** Draft for review
- **Scope:** Spec 1 of 2. This spec covers the `toolbox/differential/` algorithm refactor.
  The `panel_bst_dynamics` three-state GUI wiring is **Spec 2** (separate).

## 1. Problem / motivation

`bst_operators` is meant to be THE orchestrator for differential operators — resolve the
right operator node(s) and route them to the engines — mirroring the `bst_eigen ↔ tess_eigen`
pattern. For the Helmholtz path it is bypassed, and `bst_helmholtz` has absorbed
responsibilities that don't belong to it. Concrete findings (from reading the module):

1. **Orchestrator bypassed / resolution duplicated.** `view_helmholtz`,
   `process_helmholtz_events`, and `process_vortex_track` each *independently* resolve
   `Cov = tess_operators(...)`, `LBO = tess_operators(...)`, `Mani = tess_manifold(...)`,
   `Surf = in_tess_bst(...)` and call `bst_helmholtz('Prepare'/'Frame')` directly. The
   resolution logic is copy-pasted in every caller.

2. **`bst_helmholtz` is overloaded.** It bundles: operator assembly (`i_prepare_vertex`
   rebuilds `Wfv`, `Fvf`, the rotated gradient `Sx/Sy/Sz`, the pinned Cholesky of `K`, the
   1-ring adjacency); the div/curl + Hodge-potential algorithm (`Frame`); reconstruction
   (`Virr/Vsol/Vharm/HarmFrac`); vortex **core detection** (`FindCoresOp`, `FindCores`,
   sub-vertex quadratic fit — being retired, see §3); a **face** domain with no live consumer;
   a whole-series `Decompose` loop; and stateful `Prepare`/`Frame` verbs. It is not "primarily
   the algorithm that provides divergence and curl."

3. **Dependency inverted.** The ambient (3-D source field) div/curl math actually lives
   *inside* `bst_helmholtz`; `bst_divergence`/`bst_curl`'s `'Ambient'` branches are thin
   **re-exports** that call back into `bst_helmholtz('Decompose')`. This is backwards from
   the intended "helmholtz *sequences* `bst_divergence`/`bst_curl`."

4. **Dead operand + stale contract.** Every ambient/Helmholtz caller threads an `LBO`
   operand through as `OperatorNode{2}`, but the flat-covariant `i_prepare_vertex` reads
   only `OperatorNode{1}` (the Covariant node); `K`/`M` come from `Cov.Operator`/`Cov.Mass`.
   `LBO` is unused. The header still documents the old `{DiracNode, LBONode}` contract.

5. **Face-domain Helmholtz has no consumer.** All live callers use `Domain='vertex'`.

## 2. Goals

- `bst_operators` is the single **resolver/router** (the `bst_timefreq` analog) for the
  differential file types; callers stop resolving operators themselves.
- The Helmholtz decomposition becomes a **pure, stateless algorithm** reached through the
  **standard** Brainstorm process hook `process_helmholtz('Compute', …)`. **No invented
  hooks** (no `'Prepare'`, no `'Frame'`). Only `GetDescription`/`Run`/`Compute`/`FormatComment`.
- `bst_divergence` and `bst_curl` are the **real homes** of ambient divergence and vorticity.
- The three-state consumption (navigate / detect / save) all call the **same** stateless
  `Compute`, mirroring the `process_bandpass('Compute')` + `VisualizationFilters` precedent.
- `bst_helmholtz` and `view_helmholtz` are **retired**; **cores leave** to `bst_vortex_*`.

## 3. Non-goals (this spec)

- The `panel_bst_dynamics` three-state GUI wiring (`GlobalData.Dynamics` ephemeral options
  struct + a `FilterLoadedData`-style hook) — **Spec 2**.
- **Vortex / core detection.** The current detection stack (`bst_vortex_persistence`,
  `bst_vortex_track`, `bst_vortex_link_step`, `process_vortex_track`, and the core-detection
  code inside `bst_helmholtz`) is being **retired to `dev/experimental`** and rebuilt from
  scratch, better integrated into this architecture, in a later effort. This spec **does not
  use `bst_vortex_persistence`** and does not reimplement detection — it only produces the
  `Div`/`Curl` fields that a future detector will consume. Core detection is simply **removed**
  from the algorithm, not relocated.
- **Face-domain** Helmholtz — no live consumer; not ported. May return later as a separate
  `Domain` option if a consumer appears.
- `tess_operators` find-or-create (already implemented) and the Covariant operator math
  (unchanged).

## 4. Established-pattern anchors

- **Online filter:** `process_bandpass('Compute', F, …)` is the pure algorithm; `bst_memory`'s
  `FilterLoadedData` applies it *ephemerally* per loaded block via
  `GlobalData.VisualizationFilters`; the `process_bandpass` `Run` saves. `bst_timefreq` is
  **never** called by the filter path. → the ephemeral hook calls the pure algorithm directly,
  parameterized by a small `GlobalData` options struct.
- **`tess_cholesky`:** lazy, cached factor of an operator node (resolution: factor on node →
  in-session `persistent MEM` → compute+cache), already consumed by `bst_poisson`. Gives
  "factor once, reuse per displayed frame" with **no `Prepare` hook** — the cache is owned by
  the operator layer, keyed `Variant|surface|hh|pin`.

## 5. Target architecture

```
bst_operators        on-file orchestrator + operator RESOLVER  (the bst_timefreq analog)
   │  resolves Covariant via tess_operators / in_bst_operator; routes to engines; saves files
   ▼
process_helmholtz('Compute', J, Op)      PURE per-frame algorithm (stateless, I/O-free):
   ├── bst_divergence(J, …, 'Ambient', Cov)   ambient divergence  ← MOVED in from bst_helmholtz
   ├── bst_curl(J,        …, 'Ambient', Cov)   ambient vorticity   ← MOVED in from bst_helmholtz
   └── bst_poisson(Cov, rhs) → tess_cholesky   scalar potentials φ/ψ (cached factor, no Prepare)
```

Three states, all calling the same stateless `Compute`:

| State | Mechanism | Filter analog |
|---|---|---|
| **Navigate** | raw source frame, no transform | filter disabled |
| **Detect** (ephemeral, Spec 2) | `panel_bst_dynamics` hook → `process_helmholtz('Compute', Jframe, Op)` per displayed frame; `Op` resolved once via `bst_operators`, held in `GlobalData` | `FilterLoadedData` + `VisualizationFilters` |
| **Save** (on-file) | `process_helmholtz('Run')` loops `Compute`, writes result files | `process_bandpass` Run |

## 6. Component-level design

### 6.1 `bst_divergence` — real home of ambient divergence
In the `'Ambient'` branch, replace the delegation to `bst_helmholtz('Decompose')` with the
direct flat-covariant computation, per hemisphere, from `C = Cov.Covariant{hh}`:
- split `C.ScalarGrad` → `Gx,Gy,Gz` `[nFh x nVh]`; `Nf=C.FaceNormal`, `Af=C.FaceArea`, `Faces=C.Faces`;
- build the area-weighted face→vertex map `Wfv` (as in the current `i_prepare_vertex`);
- per-face surface divergence `divF = Gx*Jx + Gy*Jy + Gz*Jz`; `Div(vH) = s * (Wfv * divF)`, `s=+1`
  (calibrated: the flat-covariant divergence is the true divergence, includes the −2H(J·N)
  mean-curvature coupling — so the existing curvature note in the header stays valid; **no
  separate curvature term**).

### 6.2 `bst_curl` — real home of ambient vorticity
In the `'Ambient'` branch, replace the delegation with:
`curlV = [Gy*Jz−Gz*Jy, Gz*Jx−Gx*Jz, Gx*Jy−Gy*Jx]`; `omF = sum(curlV .* Nf, 2)`;
`Curl(vH) = s * (Wfv * omF)`.

### 6.3 Signatures — drop the dead `LBO`
`bst_divergence`/`bst_curl` ambient signatures become
`bst_<op>(J, ManifoldMat, 'Ambient', Surf, Cov)` (the `LBO` operand is removed everywhere it
was threaded). Update the three callers + `bst_operators`. *(Decision D1 — confirm.)*

### 6.4 `bst_poisson` — potentials home
`process_helmholtz('Compute')` recovers the scalar potentials φ (from the divergence source)
and ψ (from the vorticity source) via the weak Hodge route. The weak-source assembly (P1
sampling `Fvf`, rotated gradient `Sx/Sy/Sz`, `G'·W·Jf`) is Hodge-specific and stays in
`process_helmholtz('Compute')`; it feeds an assembled RHS to **`bst_poisson`**, which owns the
pinned solve via `tess_cholesky`. `bst_poisson` currently guards `Variant=='Laplace-Beltrami'`;
extend it to also accept the **`'Covariant'`** node (whose `Operator` is the SPD cotan
stiffness `g'Wg = cotanL`), keeping one Poisson home. *(Decision D2 — confirm vs. solving via
`tess_cholesky` inline.)*

### 6.5 `process_helmholtz` (new) — the algorithm + the save path
- **`Compute(J, Op[, OPTIONS])`** *(static, pure, I/O-free)* — input `J = [3nV x 1]` (or
  `[3nV x nT]`) ambient source; `Op` = resolved Covariant node. Sequences
  `bst_divergence` + `bst_curl` + `bst_poisson` + reconstruction. Returns a struct with
  **`Div`, `Curl`** (primary) and `Phi, Psi, Virr, Vsol, Vharm, HarmFrac, Fmag, Jnormal`
  (full vertex Hodge decomposition). **No cores.** Preserves the small-amplitude `HarmFrac`
  guard (`harmDen > 0`, *not* `max(harmDen,eps)` — the ~1e-22 denominator note).
- **`Run(sProcess, sInputs)`** — the **Save** state. Resolve `Surf`/`Cov` via `bst_operators`;
  loop `Compute` over time (or event times); write `results_` file(s) (Curl primary;
  Div/Phi/Psi optional per OPTIONS). Folds in the `process_helmholtz_events` behavior as a
  Run option.
- **`GetDescription`** (category Differential; options: FieldType, Gauge, output selection,
  event-time vs whole-series) and **`FormatComment`**.

### 6.6 `bst_operators` — resolver/router
The `'helmholtz'` Method delegates to `process_helmholtz('Compute', F, Cov)` instead of
`bst_helmholtz('Decompose')`. Operator resolution stays in `bst_operators`
(`Cov = tess_operators(SurfaceFile,'Covariant')`). `'divergence'`/`'curl'` ambient Methods drop
the dead `LBO`.

### 6.7 Core detection — removed, not relocated
The core-detection code inside `bst_helmholtz` (`FindCoresOp`, `FindCores`, `i_subvertex`,
`i_make_core`, `i_empty_cores_vertex`, and the face-core variants) is **deleted with
`bst_helmholtz`** and **not reimplemented**. `process_helmholtz('Compute')` produces the
`Div`/`Curl` (vorticity / source-sink) fields only; detecting cores *from* those fields is the
job of the new detector built in a later effort. This spec **must not** call
`bst_vortex_persistence` (which is itself being retired to `dev/experimental`).

### 6.8 Retirements
- **`bst_helmholtz.m` deleted.**
- **`view_helmholtz.m` — compute repointed, NOT deleted in Spec 1.** It is already wired into
  the dynamics GUI (`view_dynamics.m:78` opens it; `panel_bst_dynamics` drives its
  `SetComponent`/`SetSmoothing`/`UpdateFrame` verbs; `panel_helmholtz` is its options panel),
  so deleting it now would break `panel_bst_dynamics`/`view_dynamics` — that rewiring is Spec 2.
  Spec 1 only severs its `bst_helmholtz` dependency: resolve `Cov` once via `bst_operators`
  (kept in figure app-data as today), and per displayed frame call
  `process_helmholtz('Compute', Jframe, Cov)`. Display verbs untouched. `view_helmholtz` +
  `panel_helmholtz` **deletion is deferred to Spec 2**.
- **`process_helmholtz_events`** re-implemented on `process_helmholtz('Compute')` (or folded
  in as a Run option). This is a Helmholtz *maps* process (J / |J| / Phi / Psi at event
  times) — it does **no** core detection, so it migrates normally.
- **Vortex-detection stack out of scope.** `process_vortex_track`, `bst_vortex_track`,
  `bst_vortex_link_step`, and `bst_vortex_persistence` are part of the detection rebuild and
  are **retired to `dev/experimental`** — they are *not* migrated to `process_helmholtz` here.
  Because `process_vortex_track` calls the (deleted) `bst_helmholtz`, it must be **moved out at
  the same time** so nothing in `toolbox/` references the deleted engine.

## 7. Error handling

- `Compute` validates `size(J,1) == 3*nV`; clean error otherwise (mirror `bst_operators` guards).
- Missing Covariant node → `bst_operators`/`tess_operators` find-or-create (already in place);
  missing Structures atlas → existing `tess_operators` guard.
- `tess_cholesky` errors if the pinned `K` block is not SPD (existing).

## 8. Testing

- **Div/Curl parity:** `process_helmholtz('Compute', J, Cov)` reproduces the *old*
  `bst_helmholtz('Frame')` `Div/Curl/Phi/Psi/HarmFrac` to ~1e-12 on the canonical cortex +
  a synthetic source field. Basis: `dev/test_helmholtz_covariant.m`,
  `dev/test_ambient_divcurl.m` (update to the new API; assert against a captured baseline
  taken *before* the refactor).
- **Engine parity:** `bst_divergence`/`bst_curl` ambient outputs equal the old `H.Div`/`H.Curl`.
- **Maps end-to-end:** `process_helmholtz_events` (re-implemented on `Compute`) produces the
  same J / |J| / Phi / Psi maps as before on the alpha test block.
- *(No cores parity test — core detection is removed; see §3, §6.7.)*

## 9. Migration / sequencing (for the implementation plan)

1. Move ambient **div** into `bst_divergence`; ambient **curl** into `bst_curl`. Temporarily
   leave `bst_helmholtz` delegating to them to prove bit-parity.
2. Create `process_helmholtz` (`Compute` sequences div+curl+poisson+reconstruction; `Run` saves).
3. Repoint `bst_operators` `'helmholtz'` Method and `process_helmholtz_events` to
   `process_helmholtz('Compute')`; drop the dead `LBO` operand.
4. Move the vortex-detection stack (`process_vortex_track`, `bst_vortex_track`,
   `bst_vortex_link_step`, `bst_vortex_persistence`) to `dev/experimental` (out of scope; see §3).
5. Repoint `view_helmholtz`'s compute to `process_helmholtz('Compute')` (do **not** delete it).
   Delete `bst_helmholtz`. Verify nothing in `toolbox/` still references the deleted engine.
6. Update dev tests to the new API.

## 10. Open decisions to confirm

- **D1:** Drop the dead `LBO` operand from ambient div/curl signatures (recommended) vs keep
  for signature stability.
- **D2:** Extend `bst_poisson` to accept the `'Covariant'` node (recommended, one Poisson home)
  vs solve via `tess_cholesky` inline in `process_helmholtz`.
- **D3:** Face-domain Helmholtz dropped (no consumer) — confirm.
- **D4:** `panel_helmholtz.m` fate (retired with `view_helmholtz` vs repurposed) — deferred to
  Spec 2 but flagged here for the openers audit.
