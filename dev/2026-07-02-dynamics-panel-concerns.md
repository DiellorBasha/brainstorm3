# Dynamics panel — concerns, issues & optimization opportunities (living doc)

> **Purpose:** a running log of concerns, design smells, bugs, and optimization opportunities noticed
> while walking through the Dynamics-panel workflow. This is the **source document** for a future
> development-session item list — entries here graduate into concrete tasks.
>
> **Companion:** [`dev/2026-07-02-dynamics-gui-workflow-and-tests.md`](2026-07-02-dynamics-gui-workflow-and-tests.md)
> — the validated GUI workflow + test targets. Each concern below cross-references the workflow doc's
> heading (§) and, where relevant, the code line.
>
> **Legend:** severity `⛔ bug` · `🐢 perf` · `🏗️ design` · `✨ opportunity` — priority `🔴 HIGH` · `🟠 MED` ·
> `🟡 LOW` — status `OPEN` · `SCOPED` · `DONE`.
> **How to add:** append a `### C<n>` entry; keep the cross-refs; note effort/risk; link related entries.

---

## Index
| ID | Title | Sev | Prio | Status | Workflow ref |
|----|-------|-----|------|--------|--------------|
| **N1** | **North-star: pipelines as composed operators (sensors → measurement, reused over frames)** | 🏗️✨ | ⭐ FRAMING | OPEN | (all) |
| **C3** | **Extract atom design to an eigen-node-scoped designer; Dynamics becomes a consumer** | 🏗️✨ | 🔴 HIGH | OPEN | §0; §1; §3-A |
| **C6** | **Atom seed-pick abuses the geodesic Region tool (heavy, clobbers quivers, scroll-resizes, dies)** | ⛔🏗️ | 🔴 HIGH | OPEN | §3-A |
| C1 | Covariant operator is eagerly loaded into figure appdata at panel open | 🐢🏗️ | 🟠 MED | OPEN | §0 Open the panel |
| C2 | No unified load/unload lifecycle for operators (two ad-hoc caches) | 🏗️✨ | 🟠 MED | OPEN | §0; §3-F/§3-G |
| C4 | Modularize Atom-filter → Helmholtz as composable process stages (field is the interface) | 🏗️✨ | 🟠 MED | OPEN | §3-F; §3-G |
| C5 | Differential mode kernel: fold Helmholtz∘realize into a reusable `[nV×K]` operator | 🏗️🐢✨ | 🟠 MED | OPEN | §3-F; §3-G |
| C7 | Atom properties hidden; restructure panel into a resizable list/detail split | 🏗️ | 🟠 MED | OPEN | §1; §3-A |
| C8 | Consolidate Measure / Helmholtz / Threshold into one Kernel-selection menu (registry front-end) | 🏗️✨ | 🟠 MED | OPEN | §1; §3-F/§3-G |

---

### N1 — North-star: pipelines as composed operators (sensors → measurement, reused over frames)
**Type:** 🏗️ architecture framing · **Status:** OPEN (guiding principle; C1–C5 are instances)

**The principle.** Most analysis pipelines here are chains of **linear maps between per-frame
representations** — sensor space → coefficient space → field space → differential space. Because linear
maps compose, a chain collapses into **one operator** that maps sensor data directly to the measurement,
**precomputed once and applied per frame on the fly** — never materializing the intermediates. This is
Brainstorm's `ImagingKernel` philosophy (`J = K·d`, no dense source maps) generalized to *every* stage.

```
d  --A-->  c   --g-->  cf  --R-->  J  --H_O-->  O          collapses to   O = (H_O·R·g·A)·d = K_O·d
sensor    coeff     filtered    field       measurement
```

**Instances already in the codebase (this is emergent, not hypothetical):**
| Operator | Maps | Folds away |
|---|---|---|
| `ImagingKernel` (core) | sensors → source field | the forward/inverse solve |
| **Mode kernel** `A = Φᵀ M P K` (`dirac-connectome`) | sensors → eigen-coefficients | reconstruct-then-project (~500×) |
| **`ScalogramEnergy` Gram** `E = cᵀ G c` (Item 2, shipped) | coefficients → per-band energy | the `[nV×nT×M]` field (~97 GB→~300 MB) |
| **Differential mode kernel** `B_O = H_O·R` (**C5**, proposed) | coefficients → Div/Curl/Φ/Ψ | the per-frame Poisson solve |

**Three tiers of the algebra (the boundary conditions — where it holds / breaks):**
1. **Linear** → a single matrix `K_O` (source map, projection, realize, divergence/curl/potential, band
   projection). Compose freely; precompute the product.
2. **Bilinear / quadratic** → a **Gram** operator `cᵀ G c` (energies, power, variance). `ScalogramEnergy`
   is this tier — quadratic measures fold into a `[K×K]` Gram, not a single matrix, but still no field.
3. **Nonlinear** → the **materialization boundaries**: magnitude `|J|`, ratios (`HarmFrac`), thresholds /
   level-sets (Localize), argmax / peak-finding, data-adaptive normalization. These must run on a
   realized output; they mark where a composed operator must *stop* and hand off a concrete field.

**Extensions:**
- **Time axis.** Temporal / dynamic (ts/js) stages (bandpass, wavelet-in-time) are operators *along time*
  (convolutions), not per-frame — they compose in the time dimension instead of collapsing across it.
- **Data-adaptive operators.** Beamformers / data-covariance whitening depend on the data, but recompute
  the operator **once per block** and it stays linear per frame — still fits, just re-derived per block.

**The opportunity (what this is becoming).** A small **operator-composition framework**: each stage
declares its `(input space, output space, tier)`; the framework **composes and caches** the product as a
reusable kernel keyed by `(surface, eigenbasis, filter)`. Then "measurements" — source map, scalogram,
Helmholtz, connectivity, … — become **composed operators from sensors**, and the panel/process layer just
*selects and applies* them. This subsumes C1/C2 (operators are the cached, load/unload-able units), C4
(stages compose; the interface can be coefficients, not a materialized field), and C5 (one such composition).
**Related:** nxr `operatormat.Registry`/`bst_nxr_registry`, `dev/2026-07-02-frame-based-inverse-roadmap.md`,
memories `eigen-module-reorg`, `operator-resolution-refactor`, `eigen-spectral-consolidation`.

**Effort/Risk:** large / medium (a framework, built incrementally — each C-item lands a piece). Not a
single task; a **direction** that the numbered concerns realize one operator at a time.

#### N1.a — Naming & type taxonomy (⚑ don't blanket-prefix everything `IK_`)
The registry has **two layers on different axes**; only the sensor-anchored composites are "imaging kernels".
Source-side primitives are **inverse-independent** (reusable across any source/subject/method) — that's their
value, and an `IK_` prefix would falsely claim a sensor dependence they don't have.

**Layer 1 — generalized imaging kernels `ik_<measurement>` (sensor → measurement; fold in `K`):**
| Name | Maps | Recipe |
|---|---|---|
| `ik_source` | sensor → field | `K` (classic ImagingKernel) |
| `ik_coeff` | sensor → coeff | `manifold_ft ∘ embed ∘ K` (= the mode kernel `A`) |
| `ik_field` | sensor → filtered field | `manifold_ift ∘ g(λ) ∘ ik_coeff` |
| `ik_div`/`ik_curl`/`ik_helmholtz` | sensor → differential map | `H_O ∘ manifold_ift ∘ g(λ) ∘ ik_coeff` (C5) |
| `ik_energy` | sensor → band energy | `Gram ∘ g(λ) ∘ ik_coeff` (ScalogramEnergy from sensors) |

**Layer 2 — source-side primitives (NO `ik_`; canonical names already exist — catalogue, don't rename):**
| Vague alias used above | Real primitive (function) | Family |
|---|---|---|
| `Project` | **`manifold_ft`** (`ΦᵀM`, field→coeff) | `spec_` |
| `Realize` | **`manifold_ift`** (`Φ`, coeff→field; +`manifold_quat_imag` for Dirac) | `spec_` |
| spectral filter | **`bst_eigfilter_kernel`** `g(λ)` (the existing eigfilter registry) | `spec_` |
| `Embed` | **`manifold_quat_imag`** + `bst_eigenfilter('RowMap')` | (algebraic) |
| gradient | **`ScalarGrad`** = `[Gx;Gy;Gz]` (Covariant node) | `hodge_` |
| divergence / curl | **`bst_divergence`** / **`bst_curl`** | `hodge_` |
| Poisson solve | **`bst_poisson`** | `hodge_` |
| solver (factor+back-sub) | **`tess_cholesky`** (pinned) | `hodge_` |
| face↔vertex | **`bst_face2vertex`** | `hodge_` |
| eigenbasis Φ/Λ | **`bst_eigen`/`tess_eigen`** node | `spec_` |
| the inverse | **`ImagingKernel`** (Results) | `inv_` |

**⚑ Derivation rule:** `ik_` is a **property of the type signature**, not a hand-assigned prefix — a kernel
is an imaging kernel **iff its input space includes `sensor`**. The mode kernel is `ik_coeff` *because*
`manifold_ft ∘ embed ∘ K` contains `K`; `manifold_ft` alone is not, *because* it starts from `field`. Source→
measurement composites (no sensors, e.g. `B_O = H_O·R` on coefficients) get `<op>_kernel`
(`helmholtz_kernel`, `energy_gram`), not `ik_`.

**Registry entry shape:** `{ name, inSpace→outSpace, tier(linear|bilinear|nonlinear), primitives:[…],
deps:[surface, eigen-node, cov-op, …], sensorAnchored(derived) }`. `manifold_ft`/`manifold_ift` are the
canonical `Project`/`Realize` primitives; the mode kernel is the **first recipe** (`manifold_ft ∘ embed ∘ K`,
literally built by calling `manifold_ft` on the embedded ImagingKernel in `i_vector_modekernel`,
`panel_bst_dynamics.m:1256–1257`) — so the registry's first primitives and first recipe already exist in code.
⚑ dep nuance: `Project`=`manifold_ft` needs `(Φ, mass M)`; `Realize`=`manifold_ift` needs only `Φ` (synthesis
carries no metric) — different invalidation keys.

---

### C3 — Extract atom design to an eigen-node-scoped designer; Dynamics becomes a consumer
**Severity:** 🏗️ design + ✨ opportunity · **Priority:** 🔴 HIGH · **Status:** OPEN
**Workflow ref:** §0 *Open the panel* (Dynamics currently owns both design AND real-source); §1 *Panel
anatomy* (the Atom + Frame sections are the design half); §3-A *Design an atom's impulse response*.
**Code:** duplicated design surfaces — `view_atom_designer.m` + `panel_atom_designer.m` (surface-node
"Design atom...", `tree_callbacks.m:1198`) **and** the design half of `panel_bst_dynamics.m` (Atom/Frame
sections, `OnCreateAtom`/`OnSetOperator`/`OnDesignFrame`/`i_atom_realise*`). Both call the same core
`bst_eigenfilter('Atom')` / `bst_eigenwavelet` but have **diverged** at the panel level.

**Observation (the design is data-free — proven by trace).** An atom's Design/impulse is a pure function
of `(operator eigenbasis Φ/Λ, kernel g(λ), seed vertex, direction)`. It reads **neither** the imaging
kernel **nor** sensor data — those enter only at Apply/Analyze/Localize/Helmholtz (the "real-source"
verbs). The stored atom `G` holds only `{Operator name, KernelName, KernelParams, vertices, SurfaceFile}`
(`i_default_atom`, `panel_bst_dynamics.m:484`; the in-code comment at `:510` states it "does not depend
on the dSPM"). The recording contributes only its **sample rate** (a nominal time axis for dynamic
kernels), not its data. So "design" is cleanly separable from the inverse+data context.

**The problem is duplication + drift, not a missing separation.** A standalone designer already exists
(`view_atom_designer`, data-free, surface-node, nominal 1 s @ 100 Hz), but the Dynamics panel re-grew a
richer design half (full operator set incl. **Dirac / Dirac-Connectome**, the **frame designer**) that
was never back-ported. Two half-separations, one stale:

| | `view_atom_designer` (surface node) | design-half of `panel_bst_dynamics` (source node) |
|---|---|---|
| Data context | none (data-free — correct) | inherits inverse kernel + recording |
| Operators | scalar only (LB-Connectome, Laplace-Beltrami) | **full set incl. Dirac, Dirac-Connectome** |
| Time base | nominal 1 s @ 100 Hz | recording Fs / 4 s window |
| Frame designer | no | yes (itersine tight frame) |
| Realiser | `bst_eigenfilter('Atom')` | `bst_eigenfilter('Atom')` — same core, **duplicated UI** |

**Proposed direction.** Make design a **surface/eigen-node concern** and reduce Dynamics to a **consumer**:
- **Bind the designer to the EIGEN NODE, not just the surface.** An atom = `g(λ)` on a spectrum + a seed;
  the spectrum *is* the eigen node (inherits Φ/Λ, λ-range, mode count directly; seed resolves via
  eigen-node → operator → `ParentSurface`).
- **Dynamics inherits the inverse kernel + time window** and loads an atom bank keyed to an eigen node,
  owning only the real-source verbs.

**⚑ The correctness reason for eigen-node granularity (not surface-level).** Dynamics applies filters on
the **inverse's own `DiracEigenFile` basis**, not a fresh canonical one. A *surface*-level designer builds
a canonical basis that **won't match** what Dynamics applies → the impulse you designed ≠ the filter that
runs. Eigen-node binding makes **design-basis == apply-basis by construction**. This is the crux — it also
resolves the "which basis?" ambiguity behind C1/C2.

**Tensions to resolve (contracts):**
1. **Cross-node application guard.** An atom bank is valid *for a given eigen node*; if Dynamics applies it
   on a source whose inverse used a different node, **warn/guard**, don't silently mismatch. (Main risk.)
2. **Dynamic-kernel time base.** Keep atom params in **physical units** (mm via `2π/√λ`, s, m/s — already
   how `view_atom_designer` parametrizes) so they're Fs-independent and re-realise at any rate; contract:
   "params physical, time axis context-supplied."
3. **UX — separate concern, not necessarily separate pixels.** Make `panel_atom_designer` a **shared
   docked component** summonable in both the eigen-node view and inside a Dynamics session (design in
   context of the real source without a two-window round-trip).
4. **No third divergence.** Both entry points must end up driving **one** design implementation
   (`panel_atom_designer` + the `bst_eigenfilter`/`bst_eigenwavelet` core).

**Staged plan (each independently shippable):**
1. **De-duplicate:** port Dynamics' richer design UI (Dirac/Dirac-Connectome operators, frame designer)
   into `panel_atom_designer`; both panels drive one implementation. No behavior change; removes drift.
2. **Re-anchor to the eigen node:** add "Design atoms" on the eigen_ node (inherits Φ/Λ); keep the
   surface-node entry as "pick/attach an eigen node".
3. **Slim Dynamics to a consumer:** load an atom bank keyed to an eigen node; own only Apply/Analyze/
   Localize/Helmholtz, inheriting the inverse kernel + time window. Plumbing partly exists
   (`bst_dynamics('Save'/'Load')`, `process_source_atoms`).

**Effort/Risk:** medium-large / medium (multi-panel refactor + a basis-consistency contract), but each
stage ships independently and stage 1 is pure de-duplication. **Subsumes** the operator-lifecycle
questions in **C1/C2** (the eigen-node designer owns its basis/operator). **Related memories:**
`atom-designer-architecture`, `dynamics-atoms-system`, `dynamics-design-preview-modes`,
`atom-operator-applicability`, `eigen-module-reorg`.

---

### C1 — Covariant operator is eagerly built/loaded into figure appdata at panel open
**Severity:** 🐢 perf + 🏗️ design · **Status:** OPEN
**Workflow ref:** §0 *Open the panel* — the statement *"The cortex figure carries a `DynamicsOverlay`
appdata with the source link + the surface **Covariant** operator (`D.Cov`, for Helmholtz)."*
**Used by (workflow):** §3-F *Helmholtz on the filtered field* (step 2) and §3-G *Raw-source differential
maps* (step 1–2); listed in the §1 toolbar table (Measure, Helmholtz rows).
**Code:** load site `view_dynamics.m:371` (`Cov = tess_operators(SurfaceFile,'Covariant')`, behind the
`"Loading the covariant operator..."` progress bar); stashed at `view_dynamics.m:390–393`; consumed
per-frame at `view_dynamics.m:431` (`i_dynamics_overlay` → `process_helmholtz('Compute', Jt, D.Cov)`).

**Observation.** `i_open_source_figure` builds `Cov` **unconditionally** at open, even though the entire
atom workflow (Design / Apply / Analyze / Localize) never uses it. `Cov` is needed **only** for the
Measure differential maps and Helmholtz-filtered. `tess_operators(...,'Covariant')` is find-or-create,
so on a surface that has never had it, this **builds** it (nxr cotan-Laplacian + galerkin mass + vertex
frames + the `g'Wg == cotanL` consistency check) — a multi-second cost imposed on every session open,
for a feature most sessions don't touch.

**Why it's in appdata at all (keep this part).** Storing the *handle* in `DynamicsOverlay` is correct —
it is figure-scoped render state, not a one-shot computation:
1. **Hot-path residency** — once a Measure is active, `Cov` is hit on every time-cursor scrub; it must be
   memory-resident, not re-read per frame.
2. **Co-located cache** — `D.Cache` (per-frame Hodge results) lives beside `D.Cov`; the repaint hook reads both.
3. **Lifetime coupling** — the overlay (and `Cov`) is owned by `hFig`, freed on figure close (idiomatic,
   like `Surface`/`TessInfo` appdata).

**Proposed direction (lazy-load *into* appdata).** Keep the appdata home; defer the load to first
differential use:
```matlab
% i_open_source_figure: D.Cov = [] (do NOT build here)
function D = i_ensure_cov(hFig)             % called by OnMeasurement / i_dynamics_overlay / OnHelmholtzFiltered
    D = getappdata(hFig,'DynamicsOverlay');
    if isempty(D.Cov)
        D.Cov = tess_operators(D.SurfaceFile,'Covariant');   % find-or-create, once
        setappdata(hFig,'DynamicsOverlay', D);
    end
end
```
Same hot-path residency + caching; the atom-only workflow pays nothing.

**Trade-off.** Eager load buys *fail-fast* (operator errors, e.g. `CovariantInconsistent`, surface at
open, not mid-scrub). Preserve it cheaply by keeping a **cheap** capability check (unconstrained /
nComponents==3 / surface has faces) at open, and deferring only the **expensive** build to first use.

**Effort/Risk:** small / low. Good companion to the `OnLocalizeBands` fix (both = "do expensive work
only when the user asks"). **Leads into → C2.**

---

### C2 — No unified load/unload lifecycle for operators (two ad-hoc caches)
**Severity:** 🏗️ design + ✨ opportunity · **Status:** OPEN
**Workflow ref:** §0 *Open the panel* (Covariant) and §3-A..G (atom operators loaded via `i_atom_axes`).
**Code:** `Cov` cached in figure appdata (`view_dynamics.m:390`); atom eigen-axes cached in a root-`0`
`containers.Map` `DynamicsAtomAx` (`panel_bst_dynamics.m:553/576/593`).

**Observation.** Operators are find-or-create **pure artifacts** (determined by `surface|variant|tau`,
persisted as `operator_*.mat`, idempotent to reload/rebuild) — so unload-and-reload-on-demand is always
safe. Yet the panel caches them **two different ways** with **no eviction**: figure-appdata for `Cov`,
a session-global map for atom axes. Both grow monotonically until the figure/session closes. On the
17 GB Mac Mini this is exactly the shape that produced the earlier OOM (operators + eigenbases are
tens–hundreds of MB each; see the mode-kernel history in the `dirac-connectome-operator` memory).
⚑ **Live symptom (walkthrough):** the **first atom** has a long render delay — `OnCreateAtom` →
`i_atom_axes` → `bst_eigen('Axes')` builds/loads the 60-mode eigenbasis on first use (`panel_bst_dynamics.m:591`),
then realises. A cached/precomputed eigenbasis (or clearer progress) would remove the first-create stall;
it's the same operator-lifecycle gap, felt as latency rather than memory.

**Proposed direction (managed session operator cache).** One registry keyed `surface|variant|tau` with:
- **pin-the-active-view** — never evict the operator(s) backing what's currently on screen (the
  differential overlay's `Cov`, the selected atom's eigen-axes);
- **LRU / memory-budget eviction** of idle entries;
- **prefer persisted nodes** (reload is a cheap `.mat` read; only `NoSave` entries force a rebuild);
- **invalidation** hooked to operator/surface recompute (reuse the cascade-delete already in `tess_operators`).
This subsumes both current caches (C1's appdata `Cov` becomes a pinned entry while a Measure is active).

**Reuse existing infra, don't invent a third.** Align with (a) Brainstorm's `bst_memory` load/share/
**unload** lifecycle (the established pattern for datasets), and (b) the nxr `operatormat.Registry` /
`bst_nxr_registry` (v0.2.0 operator-registry). Investigate whether one of these can host the cache.

**Open questions:**
- Where does the cache live — a `bst_memory`-style global, or a dedicated `bst_operator_cache`?
- Budget policy: fixed MB, fraction of free RAM, or count-based?
- Interaction with the existing `DynamicsAtomAx` map and any figure-appdata pins.

**Effort/Risk:** medium / medium (touches caching + memory lifecycle; needs a benchmark on the Mac Mini).
Depends on C1 landing first (establishes the lazy pattern). **Related:** `dirac-connectome-operator`
(mode-kernel / OOM history), `operator-resolution-refactor` (find-or-create), `nxr-compute-plugin`
(operator registry).

---

### C4 — Modularize Atom-filter → Helmholtz as composable process stages (field is the interface)
**Severity:** 🏗️ design + ✨ opportunity · **Priority:** 🟠 MED · **Status:** OPEN
**Workflow ref:** §3-F *Helmholtz on the filtered field*, §3-G *Raw-source differential maps*.
**Code:** `process_helmholtz('Compute', J, Cov)` (`process_helmholtz.m:39`) = the modular decompose unit;
`process_helmholtz_events.m` = the pipe **already existing** with a *temporal-bandpass* pre-filter
(`J = ImagingKernel*bandpass(sensors); Ht = process_helmholtz('Compute', J, Cov)`, `:168–184`);
`OnHelmholtzFiltered` = the interactive one-frame version with the *atom* filter as the pre-filter;
realization via `i_vector_coeffs`/`i_dirac_recon_display` (`panel_bst_dynamics.m`).

**Observation.** Atom/Frame analysis and Helmholtz are **two orthogonal transforms that compose through a
serializable interface — a cortical vector field `[3nV × nT]`** — not a necessary fusion:
- `process_helmholtz('Compute', J, Cov)` is agnostic to where `J` came from (raw / bandpassed / atom-filtered).
- The atom filter is agnostic to what's downstream (magnitude display / scalogram / Helmholtz / localize).
- The pipe **already exists in process form** (`process_helmholtz_events`); the atom/frame filter is simply
  **another pluggable stage-1 pre-filter** (spectral `g(λ)` instead of a temporal bandpass).
- Even inside the panel the two toggle independently (Analyze = atom, no Helmholtz; Measure = Helmholtz, no
  atom; Helmholtz-filtered = atom→Helmholtz). So the "integration" is UX, not computation.

**Why Helmholtz is necessarily the TAIL stage (the realized field is the boundary).** The efficient atom
path stays in **Dirac mode-coefficient space**: the source enters the eigenbasis via
`c = ImagingKernelMode · data` (or the cached vector mode kernel `A·data`), the atom filters there
(`cf = g(λ).·c`, `[K × nT]`, cheap), and the vector field is **realized only on demand, per frame/window**
(`Φ·cf → [nV × 3]`). Helmholtz uses the **Covariant grad/curl operator**, which acts on the *spatial*
field, **not** on eigen-coefficients (a different basis) — so it must consume the **realized** field and
therefore sits after coefficient-domain filtering AND per-frame realization. *(In principle `Φ` could be
folded into a "differential mode kernel" `H·Φ` to skip the full-field intermediate, but the outputs are
per-vertex fields anyway and Helmholtz's operator is spatial — so the realized field is the natural, correct
stage boundary.)*

**Order matters → keep the stages explicit.** `filter ∘ Helmholtz ≠ Helmholtz ∘ filter` (linear but in
different bases: Dirac source eigenbasis vs Covariant grad/curl; live `corr 0.14`). The efficient path is
inherently *filter → realize → Helmholtz*; the reverse would need a *scalar* eigenfilter on the differential
maps (a different operator). Because order is a scientific choice, expose it as composable stages, not a monolith.

**Proposed modular shape:**
- **Stage A — field producer** (process): reconstruct a vector source field over a window/events with a
  pluggable pre-filter `none | temporal-bandpass | atom-eigenfilter g(λ) | eigenwavelet-per-band` →
  standard results file. (`process_helmholtz_events` is the bandpass instance; add the atom/wavelet
  instances by reusing `i_vector_coeffs`/`i_dirac_recon`.)
- **Stage B — decompose** (process): already `process_helmholtz('Compute')` → Div/Curl/Φ/Ψ maps.
- **Optional:** per-band Helmholtz (decompose each frame member → a *differential scalogram*) — only makes
  sense as a pipe.
- Panel stays the **interactive explorer** (live per-frame `D.Cache` scrub, design-in-context); processes own
  batch/repro/compose.

**Effort/Risk:** medium / low-medium — mostly *exposing the atom-filtered field as a stage-A process*
(a sibling of `process_helmholtz_events`), then chaining into `process_helmholtz`. **Related:** C3 (design
extraction — atoms become the reusable filterbank this pipe consumes), C1/C2 (each stage owns/loads its own
operator; the field is the only shared intermediate). **Memories:** `helmholtz-events-process`,
`bst-helmholtz-consolidation`, `dirac-connectome-operator` (mode kernel).

### C5 — Differential mode kernel: fold Helmholtz∘realize into a reusable `[nV×K]` operator
**Severity:** 🏗️ design + 🐢 perf + ✨ opportunity · **Priority:** 🟠 MED · **Status:** OPEN
**Workflow ref:** §3-F *Helmholtz on the filtered field*, §3-G *Raw-source differential maps*.
**Code:** `process_helmholtz('Compute')` — `Div/Curl` already batched (`:47–48`), but `Φ/Ψ/Virr/Vsol`
carry a **per-frame Poisson solve** in the `for t=1:nT` loop (`:63–82`); `R` = realize
(`i_dirac_recon`/`i_dirac_recon_display`); `A` = mode kernel (`i_mode_coeffs`/`i_vector_modekernel`).
**This is the kernelized form of C4's stage B, and an instance of N1.**

**The algebra.** Every Hodge output is linear in the field `J = R·cf`: `O = H_O·J`. So
`O = H_O·R·g·A·d`, and the products associate:
- `B_O = H_O·R`  `[nV×K]` — **Helmholtz-on-eigenmodes**, atom- and data-independent (pure surface+eigenbasis).
  The reusable core; cache as a derived operator node keyed `(surface, eigen-node, output)`.
- `K_O = B_O·g·A`  `[nV×nCh]` — full fold to sensors (**"differential imaging kernel"**): `Div(t)=K_div·d(t)`,
  one matmul, no field, no coefficients. Atom-specific via `g`; for the raw map (`g=I`) it's `H_O·ImagingKernel`.

**The win: eliminates the per-frame Poisson solve.** Build `B_φ = H_φ·R` by applying `H_φ` to each of the
`K` columns of `R` → **`K` solves once** (cached-Cholesky back-subs, sub-second for `K≈180`). Then `Φ` over
the whole recording is `B_φ·(g·C)` — a single `[nV×K]·[K×nT]` matmul, **zero per-frame solves** (on the test
recording: 74381 solves → 180). Same structural win as the mode kernel.

**Caveats (honest boundaries):**
1. **Linear outputs only** — `Div, Curl, Φ, Ψ, Virr, Vsol`. `Fmag=|J|`, `Hmag`, `HarmFrac` are nonlinear
   (N1 tier-3) → compute from the realized frame when needed.
2. **Amortized, not one-shot** — a single interactive frame is already cheap with the sparse solve, and
   `B_O` costs `K` solves to build; folding wins across many frames (whole-window analysis, differential
   scalogram, scrub-heavy overlays). *(This corrects an earlier "not worth folding" — true only for one frame.)*
3. **Density/storage** — `B_O` dense `[nV×K]` (~14 MB/hemi/output), one per output; modest, cacheable.
4. **Lifecycle** — depends on `(surface, eigenbasis, cached Cholesky)`; invalidate on recompute (C1/C2).

**Payoff for C4.** Stage B need not be a per-frame `process_helmholtz` call — it can *be* `B_O`, so the
atom↔Helmholtz interface becomes **coefficients** `cf` (not a materialized field) — strictly more in the
kernel-on-the-fly spirit. **Related:** N1 (the general pattern), C4 (staging), C1/C2 (operator lifecycle);
memories `helmholtz-events-process`, `bst-helmholtz-consolidation`, `dirac-connectome-operator`.

---

### C6 — Atom seed-pick abuses the geodesic Region tool (heavy, clobbers the vector display, scroll-resizes, and dies)
**Severity:** ⛔ bug + 🏗️ design · **Priority:** 🔴 HIGH · **Status:** OPEN (live-observed during walkthrough)
**Workflow ref:** §3-A *Design an atom's impulse response* (the Localize "SEL" step).
**Code:** `OnLocalize` (`panel_bst_dynamics.m:649`) → `bst_geodesic_tool('Toggle')`; `SyncSource` (`:298`)
uses only `gs.seed`; `bst_geodesic_tool.m` (a general **Region tool** = seed + geodesic disk + scroll-grow).
**Contrast:** `view_atom_designer` already seed-picks the *right* way via the one-shot **`WaveletDesignerPick`**
hook (`panel_bst_dynamics.m:766–819` has the same hook for seed *direction*, unused for the seed vertex).

**Observed (all live).** The Localize toggle repurposes a heavyweight Region/disk tool, but the atom only
needs *which vertex was clicked* — so every extra thing it does is wasted or harmful:
- **Delay + wasted work:** `Toggle('on')` runs `tess_scout_area('prewarm')` → **pre-factorizes a heat
  operator** (`bst_geodesic_tool.m:67`) just to arm a click; each click runs a **geodesic-distance solve**
  (`Seed`, `:100`).
- **Clobbers the Design display:** each click `Draw`s a **3 mm geodesic disk as a scalar region** → this is
  what replaces the atom's **vector quivers** with a scalar colormap (before/without any intent to).
- **"Two jobs" / scroll-resize:** `OnScroll → Grow → Draw` (`:141–147`) grows/shrinks the disk radius — the
  geodesic **area** tool bleeding into what should be a pure seed pick. It **does not belong in atom design**.
- **No armed-mode affordance:** only the cursor changes to a crosshair (`:62`); no on-cortex indication.
- **⛔ Dies after a few clicks** (repro pending): the pick depends on `isDynamicsGeodesicPick` figure-appdata
  + figure_3d routing to `OnClick`; a repaint that resets it, or a throw in `Seed`/`Draw`/`i_atom_preview`
  (e.g. the "not realisable" guard, or `tess_scout_area` failing on a vertex), leaves it armed-but-dead.

**Proposed direction.**
1. Replace the seed-pick with the **one-shot `WaveletDesignerPick`** hook (native, occlusion-correct, no
   geodesic solve, no disk) — the same mechanism `view_atom_designer` uses. This alone fixes the delay, the
   quiver-clobber, and very likely the "dies after N clicks" bug. **(Another facet of C3 — the panel
   duplicated the designer *worse*; consolidating onto the designer's pick removes this.)**
2. Add a clear **armed-mode affordance** on the cortex (not just the cursor).
3. **Extract the geodesic Region/area tool to its own control** (or out of atom design entirely) — it's a
   separate concern; a seed pick must not paint or resize a region.
4. Before/while fixing: **instrument the pick path** (log `seed, tool-state, figure-appdata, caught error`
   per click) to capture the exact death mechanism of the "stops after N clicks" bug.

**Effort/Risk:** small-medium / low (the clean pick already exists in `view_atom_designer`). **Related:** C3
(consolidate onto the designer's pick), C7 (the design UI restructure), memory `atom-designer-architecture`
(WaveletDesignerPick hook).

---

### C7 — Atom properties are hidden; restructure the panel into a resizable list/detail split
**Severity:** 🏗️ design (UX) · **Priority:** 🟠 MED · **Status:** OPEN (live-observed during walkthrough)
**Workflow ref:** §1 *Panel anatomy* (the Atom section), §3-A.
**Code:** the only atom readout is the single `jAtomInfo` label (`panel_bst_dynamics.m:173`), tucked in the
Atom section next to the differential **Measure** menu — so they read as one thing and the atom's properties
are effectively hidden under the Measure toggles.

**Observed.** No clear place shows the selected atom's properties; they're spread across the Atom-section
controls + one info label, adjacent to the unrelated Measure menu.

**Proposed direction.** Restructure the panel like the **Record panel's Events** component: a **resizable
split** — **atoms list on top**, a **live selected-atom detail pane below** that updates as the atom changes.
Detail shows the key fields: **Operator · vertex · time · scale · frequency**. (Swing `JSplitPane`, list-over-
detail; "horizontal split" in the sense of a top/bottom stack with a draggable divider.) Move the differential
**Measure** menu out of the atom readout's neighborhood so it no longer reads as an atom property. Build this
inside **`panel_atom_designer`** once C3 moves design there (so it's not built twice).

**Effort/Risk:** medium / low (GUI restructure, `gui/panel_bst_dynamics.m` `CreatePanel`; or `panel_atom_designer`
per C3). **Related:** C3 (home for the design UI), C6 (the pick affordance lives in the same panel).

---

### C8 — Consolidate Measure / Helmholtz / Threshold into one Kernel-selection menu (the registry front-end)
**Severity:** 🏗️ design (UX) + ✨ opportunity · **Priority:** 🟠 MED · **Status:** OPEN (walkthrough decision)
**Workflow ref:** §1 *Panel anatomy* (toolbar), §3-F/§3-G (Measure / Helmholtz), §3-A (Threshold).
**Code:** `OnMeasureMenu`/`OnMeasurement` (`:208/224`, sets `D.Op` → raw differential overlay);
`OnHelmholtzFiltered` (`:1251`, same kernels on the atom-filtered field); `OnThresholdMenu` (`:1704`, sets
`st.atomThreshold` for the Scout+Event level-set export).

**Observation.** Three separate toolbar actions are really **one selector plus context**: the differential
maps *are* kernels (C5/N1), so "which map" should be a **menu pick**, not three hardcoded buttons. Today
`Measure` = differential kernel on the **raw** source; `Helmholtz (filtered)` = the **same** kernels on the
**atom-filtered** source; `Threshold` = a level-set param for the marker export.

**Proposed direction — one `Kernel` / `Measurement` menu (two axes collapse + one odd-one-out):**
1. **What kernel** (the menu, grouped by N1 tier): `Magnitude` (tier-3 reduce, the default source display) ·
   `Divergence`/`Curl`/`Potential`/`Stream` (tier-1 linear — C5's `ik_div`/`ik_curl`/…) · `Energy` (tier-2
   Gram, the scalogram). **These entries are the registry's registered output kernels** — the menu is the
   N1 front-end.
2. **Which input** (raw vs atom-filtered): stops being two buttons. `Measure` and `Helmholtz (filtered)`
   collapse into **the same menu pick under different Apply/atom context** (no atom / Apply-off → raw;
   atom applied → filtered). Retire the standalone `Helmholtz (filtered)` button.
3. **Threshold is the tier-3 marker output**, not a linear kernel: keep it in the menu as the "**level-set →
   Scout+Event marker**" output mode, carrying its own 0..1 param and producing a marker (not a scalar
   overlay). Consolidated into one menu, but structurally the nonlinear tail (N1 tier-3 boundary).

**Why this is the right shape.** It makes the panel a **kernel selector + applier** — exactly N1's endgame
("the panel/process layer just selects and applies composed operators"). Selecting `Curl` = compose
`ik_curl ∘ [atom filter?] ∘ source` and display; the menu is populated from the operator/kernel registry,
so new measurements appear as new entries with no new buttons.

**Open questions:** (a) is raw-vs-filtered purely inferred from Apply state, or an explicit sub-toggle?
(b) where do bilinear (energy) and tier-3 (magnitude, marker) outputs render vs. the linear differential
overlays (colormap: `source` vs `stat2`)? (c) does the menu live on the cortex figure (per-figure `D.Op`) or
the panel?

**Effort/Risk:** medium / low-medium (UI consolidation now; grows into the registry front-end as N1 lands).
**Related:** N1 (this is its UI), C5 (differential kernels are the menu entries), C4 (raw/filtered = the
pre-filter axis), C7 (kernel menu vs. atom-detail pane must not be conflated — the C7 complaint was Measure
sitting *in* the atom readout).

---

## Backlog seed (promote to a dev-session item list)
> **N1 is the framing** — the items below each land one composed operator toward it.
- [ ] **C3 🔴** eigen-node-scoped atom designer; Dynamics → consumer. Stage 1 = de-duplicate design UI
  into `panel_atom_designer`; stage 2 = re-anchor to eigen node; stage 3 = slim Dynamics. *(med-large/med;
  subsumes C1/C2)*
- [ ] **C1 🟠** lazy-load `Cov` into appdata + cheap fail-fast check at open. *(small/low)*
- [ ] **C2 🟠** unified pinned+LRU operator cache; reuse `bst_memory`/nxr registry. *(medium/medium, after C1)*
- [ ] **C4 🟠** expose atom-filtered field as a stage-A process (sibling of `process_helmholtz_events`) →
  chain into `process_helmholtz`; panel stays the interactive explorer. *(medium/low-med; pairs with C3)*
- [ ] **C5 🟠** differential mode kernel `B_O = H_O·R` (kernelizes C4 stage B; kills the per-frame Poisson
  solve; linear outputs only). *(medium/medium; after C4)*
- [ ] **N1 ⭐** operator-composition framework: stages declare `(in-space, out-space, tier)`; compose+cache
  the product as a reusable sensors→measurement kernel. *(large; built incrementally by C1–C5)*
- [ ] **C6 🔴** replace geodesic-tool seed-pick with one-shot `WaveletDesignerPick`; armed-mode affordance;
  extract the geodesic Region/area tool out of atom design; instrument the "dies after N clicks" bug.
  *(small-med/low; pick already exists in `view_atom_designer`)*
- [ ] **C7 🟠** atom list/detail split-pane (Operator·vertex·time·scale·frequency, live), à la Record's Events;
  move Measure menu away from the atom readout. *(medium/low; build in `panel_atom_designer` per C3)*
- [ ] **C8 🟠** one Kernel-selection menu = Measure + Helmholtz + Threshold; raw/filtered from Apply context;
  entries populated from the kernel registry (N1 front-end). *(medium/low-med)*
