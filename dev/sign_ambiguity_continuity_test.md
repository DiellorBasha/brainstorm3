# Connection-eigenmode current representation: sulcal-wall continuity test

**Status:** design (brainstorming output) — 2026-06-05
**Branch (proposed):** `feature/sign-ambiguity-continuity`
**Depends on:** `ConnEigenmodes` axis (M2), `bst_conn_phase` / `bst_tangent_face2vertex` (M3 Plan B), nxr `vertexFrame` (M3 Plan A)

---

## 1. Claim under test

The connection-Laplacian (Fiedler) frame does **not** correct the constrained source sign — it **replaces** a representation with an inherent sign ambiguity (the signed normal projection `s = J·n̂`) with one that encodes current direction as a **continuous, gauge-consistent phase** on the tangent bundle (`z = J·ê₁ + i·J·ê₂`). The payoff is *continuity across cortical folds*, and it follows from the **smoothness of the Fiedler frame** — specifically that the smoothest tangent field co-rotates with curvature and is therefore *decoupled from the folding* that makes the constrained normal flip.

This is **not** a claim that the connection eigenmodes fix the unconstrained MNE's sign: the unconstrained estimate has no sign ambiguity. The standard unconstrained MNE lives in the **global xyz** frame and is displayed as a directionless norm `|J|`. The connection-eigenmode representation's advantage over that standard is an **anatomically-referenced, gauge-consistent, continuous phase** where the alternatives give either a directionless magnitude or direction in arbitrary world axes.

| Representation | Sign ambiguity? | Direction info? | Reference frame |
|---|---|---|---|
| Constrained `s = J·n̂` | **Yes** (flips across folds) | in/out only | normal `n̂` |
| Unconstrained, norm `\|J\|` | No | **none** (discarded) | — |
| Unconstrained, xyz vector | No | yes, but arbitrary world axes | global xyz |
| Connection-eigenmode `z` | No | yes, **continuous phase** | Fiedler frame |

The test is built in two components. **Component 1 is data-free** (geometry only) and is the conceptual core; **Component 2 uses the real MEG vector field** and is the applied confirmation.

---

## 2. Component 1 — data-free frame-smoothness, partitioned by SulciMap

**Purpose.** Show that the Fiedler tangent field's intrinsic variation is *decoupled from cortical folding*: the quantity that flips the constrained sign (the normal `n̂`) turns sharply at sulci, while the Fiedler field's covariant variation is flat across the sulcus/crown boundary, with its energy concentrated only at the ~2 topological singularities — nowhere near the sulci. This component needs only a surface carrying a `ConnEigenmodes` axis; no kernel or recordings.

**Inputs.** A DB-registered cortex surface (the test uses the real 20484V cortex; SKIP if absent, matching the benchmark idiom).

**Geometry & fields.**
- `TessMat = in_tess_bst(SurfaceFile)`; `Vtx`, `Fcs = double(Faces)`, `Nv = VertNormals`, `VertConn`.
- `SulciMap`: use `TessMat.SulciMap` if present, else compute via `tess_sulcimap(TessMat)`. Binary per-vertex label (1 = sulcal/negative-curvature vertex, 0 = gyral crown).
- `Curvature = tess_curvature(Vtx, VertConn, Nv, 0.1)` (folding magnitude reference).
- Fiedler field: `ConnEig = bst_conn_eigenmodes_ensure(SurfaceFile)`; `mctx = nxr.manifold.context(Vtx, Fcs)`; `vFrame = nxr.manifold.measure.vertexFrame(mctx)`; `R = bst_conn_phase(ConnEig, vFrame, 'Rank', 1, 'nSing', 2)`. `R.Field` is the gauge-independent 3D Fiedler tangent field `[nV×3]`; define `ê₁ = R.Field ./ |R.Field|`, `ê₂ = n̂ × ê₁` (the Fiedler frame). `R.Singularities` gives the singularity vertices per component.

**Per-edge variation metrics** (over mesh edges `(i,j)` from `VertConn`, within a connected component):
- **Normal angular variation** `δn(i,j) = acos(clamp(n̂ᵢ·n̂ⱼ))` — geometric folding density; high at sulci.
- **Fiedler covariant angular variation** `δf(i,j)` — the angle between `ê₁(i)` and `ê₁(j)` *after parallel transport along the edge via the connection*, so the comparison is gauge-correct (a raw frame-to-frame angle would conflate the arbitrary gauge with real rotation). Parallel transport uses the same nxr connection that defines the eigenmodes; obtain the per-edge transport/rotation from the connection operator used in `tess_connection_laplacian` (the off-diagonal `K(i,j)` phase encodes the Levi-Civita transport). If a clean per-edge transport accessor is not directly available, derive it from the connection Laplacian's off-diagonal complex argument — document the source explicitly in the code.

**The contrast (the deliverable statistic).** Partition edges by whether they touch a sulcal vertex (`SulciMap`) and compare distributions:
- `δn` (and the constrained scalar's spatial variation, see below) is **elevated on sulcal edges** vs. crown edges.
- `δf` is **statistically flat** across the sulcal/crown partition (no elevation at sulci), and its largest values localize at `R.Singularities`, not at sulci.

Report: median/IQR of `δn` and `δf` by partition; correlation of `δf` vs. `δn` (expected near zero) vs. correlation of `δn` vs. `Curvature` (expected high); fraction of total `δf` energy within a small geodesic radius of the singularities.

**Illustrative companion — sulcal-crossing profiles.** Using `SulciMap` connected components (label via `VertConn`), pick ~3 representative sulci. For each, trace a short geodesic path crossing the sulcus (crown → wall → fundus → wall → crown) and plot along arc length: (a) the normal's turning angle (rotates ~180°), (b) the Fiedler frame `ê₁` turning angle (smooth, no jump), and — once Component 2's `J` is available — (c) the constrained scalar `s` (flips sign at the fundus). Component 1's figure uses (a)+(b) only; (c) is added when run together with Component 2.

---

## 3. Component 2 — real-data sulcal-wall current continuity

**Purpose.** On the actual unconstrained MEG current at the auditory M100, show that the constrained scalar flips sign across facing sulcal walls while the Fiedler-frame phase stays continuous — conditioned on the physical 3D current genuinely being continuous there — and that the Fiedler frame does this where the two competing frames (global-xyz, tess_tangents) do not.

**Vector field.**
- Unconstrained kernel results file (auditory protocol, Subject01 deviant study). Apply via Brainstorm's standard path so channel selection is correct: `Results = in_bst_results(ResultsFile, 1)` → `Results.ImageGridAmp` is `[3·nV × nTime]`, `nComponents = 3`.
- Reshape to `J(v,:,t)` `[nV×3×nTime]` (component-major per Brainstorm's `[x1 y1 z1 x2 y2 z2 …]` interleave — verify against `Results.nComponents` and the `GridLoc`/`Atlas` ordering before trusting the reshape).
- Peak selection: compute global field power `sqrt(sum(|J|², over all v,components))` per time sample within the auditory window (e.g. 0.05–0.15 s) and take the argmax sample `t*`. **Not** a hardcoded 100 ms. `J* = J(:,:,t*)` `[nV×3]`.

**Frames at each vertex** (all `[nV×3]`):
- **Fiedler:** `ê₁ = R.Field./|R.Field|`, `ê₂ = n̂ × ê₁` (from Component 1).
- **tess_tangents (incumbent):** `[Uf,~] = tess_tangents(SurfaceFile,'NoSave',1)`; `[Uv,Vv] = bst_tangent_face2vertex(Fcs, Uf, Nv)`; frame `(Uv, Vv)`.
- **global-xyz (= standard unconstrained MNE):** fixed `[1 0 0]`, `[0 1 0]` projected to the tangent plane per vertex (or report phase in raw xyz — document choice; the point is that it is not surface-adapted).
- Constrained scalar (artifact exhibit): `s = sum(J* .* n̂, 2)`.

For each frame, the tangential complex coordinate is `z = (J*·frame.e1) + i·(J*·frame.e2)`; phase `θ = atan2(J*·e2, J*·e1)`; tangential magnitude `|z|`.

**Sulcal-wall pair detection.**
1. Restrict candidate vertices to `SulciMap == 1`.
2. Among those, find pairs that are (a) close in 3D — Euclidean distance < 3 mm; (b) anti-aligned normals — `n̂ᵢ·n̂ⱼ < −0.7`; (c) **not** mesh-adjacent — `j` outside `i`'s N-ring (`N = 3`) via `VertConn` powers, so they are genuinely across a fold and not neighbors on the same wall.
   - Implementation: `knnsearch` / `rangesearch` on the restricted vertex coordinates for the 3 mm radius, then filter by the normal and N-ring criteria.

**Conditioning on physical continuity (rigor step).** Keep only wall pairs where the underlying 3D current is genuinely continuous: `J*ᵢ·J*ⱼ > 0` **and** magnitudes comparable (`min(|J*ᵢ|,|J*ⱼ|)/max(…) > 0.5`). This restricts the test to pairs where the physics is continuous, so we only ask which representation *preserves* that continuity — guarding against the "you just smoothed real structure away" critique.

**Metrics over the conditioned pairs.**
- **Constrained sign-flip rate:** fraction with `sign(sᵢ) ≠ sign(sⱼ)`. Expected high (the artifact).
- **Tangential magnitude discontinuity:** relative difference `||z|ᵢ − |z|ⱼ| / max(|z|ᵢ,|z|ⱼ)` (frame-independent; same for all frames since `|z| = |J_tangential|`). Expected low — establishes the magnitude is continuous.
- **Phase discontinuity per frame:** circular distance `wrap(|θᵢ − θⱼ|)` to `[0, π]`, computed separately under Fiedler, tess_tangents, and global-xyz. Expected: Fiedler lowest, global-xyz highest/meaningless, tess_tangents intermediate (degraded near its FreeSurfer-pole singularities).

---

## 4. File structure

All under `dev/benchmarks/sign_ambiguity/` (mirrors `dev/benchmarks/` script+report idiom: a driver, focused helpers, a figures/stats output dir).

| File | Responsibility |
|---|---|
| `sign_ambiguity_run.m` | Driver. Resolves surface + kernel, runs Component 1 always, Component 2 if a usable unconstrained kernel is found; writes figures + a summary `stats.md`/`stats.csv`. SKIPs cleanly (prints `SKIP:` and returns) when the surface or kernel is absent. |
| `sa_frames.m` | `[Fiedler, FsFrame, R, vFrame] = sa_frames(SurfaceFile)` — builds the three (Fiedler / tess_tangents / global-xyz handled inline) frames + Fiedler decode; shared by both components. |
| `sa_sulcal_walls.m` | `pairs = sa_sulcal_walls(Vtx, Nv, SulciMap, VertConn, opts)` — sulcal-restricted facing-wall pair detection (3 mm, anti-normal, N-ring exclusion). |
| `sa_smoothness.m` | Component 1 metrics: per-edge `δn`, `δf` (gauge-correct covariant), SulciMap partition stats, singularity-energy fraction. |
| `sa_crossing_profile.m` | Component 1 illustrative profiles: surface shortest path between facing-wall pairs (descends through the fundus); cumulative normal vs. Fiedler-frame turning along arc length, plus the constrained scalar overlay when current is supplied. |
| `sa_continuity.m` | Component 2 metrics: applies kernel, peak selection, per-frame phase, conditioning, pair metrics. |
| `sa_figures.m` | Renders: (C1) `δn` vs `δf` by SulciMap partition; (C1) sulcal-crossing profiles; (C2) per-frame phase-discontinuity histogram + summary bars. |

**Reused, not reimplemented:** `bst_conn_eigenmodes_ensure`, `bst_conn_phase`, `bst_tangent_face2vertex`, `tess_tangents`, `tess_sulcimap`, `tess_curvature`, `tess_normals`, `in_bst_results`, `in_tess_bst`, `nxr.manifold.context` / `vertexFrame`.

---

## 5. Deliverables

- `dev/benchmarks/sign_ambiguity/sign_ambiguity_run` output dir: `figures/` (PNG + FIG) and `stats.md` + `stats.csv`.
- **Component 1 result statement:** Fiedler covariant variation is uncorrelated with SulciMap (flat across the partition) while normal/curvature variation is elevated at sulci; Fiedler energy concentrates at singularities.
- **Component 2 result statement:** across physically-continuous facing-wall pairs, the constrained map flips sign at a high rate while tangential magnitude is continuous and the Fiedler-frame phase discontinuity is lower than tess_tangents and far lower than global-xyz.

---

## 6. Testing

This is an analysis/benchmark, not a unit-tested library addition, but each helper gets a focused smoke test under `dev/tests/`:
- `test_sa_sulcal_walls.m` — on the real cortex, returns a non-empty set of pairs; every pair satisfies the 3 mm / anti-normal / non-adjacent / both-sulcal criteria. SKIP if no 20484V cortex.
- `test_sa_smoothness.m` — `δf` median on sulcal edges is not significantly greater than on crown edges (the decoupling claim); `δn` median on sulcal edges *is* greater. SKIP if no cortex.
- `test_sa_continuity.m` — with the unconstrained kernel present, the pipeline runs end to end and the conditioned pair set is non-empty; constrained sign-flip rate > Fiedler-phase-discontinuity-normalized rate. SKIP if no kernel.
- Driver `sign_ambiguity_run.m` is exercised by running it (produces figures) — verified manually, not asserted.

All tests follow the existing idiom: `brainstorm nogui` if needed, install/load `nxr-compute`, locate the real cortex, `fprintf('SKIP: …')` and return when prerequisites are absent.

---

## 7. Out of scope (explicitly not built)

- Approach 2 (synthetic ground-truth forward simulation) — deferred; would be the rigor backstop if reviewers challenge "continuity ≠ over-smoothing."
- Approach 3 (spectral compactness curves) — deferred.
- Any change to the inverse modeling, the kernels, or the `ConnEigenmodes` axis itself.
- Between-subject / group analysis (this is within-subject, intrinsic).
- A GUI viewer for these results (the existing `view_connection_phase` already shows the phase interactively; this is a quantitative offline analysis).

---

## 8. Open implementation risk to resolve during build

The **gauge-correct per-edge parallel transport** for `δf` (Component 1) is the one piece without a ready-made accessor. The plan is to read the transport rotation from the off-diagonal complex argument of the connection Laplacian `K` (built by `tess_connection_laplacian`), which *is* the Levi-Civita connection nxr uses. The implementer must confirm the sign/orientation convention of `arg(K(i,j))` against `vertexFrame` before trusting `δf`; a quick check is that a flat field (constant in the transported frame) yields `δf ≈ 0` everywhere. If the convention can't be pinned down cleanly, fall back to measuring `δf` as the angle between `R.Field(i)` and `R.Field(j)` projected into a common tangent plane (extrinsic, slightly conflates curvature) and document the approximation.
