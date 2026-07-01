# Dynamics panel — analyze / reconstruct (sub-project C) — design

**Date:** 2026-06-30
**Status:** design (approved; to be implemented in a fresh session)
**Depends on:** A (panel cleanup) + B (frame core), both shipped + pushed to `origin/development`.
**Part of:** the Dynamics-portal program (A→B→C→D). See [[dynamics-portal-optimization]].

---

## 1. Motivation

A and B make the atom bank a designed, robust **frame**. C is the payoff: actually **run source maps
through the frame** and read the result efficiently. Analyzing a windowed source `F[nV×nT]` through the
M-member frame decomposes it into M spatial-scale bands over time; the compact, headline representation
is a **spatial scalogram** — per-member energy `e_m(t)=Σ_v W(v,t,m)²`, i.e. which spatial scales carry
the source energy, when. A **synthesis residual** measures how completely the frame captures the source
(≈0 for a tight itersine frame — a tightness-on-real-data check). C also localizes each band into a
marker atom (`JTVAtoms`).

**The reuse payoff:** B's cached projection `C = Phiᵀ·B·F` (per hemisphere, per window) is exactly what a
multi-member frame analysis needs — applying all M members is `W(:,:,m)=Phi·(g_m(λ)·C)`, so the scalogram
comes almost free from B's cache. That was the stated purpose of B's cache.

## 2. Decisions (from brainstorming)

- **Primary output = spatial scalogram** (per-member energy over time) stored as a Brainstorm
  `TimefreqMat` and viewed with the **existing timefreq viewer**, plus a **synthesis-residual** metric.
- **Window-centric compute.** The time **window** is the unit of computation. The panel window-preview is
  the default; a whole-series `process_` exists but runs **only when explicitly asked**, and even then
  iterates window-by-window (never loads the whole series).
- **Shared core reused by panel + process** (DRY), reusing B's cached projection.
- **Scope = static spatial scalogram + JTVAtoms localization.** The scalogram is the STATIC spatial-frame
  analysis (`bst_eigenwavelet('Analysis')`); dynamic ts/js atoms are excluded (consistent with B's spatial
  coverage). The unified static+dynamic (JTV spatiotemporal) scalogram is deferred.
- **Scalogram granularity = 3 rows: {Global, LH, RH}.** Compact; shows hemispheric asymmetry of scale
  content. A full per-source scalogram `[nV×nT×M]` is impractical (~460 MB/window) and is NOT stored — `W`
  is transient (in-window) for energy + JTVAtoms.
- **JTVAtoms output → a SEPARATE `dynamicsmat`** (a `dynamics_*.mat` distinct from the generator bank),
  projectable to Scouts+Events later. Keeps generators and localized markers cleanly apart.

## 3. Shared core: `bst_eigenwavelet('Scalogram', ax, frame, C)`

A new verb on `bst_eigenwavelet` (natural home — it already owns `Analysis`/`Synthesis`). Inputs: the
axes `ax` (from `bst_eigen('Axes')`: per-hemi `Phi`, `Lambda`, `Mass`, `GlobalVertices`); the `frame`
(the current bank's static members as `{g_m(λ)}` handles — reuse the ad-hoc frame B already builds from
the atoms' kernels, static-only); and the per-hemisphere projection `C` (cell, `C{h}=manifold_ft(Phi{h},
Mass{h}, Fr(gv,:))` — from B's `i_apply_projection` cache, or computed per window by the process).

Returns a struct:
- `energy` `[3 × nT × M]` — `Σ_v W(v,t,m)²` for rows {Global, LH, RH}. Per hemisphere, `W_h(:,:,m) =
  manifold_ift(Phi{h}, g_m(Lam_h).*C{h})`; global = LH+RH sum.
- `residual` `[1 × nT]` (and a scalar summary) — relative reconstruction error `‖F_modal − Frec‖ / ‖F_modal‖`
  per time sample, using the **canonical dual** per mode: `Frec_h = manifold_ift(Phi{h}, dual .* C{h})`
  with `dual(λ) = Σ_m g_m(λ)² / max(Σ_m g_m(λ)², tol)` (≈1 where the frame covers, ≈0 in a spectral gap).
  So `Frec = F` on covered modes and the residual is exactly the source energy in the frame's **gaps**:
  ≈0 for a gap-free (tight) itersine frame, >0 for a loose bank. (Normalizing by the scalar bound `A=min Σg²`
  would over-amplify interior modes and is NOT used; the per-mode dual is correct.)
- `centers` `[1 × M]` — each member's characteristic spatial scale (`√λ` gain-weighted centroid; also
  provide mm via `2π/√λ·1000`), for the TimefreqMat `Freqs` axis.
- `W` (optional out) — the transient `[nV × nT × M]` coefficient field, for JTVAtoms.

Static-only: members whose kernel domain is ts/js are excluded from `frame` (the panel's frame builder
already gathers static members only). Scalar operators (LB/LB-Connectome) only, matching Apply scope.

## 4. Panel "Analyze (window)"

A toolbar button (`OnAnalyzeWindow`) in `panel_bst_dynamics.m`. On click:
1. Resolve the current operator's `ax` and the static frame members (reuse B's `i_frame_response` gather).
2. Obtain `C` for the cursor's 4 s window — reuse B's `i_apply_projection` cache (compute if absent).
3. `scal = bst_eigenwavelet('Scalogram', ax, frame, C)`.
4. Build a **transient** scalogram `TimefreqMat` (§6) and open it with `view_timefreq` (spatial-scale ×
   time spectrogram, 3 rows). Not saved to the DB.
5. Set a residual readout in the Frame section (e.g. `residual 3.2%`).

Guarded like Apply: a real source must be linked; scalar operators only (else the "scalar-only for now"
message). Default, instant (reuses the cache).

## 5. `process_source_frame` (opt-in, whole series)

A new `toolbox/process/functions/process_source_frame.m`. Runs **only when explicitly invoked** (a
"Run on whole series…" affordance on the panel launches it, or it is run from the Process tab). It:
1. Takes a source result + the frame (passed via the dynamics table / options).
2. Iterates the series in **contiguous, non-overlapping windows** (reusing the panel's window size +
   `i_cursor_window` logic), computing `C` and `Scalogram` per window and concatenating along time.
3. Concatenates the per-window `energy` along time → one saved scalogram `TimefreqMat` (3 rows, full time
   axis) + a residual time series, registered in the DB via `db_add`.
Never loads the whole series into memory at once (per-window paging, mirroring RAW block loading).

## 6. Scalogram representation (`TimefreqMat`)

Mirrors how `bst_eigen` stores an eigen-spectrum as a TimefreqMat (Freqs = `√λ`):
- `TF` `[3 × nT × M]` (rows × time × "frequency"), `Measure='power'`.
- `RowNames = {'Global','LH','RH'}`.
- `Freqs` = the M member scale centers (store `√λ` centroid; the viewer axis reads as spatial frequency,
  labeled mm where the panel controls it).
- `Time` = the window's time vector (panel) or the full series (process).
- `Comment` = e.g. `Frame scalogram (itersine×6, Laplace-Beltrami)`.
- `DataFile`/provenance point to the analyzed source result.
Viewed via the existing `view_timefreq` (no new viewer).

## 7. JTVAtoms localization → separate `dynamicsmat`

A "Localize bands" action (`OnLocalizeBands`): `T = bst_eigenwavelet('JTVAtoms', W, ax, thr)` localizes
each band `W(:,:,m)` into an atom — peak-energy time `iRef`, seed = peak vertex at `iRef`, time window =
where band energy ≥ `thr·max`, region = cortical level set at `iRef`. Returns a `dynamicsmat`. Save it as
a **separate** `dynamics_*.mat` (not the generator bank) via `bst_dynamics('Save')`; openable via
`view_dynamics`; projectable to Scouts+Events later. `thr` reuses the panel's level-set threshold
(`st.atomThreshold`, default 0.5). `W` is the transient field from step 3 (recompute if not held).

## 8. Components & files
- **New verb:** `bst_eigenwavelet('Scalogram', …)` (+ tests).
- **New process:** `toolbox/process/functions/process_source_frame.m`.
- **Panel (`panel_bst_dynamics.m`):** `OnAnalyzeWindow`, `OnLocalizeBands`, a transient-TimefreqMat
  builder helper, a residual readout in the Frame section; an "Analyze (window)" toolbar button + a
  "Localize bands" action + a "Run on whole series…" affordance.
- **Reuse:** B's `i_apply_projection` cache + `i_frame_response` frame gather; `bst_eigenwavelet`
  `Analysis`/`Synthesis`; `view_timefreq`; `bst_dynamics`/`view_dynamics`; `db_add`.

## 9. Testing
- **Headless (pure, controller-run):** `Scalogram` on a synthetic field + a known tight frame → residual
  ≈0 and Parseval energy conservation (`Σ_m energy ≈ ‖F‖²` up to the frame constant); on a deliberately
  non-tight bank → residual >0. `JTVAtoms` on a synthetic W with a known per-band peak → one atom/band at
  the expected seed vertex + time window.
- **Live (controller, MCP):** on the sub-MTL0002 Dirac→LB itersine frame — Analyze (window) opens the
  3-row scale×time spectrogram, residual ≈0 (tight); Localize bands → a dynamics_*.mat with M atoms,
  openable; screenshots. Whole-series process on a short window range → saved scalogram TimefreqMat.

## 10. Out of scope (→ D / later)
- Unified static+dynamic (JTV spatiotemporal) scalogram; per-source or per-scout scalogram granularity.
- Dirac/vector real-source analysis + the filtered-**sensor** view (D).
- Projecting the localized `dynamicsmat` to Scouts+Events (a later step; C just saves the table).
