# Dynamics — Remaining Work RESUME STATE (2026-07-02)

> **✅ STATUS 2026-07-02: BOTH ITEMS DONE + LIVE-VERIFIED, committed `c9b7b5fa` on `development`.**
> - **Item 1 (Task 6 — Helmholtz filtered):** `OnHelmholtzFiltered` + toolbar button in
>   `panel_bst_dynamics.m`. Live: filtered Curl differs from raw (corr 0.14, |Curl| 384 vs 195),
>   Div/Curl finite, HarmFrac 0.089.
> - **Item 2 (windowed scalogram):** user chose **A (Morlet full-time)**. `bst_eigenwavelet('ScalogramEnergy')`
>   (coeff-space, per-hemi Euclidean Gram, no `W`; headless equivalence to `Scalogram` = 4e-16) +
>   `OnAnalyzeWindow` now full-recording coeff-space (`i_project_fulltime` + `i_hemi_membership`). Live:
>   62 s in 0.16 s, Time == recording (co-displays, no conflict), single-block LH/RH split fixed
>   (LH 3.56e3 / RH 4.49e3, Global=LH+RH to 1e-10).
> - Cosmetic follow-up (not a defect): the timefreq viewer labels the spatial scale-center as "Hz" and a
>   1-member frame sits off the auto-scaled axis → design a multi-member frame for a filled scalogram.
>
> ---
>
> **Purpose:** hand off two remaining implementation items to a FRESH session. Read this top-to-bottom
> before touching code. Baseline is fully committed on `development`; both items are additive.
>
> **The two items:**
> 1. **Task 6 — Helmholtz/differential on the Dirac-Connectome-filtered field** (last item of the
>    Dirac-Connectome plan `dev/2026-07-02-dirac-connectome-plan.md`).
> 2. **Windowed-scalogram display** — resolve the time-base conflict; compute energies in coefficient
>    space (Morlet-style full-time vs PSD-style spectrum — a UX decision for the user).

---

## 0. How to resume (verified recipe, this exact session)

Environment: **Mac Mini, 17 GB RAM**, protocol **`preventad`**, dev repo
`~/workspace/research/code/brainstorm3`. MATLAB MCP session with GUI. **Never use `clear`** (wipes
GlobalData, hangs the session — use `rehash`; edited `.m` auto-reload).

Boot + open a live Dirac-dSPM Dynamics session on the known test source:

```matlab
if ~brainstorm('status'), brainstorm; end          % full GUI (nogui can't host the Dynamics panel)
srcFile = ['link|sub-MTL0002/@rawsub-MTL0002_ses-02_task-rest_run-01_meg_resample_notch_high/results_DiracEig_KERNEL_260625_2215.mat|' ...
           'sub-MTL0002/@rawsub-MTL0002_ses-02_task-rest_run-01_meg_resample_notch_high/data_0raw_sub-MTL0002_ses-02_task-rest_run-01_meg_resample_notch_high.mat'];
hFig = view_dynamics('FromResult', srcFile);        % opens cortex + Dynamics panel; sets DynamicsTarget/Overlay
% Build a Dirac-Connectome atom with a STATIC kernel (the default 'diffusion' is DYNAMIC and Apply bails):
panel_bst_dynamics('OnCreateAtom');
panel_bst_dynamics('OnSetOperator','Dirac-Connectome');
ctrl = bst_get('PanelControls','Dynamics');  ctrl.jKernel.setSelectedIndex(3);  % index 3 = 'heat' (static)
panel_bst_dynamics('OnKernelChange');
% Place a left-hemi seed:
st = getappdata(0,'DynamicsTarget');  Surf = in_tess_bst(st.T.SurfaceFile,0);  [ir,il] = tess_hemisplit(Surf);
seed = il(round(numel(il)/2));  ia = st.curAtom;
st.T.Groups(ia).vertices = seed;  st.atomSeed = seed;  setappdata(0,'DynamicsTarget',st);
panel_bst_dynamics('i_atom_apply');                 % should paint fiber-spread cortex + quivers, no sensor
```

**Gotchas (all learned the hard way this session):**
- The static-kernel list (no `domain` field in `bst_eigfilter_kernel('info',k)`): heat(3), diffgauss(1),
  flat(2), ideal(4), inverse_heat(5), log(6), matern(7), mexhat(8), power(9), tikhonov(10). Dynamic
  (rejected by Apply): diffusion(0), dampedwave(11), kleingordon(12), wave(13), gabor(14), resonator(15),
  stmatern(16), travwave(17).
- **`i_cursor_window` uses the GLOBAL time cursor.** Fresh session parks it at t=0, so windows read
  `[0,4]s`. To analyze elsewhere, set the cursor first (`panel_time('SetCurrentTime', 30)` — but note this
  fires a figure repaint; if the figure is half-torn it can error; reopen the session if so).
- The Dynamics session tears down (`DynamicsTarget` empties) if an action throws hard. Just re-run the
  reopen recipe.
- **"Disk full, disconnected or read-only" is NEVER disk** — it's `bst_save.m:83`'s hardcoded catch-all
  masking the real error (was OOM before the mode-kernel fix). 584 GB free; the volume writes fine.

---

## 1. Baseline already committed (do NOT redo)

Dirac-Connectome operator (Tasks 1–5) + mode-kernel perf, on `development`:

| commit | what |
|---|---|
| `2a0dbe47` | `bst_lift_connectome_dirac` (lift scalar eigenbasis → quaternion) |
| `16a37bbe` | Dirac-Connectome axes (`i_atom_axes`) + `bst_eigenfilter('RowMap')` |
| `870ecacb` | guard empty GlobalVertices block (whole-brain single-block ax) |
| `52b7c66e` | 5th operator "Dirac (connectome)" + 5-way `i_gate_mask` |
| `21f7ddb3` | Dirac-Connectome Apply (fiber cortex + quivers, no sensor) |
| `68aae27a` | scalogram + localize admit Dirac-Connectome |
| `d6380703` | `i_apply_projection` dispatch on `ax.Variant` (bug from live verify) |
| `7ee94679` | **mode kernel** `A=Φᵀ M P K` (12.6 GB→~1 GB, 500×); `i_dirac_recon_display` (chunked); save-robustness |
| `1514ff26` | frame-based inverse research roadmap (separate, strategic) |

**Live-verified:** fiber-spread atom crosses hemispheres (ipsi 76% / contra 24% vs surface-Dirac 100%/0%);
Apply paints cortex+quivers, no sensor; mode-kernel coeffs identical to reconstruct-then-project to 1.8e-14;
Apply peak RAM ~1.2 GB; Analyze computes + **saves** (session no longer OOM-crashes).

**Key functions (panel_bst_dynamics.m) for the remaining work:**
- `i_vector_modekernel(st,ax,D)` → `{A_h}` `[K×nCh]` cached; `i_vector_coeffs(st,ax,D,iWin)` → `cCell` via
  `A·data` (reconstruct-then-project fallback for kernel-free sources).
- `i_dirac_recon_display(ax, cCell, midCol)` → `[mag(nV×nT), V3mid(nV×3)]` (time-blocked; mid frame only).
- `i_apply_projection(st,ax,D,iWin,nV)` — dispatches on `ax.Variant`; caches by `(srcResult,iWin,Variant)`.
- `i_atom_apply` — vector branch (Dirac + Dirac-Connectome); `OnAnalyzeWindow` / `OnLocalizeBands`.
- `D = getappdata(hFig,'DynamicsOverlay')` carries **`D.Cov`** (surface Covariant operator, from
  `view_dynamics FromResult`) — needed for Helmholtz.

---

## 2. ITEM 1 — Task 6: Helmholtz on the Dirac-Connectome-filtered field

**Goal:** from a Dirac-Connectome (or Dirac) atom, take the filtered current field at the **cursor frame**,
run the Hodge–Helmholtz decomposition, and paint/report divergence & curl.

**Spec (from `dev/2026-07-02-dirac-connectome-plan.md` Task 6), UPDATED for the mode kernel:**

Add `OnHelmholtzFiltered()` in `panel_bst_dynamics.m`:
1. `[ctrl, st] = i_cs()`; `D = getappdata(st.hFig,'DynamicsOverlay')`; `variant = i_atom_op(st)`.
   Gate: only `Dirac` / `Dirac-Connectome` (both have a 3-vector field). Bail with an info message otherwise.
2. `ax = i_atom_axes(st, variant)`; `[kernel,kp] = i_selected_generator(st,ctrl,max(ax.Lambda{1}(:)))`;
   `g = bst_eigfilter_kernel(kernel,kp)` — require a **static** kernel (mirror the Apply domain check).
3. `iWin = i_cursor_window(D.srcDS, D.srcResult, 4)`; compute the cursor's column within the window
   (`midCol`) — reuse whatever mapping the overlay uses (`i_dynamics_overlay` maps the global cursor →
   window column; factor that out or replicate).
4. **Cheap path (use the mode kernel, one frame only):**
   ```matlab
   cCell = i_vector_coeffs(st, ax, D, iWin);                 % A·data, tiny
   cf = cCell; for h=1:numel(cf), if ~isempty(cf{h}), cf{h} = g(ax.Lambda{h}(:)).*cf{h}; end, end
   [~, V3mid] = i_dirac_recon_display(ax, cf, midCol);       % [nV×3], ONE frame
   V3col = reshape(V3mid.', [], 1);                          % [3nV×1] interleaved x,y,z
   Ht = process_helmholtz('Compute', V3col, D.Cov);         % CONFIRM signature + return fields
   ```
5. Paint/report. `view_dynamics` already has a `PickScalar` verb with cases
   `Divergence→Ht.Div`, `Curl→Ht.Curl`, `Potential→Ht.Phi`, `Stream→Ht.Psi` (see `view_dynamics.m` ~L266).
   Reuse that to paint the chosen scalar on the cortex; report `Ht.Div`/`Ht.Curl` finite in the info line.
6. **Wire a toolbar/menu entry** "Helmholtz (filtered)" near the Analyze/Localize buttons (built in
   `CreatePanel`). Follow the existing `gui_component('button',…,@(h,e)bst_call(@OnHelmholtzFiltered))` idiom.

**Before coding, CONFIRM:**
- `process_helmholtz('Compute', V3col, Cov)` exact input shape (`[3nV×1]` interleaved vs `[nV×3]`) and its
  return struct field names (`Div/Curl/Phi/Psi/Vsol`?). Check `process_helmholtz.m` + the
  `helmholtz-events-process` / `bst-helmholtz-consolidation` memories (bst_helmholtz was consolidated INTO
  `process_helmholtz('Compute')` + `bst_divergence`/`bst_curl`).
- `D.Cov` is populated (it is, from `view_dynamics FromResult`).
- Manifold normals gotcha: manifold face normals are inward-signed; Div/Phi/source-sink flip sign under
  orientation, Curl/Psi/Vsol are sign-invariant (see `manifold-face-normals-inward` memory) — don't be
  alarmed by a sign, but note it.

**Verify (live):** Dirac-Connectome atom → Helmholtz (filtered) → a Div/Curl map that DIFFERS from the raw
-source Helmholtz (it's the fiber-spread filtered field); `Ht.Div`/`Ht.Curl` finite; screenshot. **Commit:**
`feat(dynamics): Helmholtz/differential on the Dirac-Connectome-filtered field`.

**Risk:** low-medium. The field build is cheap (mode kernel + 1 frame). Main unknowns are the
`process_helmholtz` signature and the toolbar wiring — both are mechanical once confirmed.

---

## 3. ITEM 2 — Windowed-scalogram display (resolve the time-base conflict)

**The problem (fully diagnosed this session):** `OnAnalyzeWindow` saves a `timefreq` whose `Time` is the
**actual** window times (e.g. `[30,34]s` at cursor 30 — NOT rebased; verified). Brainstorm enforces ONE
global time base, so a 4 s-window timefreq can't be co-displayed with the loaded 62 s recording →
`view_timefreq` throws "Time definition … not compatible". Removing `DataFile` did NOT fix it (it's the
global-time check, not the link). The **save itself works**; only the auto-display conflicts.

**What we currently save:** `FileMat.TF = scal.energy [3 × nT × M]` (Global/LH/RH × time × scale) — already
just the energies, not the field. Good.

**How Brainstorm handles its own windowed transforms (`bst_timefreq.m`):**
- Source TF is **kernel-based** — feeds `ImagingKernel` into the transform, keeps TF in mode/kernel space,
  reconstructs per-vertex only on display (never materializes `[nVertex×nTime×nFreq]`). Same trick as our
  mode kernel.
- **PSD collapses time** to `[first last]` (2 points → a spectrum, no time axis to conflict); **Morlet keeps
  full time** (`FileMat.Time = OPTIONS.TimeVector`).

**Cost of doing the scalogram "like bst_timefreq" (benchmarked, exact sizes, this session):**
- Coeffs over the WHOLE 62 s recording `c=A·data [180×74381]`: **0.27 s / 107 MB**.
- Energies coeff-space `[3×N×8]`: **0.08 s / 14 MB**. → **total ~300 MB, a few seconds.**
- vs per-vertex over full recording `[nV×N×M]` = **97.5 GB** (12.2 GB for 1 scale) — prohibitive; avoid,
  reconstruct on demand.

**DECISION FOR THE USER (present both, then implement the chosen one):**
- **(A) Morlet-style — full-time-resolved energies over the recording's own time base.** Compute
  `E[3×N×M]` in coefficient space spanning the SAME time as what's loaded → `Time` matches → **no conflict**,
  keeps full time localization. Recommended if time-resolution matters (it does, for the alpha-vortex work).
- **(B) PSD-style — per-scale spectrum.** Collapse time → `[3×2×M]` (window bounds) → a spectrum,
  trivially cheap, auto-displays, but loses the time course.
- **(C) Save-only + open via Navigator.** Keep the window scalogram; don't auto-open (kill the popup); user
  opens it deliberately. Smallest change; keeps info; worse UX.

**Implementation notes (for A, the recommended path):**
- Replace `bst_eigenwavelet('Scalogram')`'s per-vertex reconstruction with **coefficient-space energy**:
  total energy per scale/time = `‖g_m(λ)·c(:,t)‖²` (Parseval, orthonormal Φ). For the **Global/LH/RH** rows
  precompute per-hemi **Gram** matrices `G_h = Φ(hemi-rows)ᵀ M(hemi-rows) Φ(hemi-rows)` `[K×K]` once, then
  `E_h(t) = cf(:,t)ᵀ G_h cf(:,t)` (whole-brain single block → split the quaternion rows by LH/RH vertex sets
  via `tess_hemisplit`). This removes the last per-vertex materialization from Analyze entirely.
- Decide the analysis extent: cursor window vs the full loaded recording vs a time selection. If you want
  co-display with the recording, `Time` must equal the loaded time — so either analyze the full loaded
  extent, or set the global window to the analysis window before display.
- Keep `scal.W` (per-vertex `[nV×nT]`) ONLY for Localize (`OnLocalizeBands`/`JTVAtoms`), and reconstruct it
  time-blocked (via `i_dirac_recon_display`-style chunking) or on demand — never store it full over long
  recordings.
- `i_scalogram_timefreq` already saves standalone (no DataFile) — keep that.

**Verify:** Analyze on the alpha segment → energies display without the popup; time axis is the true window
times; memory flat. **Commit:** `feat(dynamics): coefficient-space windowed scalogram (<chosen presentation>)`.

**Risk:** medium — touches `bst_eigenwavelet('Scalogram')` (or adds a coeff-space variant) and the
LH/RH Gram; the UX decision (A/B/C) must be settled with the user first.

---

## 4. Test targets & validation assets

- **Dirac-dSPM source:** the `sub-MTL0002` link above (unconstrained, nComp=3, ImagingKernelMode + eigen
  file). 276 such results exist in `preventad`.
- **Alpha vortex segment** (for item 2 / roadmap Phase 5): `S01_AEF_01_notch`, 20–25 s, 7–13 Hz, ~10.55 Hz
  right parieto-occipital burst ~22.6 s (different protocol/subject — `alpha-example-block` memory).
- **Synthetic-on-real-cortex harness:** `dev/benchmarks`.

## 5. Related docs
- `dev/2026-07-02-dirac-connectome-plan.md` — the original 6-task plan (Task 6 spec).
- `dev/2026-07-02-frame-based-inverse-roadmap.md` — the strategic inverse roadmap (separate, longer horizon).
- Memory: `dirac-connectome-operator` (has the mode-kernel + "disk-full was never disk" + open items).
