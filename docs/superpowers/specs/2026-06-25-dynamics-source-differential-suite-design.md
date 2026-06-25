# Dynamics "Source" differential suite — design (Spec 2)

- **Date:** 2026-06-25
- **Status:** Draft for review
- **Scope:** Spec 2 of the differential-Helmholtz line. Spec 1 (merged) made
  `process_helmholtz('Compute')` the stateless flat-covariant algorithm. This spec wires the
  **"Source" axis** of the dynamics panel: an ephemeral, per-frame differential overlay on the
  cortex (Divergence / Curl / Potential / Stream), folded into `panel_bst_dynamics`, retiring
  `view_helmholtz` / `panel_helmholtz`.

## 1. Problem / motivation

The Helmholtz visualization is trapped in a standalone `view_helmholtz` figure with its own
controls panel (`panel_helmholtz`), separate from the dynamics panel that owns the
navigate / detect / save workflow and the atom model. `panel_bst_dynamics` currently *drives*
`view_helmholtz` through verbs (`SetComponent`/`SetSmoothing`/`UpdateFrame`) and reads its
per-frame cache — an awkward split. We want the differential maps to be a **first-class,
dynamics-native overlay**: the "Source" axis of the atoms/dynamics system.

The mechanism for this already exists and is idiomatic: `figure_3d` maps source results to the
cortex over time (free time-nav from the source-mapping kernel) and renders a scalar colormap +
raw source-vector quiver natively; `panel_surface('FireCustomOverlay')` fires a per-figure
**`CustomOverlayFcn`** at the end of every per-frame update — the canonical insertion point for
an **intermediate, ephemeral transform** that overwrites the displayed scalar without saving.
The differential step is exactly that transform.

## 2. Goals

- A **dynamics-owned `CustomOverlayFcn`** on a standard `figure_3d` source figure computes the
  per-frame differential overlay via `process_helmholtz('Compute')`; nothing is saved unless the
  panel commits an atom.
- A **generalized differential-operator selector** in `panel_bst_dynamics` with four entries —
  **Divergence, Curl, Potential, Stream** — each a field of one `Compute` result.
- The dynamics atom "Source" recording (`OnRecord`/`OnSaveCursor`) reads the selected
  differential field from the overlay cache.
- **`view_helmholtz` and `panel_helmholtz` are deleted**; the tree menu redirects to
  `view_dynamics`.

## 3. Non-goals (explicitly deferred)

- **Eigenmode smoothing.** The `view_helmholtz` "Smooth" filters are NOT carried over. They
  become the **eigenvalue axis** of the dynamics panel + atoms, which needs its own
  navigate/detect/save state analysis — a separate later step.
- **A purpose-built vortex / core detector.** The retired stack (`dev/experimental/`) is being
  rebuilt; the existing `bst_dynamics('Extrema')` recording continues to work, unchanged.
- **`process_poisson`.** Not built. `Compute` already returns `Phi`/`Psi`; a composable
  Poisson-on-arbitrary-fields operator (and the `bst_operators` `'poisson'` Method seam) is for
  when the selector broadens beyond Helmholtz.
- **Isolines and other descriptors.** Later additions to the dynamics panel.
- **Porting anything from `view_helmholtz`** — no irrotational/solenoidal model, no `Vtot`
  quiver override. The overlay is written fresh; the raw source-vector quiver stays native to
  `figure_3d`.

## 4. Established-pattern anchors (verified)

- `figure_3d` time pipeline: `bst_figures('FireCurrentTimeChanged')` → `panel_surface(
  'UpdateSurfaceData')` (writes the RMS scalar to `TessInfo.Data`) → `UpdateSurfaceColormap` →
  `figure_3d('UpdateSurfaceColor')` → `PlotSourceVectors` → `panel_surface('FireCustomOverlay')`.
- `CustomOverlayFcn` (`panel_surface.m:2082/2088`): `setappdata(hFig,'CustomOverlayFcn',@(h)f(h))`;
  fires per frame after the standard scalar+colormap; the place to overwrite `TessInfo.Data`.
- Unconstrained (`nComponents==3`) results render the RMS-norm scalar + raw 3-vector quiver
  automatically; the overlay substitutes its own scalar.
- `process_helmholtz('Compute', J, Cov)` returns `Div, Curl, Phi, Psi` (+ Fmag/Hmag/Virr/Vsol/
  Vtot/Hresid/HarmFrac) from one coupled solve.

## 5. Architecture

```
view_dynamics(resultsFile)                          entry point (tree menu redirects here)
  ├─ open a standard figure_3d SOURCE figure on the unconstrained (nComponents==3) result
  │     (view_surface_data) — native time-nav, RMS scalar, raw source-vector quiver
  ├─ Cov = tess_operators(Surf,'Covariant')          (find-or-create)
  ├─ setappdata(hFig,'DynamicsOverlay', struct(Cov, Op='Divergence', Cache=Map, srcDS, srcResult, iTess, nV))
  ├─ setappdata(hFig,'CustomOverlayFcn', @(h) view_dynamics('Overlay', h))
  └─ panel_bst_dynamics('SetTarget', hFig, T)

panel_bst_dynamics                                   controls + 3-state workflow (exists)
  └─ Source operator selector: Divergence | Curl | Potential | Stream
        OnMeasurement: write Op into the figure's 'DynamicsOverlay'; view_dynamics('RefreshOverlay', hFig)

view_dynamics('Overlay', hFig)   = i_dynamics_overlay   (the CustomOverlayFcn; per frame)
  iT  = current time index;  D = getappdata(hFig,'DynamicsOverlay')
  if ~isKey(D.Cache, iT):
      Jt = bst_memory('GetResultsValues', D.srcDS, D.srcResult, [], 'CurrentTimeIndex', 0)   % raw 3-vector
      if size(Jt,1) ~= 3*D.nV, return, end
      D.Cache(iT) = process_helmholtz('Compute', Jt, D.Cov)     % one call, all fields
  Ht  = D.Cache(iT)
  scal = i_pick_scalar(Ht, D.Op)        % Divergence->Div, Curl->Curl, Potential->Phi, Stream->Psi
  TessInfo(D.iTess).Data = scal; set ColormapType ('stat2' signed) + DataMinMax; UpdateSurfaceColormap
```

`Ht` is cached **per time index**, so switching operator is a free re-select; only a time change
recomputes. The figure is self-contained — `CustomOverlayFcn` reads everything from its own
`'DynamicsOverlay'` appdata, so it updates correctly even when the panel is not focused.

**Math note (all four are scalars).** On a 2-manifold the rotational invariant of a vector field
is the *scalar* surface vorticity `ω = (∇×J)·n` — the normal projection of the 3-D curl — not a
vector. The four operators are two conjugate pairs: **Divergence** `div J = Δφ` ↔ **Potential**
`φ`; **Curl** `ω = (∇×J)·n = Δψ` ↔ **Stream** `ψ`. `process_helmholtz('Compute')` forms the full
3-D curl per face then projects onto `n` (`Ht.Curl`); the tangential curl components describe
out-of-surface twisting and are carried by the separate normal DOF (`Jn·n`), so the scalar
vorticity is the complete in-surface rotational signal (sign = handedness, `ψ` extrema = cores).

## 6. Components

| Component | Change |
|---|---|
| **`view_dynamics.m`** | Replace the `view_helmholtz(...)` open with: open a `figure_3d` source figure (`view_surface_data` on the unconstrained result), guard `nComponents==3`, resolve `Cov`, install `'DynamicsOverlay'` + `CustomOverlayFcn`, then `SetTarget`. Add verbs: `'Overlay'` (the per-frame `i_dynamics_overlay`) and `'RefreshOverlay'` (re-select from cache + repaint, for operator switches). Add the pure helper `i_pick_scalar(Ht, Op)`. |
| **`panel_bst_dynamics.m`** | Generalize the Measurement row from `jMeasPot`/`jMeasStr` to a four-way **Divergence / Curl / Potential / Stream** selector. `OnMeasurement` sets `st.curOp` + writes `Op` into the figure's `'DynamicsOverlay'` + `view_dynamics('RefreshOverlay', hFig)` (replaces `view_helmholtz('SetComponent')`). `OnRecord`/`OnSaveCursor` read `Ht` from `'DynamicsOverlay'.Cache` and switch on the four operators → `Ht.Div/Curl/Phi/Psi`; the atom `Function` is the operator name. `SetTarget` drops the `HelmholtzState.Lambda` bootstrap (smoothing deferred). |
| **`i_dynamics_overlay` / `i_pick_scalar`** (new, in `view_dynamics`) | The `CustomOverlayFcn` (~40 lines) + the pure scalar selector. Written fresh, not ported. |
| **`view_helmholtz.m`** | **Deleted.** |
| **`panel_helmholtz.m`** | **Deleted.** |
| **`tree_callbacks.m`** (~:1885) | The "Helmholtz / vorticity (Dirac)" source-results menu entry redirects to `view_dynamics` (merge with the existing dynamics entry if present). |

## 7. Data flow (three states)

- **Navigate** — scrub time / panel axis → `figure_3d` native pipeline → `i_dynamics_overlay`
  recomputes (cache miss) and paints the selected differential scalar. Operator switch →
  `RefreshOverlay` → cache hit → repaint. Nothing saved.
- **Detect (ephemeral)** — the differential overlay is the feedback; existing extrema recording
  (`OnRecord` → `bst_dynamics('Extrema')`) runs over the selected field.
- **Save** — `OnSaveCursor`/`OnRecord`/`OnSaveDetection` → `bst_dynamics('Save')`; the "Source"
  atom records `strength`/`charge` from `Ht.(Op)`, `Function` = operator.

## 8. Error handling

- Open: clean `bst_error` if the result is not unconstrained (`nComponents~=3`). `Cov` via
  `tess_operators` (existing Structures-atlas guard).
- `i_dynamics_overlay` runs inside `panel_surface`'s `try/catch` `FireCustomOverlay` (never
  crashes the figure); guards `size(Jt,1)~=3*nV` (skip frame). Cache keyed by time index;
  deterministic, no invalidation (smoothing deferred).

## 9. Testing

GUI code; Brainstorm runs nogui (headless-limited). Strategy:
1. **Pure helper unit test** — `i_pick_scalar(Ht, Op)` for all four operators → `Div/Curl/Phi/Psi`.
2. **Atom-recording test** (extend `test_detect_save` / `test_dynamics_atoms`) — prime a figure's
   `'DynamicsOverlay'.Cache`, run `OnSaveCursor`/`OnRecord`, assert the atom `Function` and
   `strength`/`charge` come from the selected differential field.
3. **Open smoke test** — `view_dynamics(unconstrainedResult)` opens a figure with
   `CustomOverlayFcn` installed and the first frame's `TessInfo.Data` equals the differential
   scalar (not the RMS norm).
4. **Cleanliness + regression** — zero `view_helmholtz`/`panel_helmholtz` references in
   `toolbox/`; existing dynamics regression tests pass after the repoint.

## 10. Decisions (confirmed)

- **D1 — selector widget: combobox** (for now; folds into the Atoms menu in a later step).
- **D2 — operator labels: "Divergence / Curl / Potential / Stream"** (no irrotational/solenoidal).
- **D3 — colormap: `'stat2'`** (diverging) — all four are signed scalars.
- **Curl is the scalar surface vorticity `(∇×J)·n`** (see the Math note in §5); the code is
  mathematically correct as written.
