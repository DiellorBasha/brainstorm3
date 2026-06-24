# Geodesic Area tool → dynamics suite — design

**Goal:** Move the custom geodesic Area tool (heat-distance cortical disk) out of `panel_scout` and into
the dynamics suite as a dynamics-owned tool that emits a **transient disk** (`seed` center + `radius`
extent + vertices) — the native **source-axis point+extent primitive** of the atom tensor — and rewire
the atom panel's Capture to use it instead of a Scout.

**Status:** design, ready to plan. Phase 2 of the atom-tensor architecture
(`docs/superpowers/specs/2026-06-24-atom-tensor-architecture-analysis.md`, §2.8 / §6 step 2). Prerequisite
for the source-axis navigator block (a later phase).

**Author:** Diellor Basha, 2026-06-24

---

## 1. Motivation

The geodesic Area tool was added as a **custom extension to `panel_scout`** (not standard Brainstorm
scout behaviour): a click seeds a heat-distance disk, scroll grows it, and it is stored as a *scout* in
the atlas. In the atom-tensor model a seed (`center`) + radius (`extent`) is exactly the **source-axis
Localization**, so the tool belongs in the dynamics suite, where it can emit a *transient* disk read
directly by the atom panel — with no scout, hence no atlas-pollution / bulk-scout-delete hazard.

This move does **not** build the source-axis navigator block (a later phase); it delivers the
dynamics-owned tool, the figure-interaction plumbing, and the Capture rewire, so the source primitive is
ready.

## 2. Current footprint (what exists today)

- **Engine** `toolbox/anatomy/tess_scout_area.m` (heat-distance disk; pure geometry). Called by
  `panel_scout` only — but by **both** the Area tool **and** the geodesic **Line** tool
  (`'prewarm'`/`'path'` verbs). **Stays in anatomy** (shared geometry primitive).
- **Area tool UI/interaction** in `panel_scout.m`: `jToggleArea`, `IsAreaToolActive`, `AreaToolToggle`,
  `CreateScoutArea`, `AreaToolScroll`, `AreaCache`; the `if IsAreaToolActive()` dispatch inside
  `CreateScoutMouse`. Output = a `CreateScout` in the atlas.
- **Geodesic Line tool** in `panel_scout.m`: `jToggleGeo`, `GeodesicToolToggle`, `CreateGeodesicMouse`,
  `i_geodesic_edgeparam`, and the shared helper `i_get_scout_surface`. **Stays in scouts.**
- **figure_3d hooks**: cortical click → `panel_scout('CreateScoutMouse', hFig)` when
  `getappdata(hFig,'isSelectingCorticalSpot')` (line ~624); scroll → `panel_scout('AreaToolScroll', …)`
  (line ~944). Clicked surface vertex is resolved by `select3d(TessInfo(iTess).hPatch)` → `vi`.
- **Capture coupling**: `panel_bst_dynamics.m:394` reads `panel_scout('GetSelectedScouts')`.

## 3. Decisions (locked in brainstorming)

1. **Move the Area tool's UI + interaction into the dynamics suite**, emitting a transient disk (not a
   scout); **remove** it from the Scout panel. Keep `tess_scout_area` in `toolbox/anatomy`. Leave the
   geodesic **Line** tool in scouts.
2. **figure_3d gets two small parallel branches** routing click/scroll to the dynamics tool when its pick
   flag is set — no shared-dispatcher refactor.
3. The atom panel's **Capture** reads the dynamics tool's state directly (no `GetSelectedScouts`).

## 4. Architecture

### 4.1 New module `toolbox/dynamics/bst_geodesic_tool.m`

Verb-dispatched (`eval(macro_method)`), in the dynamics suite (auto on path). Holds the transient
heat-disk tool state and reuses `tess_scout_area`. Near-lift of `CreateScoutArea` + `AreaToolScroll` +
`AreaCache`, with scout-creation swapped for transient disk state + overlay.

```
bst_geodesic_tool('Toggle', onoff)   enter/exit dynamics cortical-pick mode: set appdata
                                       'isDynamicsGeodesicPick' on the scout-capable 3-D figures,
                                       cross cursor, prewarm tess_scout_area('prewarm', Surf);
                                       on exit restore the pointer + clear the flag + Clear the overlay.
bst_geodesic_tool('IsActive')        -> logical (the toggle state)
bst_geodesic_tool('OnClick', hFig)   resolve vi=select3d(TessInfo(iTess).hPatch); Seed(hFig, vi)
bst_geodesic_tool('Seed', hFig, vi)  [verts, phi] = tess_scout_area(Surf, vi, R0=0.003); cache
                                       struct(seed=vi, phi, radius=R0, vertices=verts, SurfaceFile, hFig);
                                       Draw. (programmatic/test entry; bypasses select3d)
handled = bst_geodesic_tool('OnScroll', scrollCount)   if active+seeded: R = max(STEP, R - scrollCount*STEP);
                                       verts = tess_scout_area(Surf, seed, R, phi) (cached phi, no re-solve);
                                       update cache; Draw; return true. Else return false.
st = bst_geodesic_tool('GetState')   -> struct(seed, pos, radius, vertices, SurfaceFile) or [] if unseeded
bst_geodesic_tool('Draw', hFig)      transient disk overlay: patch of faces with all 3 verts in the disk,
                                       FaceAlpha ~0.3, EdgeColor none, tag 'GeodesicToolDisk' (low-level,
                                       NextPlot='add' to avoid the axes-reset trap). Clears prior overlay first.
bst_geodesic_tool('Clear', hFig)     delete the 'GeodesicToolDisk' overlay.
```

State cache is `persistent` (one active disk at a time), mirroring `AreaCache`. `pos` = the seed's 3-D
coordinates from the surface vertices (for the source Localization `center`).

`R0 = 0.003` (3 mm) and `STEP = 0.003` mirror the current Area tool constants.

### 4.2 `toolbox/gui/figure_3d.m` — two parallel branches

Click handler (~line 623), dynamics pick takes precedence when active:

```matlab
if getappdata(hFig, 'isDynamicsGeodesicPick')
    bst_geodesic_tool('OnClick', hFig);
elseif isSelectingCorticalSpot
    panel_scout('CreateScoutMouse', hFig);
end
```

Scroll handler (~line 944). The current `panel_scout('AreaToolScroll', …)` call is **replaced** by the
dynamics tool's scroll (the Area tool is removed from scouts, 4.3), ahead of the default zoom:

```matlab
% was: if panel_scout('AreaToolScroll', double(event.VerticalScrollCount)) ...
if bst_geodesic_tool('OnScroll', double(event.VerticalScrollCount))
    % consumed: grew the dynamics disk
else
    % ... default zoom (unchanged)
end
```

`OnClick`/`OnScroll` return early / `false` when the tool is inactive, so scout selection and default zoom
are unchanged when the dynamics tool is off.

### 4.3 `toolbox/gui/panel_scout.m` — remove the Area tool only

- **Remove:** `jToggleArea` (button + its row), `AreaToolToggle`, `CreateScoutArea`, `AreaToolScroll`,
  `AreaCache`, `IsAreaToolActive`, and the `if ~isVolumeAtlas && IsAreaToolActive()` dispatch block inside
  `CreateScoutMouse`.
- **Keep:** the geodesic **Line** tool (`jToggleGeo`, `GeodesicToolToggle`, `CreateGeodesicMouse`,
  `i_geodesic_edgeparam`), `i_get_scout_surface` (the Line tool uses it), and the `IsGeodesicToolActive()`
  dispatch in `CreateScoutMouse`.
- The `jToggleArea` reference in the tool-release block (`AreaToolToggle`/`GeodesicToolToggle` were
  mutually exclusive) is updated: the Line tool no longer needs to deselect a (now-removed) Area toggle.

### 4.4 `toolbox/gui/panel_bst_dynamics.m` — rewire Capture + add the Region-tool toggle

- `OnCaptureRegion` (`:394`): replace `[sScout,~,sSurf] = panel_scout('GetSelectedScouts')` with
  `st = bst_geodesic_tool('GetState')`. Guards: if `isempty(st)` → warn "Seed a region with the Region
  tool first." If `st.SurfaceFile` differs from the atom table's surface (via `file_compare`, with the
  group-level fallback) → warn "different surface." Then `seed=st.seed; pos=st.pos; hemi=1+(pos(2)<0);`
  and `AttachRegion(G, o, st.vertices, seed, pos, hemi)` — identical downstream to today.
- Add a **Region tool** toggle in the Record/Capture row that calls `bst_geodesic_tool('Toggle', …)`, so
  the user activates the heat-disk tool from the panel (click a vertex, scroll to grow), then **Capture**.
- The "Capture region → active atom" button/menu item is unchanged; only its data source changes.

## 5. Data flow (the new Capture path)

```
user: Region tool ON  → bst_geodesic_tool('Toggle', 1)  (figures enter dynamics pick mode)
user: click cortex    → figure_3d → bst_geodesic_tool('OnClick') → seed disk (transient overlay)
user: scroll          → figure_3d → bst_geodesic_tool('OnScroll') → grow/shrink disk
user: select an atom occurrence in the list  (the active atom)
user: Capture         → OnCaptureRegion → st=bst_geodesic_tool('GetState')
                         → bst_dynamics('AttachRegion', G, o, st.vertices, st.seed, st.pos, hemi)
                         → i_apply (redraw) + auto-save
```

No scout is created at any point; the disk is a transient overlay snapshotted (a copy) into the atom.

## 6. Testing

- **`dev/test_bst_geodesic_tool.m`** (new, headless; Brainstorm live, a cortex figure open or a surface
  loaded):
  - `Seed(hFig, vi)` at an explicit vertex → `GetState` returns `seed==vi`, a non-empty `vertices` disk
    around `vi`, and `pos` = that vertex's coordinates; `SurfaceFile` set.
  - Growing (drive the radius path via `OnScroll(-1)` or a `Seed` at larger R) enlarges the vertex set
    monotonically; shrinking reduces it; radius floors at one `STEP`.
  - `GetState` is `[]` before any seed; `Clear` removes the `'GeodesicToolDisk'` overlay.
  - (Interactive `OnClick`/`select3d` is not unit-tested; `Seed` covers the engine path.)
- **`dev/test_dynamics_atoms.m`** — T6/T7 build regions via `tess_scout_area` + `AttachRegion` directly
  (not via Scouts), so they are unaffected and must stay **8/8**. Add a check that `OnCaptureRegion` reads
  `bst_geodesic_tool('GetState')`: seed the tool programmatically, select an occurrence, call
  `OnCaptureRegion`, assert `region`/`vertices` written to that occurrence — proving the Scout decoupling
  end-to-end.
- **Regression**: a Scout-panel smoke check that the Scout panel still loads and the geodesic **Line**
  tool still works (Area toggle gone, no errors).

## 7. Out of scope

- The source-axis **navigator block** (Time/Freq/Source/Scale uniform blocks) — a later phase; this phase
  only delivers the tool + plumbing + Capture rewire.
- Renaming `tess_scout_area` (kept; both the dynamics tool and the scout Line tool use it).
- Moving the geodesic **Line** tool (stays in scouts).
- Soft/weighted (untruncated heat-kernel) disk — the wavelet form is future (§2.6 of the analysis); this
  tool emits the hard disk.

## 8. Risks / notes

- **figure_3d is shared plumbing**: the parallel branches must be strictly gated on
  `isDynamicsGeodesicPick` / the tool's active+seeded state so scout selection and default zoom are
  untouched when the dynamics tool is off. The `OnScroll` early-`false` return is the guard.
- **Mutual exclusion**: activating the dynamics Region tool should deselect the scout New/Area/Line modes
  (and vice versa) to avoid two cortical-pick modes at once — the dynamics `Toggle` clears
  `isSelectingCorticalSpot` on the figures, and the scout tools already clear their own on activation.
- **Removing `AreaToolScroll`** from `panel_scout` requires the figure_3d scroll edit in the **same**
  change, or the deleted call errors. Plan them as one task.
