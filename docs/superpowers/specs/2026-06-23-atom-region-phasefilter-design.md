# Atom geodesic regions + phase filtering — design

**Goal:** Make the Dynamics atom a *true* spatiotemporal atom by (1) letting the user attach a
geodesic cortical region to the active atom (turning a time-only phase marker into a localizable
time+space atom), and (2) filtering the atom list/cortex by oscillation phase so analysis can focus
on one phase or a combination.

These are two independent features sharing one panel. Both are display/recording layer changes over
the existing `panel_bst_dynamics` + `view_dynamics` + `bst_dynamics` stack; no detector or solver
changes.

## Context (current state)

- After **Detect windows**, the table holds a band-window extended stack `alpha (8-13 Hz)` plus 4
  phase-marker children `alpha_peak/_trough/_rising/_falling`. Phase markers are **time only**
  (`times` filled, `vertices/pos` empty) and draw no cortex markers.
- The atom model (`db_template('atomgroup')`) stores **per-occurrence** `vertices`/`pos`/`hemi` — a
  single seed per occurrence. There is **no field for a region** of vertices.
- The geodesic Area tool lives in `panel_scout`: `tess_scout_area(SurfaceFile, Seed, Radius, [phi])`
  grows the heat-distance disk; the figure's click (`isSelectingCorticalSpot` → `panel_scout('CreateScoutMouse')`)
  and scroll (`panel_scout('AreaToolScroll')`) plumbing is hardwired to `panel_scout`.
- `panel_scout('GetSelectedScouts')` is a callable verb returning `[sScouts, iScouts, sSurf, iSurf]`.
- `gui_component('checkboxmenuitem', ...)` exists (a `JCheckBoxMenuItem`).

## Decisions (locked during brainstorming)

1. **Region tool = reuse the Scout Area tool.** The user grows a disk with the existing Scout Area
   tool (familiar scroll-to-grow), then a Dynamics **Capture region → active atom** action snapshots
   the selected scout's vertices into the atom. No new figure plumbing; the Scout atlas is touched
   only transiently (we copy vertices out; we never persist a scout into the atom).
2. **Storage = per-occurrence `region` field** on `atomgroup` (cell array). Empty ⇒ point atom
   (back-compatible).
3. **Active atom = the selected list occurrence.** Capture writes `region{i}` and sets the seed
   `vertices(i)`/`pos(i)` from the scout seed, localizing a previously time-only marker.
4. **Phase filter = a "Show phases" submenu** under the Atoms menu (checkable items), affecting both
   the flat list and the cortex markers/regions. Non-destructive.

## Architecture & data model

One new field on `db_template('atomgroup')`:

```
region   % cell {1×N}: region{i} = [v1 v2 … vk] cortex vertex indices for occurrence i; {} = point atom
```

`vertices(i)`/`pos(i)` remain the **seed**; `region{i}` is the geodesic disk grown around it. The
field rides along in the `.mat` and defaults to `{}` via the template, so existing tables and the
`Record`/`Detect` paths keep working unchanged.

The schema-equality regression (`isequal(fieldnames(T2.Groups), fieldnames(db_template('atomgroup')))`)
continues to hold because the field is added to the template itself.

**New orchestrator verb** (keeps occurrence arrays aligned; panel stays thin; logic is unit-testable):

```matlab
G = bst_dynamics('AttachRegion', G, o, regionVerts, seed, pos, hemi)
%   Pads vertices/pos/hemi/strength/charge to full length N = size(G.times,2) with NaN, and
%   region to length N with []. Then sets occurrence o: region{o}=regionVerts(:)',
%   vertices(o)=seed, pos(o,:)=pos, hemi(o)=hemi. Returns the updated G.
%   A time-only marker (vertices=[] pos=[]) becomes partially localized: occurrence o has a
%   finite seed + a region; the other occurrences stay NaN/[] (still time-only).
```

### Files touched

| File | Change |
|------|--------|
| `toolbox/db/db_template.m` | add `region` field to the `atomgroup` template |
| `toolbox/dynamics/bst_dynamics.m` | add the `AttachRegion` verb |
| `toolbox/gui/panel_bst_dynamics.m` | Capture action (menu item + Record-section button) + `OnCaptureRegion`; `Show phases` submenu + `st.showPhase` state + `OnTogglePhase`; phase filtering in `i_window_atoms`; pass `showPhase` into Redraw |
| `toolbox/gui/view_dynamics.m` | `Redraw`: draw region patches; honor `showPhase`; tolerate partial localization (NaN pos rows) |
| `dev/test_dynamics_atoms.m` | T6 (Capture localizes a phase marker) + T7 (phase filter) |

## Feature A — Capture region → active atom

```
1. User selects a phase-marker row in the right list, e.g. "22.600s  peak".
   → ACTIVE atom = that occurrence (st.occMap row → group g, occurrence o).
2. User activates the existing Scout "Area" tool, clicks a vertex, scrolls to grow
   the geodesic disk (heat distance; panel_scout owns this, unchanged).
3. User clicks "Capture region → active atom":
      [sScout,~,sSurf] = panel_scout('GetSelectedScouts');
      validate: an occurrence is selected; a scout is selected; sSurf matches the atom SurfaceFile.
      seed = sScout.Seed (or sScout.Vertices(1) if no seed);
      pos  = Surf.Vertices(seed,:);  hemi = 1 + (pos(2) < 0);   % SCS Y>0 = left
      G = bst_dynamics('AttachRegion', G, o, sScout.Vertices, seed, pos, hemi);
      st.T.Groups(g) = G;  i_apply(st);  auto-save.
   The "22.600s peak" atom now has a seed + region → it draws on the cortex and is highlightable.
```

Capture is exposed in **two** places (both call `OnCaptureRegion`): an item in the Atoms menu bar
(beside *Record at cursor*) and a button in the **Record** section (it is the spatial sibling of
*Record at cursor*).

**Guard rails (warn-and-abort, no silent failure):**
- No occurrence selected → "Select an atom in the list first."
- No scout selected (`GetSelectedScouts` empty) → "Grow a region with the Scout Area tool first."
- Scout surface ≠ atom SurfaceFile → "The selected region is on a different surface." abort.
- The active occurrence must be **simple** (a phase child or recorded group); a window (extended)
  leaf is not a capture target → "Select a single atom, not a time window."

The captured vertices are a **snapshot copy**; the transient scout is left as the user's selection
and never referenced by the atom (keeps atoms self-contained and avoids the bulk-scout-delete atlas
hazard).

## Feature B — Show-phases filter

```
st.showPhase = [1 1 1 1]    % peak trough rising falling; all on by default
Atoms ▾ → Show phases ▸ → ☑peak ☑trough ☑rising ☑falling   (checkboxmenuitem each)
toggling phase k → flips st.showPhase(k) → i_apply (rebuild tree list + redraw)
```

Applied consistently in two places:
- `i_window_atoms` — skip occurrences whose `phase` is toggled off (flat list shrinks).
- `view_dynamics('Redraw')` — skip markers/regions for phase-child groups whose `phase` is toggled
  off (cortex follows the list).

Pure display filter — never deletes atoms. Groups with no `phase` (recorded Function groups) are
always shown. The phase→index mapping is `{'peak','trough','rising','falling'}`.

## Rendering (view_dynamics Redraw)

For each occurrence `o` of a group `g` that is shown (phase not filtered out) and has a non-empty
`region{o}`, draw a translucent cortex patch of the faces fully inside the region:

```matlab
inReg = false(nVert,1);  inReg(G.region{o}) = true;
fIn   = all(inReg(Surf.Faces), 2);                      % faces with all 3 verts in the region
patch('Faces',Surf.Faces(fIn,:), 'Vertices',Surf.Vertices, 'Parent',hAxes, ...
      'FaceColor',G.color, 'FaceAlpha',0.35, 'EdgeColor','none', ...
      'Tag', sprintf('AtomRegion%d_%d', g, o));
```

- Region patches are cleared each redraw alongside the existing wipe, by extending the delete to the
  `AtomRegion` tag prefix (`delete(findobj(hAxes,'-regexp','Tag','^AtomMarker'))` gains a sibling for
  `^AtomRegion`).
- The seed still draws as the group's `AtomMarker%d` point.
- `GroupsPosOff{g}` stays occurrence-aligned: it is built from the full-length `pos` (NaN rows for
  unlocalized occurrences). `line()` skips NaN points, so unlocalized occurrences draw no marker and
  their highlight stays hidden. Partial localization therefore needs no special-casing in the
  occurrence→position map.
- `DynamicsSurf` already caches `Surf.Faces`/`Surf.Vertices`, so no surface reload.

## Testing (dev/test_dynamics_atoms.m, headless)

Extends the existing suite (currently T1–T5). Both new tests reuse the kernel-link fixture; they SKIP
with the suite if no unconstrained kernel link is available.

- **T6 — Capture localizes a phase marker.** After Detect: pick a phase-child group and occurrence
  `o`; synthesize a region via `tess_scout_area(SurfaceFile, seed, 0.005)` (5 mm); call the capture
  path; assert `region{o}` non-empty, `vertices(o)`/`pos(o)` finite, the group now draws an
  `AtomMarker%d` line and an `AtomRegion%d_%d` patch on the figure, and a Save→Load round-trip
  preserves `region` (schema-equality holds).
- **T7 — Phase filter.** Toggle `peak` off; assert the selected window's flat list has no peak rows
  and the figure has no peak markers/regions; toggle `peak` back on → rows and markers restored.

Suite target: **7/7** (T2/T3/T6/T7 SKIP together when no kernel link is present).

## Out of scope (later)

- Growing the region from *within* the Dynamics panel (own toggle + figure hooks). Deferred in favor
  of reusing the Scout Area tool.
- Region-level descriptors (area, mean strength over the region, flux through its boundary).
- Phase tagging / multi-select combination logic beyond per-phase show/hide.
- Editing or re-growing a captured region in place (re-capture overwrites `region{o}`).
