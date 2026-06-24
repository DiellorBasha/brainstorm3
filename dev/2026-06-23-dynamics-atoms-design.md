# Dynamics Atoms — Spatiotemporal Sparse Markers (Phase 1 Design)

**Date:** 2026-06-23
**Author:** Diellor Basha

## Motivation

Brainstorm has **Events** (sparse in *time*, attached to recordings) and **Scouts**
(sparse in *space*, attached to surfaces) but no structure that is sparse in
*both*. Dipoles are the only spatiotemporal sparse object, but they are rigid
(built for ECD fitting) and not extensible with arbitrary descriptors. We need a
new, extensible **spatiotemporal sparse marker** — fusing the time fields of
Events with the space fields of Scouts, plus frequency, scale, and an open
descriptor bag — to hold detected source-domain singular points and, later,
trajectories and motifs for statistical analysis.

## Core principle: reference, not copy

An **Atom** stores *coordinates and a pointer*, never the underlying data. Its
(time, space, band, scale) coordinates + provenance are sufficient to
re-derive, visualize, and analyze the referenced data on demand (re-apply the
kernel to the band-passed sensors at that time, decompose, inspect that vertex).
This keeps the table tiny and infinitely re-analyzable.

## Schema

### Atom (one spatiotemporal marker) — `db_template('atom')`

```
% --- identity ---
label        char      % free type name, e.g. 'vortex'
category     char      % singular-point type 'source'|'sink'|'vortex'|'antivortex' (derived)
color        [r g b]
% --- TIME coordinate (ref into recording; cf. Events) ---
time         double    % seconds
sample       double    % sample index (convenience)
phase        char      % oscillation phase 'peak'|'trough'|'rising'|'falling'  [COLUMN]
sourceEvent  char      % originating event group, e.g. 'alpha_peak'
% --- SPACE coordinate (ref into surface; cf. Scouts) ---
vertex       double    % seed vertex index on SurfaceFile          [COLUMN: Scout]
pos          [x y z]   % position (display without loading the surface)
hemi         uint8     % 1=L, 2=R                                   [COLUMN]
region       char      % atlas region label (optional)
% --- FREQUENCY coordinate ---
band         [fLo fHi] % Hz
bandName     char      % 'alpha'|'beta'|'gamma'|...                 [COLUMN: Frequency]
% --- SCALE coordinate (dominant bst_eigen spectrum band; reserved, [] for now) ---
scale        [k1 k2]   % eigenmode index range, or []
scaleName    char      % dominant eigenspectrum band (eigen-populated later)  [COLUMN: Scale]
% --- FUNCTION (which scalar field the atom is an extremum of; cf. Scout.Function) ---
Function     char      % 'potential'(Phi)|'stream'(Psi)|'magnitude'(|J|)|...  [COLUMN]
% --- descriptors (measured at the extremum) ---
strength     double    % |J|, omega, or |Phi|/|Psi| extremum value
charge       int8      % +1/-1 topological charge (source/sink sign)
chirality    int8      % +1/-1 vortex rotation sense
persistence  double    % topological persistence
descriptors  struct    % OPEN-ENDED extensible bag (future axes: velocity, motifId, eigencoeffs, ...)
% --- provenance (refs to re-derive the data) ---
DataFile     char      % recording (time source)
ResultsFile  char      % kernel/source link (the inverse used)
SurfaceFile  char      % cortex (space source)
```

Panel columns (panel_dynamics): **Time | Phase | Freq | Scale | Function | Hemisphere | Vertex**.

### Atoms are a table, not "built from" Events/Scouts

The atoms table is the primary object; Events and Scouts are *projections* of it,
not its parents. A planned `bst_dynamics` projector will split a table into a
standard Events group (by `time`/`phase`) and/or a Scouts atlas (by `vertex`) so
Brainstorm's native event/scout navigation can be reused on demand. `Function`
(cf. `Scout.Function`) records which scalar field the atom is an extremum of.

### Grouped, nested, simple/extended (mirrors Events)

The flat per-atom array was replaced by **atom GROUPS** that reuse the Events
grouping system. Each group is an extended Event group + atom extensions:
- `times` is **[1×N] simple** (point) or **[2×N] extended** (window), exactly
  like Events; `type` tracks which.
- per-occurrence **parallel arrays** (`vertices`/`pos`/`hemi`/`strength`/...),
  group-level coordinates (`phase`/`band`/`scale`/`Function`), and a **`parent`**
  label for **nesting**.

The populate builds the refphase hierarchy:
```text
alpha (8-13 Hz)   [extended window]   parent=''
├─ alpha_peak     [simple]            parent='alpha (8-13 Hz)'  Function=magnitude
├─ alpha_trough   [simple]
├─ alpha_rising   [simple]
└─ alpha_falling  [simple]
```
`panel_dynamics` is now a **JTree**: top-level groups → nested child groups →
occurrence leaves; selecting a leaf highlights its cortex marker and jumps the
time. `view_dynamics` draws one marker set per spatial group; windows
(temporal-only, no vertex) appear in the tree but draw no markers.

### Table (collection) — `db_template('dynamicsmat')`

```
Comment      char           % node label
Atoms        struct array   % the markers (db_template('atom'))
nAtoms       double
DataFile     char           % default recording (if homogeneous)
SurfaceFile  char           % default surface
Options      struct         % detection parameters (band, gates, ...)
History      cell
```

Stored as `dynamics_*.mat` in the study folder. (Phase 1: loadable by path; a
tree node + viewer come in Phase 2 via the bst-java fork — no GUI here.)

## Milestone 1 — IMPLEMENTED (data model + GUI skeleton)

Per the refined scope: build the **atom *system*** (model + GUI + integration);
defer the principled detection science. Test data comes from a trivial
source-magnitude peak detector, NOT the Helmholtz cores.

1. **`db_template.m`** — `'atom'` and `'dynamicsmat'` cases. ✅
2. **`toolbox/dynamics/bst_dynamics.m`** — the dynamics-module orchestrator
   (own `toolbox/dynamics/` folder, mirroring `toolbox/eigen` / `toolbox/timefreq`;
   the module that works on atom files). `New` / `NewGroup` / `AddGroup` / `Flatten` / `Save` / `Load`
   (Phase 1 saves/loads by path; the `dynamics_` type is not yet registered with
   file_gettype — that comes with the Phase-2 tree node). ✅
3. **`process_source_atoms.m`** ("Detect source atoms (peaks)") — TEST POPULATE.
   Input: unconstrained kernel link. At each event time it reconstructs `|J|`
   and takes the top-N **local maxima** on the surface graph (`VertConn`) as
   atoms (category `'peak'`, strength `|J|`, full provenance). Saves a
   `dynamics_*.mat` (not returned as a tracked OutputFile — bst_process can't
   resolve the unregistered type). ✅
4. **`toolbox/gui/view_dynamics.m`** — loads a table, opens the cortex
   (`view_surface`), draws atoms as colored markers on `Axes3D` (offset ~2 mm
   along the vertex normal to clear the opaque surface; low-level `line()` so it
   never trips the NextPlot/newplot reset), and opens the panel. ✅
5. **`toolbox/gui/panel_dynamics.m`** — a `JList` of atoms; row selection
   highlights the atom's marker (a yellow `AtomSel` marker) and jumps the
   recording time (`panel_time('SetCurrentTime', ...)`). ✅

## Testing — `dev/test_dynamics_atoms.m` (3/3)

- T1: model + I/O round-trip (Save→Load identical, fields match template).
- T2: populate on MN-unconstrained kernel × `data_block001_02`, `alpha_peak`
  → 192 atoms (3×64), every Atom has valid vertex/pos/time/band/provenance.
- T3: `view_dynamics` smoke — figure with atom markers + `AtomSel` + panel with
  all 192 rows; selecting a row moves the selection marker.

## panel_bst_dynamics (increment 1 — bookkeeping core)

A Record-style atom panel (`toolbox/gui/panel_bst_dynamics.m`, panel name `Dynamics`),
superseding the `panel_dynamics` skeleton (retired to `dev/experimental/`):
- **File** menu: Open / Save / Save as a `dynamics_*` table.
- **Atoms** menu: Add / Rename / Delete / Set color / Sort (by name / time) group.
- **Reuses the Record Events-section components** (max UI reuse): a colored group
  **JList** (`BstColorListRenderer` + `BstListItem`/`setColor`; child phase groups
  *indented* under their window) on the left, the selected group's **occurrence
  JList** (`BstStringListRenderer`) on the right, in a `JSplitPane`. Select a group
  → its occurrences; select an occurrence → cortex marker highlight +
  `panel_time('SetCurrentTime', …)`.
- **Docked** as a tools tab (`gui_show(bstPanel,'BrainstormTab','tools')`), like
  Record/Helmholtz — not a floating window.
- Edits modify the in-memory table (`getappdata(0,'DynamicsTarget').T`) then
  `view_dynamics('Redraw', hFig, T)` (factored marker draw, per-group `AtomMarker<g>`
  tags, cached surface) + rebuild the tree.

The four analysis axes — **Time** = refphase, **Space** = `panel_helmholtz` (the
compact one), **Frequency** = band, **Scale** = `bst_eigen` — fold into this panel
in later increments.

## Out of scope (later phases)

- Phase 2: a DB tree node for `dynamics_*` (bst-java fork) + file_gettype/
  file_fullpath wiring so tables appear in the tree and open on double-click.
- Phase 3: principled detection (Helmholtz Sources/Cores → source/sink/vortex
  with the source-amplitude + persistence gates), replacing the test peak
  populate; then trajectory/motif grouping + statistics.
- Scale axis population (bst_eigen) — slot reserved, not filled yet.
