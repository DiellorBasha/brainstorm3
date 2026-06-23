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
label        char      % type, e.g. 'vortex' / 'source'
category     char      % 'source' | 'sink' | 'vortex' | 'antivortex'
color        [r g b]
% --- TIME coordinate (ref into recording; cf. Events) ---
time         double    % seconds
sample       double    % sample index (convenience)
sourceEvent  char      % originating event group, e.g. 'alpha_peak'
% --- SPACE coordinate (ref into surface; cf. Scouts) ---
vertex       double    % seed vertex index on SurfaceFile
pos          [x y z]   % position (display without loading the surface)
hemi         uint8     % 1=L, 2=R
region       char      % atlas region label (optional)
% --- FREQUENCY coordinate ---
band         [fLo fHi] % Hz
bandName     char      % 'alpha'
% --- SCALE coordinate (eigenmode band via bst_eigen; reserved, [] for now) ---
scale        [k1 k2]   % eigenmode index range, or []
scaleName    char      % or ''
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
2. **`toolbox/io/bst_dynamics.m`** — `New` / `NewAtom` / `Add` / `Save` / `Load`
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

## Out of scope (later phases)

- Phase 2: a DB tree node for `dynamics_*` (bst-java fork) + file_gettype/
  file_fullpath wiring so tables appear in the tree and open on double-click.
- Phase 3: principled detection (Helmholtz Sources/Cores → source/sink/vortex
  with the source-amplitude + persistence gates), replacing the test peak
  populate; then trajectory/motif grouping + statistics.
- Scale axis population (bst_eigen) — slot reserved, not filled yet.
