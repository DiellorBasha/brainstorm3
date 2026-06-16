# Cross-frame vortex tracking (Phase 2)

**Date:** 2026-06-16
**Status:** approved (design)
**Builds on:** Phase 1 (persistence-based core detection, commit 731d5a3e).

## Problem

Phase 1 detects persistence-ranked vortex cores (and sources/sinks) per frame. The
GUI is active-frame only, so it cannot link cores across time. Phase 2 adds a **batch
process** that runs detection over a time window and links cores into trajectories,
stored in a format Brainstorm can already visualize.

## Decisions (settled with user)

- **Output = Brainstorm dipoles file** (reuse `view_dipoles`/`panel_dipoles`: time
  slider, per-group selection, DB integration). A dedicated trajectory-polyline viewer
  is deferred to a later phase if the dipole representation proves limiting.
- **Linking = greedy nearest + gating**, chirality-consistent, within a max-jump radius.
- **Scope = both vortices (ψ cores) and sources/sinks (φ cores)** in one run, one dipoles
  file, groups distinguished by name.

## Components

1. **`bst_vortex_track.m`** (new, `toolbox/math/`) — pure linker, geometry-light.
   - Input: per-frame core lists (each carrying `pos` [1x3], `chirality` (+/-1),
     `persistence`, `iVertex`) + options `MinPersistence` (default 0), `MaxJump` (m).
   - Greedy chirality-consistent nearest-neighbor linking across consecutive frames:
     pre-filter cores by persistence; match each active track head to the nearest
     same-chirality core in the next frame within `MaxJump` (candidates sorted by
     distance, assigned one-to-one); unmatched next-frame core = birth, unmatched head
     = death.
   - Output: struct array `Tracks`, one per trajectory:
     `.frames` [1xL], `.iVertex` [1xL], `.pos` [Lx3], `.persistence` [1xL],
     `.chirality` (scalar +/-1), `.birthFrame`, `.deathFrame`.
   - Unit-testable in isolation (no Brainstorm dependency).

2. **`bst_dirac_helmholtz.m`** `Decompose` — also emit `H.Sources{t}` (currently only
   `H.Cores{t}`). One-line addition; `Frame` already returns `Ht.Sources`.

3. **`process_vortex_track.m`** (new, `toolbox/process/functions/`) — `results` -> `dipoles`.
   - Precondition: unconstrained (3-component) Dirac source file (same check as
     `view_helmholtz`); else `bst_report('Error', ...)`.
   - Load operators (Dirac, Laplace-Beltrami) for the surface (auto-build via
     `tess_operators` if absent); `Prepare` once.
   - Over the selected time window: `J = kernel*data` per column, `Frame` -> `Ht.Cores`
     and `Ht.Sources` per frame.
   - Link vortices and sources separately via `bst_vortex_track`.
   - Assemble one `DipolesMat` (see encoding) and save + `db_add`.

## Dipole encoding

- Each trajectory -> one dipole `.Index` (group). `DipoleNames{Index}` =
  `'Vortex+ #k'` / `'Vortex- #k'` / `'Source #k'` / `'Sink #k'`.
- Per frame in a track:
  - `.Index` = track id (group), `.Time` = frame time,
  - `.Loc` = sub-vertex `pos'` (3x1, SCS meters, same CS as surface Vertices),
  - `.Amplitude` = `chirality * (persistence/maxPers) * normal(core)` (3x1) — arrow
    points along the spin axis by handedness, length proportional to persistence,
  - `.Goodness` and `.Perform` = persistence (the viewer's amplitude-threshold slider
    then behaves as a persistence gate, consistent with Phase 1),
  - remaining dipole fields (`Origin`, `Errors`, `Noise`, ...) set to the neutral
    defaults used by `process_dipole_scanning`.
- `DipolesMat.Time` = unique frame times; `.DataFile` = input results file; `.Subset`
  and `.DipoleNames` filled per group as in `process_dipole_scanning`.
- Normals from `Surf.VertNormals(iVertex,:)`; `maxPers` = max finite persistence across
  all tracks (Inf-persistence globals use `maxPers`).

## Process options

- `timewindow` — time range to process.
- `MinPersistence` (default 0) — persistence pre-filter (units of psi/phi).
- `MaxJump` (default 0.010 m) — max core displacement between consecutive frames.
- `trackSources` (default true) — also track phi sources/sinks.

## Error handling

- Non-3-component input -> `bst_report('Error', ...)`, return `{}`.
- Empty window or zero cores -> empty dipoles file + `bst_report('Warning', ...)`.
- Missing Dirac/LBO operators -> auto-build via `tess_operators`.

## Testing (TDD)

- `bst_vortex_track` unit tests (hand-built synthetic core sequences):
  one core moving along a line over 5 frames -> one length-5 track; two non-crossing
  cores -> two tracks; one birth + one death detected; chirality mismatch blocks a link;
  a jump > `MaxJump` splits into two tracks.
- `Decompose` regression: `H.Sources` present and equals per-frame `Frame` output.
- Integration: run `process_vortex_track` over the S01 alpha window -> a dipoles file
  whose dominant vortex track persists across consecutive frames; `view_dipoles` opens it.

## Non-goals

Dedicated trajectory-polyline viewer; merge/split topology beyond birth/death; optimal
(Hungarian) assignment; null-based MinPersistence calibration. Deferred.
