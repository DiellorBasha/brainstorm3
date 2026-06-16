# Vortex-core detection by topological persistence (Phase 1)

**Date:** 2026-06-16
**Status:** approved (design)
**Scope:** Phase 1 of 2. This phase = principled per-frame core detection. Phase 2
(cross-frame tracking, a batch `process_`) is a separate later cycle that builds on this.

## Problem

`bst_dirac_helmholtz>FindCores` flags every 1-ring local max/min of the solenoidal
stream function ψ and returns all of them (34 on a real frame), ranked by nothing.
The only filter is an *amplitude* gate in `view_helmholtz` (`GateFrac × max|ω|`), which
cannot distinguish a shallow noise bump from a deep, well-separated vortex. We replace
"all extrema + amplitude gate" with "extrema ranked by **topological persistence**".

Persistence = peak-minus-saddle contrast from a single union-find pass over the
superlevel-set filtration (the elder rule). It gives every core a principled
significance score; shallow bumps die at a nearby saddle (low persistence), genuine
vortices sit far from the diagonal.

## Decisions (settled with user)

- **Per-hemisphere** detection (fixes the example's single-global bug on a disconnected
  cortex: each hemisphere yields its own Inf-persistence global core).
- Gate becomes **persistence (fraction of frame max)** — drop-in for the amplitude gate,
  no new UI controls.
- Phase 1 = detection **+ sub-vertex localization**. Cross-frame tracking is Phase 2.

## Components

1. **`bst_vortex_persistence.m`** (new, `toolbox/math/`) — pure topological detection on a
   scalar field over a mesh 1-ring. Brainstorm-ified port of the reference: classic input
   parsing (no `arguments` block, older-MATLAB compatible), geometry-free, generic (used
   for ψ cores AND φ sources). Returns a columnar struct
   `{peak, persistence, birth, death, chirality, isGlobal}`, sorted by persistence desc.
   Detects `+field` (chirality +1) and `−field` (chirality −1) separately.

2. **`bst_dirac_helmholtz.m`**
   - `Prepare`: precompute & cache **per-hemisphere 1-ring neighbor lists** (from
     `Op.VertConn` restricted to each `Op.vH{hh}` — exact, hemispheres are disconnected),
     and stash `Op.Vtx`. Reuses adjacency across frames (no per-frame rebuild).
   - `FindCores` (rewritten, generic over field+sign-source): loop hemispheres → call
     `bst_vortex_persistence` on that hemisphere's field values → map peaks back to global
     vertex indices → sub-vertex localize → assemble. Returns ALL cores (no threshold here;
     the GUI gate prunes at display so the slider stays live).

3. **`view_helmholtz.m`** — `UpdateFrame` gate switches from `|ω|` to persistence:
   `keep mk where persistence ≥ GateFrac · max(finite persistence)`; `isGlobal` cores always
   kept. Marker plotting may use sub-vertex `pos` when present. Readout: "N cores, top
   persistence X".

4. **`panel_helmholtz.m`** — relabel the gate slider/tooltip to "persistence" (mechanics
   unchanged: 0–100 → fraction).

## Output schema (struct array, back-compatible)

Keep `iVertex`, `charge` (= chirality sign), `omega` (vorticity for cores / divergence for
sources, retained for display continuity). **Add** `persistence`, `isGlobal`, `birth`,
`death`, `pos` (1×3 sub-vertex xyz). Sorted by persistence descending.

## Sub-vertex localization

Quadratic fit of the field over the core's 1-ring in a local tangent chart (project
neighbor offsets onto the tangent plane via the vertex normal); extremum at `−H⁻¹[b;c]`,
mapped back to 3D and clamped to the 1-ring radius. Degenerate/flat Hessian → fall back to
the vertex position. Refines `pos` only, not persistence.

## Testing (TDD)

- `bst_vortex_persistence` unit tests on tiny synthetic meshes with hand-computed answers:
  two unequal bumps → one finite-persistence feature (peak₂ − saddle) + one Inf global; the
  `+ψ`/`−ψ` split yields correct chirality; flat field → a single global, no spurious cores.
- Per-hemisphere: synthetic ψ with one bump per hemisphere on the canonical cortex → two
  Inf-globals (proves the disconnected-mesh fix).
- Integration: real S01 alpha frame — the persistence gate cleanly reduces the 34 cores and
  the strongest survivors sit at the known superior-parietal vortex; `view_helmholtz` opens
  without error; `Frame` still returns `Cores`/`Sources` consumed by the GUI.

## Non-goals (Phase 2)

Cross-frame trajectory linking (birth/death/merge over time), null-based MinPersistence
calibration, medial-wall collar masking. All deferred.
