# Design: Connection-Laplacian Phase Readout + Unified Cortex Viewer (M3)

- **Date:** 2026-06-04
- **Author:** Diellor Basha (with Claude)
- **Status:** Design approved — ready for implementation plan
- **Branch:** feature/connection-phase-readout (off `development`)
- **Builds on:** `tess_connection_laplacian` (M1), `ConnEigenmodes` axis (M2). See
  `dev/connection_laplacian_integration.md`, `dev/connection_eigenmodes_integration.md`,
  and project memory `connection-laplacian-phase`.

## 1. Goal

A single co-designed milestone with two intertwined deliverables:

1. **Phase-readout engine** — turn a stored connection-Laplacian eigenmode into a
   location-marking phase on the cortex (the project's north star: the first
   non-trivial eigenmode's phase as an angular location coordinate).
2. **Unified layered cortex viewer** — consolidate the three current cortex
   visualizations (tangent-frame `view_tangents`, scalar-eigenmode `view_eigenmodes`,
   and the new connection-Laplacian phase/frames) into one layered surface viewer.

They are co-designed because two viewer layers (intrinsic frames, phase) cannot be
drawn until the readout math exists. The 2D analysis plot `view_eigenmode_spectrum`
is **out of scope** (it is a plot, not a surface view).

## 2. Conceptual foundation (settled in brainstorming)

A connection-Laplacian eigenvector `φ` is, per vertex, a **complex coordinate**
`z_i ∈ ℂ` expressed in nxr/geometry-central's **arbitrary per-vertex frame**
(first-outgoing-halfedge x-axis). Two distinct objects come from it:

- **The 3D tangent field** `w_i = Re(z_i)·e1_i + Im(z_i)·e2_i`, using the per-vertex
  basis `(e1_i, e2_i)`. This is **gauge-independent** — the actual eigenmode vector
  field. nxr's frame is only the *decoder*.
- **A phase** — a circle-valued (S¹) quantity. The raw `arg(z_i)` in nxr's frame is
  noise (the frame is arbitrary per vertex); a meaningful phase requires a *smooth*
  gauge or a reference field.

**The first non-trivial eigenmode (Fiedler) is itself the smooth, globally
consistent, winding field** — this is the payoff of the connection-Laplacian
spectrum. No trivial connection is needed to *manufacture* a consistent frame; the
eigenmode already is one. On a hemisphere (a topological sphere, χ = 2),
Poincaré–Hopf forces the field's singularity indices to sum to 2 — the smoothest
field winds +1 around **each of two poles** (a dipole-like circulation), with phase
undefined exactly at those two points. **That winding is the intrinsic location
handle** (azimuthal position around the poles ≈ a longitude-like coordinate).

**Two readout modes, mirroring Brainstorm's analysis modes:**

- **Within-subject (intrinsic):** location comes from the eigenmode's own winding —
  registration-free. Rendered faithfully and reference-free by the **field
  renderings** (glyphs / streamlines / stripes) and **iso-phase contours**, whose
  level sets follow the field itself (Knöppel stripe/integration), so the winding is
  intrinsic with no external frame.
- **Between-subject:** read the field's angle against the **FreeSurfer-registered
  trivial-connection frame** (the existing per-vertex `TangentFrame`), giving an
  anatomically aligned, cross-subject-comparable phase, with a global-phase
  gauge-fix (anchor at a pole).

**Technical note (winding/index) for the field-vs-reference readout.** When phase is
read as `arg(field / referenceField)`, the winding around a loop equals
`index(field) − index(reference)` inside it. So the reference field's singularities
must be placed **differently** from the eigenmode's, or the windings cancel and the
scalar phase is flat. The FreeSurfer poles satisfy this (they are anatomically
fixed, generally distinct from the eigenmode's intrinsic singularities). The
intrinsic stripe/contour rendering integrates the field directly, so it winds with
the field's own singularities and has no such cancellation.

## 3. Components (implemented in this order)

### 3.1 nxr per-vertex frame export (`nxr-compute` repo) — critical path
geometry-central computes `vertexTangentBasis` internally but the MATLAB binding
exposes only a per-**face** frame (`measure.frame`). Add a binding command (e.g.
`nxr.manifold.measure.vertexFrame`) returning `e1`, `e2`, `normal` as `[nV×3]`
(the exact frame the connection Laplacian was built in), then rebuild and repackage
the macOS plugin (new release asset; update the `PlugDesc` version in
`bst_plugin.m`). This is the only cross-repo piece and gates everything that decodes
an eigenvector to 3D.

### 3.2 Per-vertex FreeSurfer reference frame (Brainstorm)
`TangentFrame` is stored per **face** (`tess_tangents`). Transfer it to per **vertex**
(parallel-transport each incident face's direction into the vertex tangent plane,
area-weighted average, renormalize) so the smooth registration frame lives on
vertices, where the eigenmodes do. Store as a per-vertex variant alongside the face
frame (the `TangentFrame` format already reserves a `Domain` field).

### 3.3 Phase-readout engine (Brainstorm)
A function (e.g. `bst_conn_phase`) that, given a surface's `ConnEigenmodes` + the
nxr per-vertex frame + (for between-subject) the per-vertex FS frame, and a mode
index, returns:
- `Field` — the `[nV×3]` 3D tangent field `w` (gauge-independent).
- `Phase` — `[nV×1]` phase, in the requested gauge (`'intrinsic'` | `'fs'`).
- `Magnitude` — `[nV×1]` `|z|` (→ 0 at singularities; used to modulate/mask phase).
- `Singularities` — vertex indices / locations where the field vanishes.
Default mode = the first non-trivial (Fiedler); default gauge = `'intrinsic'`.

### 3.4 Unified layered cortex viewer (Brainstorm GUI)
One 3D cortex figure (opaque, wireframe, à la `view_tangents`) + a control panel.
Because Brainstorm colors a surface by **one scalar at a time**, layers split into:

- **Surface-scalar driver (choose one):** `none/plain` · scalar-LBO-eigenmode
  (lever-driven, reuses `view_eigenmodes`/`panel_eigenmodes`) · **connection-mode
  phase** (cyclic colormap, magnitude-modulated) · connection-mode magnitude ·
  stripes (`cos(k·phase)`, B/W colormap).
- **Overlays (toggle any combination):** extrinsic FS tangent frames (glyphs) ·
  intrinsic connection frames (glyphs) · connection-mode field glyphs / streamlines ·
  **iso-phase meridian contours** · singularity markers. All overlays reuse the
  `view_tangents` line-overlay machinery.

**Controls:** connection-mode selector (Fiedler default, browsable 1…k) ·
within-/between-subject toggle (intrinsic vs FS gauge) · glyph density/length/width
(existing `view_tangents` keys) · colormap/contour settings.

**Data flow (per surface):** `bst_conn_eigenmodes_ensure` (M2) +
`bst_eigenmodes_ensure` (scalar) + nxr per-vertex frame + per-vertex FS
`TangentFrame` → for the selected mode, `bst_conn_phase` → render per layer
selection.

### 3.5 Validation
On a real hemisphere: the Fiedler field winds +1 around each of its two
singularities (index check); the intrinsic phase sweeps monotonically with an
azimuthal coordinate around the poles; the FS-gauge phase is smooth and aligned to
the reg-sphere longitude. Quantify "phase ⇒ location" (e.g. correlation of phase
with reg-sphere azimuth away from the poles).

## 4. Rendering decision

A bare phase colormap is necessary-but-not-sufficient: phase is **cyclic** (needs a
cyclic colormap, else a false seam at ±π) and **undefined at singularities** (needs
magnitude modulation). The field's winding is most faithfully shown reference-free by
the field/contour renderings. Chosen **headline**: **cyclic phase colormap,
magnitude-modulated**, **plus iso-phase meridian contours** (the coordinate grid that
makes "winds ⇒ location" legible). **Complementary layers:** field glyphs /
streamlines, singularity markers, and stripes as an aesthetic alternative to
contours. Within Brainstorm's scalar-colormap path, phase / stripes / contour-bands
are all achievable by choosing the scalar transform + colormap; glyphs/streamlines
are the line-overlay layer. nxr already provides `stripe_patterns`, `isolines`, and
`streamlines` to leverage.

## 5. Testing

Function-style tests under `dev/tests/` (MATLAB MCP), real 20484-vertex cortex,
skipping cleanly if absent (existing resolver idiom):

- **nxr per-vertex frame:** exported `e1`/`e2`/`normal` are unit, orthonormal,
  `e1×e2` aligned to the vertex normal; `[nV×3]`.
- **FS face→vertex transfer:** per-vertex frame is unit/orthonormal and smooth
  (small angular change vs incident face frames away from singularities).
- **Phase engine:** the decoded 3D field is gauge-independent (matches `smoothVertex`
  for the Fiedler up to sign/scale, as a cross-check); phase winds ±1 around each
  singularity (discrete index); intrinsic vs FS phase agree up to a smooth gauge
  difference; magnitude → 0 at singularities.
- **Viewer smoke:** each surface-scalar driver and each overlay renders without
  error and toggles on/off (following `test_view_tangents`).

## 6. Files (indicative)

| Path | Change |
|---|---|
| `nxr-compute` `src/` + `bindings/mex` + `+nxr/.../measure/vertexFrame.m` | New: per-vertex frame export; rebuild/repackage plugin. |
| `toolbox/core/bst_plugin.m` | Bump nxr-compute `PlugDesc` version/URL to the rebuilt plugin. |
| `toolbox/anatomy/tess_tangents.m` (or a helper) | Per-vertex FS frame (face→vertex transfer). |
| `toolbox/math/bst_conn_phase.m` | New: phase-readout engine. |
| `toolbox/gui/view_cortex_frames.m` + a control panel | New: unified layered viewer (final name to confirm at plan time). |
| `toolbox/gui/view_tangents.m`, `view_eigenmodes.m` | Refactor glyph/colormap cores into reusable layer modules (no behavior change). |
| `dev/tests/test_*` | New tests per §5. |

## 7. Staging into sub-plans

This spans three subsystems, so implementation is split into sequential sub-plans,
each independently testable:

- **Plan A — nxr per-vertex frame export** (`nxr-compute` + plugin rebuild + the
  binding/parity test). Unblocks everything.
- **Plan B — readout math** (per-vertex FS frame + `bst_conn_phase` + phase tests).
- **Plan C — unified viewer** (layer-module refactor of the existing viewers + the
  new shell/panel + the connection phase/frame/contour layers + smoke tests).

## 8. Out of scope / deferred

- Cross-subject *statistics* on phase (group analysis) — this milestone delivers the
  aligned per-subject phase; group-level use is later.
- The 2D `view_eigenmode_spectrum` plot (stays separate).
- Process/pipeline wrappers for batch phase computation.

## 9. Open questions / risks

- **Intrinsic gauge choice** for the within-subject scalar phase (vs. relying on
  stripes/contours, which are reference-free): the exact construction
  (e.g. integrate the field à la Knöppel, or transport-from-seed) is to be chosen and
  validated in the viewer. The winding is intrinsic regardless of choice.
- **Singularity detection** on a discrete mesh (where `|z|→0`): needs a robust
  index/zero-finding rule; affects markers and phase masking.
- **nxr frame ↔ stored eigenvector alignment:** the exported per-vertex frame must be
  the exact gauge the stored `ConnEigenmodes` were computed in; verified by the
  gauge-independence cross-check (decoded Fiedler vs `smoothVertex`).
