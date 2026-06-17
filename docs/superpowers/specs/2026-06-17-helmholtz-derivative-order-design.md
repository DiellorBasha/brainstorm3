# Helmholtz view: time-derivative order selector (Field / Velocity / Acceleration)

**Date:** 2026-06-17
**Status:** approved (design)
**Builds on:** view_helmholtz Dirac/Helmholtz active-frame viewer + persistence detection + tracking.

## Motivation

Show the rate of change of the cortical source field, not just the instantaneous field.
"Velocity" = the per-time-step time derivative of the source vectors (backward difference
over dt); "Acceleration" = the second time derivative. Displayed the same way as the field
today: norm as the cortex colormap, direction as the quiver. This reveals how fast (and,
with acceleration, how non-uniformly) the source vectors are changing.

## Decisions (settled with user)

- Derivative of the **total source field** J (the field already driving the quiver) — needs
  no decomposition to form; just fetch consecutive frames and finite-difference.
- **Detect on the derivative field**: the Helmholtz decomposition, component colormap,
  persistence markers, and tracker all operate on Dᵏ(J), so vortices/sources are those of
  the rate-of-change when an order ≥ 1 is selected.
- Ship a **3-way selector** (Field=0 / Velocity=1 / Acceleration=2) now, not just a toggle.

## Core idea (why it's small)

The only change to the data path is replacing the input field `Jt` (currently
`GetResultsValues` at the current frame) with its k-th backward time derivative `Dᵏ J(t)`.
Every downstream step — `bst_dirac_helmholtz('Frame', Op, Jt, needCores)`, the component
scalar -> colormap, `Ht.Vtot` -> quiver, the per-hemisphere persistence gate, markers, and
the trajectory tracker — runs unchanged on that field. So "detect on the velocity field"
and "norm/quiver display" are inherited, not re-implemented.

Finite differences (consecutive frames, uniform dt):
- k=0:  J(t)
- k=1:  (J(t) - J(t-1)) / dt
- k=2:  (J(t) - 2 J(t-1) + J(t-2)) / dt^2

## Components / files

1. **`toolbox/math/bst_time_derivative.m`** (new, pure): `D = bst_time_derivative(F, dt, order)`
   where `F` is `[n x (order+1)]`, columns = consecutive frames oldest->newest (last col =
   current). Returns `[n x 1]`. order 0/1/2 supported. Unit-tested.

2. **`toolbox/gui/view_helmholtz.m`**:
   - `St.Deriv` (0 default) in the state struct.
   - dispatch: add `'SetDeriv'`.
   - `SetDeriv(hFig, order)`: set `St.Deriv`, clear `St.Cache`, reset in-progress trajectory
     (`St.Tracks=[]; St.LastIT=[]`), `UpdateFrame`.
   - `UpdateFrame`: replace the `Jt = GetResultsValues(...,iT,0)` + smoothing block with a
     fetch of the last `(St.Deriv+1)` consecutive frames (`GetResultsValues` at `iT, iT-1,
     ...`), assembled oldest->newest into `F`; `dt = Time(iT)-Time(iT-1)` from the results
     time vector; `Jt = bst_time_derivative(F, dt, St.Deriv)`; then apply the existing
     eigenmode smoothing to `Jt` once (smoothing is linear, so post-difference == per-frame).
     If `iT <= St.Deriv` (insufficient history) set `Jt = []` -> blank display.
   - Blank-display path: when `Jt` is empty, set cortex Data to zeros, hide quiver + markers,
     readout = e.g. 'velocity: needs >=1 earlier frame'. Return.
   - Readout: prefix the order when > 0 (e.g. the existing count/readout text, with a
     'velocity'/'acceleration' tag) so it's clear what is shown.

3. **`toolbox/gui/panel_helmholtz.m`**: a "Derivative:" selector (Field / Velocity /
   Acceleration), radio buttons matching the Component-radio style, wired to
   `OnDeriv(panelName)` -> `view_helmholtz('SetDeriv', hFig, order)`; add controls to the
   panel struct.

## dt / units / edge cases

- `dt = Time(iT) - Time(iT-1)` (uniform MEG sampling). Velocity units = source/s, acceleration
  = source/s^2; the colormap auto-scales from the data min/max (existing logic), so units need
  no special handling.
- Insufficient history (`iT <= order`): blank field + readout note (no error).
- Frame fetch: `bst_memory('GetResultsValues', St.srcDS, St.srcResult, [], iFrame, 0)` per
  needed frame; cheap (kernel*data column).

## Invalidation

Changing the order changes the displayed field, so `SetDeriv` clears `St.Cache` (the cached
`Ht` is order-specific) and resets the in-progress trajectory. The cache key stays `iT`,
representing "current order".

## Testing (TDD)

- `bst_time_derivative` unit tests: order 0 returns the current (last) column; order 1 =
  `(b-a)/dt` on `F=[a b]`; order 2 = `(c-2b+a)/dt^2` on `F=[a b c]`; on hand-computed numbers.
- Integration in `view_helmholtz`: open on the S01 Dirac source, `SetDeriv` Velocity ->
  colormap = |dJ/dt|, quiver = velocity direction, cores (if markers on) are of the velocity
  field; a frame with insufficient history blanks; Acceleration = second difference.

## Non-goals

Higher orders than 2; centered/other difference schemes (backward only); a separate
acceleration analysis pipeline; smoothing the derivative across time (only the existing
spatial eigenmode smoothing applies).
