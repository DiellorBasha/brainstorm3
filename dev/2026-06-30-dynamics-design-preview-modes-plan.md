# Dynamics Design/Preview (Apply) modes Implementation Plan

> **For agentic workers:** Execute INLINE (superpowers:executing-plans) with live checkpoints. Task 1 (filter-apply core) is headless-TDD'd; Tasks 2-3 (clean-cortex launch + Apply UI + sensor view) are built + validated live.

**Goal:** Split the Dynamics atom preview into a default **Design** view (clean cortex, atom impulse response, no real data) and a **Preview** view (Apply toggle: reconstruct the real source over a 4 s window, filter it through the selected atom with `bst_eigenfilter('Analysis')`, and show the filtered source map + the filtered sensor timeseries).

**Architecture:** `view_dynamics` opens a clean `figure_3d` (no real source/timeseries) keeping the `srcResult` link; the overlay gains an `atom-filtered` mode. Apply reconstructs the windowed source via `bst_memory('GetResultsValues')`, filters it with `bst_eigenfilter('Analysis', F, EigenMat, OperatorMat, kernel, kp)` on the atom's operator, paints the magnitude, and forwards `Ffilt` through the head-model gain to a sensor timeseries.

**Tech Stack:** MATLAB R2023b, Brainstorm dev fork. `bst_eigenfilter('Analysis')`, `bst_memory('GetResultsValues'/'GetTimeVector')`, `bst_eigen('Axes')`, head-model gain, `figure_3d`/`view_surface`, `gui_component` toolbar toggle.

## Global Constraints

- No new dependencies; Brainstorm components only; live validation in the Brainstorm session (restart `brainstorm` if it drops; `gui_hide('Dynamics')` before re-show).
- Apply operates on the **selected atom** over the **4 s window at the cursor** (the same window `i_atom_axes` uses). Whole-bank/frame + full-record stay out.
- `Analysis` needs `EigenMat = struct('Phi',ax.Phi,'Lambda',ax.Lambda,'Variant',variant)` + `OperatorMat = struct('Mass',ax.Mass)`, where `ax = i_atom_axes(st, G.Operator)`.
- Operator/source gate: read `R.nComponents`; scalar source → only Geometric/Connectomic enabled in Set-operator; vector source → all four. Required in Apply.
- Source→sensor: `data_filt = Gain × Ffilt` (head-model gain `Gain = [nCh × 3·nV]` from the result's `HeadModelFile`). Verify availability in Task 3; fall back to a warning if no head model.
- `lint` every edited `.m`; commit after each task with `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

## File structure
- Modify: `toolbox/gui/panel_bst_dynamics.m` (filter-apply, Apply toggle, gate), `toolbox/gui/view_dynamics.m` (clean-cortex launch, `atom-filtered` overlay).
- Test: `dev/tests/test_dynamics_apply.m` (Task 1).

---

### Task 1: filter-apply core (windowed source → Analysis → Ffilt) (headless)

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m`; Test `dev/tests/test_dynamics_apply.m`

**Interfaces (produces):**
- `Ffilt = i_atom_filter_field(F, ax, variant, kernel, kp)` — wrap `bst_eigenfilter('Analysis', F, struct('Phi',ax.Phi,'Lambda',ax.Lambda,'Variant',variant), struct('Mass',ax.Mass), kernel, kp)`; return `[]` on `isError`.
- `[iWin, tWin] = i_cursor_window(srcDS, srcResult, secs)` — the sample indices of a `secs`-long window starting at the global cursor (clamped to the recording), from `bst_memory('GetTimeVector')`.

- [ ] **Step 1: failing test** — `dev/tests/test_dynamics_apply.m`: build a synthetic `ax` (scalar, like the substrate tests) + a random `F=[nV×nT]`; `Ffilt = panel_bst_dynamics('i_atom_filter_field', F, ax, 'Laplace-Beltrami', 'heat', struct('lmax',max(ax.Lambda{1}),'t',0.05))`; assert `size(Ffilt)==size(F)`, finite, and that a low-pass kernel reduces high-frequency energy (`norm(Ffilt) < norm(F)` for a non-trivial heat `t`). Assert `i_cursor_window` returns a contiguous in-range sample vector of the expected length for a synthetic time vector + cursor (use a `srcComment`/seam or a stubbed time vector).
- [ ] **Step 2: run → fail.**
- [ ] **Step 3: implement** `i_atom_filter_field` (the `Analysis` wrapper + guard) and `i_cursor_window` (cursor from `GlobalData.UserTimeWindow.CurrentTime`, window `[t0, t0+secs]` clamped; map to samples via the result's time vector).
- [ ] **Step 4: run → pass.**
- [ ] **Step 5: lint + commit** `feat(dynamics): atom filter-apply core (windowed source -> Analysis -> filtered field)`.

---

### Task 2: clean-cortex Design launch + Apply toggle + filtered-source paint (live)

**Files:** Modify `toolbox/gui/view_dynamics.m`, `toolbox/gui/panel_bst_dynamics.m`. **Live-validated.**

- [ ] **Step 1: clean-cortex launch.** In `view_dynamics`'s `i_open_source_figure` (or its caller), open the cortex with `view_surface` (geometry only) instead of `view_surface_data`, and do NOT open the recording-timeseries figure. Keep building the `DynamicsOverlay` with `srcDS/srcResult/iTess/nV` (so Apply can reach the kernel/recording), but `Op='atom'` and no native source paint. Verify the impulse-response preview still works on the clean cortex.
- [ ] **Step 2: `atom-filtered` overlay.** In `view_dynamics`, add `SetFilteredField(hFig, Ffilt, gv, tWin)` (stores `D.AtomField=Ffilt` magnitude-reduced, `D.AtomGV=gv`, `D.Op='atom-filtered'`, `D.tWin`) and an `i_dynamics_overlay` branch that paints `Ffilt(:,iT)` for the cursor frame within the window (reuse the `atom` scatter-paint; the only difference is the data source + the time mapping into the window).
- [ ] **Step 3: Apply toggle.** Add an **Apply** `ToolbarToggle` to the east toolbar (`OnApply`). ON → `i_atom_apply(st)`: `ax=i_atom_axes(st,G.Operator)`; `[iWin]=i_cursor_window(srcDS,srcResult,4)`; `F=bst_memory('GetResultsValues',srcDS,srcResult,[],iWin,0)`; reduce to the operator's row layout; `Ffilt=i_atom_filter_field(F,ax,G.Operator,G.KernelName,G.KernelParams)`; `view_dynamics('SetFilteredField',...)`. OFF → back to `i_atom_preview()` (impulse). Param/operator/seed edits while Apply is ON re-run `i_atom_apply`.
- [ ] **Step 4: operator/source gate.** On select/create, read `R.nComponents` (source result behind the figure); `i_op_enable(nC)` enables the matching Set-operator radios (`nC==1` → Geometric/Connectomic only; else all). Don't default to a disabled operator.
- [ ] **Step 5: live gate.** Launch → clean cortex (no source data/timeseries); add atom → impulse response; **Apply** → the filtered real-source map appears + scrubs over the window; editing params re-filters; Apply off → impulse; incompatible operators greyed. Fix issues.
- [ ] **Step 6: commit** `feat(dynamics): Design (clean cortex / impulse) + Preview (Apply -> filtered source) modes + operator gate`.

---

### Task 3: filtered sensor timeseries (live) — SUPERSEDED, needs its own design

> **2026-06-30 update:** Tasks 1+2 shipped (commits on `development`: filter-apply core; Design/Preview
> modes + operator gate). Task 3 as written below (`Gain × Ffilt`) is **wrong** — per the user, the sensor
> view is a **Dirac-operator-only** feature: transform the imaging kernel into the **Dirac eigenbasis**
> (`nEig × nChannels`, the Dirac dSPM), filter the eigenmodes by `g(λ)`, forward back to channels
> (`D_filt = L_eig·diag(g(λ))·K_eig·D`). It requires launching from the Dirac dSPM kernel, lifting the
> Task-2 Dirac Apply guard, and reusing `bst_dirac`. See design §6 (RESOLVED). **Do this as a fresh
> brainstorm → spec → plan**, not the steps below.

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m`. **Live-validated.**

- [ ] **Step 1: forward to sensors.** In `i_atom_apply`, after `Ffilt`: load the head-model gain `Gain` from the result's `HeadModelFile` (`in_bst_headmodel`/`bst_memory`); `data_filt = Gain × Ffilt` (`[nCh × nWin]`). Guard: if no head model, skip the sensor view + note it in the readout.
- [ ] **Step 2: display.** Show `data_filt` as a sensor timeseries over the window — open/update a lightweight timeseries figure (or a Brainstorm `view_timeseries`-style display) labelled "atom-filtered sensors". Keep it minimal (one figure, refreshed on re-filter).
- [ ] **Step 3: live gate.** Apply ON → the filtered sensor timeseries shows the filter's temporal/spatial effect mapped to channels; editing the atom updates it; Apply OFF closes/clears it.
- [ ] **Step 4: commit** `feat(dynamics): Apply shows the atom-filtered sensor timeseries (Gain x Ffilt)`.

---

## Done criteria
- Default launch = clean cortex + impulse response (no real source/timeseries); **Apply** filters the selected atom's real source over a 4 s window (`bst_eigenfilter('Analysis')`) → filtered cortex map + filtered sensor timeseries; operators incompatible with the source are disabled; `test_dynamics_apply` passes; the live gate passes.

## Risks / notes
- `bst_eigenfilter('Analysis')` `RowMap` must accept the unconstrained (3·nV) / vector / Dirac (4·nV) layouts — verify per operator; guard (skip + warn) like the impulse path.
- `F = GetResultsValues(..., iWin, ...)` over 4 s × `Fs` for unconstrained is sizeable; the window bound keeps it OK.
- Head-model gain may be absent/large; Task 3 guards. The source→sensor choice (gain forward vs imaging-kernel pinv) — gain forward is physically correct; switch to `pinv(K)` only if no head model and the user prefers the imaging-kernel route.
- Big `view_dynamics` launch change — keep the impulse/Design path working throughout; Preview lives behind the Apply toggle.
