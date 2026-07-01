# Dynamics panel — frame core (sub-project B) — design

**Date:** 2026-06-30
**Status:** design (approved; to be implemented in a fresh session)
**Depends on:** sub-project A (panel substrate cleanup, `panel_bst_dynamics.m` 1861→868, on local `development`).
**Part of:** the Dynamics-portal program (A→B→C→D). See [[dynamics-portal-optimization]].

---

## 1. Motivation

The panel's atom list is currently N *independent* single-filter atoms. The frame machinery
(`bst_eigenwavelet`: `Design`/`Evaluate`/`Bounds`/`Analysis`/`Synthesis`) — including tight frames
and frame-bound robustness — exists but is **never called by the panel**. Sub-project B makes the
bank a genuine **frame**: it surfaces frame robustness (coverage + bounds), lets the user generate a
tight (itersine) frame in one click, and makes the interactive Apply/Preview **instant** by caching
the windowed-source projection. This is the conceptual heart of the program: "atom filterbanks and
robust frames to analyze source maps as efficiently as possible."

## 2. Decisions (from brainstorming)

- **Frame model = atoms-as-frame (implicit).** The atoms sharing an operator ARE the frame; no new
  data object. Frame bounds are a property of the atoms' kernel spectral responses `g_m(λ)`, which is
  **seed-independent**: `S(λ) = Σ_m |g_m(λ)|²`, `A = min_λ S(λ)`, `B = max_λ S(λ)`, tightness `= B/A`.
- **Robustness surfaced as a coverage plot + numbers:** a small `S(λ)` plot (x-axis = mm spatial
  scale via `2π/√λ`) beside the `A`/`B`/`B/A` values. Shows *where* coverage gaps are, not just that
  they exist.
- **"Design tight frame (itersine)" = replace with N members, N adjustable.** Replaces the bank with
  a clean itersine tight family of N members spanning the operator's spectrum (preserves the
  tightness guarantee `S(λ)≈const`); confirms first if the bank is non-empty. All members seeded at a
  default vertex; re-localize individually.
- **Itersine integration = first-class registry kernel (Approach 1).** Add
  `bst_eigfilter_design_itersine.m` so a generated member is an *ordinary* kernel atom
  (`KernelName='itersine'`, `KernelParams` carrying `member`/`Nf`/`lmin`/`lmax`). Bounds, impulse
  paint, Apply, and save/load all compose through the existing paths — no special-casing.
- **B/C boundary:** B ships (1) the live bounds/coverage readout, (2) Design-tight-frame generate,
  (3) cached-projection instant Apply for the **selected** atom. Multi-member frame analysis
  (scalogram, per-member energy, synthesis residual, `JTVAtoms`, `process_`) stays in **C**, and will
  reuse the projection cache B builds.

## 3. Components & files

### 3.1 New: `toolbox/eigen/eigfilter/bst_eigfilter_design_itersine.m`
An itersine tight-frame **member** kernel. Given `KernelParams` with fields `member` (1..Nf), `Nf`,
`lmin`, `lmax`, returns a handle `g(λ)` = the `member`-th itersine window. The formula is the one
already in `bst_eigenwavelet('Design','itersine')` (half-cosine, overlap 2):

```matlab
overlap = 2;
scale   = lmax / (Nf - overlap + 1) * overlap;
kf      = @(x) sin(0.5*pi*(cos(pi*x)).^2) .* (x >= -0.5 & x <= 0.5);
g       = @(l) kf(double(l(:))/scale - (member - overlap/2)/overlap) ./ sqrt(overlap) .* sqrt(2);
```

Registry integration:
- `bst_eigfilter_kernel('list')` picks it up automatically (dispatch to `bst_eigfilter_design_<name>`).
- `bst_eigfilter_kernel('info','itersine')` returns `domain='static'` and a band-pass sign class
  (diverging/`stat2` colormap for the impulse paint — an itersine member is spectrally band-pass, so a
  seeded impulse is spatially oscillatory, not one-signed).
- `bst_eigfilter_controls('Sliders'/'ToKernel')`: itersine's parameters (`member`,`Nf`,`lmin`,`lmax`)
  are set programmatically by the generate action, NOT via the contextual scale/rate sliders. When an
  itersine atom is selected, the Atom-section sliders are hidden/disabled and the readout shows
  `itersine · member ii/Nf`. (A user does not hand-tune a single tight-frame member; they regenerate
  the family with a different N.)

### 3.2 Panel: `panel_bst_dynamics.m`
- `i_frame_response(ax, st)` → `struct('lam',..,'S',..,'A',..,'B',..,'tightness',..,'nMembers',..)`.
  Gathers, for every atom whose `Operator` equals the current operator, its `g(λ)` (built via
  `bst_eigfilter_kernel(KernelName, kp)` with `kp.lmax` from `ax`), evaluates on a dense λ-grid
  spanning the operator spectrum, sums the squares. Atoms on other operators are excluded.
- `OnDesignFrame()` — the generate action. Reads N from the spinner; if the bank is non-empty,
  `java_dialog('confirm', ...)`; builds N atoms `i_default_atom('itersine', kp_ii, seed, surf, label, op)`
  with `kp_ii = struct('member',ii,'Nf',N,'lmin',lminPos,'lmax',lmax,'vals',[])` on the current
  operator, seeded at `ax.GlobalVertices{1}(1)`; `UpdateAtomList`; refresh bounds; select member 1.
- `i_frame_refresh()` — recompute `i_frame_response` and repaint the coverage plot + labels. Called
  from `OnCreateAtom`, `i_atom_writeback` (param edit), `AtomsListValueChanged_Callback` (select),
  `OnSetOperator`, `OnDesignFrame`, and atom delete.
- Cached-`C` projection in the Apply path (see §5).

### 3.3 Response view: **reuse + extend** `toolbox/gui/view_eigfilter_response.m`
No new plot code or in-panel axes. `view_eigfilter_response` is the design-system view for eigenfilter
spectral responses and **already has a bank mode**: called with a struct `g` carrying
`g.Kernels` (cell of `@(l)` handles), `g.Active` (highlighted index), and `g.OnSelect(j)` (curve
click callback), it plots every member's `g_m(λ)` with the active one highlighted and click-to-select.
It is a tagged create-or-update companion figure (`EigfilterResponse`), so repeated calls update in
place; `view_eigfilter_response('close')` closes it.

Extend it (backward-compatibly) with two optional bank-mode fields:
- `g.Coverage` (bool) → also plot `S(λ) = Σ_m |g_m(λ)|²` as a bold line, and a faint horizontal
  reference at `mean(S)` so gaps/peaks read against "flat".
- `g.Bounds` (struct `A`,`B`,`tightness`) → annotate the axes title/corner with `A / B / B/A (✓)`.

The panel drives it: `view_eigfilter_response(struct('Kernels',{gCell}, 'Active',curAtom,
'OnSelect',@(j)SetSelectedAtom(j), 'Coverage',true, 'Bounds',bnd), lambdas, 'Frame coverage')`.
The `OnSelect` reuse gives click-a-curve → select-that-atom for free.

### 3.4 GUI (in `CreatePanel`): a bordered **Frame** section
Placed in the SOUTH stack below the existing **Atom** section — **numbers + controls only; the curve
lives in the reused response view**:
- three labels: `A`, `B`, `B/A` (tightness; a ✓ suffix when `|B/A - 1| < 0.05`);
- an `N` spinner (integer, default 6, range 2..24);
- a `Design tight frame` button (`@OnDesignFrame`);
- a `Show coverage` toggle that opens/updates (or closes) the `view_eigfilter_response` figure — the
  same open/close pattern `panel_eigenfilter_options` uses for its response figure.

Layout sketch (docked panel):

```
├ Frame ────────────────────────────┤
│  A 0.98   B 1.02   B/A 1.04 ✓     │
│  N [6 ▾]  [Design tight frame]    │
│  [ ] Show coverage                │
└───────────────────────────────────┘
```

The `S(λ)` coverage curve + per-member `g_m(λ)` tiles render in the companion
`view_eigfilter_response` figure (bank mode, `Coverage=true`), updated live on every bank change.

## 4. Bounds data flow (design-time, no data)

On any bank change (`i_frame_refresh`): resolve the current operator's `ax` (via the existing
`i_atom_axes(st, variant)`, already cached per `variant|surface`); call `i_frame_response`; set the
docked A/B/tightness labels; and, if the `Show coverage` toggle is on, update the reused
`view_eigfilter_response` figure (bank mode with `Coverage=true`, `Bounds=<A,B,tightness>`,
`Active=curAtom`, `OnSelect=@SetSelectedAtom`). Pure spectral math on ≤ a few hundred λ — instant, no
progress bar. Empty bank → labels `—`, no curve; single atom → the one-member `S(λ)` with tightness `—`
(undefined for a non-covering single band).

## 5. Cached-projection instant Apply

**Problem:** `i_atom_apply` currently, on every slider settle in Apply mode, re-fetches the windowed
source (`bst_memory('GetResultsValues')`), reduces it, and re-projects it onto the eigenbasis inside
`bst_eigenfilter('Analysis')`. That is the slow, progress-barred step.

**Refactor:** split the Apply into *project-once* and *apply-gain*.
- **Project (cached):** for the key `(srcResult, iWin-window, operator)`, compute once
  `F = i_paintable_scalar(GetResultsValues(...), nV)` (reduced scalar field over the 4 s window) and its
  per-hemisphere modal coefficients `C{h} = manifold_ft(Phi{h}, B{h}, F(gv{h},:))`. Cache on
  `getappdata(0,'DynamicsApplyCache')` as `struct('key',..,'C',{C},'nV',..,'iWin',..)`.
- **Apply-gain (cheap):** `Ffilt(gv{h},:) = manifold_ift(Phi{h}, h(λ).*C{h})` with `h = g(λ)` from the
  selected atom's kernel; reduce + normalize + paint. No `GetResultsValues`, no re-projection.
- **Invalidation:** the cache key includes the source result, the window sample indices, and the
  operator variant. It does NOT include the seed or the kernel params — Apply filters the whole field,
  so those never invalidate `C`. Changing the cursor window or operator rebuilds `C`.

`OnParamSettle` in Apply mode therefore calls only the apply-gain path → instant re-paint.

This cache is the substrate sub-project C reuses: applying **all** frame members is `Phi·(h_m(λ).*C)`
for each `m` from the same `C` (the scalogram), with no extra projection.

## 6. Operator scope & edge cases

- Bounds/coverage + generate: **all** panel operators (LB, LB-Connectome, Connection, Dirac) — the
  computation is `Lambda`-only.
- Cached instant Apply: **scalar operators only** (LB, LB-Connectome), matching the current Apply
  scope. Dirac/Connection Apply keeps the existing "scalar-only for now" guard (sub-project D).
- Mixed-operator bank: `S(λ)` is computed over the current operator's atoms only; the readout notes
  `k/N atoms on <operator>` when the bank is mixed.
- Empty bank → `—` labels, empty plot. Single atom → its `S(λ)`, tightness `—`.
- Itersine `lmin`: use the smallest positive eigenvalue `lminPos = min(λ|λ>1e-9)` (λ=0 is the constant
  mode) for the frame's lower range, matching how the panel already derives mm bounds.

## 7. Testing

**Headless (pure, MCP-run by controller):**
- `bst_eigfilter_design_itersine`: build all Nf members' `g(λ)` on a dense λ-grid; assert
  `S(λ)=Σ g_m² ` is ≈constant across `[lminPos, lmax]` (tightness `B/A < 1.05`) for a few N (4, 6, 12).
- `i_frame_response`: for a known 3-kernel set on a synthetic Lambda, assert A/B/tightness/nMembers.

**Live (MCP, controller live pass):**
- On the sub-MTL0002 Dirac session: `Design tight frame` (N=6) → coverage plot reads flat, tightness
  ✓; screenshot.
- Add a mismatched hand kernel → tightness degrades visibly; screenshot.
- Apply mode on a scalar operator: drag the kernel slider → re-paint is instant (cached-C), no
  progress bar; screenshot before/after.

## 8. Out of scope (→ C / D)
- Applying the whole frame to the source → scalogram / per-member energy / synthesis residual;
  `JTVAtoms`; a `process_` for the whole series. (C — reuses B's projection cache.)
- Dirac/vector real-source Apply and the filtered-**sensor** view. (D.)
- Per-member `h_m(λ)` overlay on the coverage plot (the richer "Coverage + per-member" option).
- Hand-tuning individual itersine members (regenerate with a different N instead).
