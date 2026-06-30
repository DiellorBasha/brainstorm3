# Atom tool in panel_bst_dynamics — design

**Date:** 2026-06-30
**Status:** design (approved direction; pending spec review → plan)
**Depends on:** the eigfilter/JTV-atom stack — `bst_eigen('Axes')`, `bst_eigenfilter('Atom')`, the js kernels, `bst_dynamics('Levelset')` (all shipped to origin/development by 075157ca).

---

## 1. Motivation

`panel_bst_dynamics` and the standalone atom designer are the two ends of one pipeline that grew up apart. The panel references the eigfilter/atom engine **zero times** — it only manages markers (`bst_dynamics`), the geodesic disk (`bst_geodesic_tool`), and the differential-maps overlay (`view_dynamics`). Its "Navigator" section (Frequency / Source / Scale axis blocks) tries to *navigate* localization axes that aren't cleanly defined, and mostly does nothing.

The unifying insight: a **dynamics atom IS a thresholded localized filter**. You design an eigfilter atom (kernel + temporal/spatial support, localized at a seed), then **threshold** it into its hard indicators — the **Scout** (spatial level set) and the **Event** (temporal level set). The `scale`/`freq` axes that felt unclean are simply the kernel's spatial (`λ`) and temporal (`ω`) support — *parameters*, not things to navigate. The existing geodesic Region tool is already the degenerate case (a static `heat` kernel thresholded spatially → a Scout); the atom tool generalizes it to the whole kernel family.

This spec makes the panel an **authoring tool** for that: a filter tool you select, set its contextual parameters, localize on the cortex, and threshold-and-store.

## 2. Scope

**In:** an "atom tool" in `panel_bst_dynamics` — filters exposed as tools with contextual parameters, click-to-localize on the cortex, live preview via the eigfilter realiser, and threshold→store into a Scout+Event marker. It **replaces the Navigator section** and **subsumes the geodesic Region tool** (`heat` preset).

**Out (deferred / untouched):**
- Detection over a recording via eigfilter-Analysis (a separate future effort).
- The differential maps (Helmholtz/Hodge via `view_dynamics`) — unchanged.
- The standalone atom designer — stays separate, as the full filter-reference lab.

## 3. UI — filters as tools, contextual parameters

The `jNav` "Navigator" panel is replaced by a compact **Atom** section. Filters are surfaced **as tools**, the way `bst_geodesic_tool` is — selecting a filter makes it the active tool and shows **only that filter's parameters** (contextual, morphing per kernel — never fixed persistent Frequency/Scale blocks).

- **Filter selector** — a compact selector of the kernel "tools" (heat, gabor, travwave, resonator, diffusion, …), grouped (static / dynamic-ts / dynamic-js) like the designer's dropdown.
- **Contextual parameters** — for the selected filter, render only its parameters with their labels/units/ranges (e.g. `heat` → Scale (mm); `gabor` → Scale + Freq (Hz) + BW; `resonator` → Freq + Q; `travwave` → Speed + RidgeW). Compact rows that change with the filter. The two "main axes" the user cares about — **temporal-frequency support** and **spatial-scale support** — are always among them where the kernel has them; extra params (Q, ridge width) appear only for the kernels that use them.
- **Localize toggle** — a toggle (replacing `jRegionTool`'s surfacing) that arms click-to-seed on the cortex: toggle on → click a vertex = the atom's seed. This is the geodesic-tool gesture, generalized.
- **Threshold** + **Store** — a threshold control and a Store action (the second step).

**Live preview:** with the tool armed and a seed placed, the realised atom field shows through the panel's `DynamicsOverlay` (the same overlay the panel already drives), updating as the contextual params change — the designer experience, inside the panel.

## 4. DRY — a shared per-kernel control spec

The designer already encodes "which params a kernel exposes, their labels/units/ranges, and how slider values map to kernel params" (`i_config_sliders` + `i_phys2kernel`). The panel needs the *same* knowledge. This is the real second consumer that justifies extracting it:

- **`bst_eigfilter_controls(kernelName)`** (new, in `eigen/eigfilter/`) returns the per-kernel **control spec**: a list of `{name, label, unit, min, max, default}` plus the **`tokernel(vals, lmax)`** mapping → the kernel param struct. Pure data + a mapping; no GUI toolkit.
- The **designer** (MATLAB) renders the spec with its uicontrol sliders; the **panel** (Java/Swing) renders it with `gui_component` fields. Each GUI keeps its own widgets; the *param definitions and the value→param mapping live in one place*.
- The designer's `i_config_sliders`/`i_phys2kernel` are refactored to consume `bst_eigfilter_controls` (so the two GUIs cannot drift). This also fixes the earlier "no `otherwise`/8 static kernels untunable" gap by giving every kernel a spec entry.

## 5. Engine wiring (consolidation around eigfilter)

The panel gains the engine the designer already rides:
- `ax = bst_eigen('Axes', …)` — bound to the panel's recording (its `DataFile` Time axis); the panel already has `bst_dynamics('Axes')` which now delegates to `bst_eigen('Axes')`, so this is in place.
- Live field: `[W,gv] = bst_eigenfilter('Atom', ax, kernel, params, seed)` → painted through the existing `DynamicsOverlay`.
- Store: `LS = bst_dynamics('Levelset', W, gv, threshold)` → Scout (`LS.scoutVertices`) + Event (`LS.eventSamples` → time window). The geodesic/heat disk is one kernel among many; `bst_geodesic_tool` is retired (its click+overlay role absorbed by the atom tool, `heat` preset).

## 6. Data model

A stored atom keeps **materialized Scout + Event** (so the DB tree, scouts, and event consumers work unchanged) **plus its generator** — `{kernel, seed, params, threshold}` — recorded on the atom group as its spatial/temporal-scale metadata. The realised field is computed lazily for preview and **never stored** ("atom = sparse reference, not a data copy"). The `scale`/`freq` axes are read from `params`, which is what cleans them up. The `db_template('atomgroup')` gains generator fields (`KernelName`, `KernelParams`, `Threshold`); `times`/`band`/`vertices`/`scale`/`region` continue to hold the materialized markers.

## 7. Behavior / flow

1. Select a filter tool (e.g. `gabor`) → its contextual params appear; `Localize` arms click-to-seed.
2. Click a cortex vertex → seed; the atom previews live through the overlay; adjust contextual params (Freq/Scale/…) → preview updates.
3. Set **Threshold** → the Scout/Event level sets update in the preview.
4. **Store** → `Levelset` materializes Scout + Event; a new atom (with its generator params) is added to the panel's atom table and the cortex.

## 8. Testing

- `bst_eigfilter_controls`: per kernel, returns a spec whose `tokernel` produces the same param struct the designer's `i_phys2kernel` produced before the refactor (regression: designer output unchanged for existing kernels).
- Atom realise+threshold path (headless, controller-run live): for a kernel + seed + params on a real cortex, `bst_eigenfilter('Atom')` → `bst_dynamics('Levelset')` yields a non-empty Scout and Event; threshold monotonicity (higher threshold → subset Scout).
- Data model: a stored atom round-trips through `bst_dynamics('Save'/'Load')` carrying `KernelName/KernelParams/Threshold`; the materialized Scout/Event match a fresh realise+levelset of the generator.
- GUI (controller live, manual): the Navigator is gone; selecting filters shows contextual params; click-to-seed + threshold + Store adds a marker; the `heat` filter reproduces the old geodesic-disk result.

## 9. Risks / notes

- Per the repo's hard policy, GUI/live-MATLAB validation is the controller's job (subagents do static checks only).
- The panel is Java/Swing (`gui_component`/`JPanel`) — the contextual-param rendering is the main new GUI work; the engine/data pieces are reuse.
- Retiring `bst_geodesic_tool`: confirm no other caller depends on it before removal (it is currently surfaced only via the panel's `jRegionTool`).
- `bst_dynamics('AttachRegion')` (the old "capture region → atom") is superseded by the threshold→store path; reconcile or retire it.
- The Frequency block's prior **band-focus** (band preset ↔ recording-bandpass display, the `NotifySelection`/`i_driving` focus state shipped 2026-06-25) is **repurposed** into the kernel's frequency-support control. The recording-bandpass *display* side-effect is dropped this phase (it can return later as a recording-display option, independent of any atom). Per the user: the frequency control's job is selecting/adjusting the atom's frequency support.
