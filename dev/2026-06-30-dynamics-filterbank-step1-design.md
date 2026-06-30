# Dynamics panel → eigenwavelet filterbank (step 1) — design

**Date:** 2026-06-30
**Status:** design (approved; pending spec review → plan)
**Depends on:** the panel atom tool already in `panel_bst_dynamics` (`i_atom_realise`/`i_atom_preview`/`SetAtomField` source overlay, the `panel_eigenfilter_design` contextual sliders, the geodesic seed-pick) — all committed on `development`.

---

## 1. Motivation / conceptual shift

Reorient the Dynamics panel from a "thresholded-markers" tool toward an **eigenwavelet filterbank** interface: atoms are **filters** you build and apply to the source maps (and, later, the recordings) while exploring.

- An **atom = a filter**: `{operator, kernel, seed, params}` — a reusable eigenwavelet reference.
- The **"threshold → Scout+Event" definition is retired** as the default. (It can return later as an optional export; it no longer *is* the atom.)
- The **atom list becomes a filterbank** — the first step toward filterbanks (atom groups, like Recordings event groups) and frames.
- Creating or selecting an atom **paints its realised spatiotemporal field on the source map** (reusing the `SetAtomField` overlay). Source-map preview only this step; the sensor-timeseries projection is a later step.

## 2. UI changes to `panel_bst_dynamics`

1. **Stacked panes.** Flip the tree|list `JSplitPane` from `HORIZONTAL_SPLIT` (side-by-side) to `VERTICAL_SPLIT` (top/bottom).
   - **Top** = the **atom list** (`atom1`, `atom2`, … — the filterbank), one row per atom.
   - **Bottom** = a **read-only detail line** for the selected atom: kernel · seed vertex · key params.
2. **"+" Create-atom button.** Copy `panel_scout`'s "Create a scout" toolbar button **verbatim** — same icon and tooltip style, adapted to "Create an atom". (Lives on the panel's existing action toolbar.)
3. **East toolbar gains the atom actions; the Atom section loses them:**
   - **Localize** — arms click-to-seed (geodesic seed-pick) to re-seed the *selected* atom; on a click the atom's seed updates and the preview refreshes.
   - **Threshold** — opens a small numeric input (default 0.5); de-emphasized, for the optional level-set only (not the atom's identity).
   - **Save** (the existing toolbar Save button) — now **persists the atom table (filterbank)** to its `dynamics_*.mat` on disk.
4. **Atom section = parameters only** — Operator (Connectomic | Geometric toggles) + Filter (grouped combobox) + the contextual parameter sliders. No Localize/Threshold/Store controls here.

## 3. Create-atom flow

Clicking **+ Create atom**:
1. Builds (or reuses cached) eigen-axes for the linked surface (`i_atom_ensure_axes`).
2. Creates a default **Diffusion** atom: `KernelName='diffusion'`, default params, **seed = centroid-nearest vertex** of the first eigenbasis block (the designer's default-seed rule).
3. Appends it to the atom table as a new `atomgroup` labelled `atomN` (auto-incremented).
4. Selects the new atom → the Atom section loads its kernel/params, and the **source preview updates immediately** (`i_atom_preview`), so the spatiotemporal effect on the source map is visible at once.

Selecting an existing atom in the list does the same load+preview. Editing the Atom-section controls re-realises the selected atom and updates the preview + the bottom detail line.

## 4. Data model

- Each atom = one `atomgroup` carrying its **generator** (`KernelName`, `KernelParams`, `vertices`=seed) + provenance (`SurfaceFile`). The `Threshold`/`region`(Scout)/`times`(Event) fields are **left unset** by default (no thresholding at create time).
- The atom table (`st.T`, a `dynamicsmat`) is the filterbank; `Save` persists it. `AtomFromKernel` (realise→threshold→Scout+Event) is **no longer on the create path** — it remains available for the optional Threshold export only.
- The selected atom's index is tracked in `st` (e.g. `st.curAtom`), driving the Atom-section load and the preview.

## 5. Reuse / what changes

Reused as-is: `i_atom_realise`/`i_atom_preview` (realise + normalize + `SetAtomField`), the `panel_eigenfilter_design` contextual sliders + operator toggles + grouped filter list, `bst_geodesic_tool` seed-pick, `i_atom_ensure_axes`.

Changed in `panel_bst_dynamics`: the split orientation + pane contents; a `+`/Create-atom toolbar button + `OnCreateAtom`; `Localize`/`Threshold` move to the toolbar (the Atom section drops them); `Save` repurposed to persist the table; atom-list build/select wired to the generator atoms.

## 6. Scope

**In:** the four UI changes + the create-default-atom flow + source-map preview + Save-persists-filterbank.
**Out / untouched this step:** sensor-timeseries projection (later); filterbank *groups* + frames (later); the legacy detection/record toolbar actions (Detect/Show/Clear) stay in place untouched — reconciling them belongs with the filterbank step.

## 7. Testing

- Headless: creating a default atom appends a generator `atomgroup` (`KernelName='diffusion'`, seed set, threshold unset); the atom-list/detail strings render from a generator atom; Save round-trips the table.
- Live (controller + user): open a dynamics session → **+ Create atom** adds `atom1` and the diffusion pattern appears on the source map; selecting/editing re-previews; **Localize** re-seeds; **Save** writes the file; the split is stacked; the Atom section holds only parameters.

## 8. Risks / notes

- Swing layout (split flip, toolbar button) is built live; the `+` button is copied from `panel_scout` verbatim.
- Reconcile the existing `OnStore`/`jStore`/`jLocalize`/`jThresh` (currently in the Atom section from the prior atom-tool work) — remove them from the Atom section and route their behavior to the new toolbar buttons.
- Brainstorm session has been unstable this session; live validation may need restarts.
