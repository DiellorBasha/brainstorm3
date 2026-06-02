# Eigenmode leadfield in the "Compute head model" dialog

**Date:** 2026-06-02
**Author:** Diellor Basha (design captured with Claude)
**Status:** Design — pending review before implementation plan
**Repo:** `research/code/brainstorm3` (branch off `development`)

---

## 1. Motivation

The eigenmode (GBF) source-mapping pipeline already exists (forward composer
`bst_eigenmode_leadfield` + mode-space inverse `bst_inverse_eigenmodes` +
processes, merged to `development`). The only GUI entry point today is the batch
**Process** "Compute eigenmode leadfield". `process_*` files are batch wrappers;
the eigenmode leadfield is a **head model** that must live on the subject and be
available for downstream eigenmode applications, so it belongs in the interactive
**"Compute head model"** dialog.

**Key framing.** The eigenmode leadfield `L̃ = L·Φ` is **not new forward physics** —
it re-expresses the cortical *source space* in the Laplace–Beltrami eigenmode
basis (K modes instead of N vertices). The physics leadfield `L` is still computed
with the user's chosen method (overlapping spheres / OpenMEEG / DUNEuro). Therefore
it belongs with the **Source space** radio group, not the forward-method dropdowns
(where it would wrongly imply a new physics method alongside DUNEuro/OS).

---

## 2. Scope

### In scope
- A new **Source space** radio in `panel_headmodel`: **"Cortex surface harmonics"**,
  with a **"Number of modes"** field (0 = all).
- `ComputeHeadModel` branch that, for this source space, computes the base cortical
  leadfield with the selected physics, composes `L̃ = L·Φ` (reusing the merged
  `bst_eigenmode_leadfield`), and saves **only** the harmonic head model node.

### Out of scope
- No change to the eigenmode math, the inverse, or the eigenspectrum tooling.
- No change to `bst_headmodeler`'s physics (it still computes the base `L`).
- The batch `process_eigenmode_leadfield` stays as-is (composes on an existing head
  model); the dialog computes the base fresh. Both entry points coexist.

---

## 3. User-facing behavior

Source space group becomes:

```
Source space:
  ( ) Cortex surface
  ( ) MRI volume
  ( ) Custom source model
  (•) Cortex surface harmonics     [ Number of modes: 0 = all ]
```

- The "Number of modes" value field is **enabled only** when "Cortex surface
  harmonics" is selected.
- The MEG/EEG/ECOG/SEEG/NIRS physics dropdowns are **unchanged** and selected
  normally — they determine the base leadfield `L`.
- Being a radio, harmonics is mutually exclusive with volume/mixed (eigenmodes are
  cortex-only) — the constraint enforces itself.
- The resulting head model node is labelled **`<base method> | harmonic`**, e.g.
  `Overlapping spheres | harmonic`, `OpenMEEG BEM | harmonic`.

---

## 4. Architecture & data flow

```
  Compute head model dialog (panel_headmodel)
    Source space = "Cortex surface harmonics", Number of modes = K
    Physics = e.g. MEG: Overlapping spheres
        │
        ▼  GetPanelContents → sMethod
    sMethod.HeadModelType   = 'surface'          (base physics is a surface run)
    sMethod.SourceCompression = 'eigenmode'      (new flag)
    sMethod.nModes          = K
        │
        ▼  ComputeHeadModel(iStudies, sMethod)   — harmonic branch
    1. Validate: subject default cortex has eigenmodes (in_tess_eigenmodes). Else error.
    2. Compute BASE in-memory: bst_headmodeler with OPTIONS.HeadModelFile = ''
       → OPTIONS.HeadModelMat = base head-model struct (Gain [nCh×3·nVert], GridOrient,
         SurfaceFile, …). No base node is written.
    3. Compose: bst_eigenmode_leadfield(baseHeadModelMat, Eigenmodes, 'nModes', K)
       → composed struct (Gain = L̃ [nCh×K], isEigenmode=1, Eigenvalues, nModes, SurfaceFile).
    4. Save ONLY the harmonic node: bst_save → headmodel_eigenmode_*.mat;
       Comment = '<base method label> | harmonic'; register via db_template('HeadModel');
       set active; UpdateNode.
        │
        ▼
  headmodel_eigenmode_*.mat   (the only node produced)  →  consumed downstream by
  bst_inverse_eigenmodes / process_eigenmodes_inverse / "Eigenmode source mapping".
```

### Why in-memory base
`bst_headmodeler` writes a node only when `OPTIONS.HeadModelFile` is non-empty; with
`''` it returns the full struct in `OPTIONS.HeadModelMat`. Computing the base
in-memory lets us compose and persist **only** the harmonic node (per the decision
to not save a separate base node), while still reusing the exact physics path.

### Vertex consistency
The base is computed on the subject's **default cortex** (`sSubject.iCortex`), the
same surface that carries the eigenmodes `Φ`, so `L` and `Φ` always share the vertex
count — sidestepping the head-model/eigenmode drift seen in OMEGA.

---

## 5. Components / files

| File | Change |
|---|---|
| `toolbox/forward/panel_headmodel.m` | Add the "Cortex surface harmonics" radio + "Number of modes" field; enable/disable the field with the radio; map it in `GetPanelContents` (`SourceCompression`, `nModes`); add the compose-and-save-only-harmonic branch in `ComputeHeadModel`; label `<base> | harmonic` in `UpdateComment`. Extend the `ctrl` struct. |
| `toolbox/forward/bst_eigenmode_leadfield.m` | **Reused unchanged** (engine). |
| `toolbox/io/in_tess_eigenmodes.m` | **Reused unchanged** (prerequisite check + Φ). |

No other files change.

## 6. Error handling

- **No eigenmodes on the default cortex:** block the harmonic run with a clear
  message ("No eigenmodes on the cortex surface — compute eigenmodes first."). If
  feasible in the panel, also visually disable the radio (or its OK path) when the
  default cortex has none; at minimum, error cleanly in `ComputeHeadModel`.
- **No default cortex:** reuse the dialog's existing "No cortex surface available"
  path.
- **Number of modes > available:** clamp to available (already handled inside
  `bst_eigenmode_leadfield`).

## 7. Testing

- **Pure:** `GetPanelContents` mapping — selecting the harmonics radio yields
  `sMethod.HeadModelType='surface'`, `SourceCompression='eigenmode'`, `nModes=K`;
  the modes field gates correctly. (Panel logic is testable via the documented
  `panel_headmodel('GetPanelContents')` contract; mirror existing panel tests.)
- **e2e (OMEGA):** call `ComputeHeadModel(iStudy, sMethod)` with the harmonic
  source space on a subject that has a cortex + eigenmodes; assert **exactly one**
  new head model node, `isEigenmode=1`, `Gain` is `[nCh × nModes]`, `Eigenvalues`
  present, Comment ends `| harmonic`, and **no base node** was added.
- **Reused:** existing `bst_eigenmode_leadfield` pure tests cover the composition
  math; `bst_inverse_eigenmodes` tests cover downstream consumption.

## 8. Open implementation detail (resolved during planning)

- Exact `GetPanelContents` field names (`SourceCompression` vs a boolean
  `isEigenmodeSpace`) and how `ComputeHeadModel` forces `HeadModelType='surface'` for
  the base `bst_headmodeler` call while branching on the compression flag.
- Whether the panel can cheaply check eigenmode availability at dialog-open to
  disable the radio, or only validates on OK (validate-on-OK is acceptable).
