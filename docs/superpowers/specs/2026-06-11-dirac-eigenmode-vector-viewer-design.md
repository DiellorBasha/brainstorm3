# Dirac Eigenmode Vector-Field Viewer — Design

- **Date:** 2026-06-11
- **Status:** Draft for review
- **Author:** Diellor Basha (with Claude)
- **Related:** `view_eigenmodes.m`, `view_leadfield_vectors.m`, `bst_dirac.m`, `tess_eigen.m`,
  `panel_eigenmodes.m`, `panel_surface.m`, `bst_figures.m`, `tree_callbacks.m`, `process_eigenmodes.m`

## 1. Goal

Replace the deprecated scalar Laplace–Beltrami `view_eigenmodes` with a **standalone vector-field
viewer** for **Dirac (D2) eigenmodes**, launched from the `eigen_` DB node. The viewer loads an
`eigen_*.mat` file, renders each eigenvector as an ambient 3D **quiver field** on the cortex, and
lets the user **cycle through modes with keyboard shortcuts** — modeled on how
`view_leadfield_vectors` cycles through sensor channels.

As part of the same work, **fully retire the `panel_eigenmodes` display-lever subsystem** (the
"spatial scale (eigenmodes)" tab and all its couplings), which existed only to drive the old scalar
viewer.

## 2. Background and conceptual framing

The current `view_eigenmodes(SurfaceFile)` loads surface-stored scalar LBO modes via the legacy
`in_tess_eigenmodes`, registers a transient scalar Source result, and drives the `panel_eigenmodes`
lever to superpose a range of modes (signed `stat2` colormap). This scalar path is deprecated.

The eigenbasis we use now is stored as an `eigen_` DB node under the parent surface
(`db_template('eigenmat')`): `Phi{1,2}` per hemisphere, `Lambda{1,2}`, `GlobalVertices{1,2}`,
`Variant`, `Provenance`. For the **Dirac** variant, each `Phi{hh}(:,k)` is a `[4·nVh × 1]`
**quaternion field** — four real components `(w, i, j, k)` per hemisphere vertex. The ambient 3D
vector at each vertex is the **quaternion vector part** `(i, j, k)` → `(x, y, z)`, dropping the `w`
slot. This is exactly the extraction `bst_dirac/local_reconstruct` performs to map mode coefficients
back to a cortical 3-vector field; here we apply it to a single eigenvector.

A 3D vector anchored at a surface vertex is a section of the pullback bundle `f*Tℝ³`: a full copy of
ℝ³ over each point, **not** the 2D tangent plane. Display is therefore purely ambient and Cartesian —
the quaternion vector-part components are drawn directly as `quiver3` `U,V,W`, with no projection and
no transform. Dirac eigenvectors are **real** (the operator is a real 4×4-block representation), so no
complex handling is needed (Connection Laplacian, a future variant, is complex).

## 3. Non-goals

- **No differential analysis** (divergence, curl, flow) — display only.
- **No LBO/Connection rendering yet.** The viewer dispatches on `EigenMat.Variant`; only Dirac is
  implemented. LBO and Connection Laplacian raise a clear `bst_error` and are added later as branches.
- **No transient DB result and no `panel_eigenmodes`.** The viewer is fully self-contained
  (`view_surface` + `quiver3` + own keyboard), like `view_leadfield_vectors`.
- **No changes to `tess_eigen`, `bst_dirac`, or the inverse.** The eigenbasis is consumed as stored.
- **No removal of the legacy compute path.** `process_eigenmodes` (legacy compute),
  `panel_eigenmodes_compute` (its dialog), and the `in_tess_eigenmodes` ecosystem are kept; only the
  display *lever* is retired.

## 4. Architecture and file changes

Two independent units of work.

### A. New viewer — rewrite `toolbox/gui/view_eigenmodes.m`

- **New contract:** `hFig = view_eigenmodes(EigenFile)` where `EigenFile` is an `eigen_*.mat` node
  path (previously a surface file). Pure subfunctions remain dispatchable via `eval(macro_method)`
  for headless tests (e.g. `view_eigenmodes('ReconstructModeField', ...)`).
- Loads `EigenMat`, resolves `ParentSurface`, dispatches on `EigenMat.Variant`. Dirac implemented;
  others → `bst_error`.
- Self-contained figure: `view_surface` (translucent gray cortex) + own `DrawArrows`/`quiver3` +
  own `KeyPressFcn`. No transient DB result, no `panel_eigenmodes`, no cleanup hook.

### B. Lever retirement

- **Delete** `toolbox/gui/panel_eigenmodes.m` (the display lever).
- `gui_initialize.m:53` — remove `gui_show('panel_eigenmodes', 'BrainstormTab', 'tools');`.
- `bst_figures.m` — delete `FireModesChanged()` (~997–1040) in full (its only co-writer of
  `GlobalData.UserModes` is the deleted lever); remove the two
  `isTabVisible('EigenModes') → panel_eigenmodes('UpdatePanel')` blocks (~890, ~1163); drop any
  now-orphaned `UserModes` initialization.
- `panel_surface.m` — remove the lever hooks at ~1816–1817 (`ApplyToColumn` wrap) and ~1860–1862
  (`IsActive` DataMinMax override); `TessInfo(iTess).Data` reverts to the plain `GetResultsValues`
  result.

### Wiring (`tree_callbacks.m` + `process_eigenmodes.m`)

- `EigenView_Callback` (~4081): replace the field-dump stub body with `view_eigenmodes(filenameFull)`
  (keep the `file_exist` guard + `bst_error`).
- Remove the cortex-node "View eigenmodes" menu item (~1185).
- Keep cortex-node "Compute eigenmodes (legacy)" (~1183) and the "Compute eigenmodes"
  Dirac/LBO/Connection submenu (~1204).
- `process_eigenmodes.m:289` — remove the `view_eigenmodes(SurfaceFile)` auto-view (and its
  "Visual confirmation" comment). Legacy compute still runs + reloads the DB; it no longer auto-opens.

### Explicitly kept (untouched)

`panel_eigenmodes_compute.m`, `process_eigenmodes` (legacy compute + its menu item), and the entire
`in_tess_eigenmodes` ecosystem (~30 files).

## 5. Viewer contract and Dirac reconstruction

### 5.1 Input resolution

```
hFig = view_eigenmodes(EigenFile)
  EigenMat = load(file_fullpath(EigenFile))      % db_template('eigenmat')
  Surface  = EigenMat.ParentSurface
  TessMat  = in_tess_bst(Surface)                % .Vertices = quiver anchors
  K        = size(EigenMat.Phi{1}, 2)            % modes actually stored (not EigenMat.K)
  switch lower(EigenMat.Variant)
    case 'dirac' -> vector viewer (implemented)
    case {'laplace-beltrami','connection laplacian'}
                 -> bst_error('Vector viewer currently supports Dirac eigenmodes only.')
    otherwise    -> bst_error('Unknown eigen variant.')
  end
```

Guards (each → clear `bst_error`): missing file; empty/short `Phi`; unset `Variant`; `Phi{hh}` row
count `≠ 4·numel(GlobalVertices{hh})`.

### 5.2 Per-mode vector field (pure subfunction, unit-testable)

```matlab
function V3 = ReconstructModeField(EigenMat, k, nVert)   % -> [nVert x 3]; zeros off-support
    V3 = zeros(nVert, 3);
    for hh = 1:2
        vH  = EigenMat.GlobalVertices{hh}(:);
        col = EigenMat.Phi{hh}(:, k);        % [4*nVh x 1], real for Dirac
        V3(vH,1) = col(2:4:end);             % quaternion i  -> x   (drop w = col(1:4:end))
        V3(vH,2) = col(3:4:end);             % quaternion j  -> y
        V3(vH,3) = col(4:4:end);             % quaternion k  -> z
    end
end
```

Components are ambient SCS Cartesian → drawn directly as `quiver3` `U,V,W`. Per-vertex magnitude
`‖V3‖` shows the mode's support and relative strength. Both hemispheres' rank-`k` mode are shown
together (both lit).

### 5.3 Sign / gauge note (documented, not handled)

Each eigenvector carries an arbitrary global sign, and degenerate quaternionic multiplets carry an
arbitrary intra-multiplet basis rotation, so absolute arrow direction within a multiplet is
gauge-dependent. The viewer shows the field as stored — a known property of eigenbasis display, not
a bug.

### 5.4 Arrow magnitude

Drawn at **true magnitude** scaled by a global `quiverSize` (like `view_leadfield_vectors`), so mode
support and relative strength are visible. An `n` key toggles unit-normalization (ε-guard); default
is true magnitude.

## 6. Figure setup, interaction, legend

### 6.1 Setup (mirrors `view_leadfield_vectors`)

- `hFig = view_surface(Surface, 0.5, [0.5 0.5 0.5], 'NewFigure')` — translucent gray cortex.
- `hAxes = findobj(hFig,'-depth',1,'Tag','Axes3D')`; `figure_3d('SetStandardView', hFig, 'left')`.
- Stash and hijack `KeyPressFcn` (original called in `otherwise`).
- Bottom-left legend `uicontrol` (green-on-black).
- `set(hFig,'Name',['Eigenmodes: ' Variant ' | ' Surface])`.
- State as nested-function closure vars: `iMode=1`, `quiverSize=1`, `quiverWidth=1`,
  `thresholdAmplitude=1`, `thresholdBalance=0`, `useNormalize=false`.

### 6.2 `DrawArrows()` (nested)

1. `delete(findobj(hAxes,'-depth',1,'Tag','eigArrows'))`.
2. `V3 = ReconstructModeField(EigenMat, iMode, nVert)`; optional unit-normalize (ε-guard).
3. Amplitude threshold via the template's CDF gate (`thresholdBalance` chooses `<=` / `>`).
4. `quiver3(Vertices(:,1), Vertices(:,2), Vertices(:,3), V3(:,1), V3(:,2), V3(:,3), quiverSize,
   'Parent',hAxes, 'LineWidth',quiverWidth, 'Color',[.3 1 .3], 'Tag','eigArrows')`.
5. Update legend.

### 6.3 Keyboard map (single-mode stepping)

| Key | Action |
|---|---|
| `Left` / `Right` | mode `−1` / `+1` (wraps `1..K`) |
| `PgUp` / `PgDn` | mode `+10` / `−10` |
| `Shift+Up/Down` | quiver length `× / ÷ 1.2` |
| `Ctrl+Up/Down` | quiver width `× / ÷ 1.2` |
| `Alt+Up/Down` | amplitude threshold `± 0.01` |
| `Alt+Enter` | toggle threshold direction (`<=` / `>`) |
| `n` | toggle unit-normalize arrows |
| `0`–`9` | standard Brainstorm views (stashed callback) |
| `h` | help dialog |
| otherwise | original `KeyPressFcn` |

### 6.4 Legend

`Mode k / K   |   λL = <Lambda{1}(k)>, λR = <Lambda{2}(k)>   |   τ = <Provenance.Tau>
[arrows: <count> | H for help]`.

### 6.5 Cleanup

None required — nothing is registered in the DB; closing the figure closes it.

## 7. Edge cases and error handling

- **Non-Dirac eigen node** → `bst_error` from the variant dispatch (LBO/Connection deferred).
- **`k` out of range** (`< 1` or `> K`) → wrap on keyboard; `ReconstructModeField` errors clearly if
  called directly with a bad `k`.
- **Row-count mismatch** (`Phi{hh}` not `4·nVh`) → `bst_error` before display.
- **Near-zero vertices** under normalize → ε-guard (no arrow / zero length), no NaN.
- **Missing `Provenance.Tau`** → legend shows `τ = ?` (non-fatal).

## 8. Testing

### 8.1 Headless unit tests — `dev/tests/test_eigenmode_vector_field.m`

Via `eval(macro_method)` dispatch (no figure):
- Synthetic `EigenMat` (2 hemispheres, small `nVh`, known `Phi`): `ReconstructModeField` maps
  `(i,j,k)→(x,y,z)`, drops `w`, scatters to the right `GlobalVertices` rows, zeros off-support.
- Mode index `k` selects the correct column; out-of-range `k` errors.
- Threshold/decimation helper (if factored): threshold reduces count; show-all keeps all.

### 8.2 Live-figure integration — `dev/tests/test_eigenmode_vector_viewer.m`

Requires the running session + a Dirac `eigen_` node on the canonical cortex
(`tess_eigen(Surface,'Dirac','K',...)` if absent):
- `view_eigenmodes(EigenFile)` returns a valid `hFig`; `Axes3D` survives; an `eigArrows` quiver with
  `numel(UData) > 0` exists.
- Stepping mode (`Right`/`PgDn`) changes the drawn field (`UData` differs).
- `quiverSize` / threshold keys change arrow length / count.
- An **LBO** eigen node → the "Dirac only" `bst_error`.
- Close cleanly (no DB node left behind).

### 8.3 Retirement regression (shell or test step)

- `grep -rn "panel_eigenmodes(" toolbox/` → empty (only `panel_eigenmodes_compute` references remain).
- `grep -rn "FireModesChanged\|UserModes\|ModesChangedCallback\|EigenView" toolbox/` → empty.

## 9. Future work (out of scope here)

- LBO (real scalar) and Connection Laplacian (complex scalar) variant branches → colormap display.
- Optional arrow coloring by direction or hemisphere.
- Differential analysis / filtering of the field in the eigen-coefficient domain.
