# view_manifold_registration — Registration-Sphere View — Design

- **Date:** 2026-06-12
- **Status:** Approved (inline execution)
- **Author:** Diellor Basha (with Claude)
- **Related:** `view_manifold.m`, `view_surface_sphere.m`, `tess_hemisplit.m`, `tree_callbacks.m`

## 1. Goal

Display the FreeSurfer registration sphere for a manifold node's parent surface, with the manifold's
singularity poles overlaid as lollipop pins — the same pins `view_manifold` draws on the cortex, now
shown on the sphere. This is the first step toward validating analytic functions (heat, wave) against
their analytic counterparts on the sphere.

## 2. Approach (settled)

- **Dedicated file** `toolbox/gui/view_manifold_registration.m` (keeps `view_manifold` focused).
- **Reuse `view_surface_sphere`** to build the sphere figure, then overlay the poles.
- **Poles only** for this version; frames and analytic colormaps extend the same overlay seam later.

## 3. Components

### 3.1 `toolbox/gui/view_manifold_registration.m` (new)

`hFig = view_manifold_registration(ManifoldFile)`:
1. `file_exist` guard; `load` the node; require `ParentSurface`, `Embedded` (1×2), `Gauge` (1×2) →
   else `bst_error` + return `[]`.
2. `nVert = size(in_tess_bst(ParentSurface).Vertices,1)`; pole ids by **reusing**
   `G = view_manifold('DeriveVertexFrame', Embedded, Gauge, nVert)` → `G.Sing` (global vertex ids).
3. `hFig = view_surface_sphere(ParentSurface, 'orig')`. If empty (no `Reg.Sphere`), return `[]`
   (`view_surface_sphere` already raised its own `bst_error`).
4. Read the displayed sphere positions: `hPatch` = the patch in `Axes3D` with `nVert` vertices;
   `sphV = get(hPatch,'Vertices')`.
5. `[ir, il] = tess_hemisplit(in_tess_bst(ParentSurface))`; `[Base,Tip] = RegSingGlyphs(sphV, G.Sing, ir, il)`.
6. `hold(hAxes,'on')`; draw blue stems (`line`, Tag `manifoldRegSingStem`) + markers (`plot3 'o'`,
   Tag `manifoldRegSing`); `P` toggles them, `H` shows help, other keys → the stashed `view_surface_sphere`
   key callback. Figure name `Manifold registration: <Surface>`.

**Macro dispatch:** `view_manifold_registration('RegSingGlyphs', sphV, poleIdx, ir, il)` routes to the
pure subfunction (headless testable); otherwise the GUI entry.

### 3.2 `RegSingGlyphs(sphV, poleIdx, ir, il)` → `[Base, Tip]` (pure)

```
poleIdx = poleIdx(:);  Base = sphV(poleIdx,:);
cL = mean(sphV(il,:),1);  cR = mean(sphV(ir,:),1);
rL = mean(vecnorm(sphV(il,:)-cL,2,2));  rR = mean(vecnorm(sphV(ir,:)-cR,2,2));
for each pole p (with center c = cL if p in il else cR, radius rr = rL or rR):
    rad = sphV(p,:) - c;  rad = rad / max(norm(rad), eps);
    Tip = Base + 0.15*rr * rad
```
Lifting from each pole's **own hemisphere** sphere center (not the global midpoint between the two
offset spheres) makes every pin point radially out of its sphere.

### 3.3 `toolbox/tree/tree_callbacks.m`

Manifold popup (`case 'manifold'`): add, after "View manifold":
```matlab
gui_component('MenuItem', jPopup, [], 'View registration sphere', IconLoader.ICON_SURFACE, [], @(h,ev)bst_call(@view_manifold_registration, filenameFull));
```

## 4. Data flow

```
view_manifold_registration(ManifoldFile)
  → load node; Surface=ParentSurface; nVert
  → G = view_manifold('DeriveVertexFrame', Embedded, Gauge, nVert)   % poles = G.Sing
  → hFig = view_surface_sphere(Surface, 'orig')
  → sphV = patch Vertices; [ir,il] = tess_hemisplit
  → [Base,Tip] = RegSingGlyphs(sphV, G.Sing, ir, il); draw lollipops; P/H keys
```

## 5. Edge cases

- Missing file / invalid node → `bst_error` + `[]`.
- No `Reg.Sphere` on the parent surface → `view_surface_sphere` errors and returns `[]`; detect the
  empty handle and return without overlaying.
- Empty `G.Sing` → sphere shown, no lollipops (non-fatal).
- Multiple patches in the axes → pick the one whose `Vertices` has `nVert` rows.

## 6. Testing

- **Headless** `dev/tests/test_manifold_reg_glyphs.m` — pure `RegSingGlyphs` via macro dispatch:
  two offset unit spheres (L center +y, R center −y), `poleIdx` at L-top and R-bottom; assert
  `Base = sphV(poleIdx)`, each `Tip` lifted outward from its own hemisphere center
  (`dot(Tip−Base, sphV(p)−center) > 0`), the L/R pins point in opposite global z, and the count
  matches `numel(poleIdx)`.
- **Live** `dev/tests/test_view_manifold_registration.m` — open on the registered manifold node;
  assert a figure returns, `Axes3D` + a sphere patch exist, a `manifoldRegSing` marker object with
  `numel(G.Sing)` points is drawn, `P` toggles it off, close cleanly.

## 7. Future work (out of scope)

- Manifold tangent frames mapped onto the sphere.
- Analytic-function colormaps (heat/wave) on the sphere for validation against closed-form solutions.
- Nesting the cortex and sphere views behind one entry point.
