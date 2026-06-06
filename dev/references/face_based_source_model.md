# Face-Based Current Flux Source Model

**Status:** Reference note — not yet implemented in production  
**Implements:** Cortical sources as primal 2-forms (current flux through triangular patches)  
**Authors:** Diellor Basha, 2026

---

## 1. Physical premise

### Vertex-based model (current state)

At each vertex v the model places a **point dipole** with no spatial extent:

```
p_v = s_v · n̂_v   [A·m]
B(r) = (μ₀/4π) (p_v × (r − r_v)) / |r − r_v|³
```

`r_v` is the vertex position and `n̂_v` is an *averaged* normal (area-weighted mean of
surrounding face normals).  Both are approximations: the source has zero area, and the
normal direction is not the exact normal of any geometric element.

### Face-based model (this note)

At each triangular face f the model assigns a **current flux density** s_n(f) [A/m²].
The total dipole moment of the patch is the area integral:

```
m_f = ∫_{face f} s_n · n̂_f dA   [A·m]
     = s_n_f · A_f · n̂_f
```

The leadfield entry is the area integral of the Green's function over the patch:

```
L(sensor, f) = ∫_{face f} G(r_sensor, r') · n̂_f dA'
```

This is a **distributed source** — `s_n_f` is a 2-form (current flux per unit area, or
equivalently the coefficient of the area form n̂_f dA).  The area weighting is built into
the source definition, not applied post-hoc.

---

## 2. Physical scale and the cortical column argument

### ico5 face geometry

The 20484-vertex bilateral surface is ico5 resolution — 10242 vertices per hemisphere.
For a closed genus-0 hemisphere (topological sphere), Euler's formula gives:

```
F ≈ 2V − 4 ≈ 20480 faces per hemisphere
```

With one hemisphere ≈ 100,000 mm² total surface area:

```
A_face = 100,000 mm² / 20,480 ≈ 5 mm²
edge   = √(4 × 5 / √3) ≈ 3.4 mm    (equilateral triangle)
```

Each ico5 triangle has area **~5 mm²** and edge length **~3–4 mm**.

### Cortical column sizes

| Unit | Diameter | Cross-sectional area |
|---|---|---|
| Minicolumn (anatomical) | 50–80 μm | 0.002–0.005 mm² |
| Macrocolumn (functional) | 0.5–1 mm | 0.2–0.8 mm² |
| **ico5 face** | **~2.5 mm** (eff. diam.) | **~5 mm²** |

One ico5 face contains approximately **1000–25000 minicolumns** and **6–25 macrocolumns**.

### What this means physically

**The face source is a mesoscale average, not a single column.**
The recovered scalar s_n(f,t) is the net normal current flux through a ~5 mm² patch —
the spatially integrated dipole moment of all cortical columns within that patch, weighted
by the area integral of the Green's function.  This is a genuine mesoscale quantity.

**The 2-form interpretation is physically honest about this.**
Current flux through a ~5 mm² patch is exactly the right description of what an ico5
source element represents.  It makes no pretense of resolving individual columns.
A point dipole approximation at a vertex implies spatial localization that is not
justified at this scale — a vertex "represents" all the same tissue as the surrounding
faces but without acknowledging the area integration that went into it.

### MEG spatial resolution vs ico5 sampling

MEG spatial resolution is empirically **~5–10 mm** — roughly one to three ico5 face
edge lengths.  The ico5 triangulation is therefore sampled at approximately the **Nyquist
rate** for MEG spatial resolution:

- Full FreeSurfer surface (~150,000 vertices per hemisphere, edge ~0.6 mm) is **massively
  oversampled** relative to MEG observability.  The vast majority of those fine-scale
  degrees of freedom are in the null space of the leadfield — MEG simply cannot see them.
- ico5 matches the MEG observable scale.  Adding more vertices does not add observable
  information; it only adds degrees of freedom that are empirically invisible.

### Implication for eigenmode truncation

The observable eigenmodes correspond to spatial wavelengths Λ_k = 2π/√λ_k >> 5 mm.
The first ~200 observable LBO eigenmodes span spatial scales from the whole cortex down
to roughly 5–10 mm.  The face-based source model with ~40,000 faces is thus dramatically
overdetermined — the ~200 observable modes each cover many faces.

**This strengthens the face-based approach:** because each face is a proper
area-integrated quantity, the eigenmode expansion on face-based sources is a
*geometrically honest compression* — smooth basis functions that match the actual
resolution of the measurement.  The point-based approach creates a false impression that
vertex-level resolution is meaningful, when MEG is seeing only the smooth spatial integrals
that the eigenmodes represent.

The practical consequence: eigenmode truncation at K ≈ 200 modes is not a lossy
approximation at ico5 resolution.  It is a lossless description of the MEG-observable
source space.

---

## 3. DEC interpretation

In Discrete Exterior Calculus the correct type for a current sheet on a surface is a
**primal 2-form**:

| Object | DEC type | Space | Physical meaning |
|---|---|---|---|
| s_n_f | Primal 2-form ∈ Ω²(M) | Per face scalar | Normal current flux (A·m) |
| t̂₁_f, t̂₂_f components | Primal 1-form ∈ Ω¹(M) | Per edge scalar | Tangential current along edge |
| Scalar potential | 0-form ∈ Ω⁰(M) | Per vertex scalar | (not directly used here) |

The face-based model treats s as a 2-form, which is the *natural* DEC type for
sources on a 2-manifold — same mathematical type as the volume form.

---

## 4. BEM consistency

The standard MEG forward model is already face-based for the **head boundaries** (scalp,
skull, CSF) via BEM.  BEM solves surface integrals over triangular patches; it is
fundamentally a face method.  The cortical surface is itself a boundary (the inner surface
of the CSF), yet is treated with point sources at vertices — a geometric inconsistency.

The face-based source model removes this inconsistency.  The primary current contribution
becomes:

```
B_primary(r) = (μ₀/4π) Σ_f s_n_f ∫_{face f} n̂_f × ∇G(r,r') dA'
```

This has the same mathematical form as the BEM volume-current integrals.  The forward
model is now geometrically unified: one type of area-integral computation for all
boundary surfaces including the cortex.

---

## 5. What changes numerically

**For distant sensors** (sensor-to-source distance >> face size):

```
∫_{face f} G(r,r') dA' ≈ G(r, x_f) · A_f   (centroid rule, O(h²/d²) error)
```

where x_f = (r_i + r_j + r_k)/3.  Numerically negligible for h~3mm, d~100mm
(error ~0.1%).

**For high-curvature regions** (sulcal fundi, gyral crowns): the face normal n̂_f is
exact for each patch; vertex normals near sharp bends are smoothed by averaging.
Face-based forward entries near sulcal fundi are more accurately oriented.

**Conditioning:** The leadfield null-space structure is more precisely characterised.
Near-cancellation of opposite sulcal wall columns is exact when face normals are exactly
antiparallel; vertex normal averaging introduces a small artificial decorrelation that
biases the singular value spectrum.

---

## 6. The face Laplacian — null space issue on surfaces with boundary

**Important numerical finding (test_phase_recovery_v4):**  
The discrete face Laplacian `d₁ · ★₁⁻¹ · d₁ᵀ` is numerically degenerate  
(RCOND ≈ 5e-17) when computed on a hemisphere (manifold with boundary).

Cause: 11.8% of edges have negative cotan weights (obtuse triangles). The DEC
`★₁` is not positive definite; `★₁⁻¹` has negative entries; the product
`d₁ · diag(negative) · d₁ᵀ` is indefinite and the boundary creates a large null space.

**Practical fix:** Use vertex LBO eigenmodes (from `d₀ᵀ · ★₁ · d₀`, which IS PSD
for a triangulated surface) interpolated to faces:

```matlab
Phi_f(face, k) = (Phi_v(i,k) + Phi_v(j,k) + Phi_v(l,k)) / 3   % corner average
```

Then orthonormalize under the face mass matrix `M_f = ★₂` via Cholesky:

```matlab
G = Phi_f' * M_f * Phi_f;          % Gram matrix  [K x K]
[R,~] = chol(G);
Phi_f_orth = Phi_f / R;             % M_f-orthonormal face basis
```

This gives machine-precision M_f-orthogonality (off-diag < 1e-15).

**Why this works physically:** Vertex LBO eigenmodes are smooth; averaging to face
centroids preserves this smoothness. The resulting face-indexed modes are smooth functions
on the dual mesh, appropriate for representing slowly varying current flux patterns.

**When a true face Laplacian is needed:** The correct formulation requires a cotan-weight
clamping (set negative weights to zero or to a small positive value) before computing
`d₁ · clamp(★₁)⁻¹ · d₁ᵀ`. This changes the spectrum but restores positive definiteness.

## 7. The face Laplacian and its eigenmodes

### Operators (all from `nxr_compute('assembleDECOperators', h)`)

```
d₀ : [nE × nV]    exterior derivative on 0-forms
d₁ : [nF × nE]    exterior derivative on 1-forms
★₁ : diag [nE]    Hodge star on 1-forms (cotan weights)
★₂ : diag [nF]    Hodge star on 2-forms (face areas)
```

### Face scalar Laplacian (dual Laplacian on primal 2-forms)

```
Δ_face = d₁ · diag(1 ./ diag(★₁)) · d₁'   [nF × nF]
M_face = ★₂                                  [nF × nF, diagonal face areas]
```

Eigenvalue problem: `Δ_face · ψ = λ · M_face · ψ`

### How the spectrum differs from the vertex LBO

For a fine triangulation of a 2-manifold (nF ≈ 2·nV, nE ≈ 3·nV):

- Vertex LBO Δ_vert acts on Ω⁰ (dim nV), face Δ_face acts on Ω² (dim nF ≈ 2·nV)
- Both have one zero eigenvalue (DC / constant)
- The non-zero eigenvalues are NOT the same set (Hodge duality holds in the continuum
  but NOT in the discrete setting because nV ≠ nF)
- The face Laplacian produces ~2× more eigenmodes; eigenvalues are numerically different
- Face eigenmodes ψ_k(f) are functions on face barycenters, not vertices — they live on
  the dual mesh as 0-forms

For our 20484-vertex mesh (10242 per hemisphere, 20480 faces per hemisphere):
  vertex LBO → 10241 non-DC modes per hemisphere
  face  LBO  → 20479 non-DC modes per hemisphere

### Spatial gradient of a face-indexed phase field

Given Φ(f,t) ∈ Ω²(M) (phase at each face), the exterior derivative (codifferential) is:

```
δΦ = ★₁⁻¹ · d₁' · ★₂ · Φ   ∈ Ω¹(M)   [nE × nT]
```

This gives a **per-edge 1-form**: the phase difference between the two adjacent faces
weighted by their shared dual-edge length.  The direction of δΦ on edge e tells you
which adjacent face is ahead in the wave.

To obtain a tangent vector at each face, restrict δΦ to the three edges of that face
and solve a 2×3 least-squares problem in the face tangent plane — or use the trivial
connection face frame directly.

---

## 8. Implementation recipe

### Step A: Face-based constrained leadfield (approximation)

Requires: unconstrained Gain [nCh × 3nV], face normals n̂_f [nF × 3], face areas A_f [nF × 1]

```matlab
Gain_u = double(hm.Gain);          % [nCh × 3nV]  unconstrained
FaceNormals = compute_face_normals(Vtx, Faces);   % [nF × 3]  exact
A_f = face_areas(Vtx, Faces);                     % [nF × 1]

L_face = zeros(nCh, nF);
for f = 1:nF
    i = Faces(f,1); j = Faces(f,2); k = Faces(f,3);
    G_avg = (Gain_u(:,3*i-2:3*i) + Gain_u(:,3*j-2:3*j) + Gain_u(:,3*k-2:3*k)) / 3;
    L_face(:,f) = G_avg * FaceNormals(f,:)' * A_f(f);
end
% Error: O(h²/d²) ≈ 0.1% for h~3mm, d~100mm
```

### Step A (exact): Face-based constrained leadfield from centroid evaluation

For the os_meg overlapping spheres analytic model, re-evaluate at face centroids:

```matlab
x_f = (Vtx(Faces(:,1),:) + Vtx(Faces(:,2),:) + Vtx(Faces(:,3),:)) / 3;
% Call the os_meg analytic formula at x_f with orientation n̂_f:
L_face = bst_meg_leadfield(x_f, FaceNormals, sensor_positions, head_model_params);
% This is the exact face-centroid leadfield; no vertex normal approximation enters.
```

**When to use exact vs approximation:**
- Testing and prototyping: approximation is fine (0.1% error)
- Production face-based head model: use centroid evaluation

### Step B: Face Laplacian eigenmodes

```matlab
h = nxr_compute('create', Vtx, Faces);
dec = nxr_compute('assembleDECOperators', h);
nxr_compute('destroy', h);

L_f = dec.d1 * diag(1./diag(dec.hodge1)) * dec.d1';   % [nF × nF]
M_f = dec.hodge2;                                        % [nF × nF] diagonal

% Solve for K face eigenmodes (left hemisphere faces only)
[Psi, Lam] = eigs(L_f(lhF, lhF), M_f(lhF, lhF), K, 'smallestabs');
lam_f = real(diag(Lam));
% Psi: [nLHF × K]  face-indexed eigenmodes
```

### Step C: Eigenmode face leadfield (face-based Eigen-MNE)

```matlab
L_tilde_face = L_face(:, lhF) * Psi_face;   % [nCh × K]
% Standard regularized inverse on L_tilde_face → θ(t) [K × nT]
% Reconstruction: s_eig(f,t) = Psi_face * θ(t)  [nLHF × nT]
```

### Step D: Spatial gradient of face phase field

```matlab
% Phi_f: [nF × nT]  phase field on faces (primal 2-form)
hodge1_inv = spdiags(1./diag(dec.hodge1), 0, size(dec.hodge1,1), size(dec.hodge1,1));
dPhi = hodge1_inv * dec.d1' * dec.hodge2 * Phi_f;   % [nE × nT]  codifferential
% dPhi(e,t) = phase difference between the two faces sharing edge e,
%             scaled by dual-edge length (face-centroid to face-centroid distance)
```

---

## 9. Why this matters for wave detection

The sign of s_n(f,t) is physically unambiguous:
- Positive = net outward current flux through face f at time t
- Negative = net inward current flux
- Sign determined entirely by mesh winding convention, consistent everywhere

The phase Φ(f,t) = arg(s_n_analytic(f,t)) is a well-defined angular quantity:
- No vertex-normal approximation in the sign
- No leakage between normal and tangential components
- Spatial gradient δΦ is exactly the per-edge phase difference in the DEC sense

This is the correct physical and mathematical framework for detecting traveling waves as
propagating wavefronts of cortical current flux.

---

## 10. Open questions and future work

1. **Exact centroid leadfield for os_meg:** Implement the analytic formula at face
   centroids rather than approximating from vertex columns.  The os_meg formula is a
   closed-form expression and can be evaluated at any point.

2. **BEM-consistent cortical source:** Formally integrate the face-based cortical source
   into the BEM system as another boundary element, treating the pial surface as a
   current-carrying BEM boundary.  This would give a fully consistent forward model.

3. **DEC-based face LBO eigenmodes in production Brainstorm:** Build a face-based
   eigenmode head model analogous to bst_eigenmode_leadfield but using face LBO modes
   instead of vertex LBO modes.  The eigenvalue problem uses Δ_face and M_face from DEC.

4. **Phase gradient operator in Brainstorm:** Implement the codifferential δΦ as a
   standard Brainstorm operator for wave detection on face-indexed source maps.

5. **Comparison with vertex LBO eigenmode spectrum:** Quantify how the face LBO
   eigenvalues and eigenmodes differ from vertex LBO, and whether the face eigenmodes
   provide better sign consistency (fewer sign errors from the inverse) than vertex
   eigenmodes in the phase recovery tests.
