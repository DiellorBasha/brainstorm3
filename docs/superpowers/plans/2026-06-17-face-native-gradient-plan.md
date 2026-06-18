# Face-native dual-mesh gradient (`gradFace` / `lapFace`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (or subagent-driven-development) to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Add two first-class nxr-compute operators — a barycentric dual-mesh face gradient `gradFace` `[3F×F]` and its Galerkin face Laplacian `lapFace` `[F×F]` (`K̃ = gradFaceᵀ W_F gradFace`) — then rewrite `bst_dirac_helmholtz_face` to a fully face-native Helmholtz with an exact (HarmFrac→0) round-trip.

**Architecture:** `gradFace` is the Green–Gauss / DEC barycentric gradient, the exact dual of the vertex FEM gradient: `grad ψ̃|_f = (1/2A_f) Σ_k (ψ̃_{g_k}−ψ̃_f)(n_f × ℓ_k)` over face `f`'s three edge-neighbors `g_k` (closed mesh → always 3). It annihilates constants (`Σ_k ℓ_k = 0`) and is tangent by construction. `lapFace` is built from the same triplets so symmetry/PSD are free; its kernel is constants-only (genus-0). Everything else (skew-gradient `n_f×G̃`, div/curl source maps) is trivial MATLAB on top.

**Tech Stack:** C++17, geometry-central, Eigen, MEX + WASM (nxr-compute); MATLAB / Brainstorm (consumer). Spec: `docs/superpowers/specs/2026-06-17-face-native-gradient-design.md`.

**Preconditions:** nxr-compute on a fresh feature branch off `main` (the `diracFaceIntrinsicD` work is on `main`). Brainstorm `development` branch. TutorialAuditory protocol + toolbox/dev paths loaded in the MATLAB MCP. Build host has cmake/ninja + the geometry-central submodule.

---

## PHASE C — nxr-compute operators (C++)

### Task C0: Feature branch

- [ ] **Step 1: Branch.** In `~/workspace/research/code/nxr-compute`:
```bash
git checkout main && git pull && git checkout -b feat/face-native-gradient
```

---

### Task C1: `gradFace` builder + C++ tests

**Files:**
- Create: `nxr-compute/src/face_gradient.cpp`
- Modify: `nxr-compute/include/nxr/operators.h` (or wherever `ops::dirac` builders are declared — match `matrixFaceIntrinsic`'s declaration site; grep `matrixFaceIntrinsic` in headers)
- Modify: `nxr-compute/CMakeLists.txt` (add `src/face_gradient.cpp` to the library sources; add a `test_face_gradient` executable + `add_test`)
- Create: `nxr-compute/test/test_face_gradient.cpp`

- [ ] **Step 1: Write the failing C++ test** `test/test_face_gradient.cpp` (model the header/EXPECT/icosphere boilerplate on `test/test_dirac_face_operator.cpp`). Tests for the gradient alone:
```cpp
#include "nxr/compute.h"
#include "nxr/facets.h"
#include <Eigen/Eigenvalues>
#include <cmath>
#include <iostream>
using namespace nxr::manifold;
static int g_failures = 0;
#define EXPECT(cond,msg) do{ if(cond){std::cout<<"  [PASS] "<<msg<<"\n";} \
  else{std::cout<<"  [FAIL] "<<msg<<"\n"; ++g_failures;} }while(0)
// (copy icosphere(V,F) from test_dirac_face_operator.cpp — 12 verts, 20 faces, closed genus-0)

static void testGradFace() {
    std::cout << "\n=== gradFace: barycentric dual gradient ===\n";
    std::vector<double> V; std::vector<int32_t> F; icosphere(V,F);
    Manifold m(V.data(), 12, F.data(), 20);
    const int Fn = m.nF();
    Eigen::SparseMatrix<double> G = ops::facegrad::gradient(m);
    EXPECT(G.rows()==3*Fn && G.cols()==Fn, "gradFace is [3F x F] = [60, 20]");
    // constant precision: grad of a constant per-face field is zero
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(Fn);
    EXPECT((G*ones).norm() < 1e-10, "gradFace annihilates constants");
    // tangency: each per-face 3-vector is orthogonal to that face's normal
    m.geometry().requireFaceNormals();
    Eigen::VectorXd psi = Eigen::VectorXd::Random(Fn);
    Eigen::VectorXd g = G*psi;
    double maxdot = 0; int fi = 0;
    for (auto f : m.mesh().faces()) {
        auto n = m.geometry().faceNormals[f];
        Eigen::Vector3d gv(g[3*fi], g[3*fi+1], g[3*fi+2]);
        Eigen::Vector3d nn(n.x, n.y, n.z);
        maxdot = std::max(maxdot, std::abs(gv.dot(nn))); ++fi;
    }
    EXPECT(maxdot < 1e-9, "gradFace output is tangent (perp to face normal)");
}
int main(){ testGradFace(); std::cout<<(g_failures? "\nFAILURES\n":"\nALL PASS\n"); return g_failures?1:0; }
```

- [ ] **Step 2: Declare + implement `ops::facegrad::gradient`.** Add the declaration next to the other op builders, then create `src/face_gradient.cpp`:
```cpp
#include "nxr/compute.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include <Eigen/Sparse>
#include <array>
namespace nxr::manifold::ops::facegrad {
using namespace geometrycentral; using namespace geometrycentral::surface;

// Barycentric dual-mesh gradient of a per-face scalar: per-face ambient 3-vector [3F x F].
//   grad psi|_f = (1/2A_f) * sum_k (psi_{g_k} - psi_f) * (n_f x l_k)
// where edge k of face f (local verts (k+1)%3 -> (k+2)%3) has edge vector l_k and the
// neighbor face g_k across it. n_f x l_k is the tangent in-plane edge normal (length |l_k|);
// sum_k l_k = 0 so constants are annihilated and the self (psi_f) term vanishes. Closed-mesh
// v1: every face has exactly 3 edge-neighbors (fail loud on an open boundary).
Eigen::SparseMatrix<double> gradient(Manifold& m) {
    auto& mesh = m.mesh(); auto& geom = m.geometry();
    if (mesh.hasBoundary())
        throw core::Error(core::ErrorCode::InvalidInput,
            "facegrad::gradient: open boundary unsupported (closed-mesh v1)",
            "Pass a closed mesh (e.g. a FreeSurfer hemisphere).");
    geom.requireFaceNormals(); geom.requireFaceAreas();
    const int Fn = m.nF();
    auto vec3 = [](const Vector3& u){ return Eigen::Vector3d(u.x,u.y,u.z); };
    std::vector<Eigen::Triplet<double>> T; T.reserve((size_t)Fn*9);
    for (Face f : mesh.faces()) {
        const int fi = (int)f.getIndex();
        const double A = geom.faceAreas[f];
        if (A <= 0.0) throw core::Error(core::ErrorCode::InvalidInput,
            "facegrad::gradient: degenerate (zero-area) face",
            "Face " + std::to_string(fi) + "; fix mesh quality first.");
        Eigen::Vector3d nf = vec3(geom.faceNormals[f]);
        std::array<Eigen::Vector3d,3> P{}; std::array<int,3> vid{}; int c=0;
        for (Vertex v : f.adjacentVertices()){ P[c]=vec3(geom.inputVertexPositions[v]); vid[c]=(int)v.getIndex(); ++c; }
        // halfedges of f in order; he between local (k+1)%3 and (k+2)%3 lies opposite local k.
        // its twin's face is the neighbor across that edge.
        int k = 0;
        for (Halfedge he : f.adjacentHalfedges()) {
            // he goes P[k] -> P[(k+1)%3]; this is the boundary edge between those two verts.
            Eigen::Vector3d l = P[(k+1)%3] - P[k];          // boundary edge vector (CCW)
            Eigen::Vector3d w = nf.cross(l) / (2.0*A);      // tangent in-plane edge normal / 2A
            Face g = he.twin().face();                       // neighbor across this edge
            const int gi = (int)g.getIndex();
            // (psi_g - psi_f) * w  ->  +w on column g, -w on column f
            for (int d=0; d<3; ++d) {
                if (w[d]!=0.0) { T.emplace_back(3*fi+d, gi, w[d]); T.emplace_back(3*fi+d, fi, -w[d]); }
            }
            ++k;
        }
    }
    Eigen::SparseMatrix<double> G(3*Fn, Fn); G.setFromTriplets(T.begin(), T.end()); G.makeCompressed();
    return G;
}
} // namespace
```
*(Orientation note: if the round-trip test in Task D yields `corr(ψ̃,ψ0) ≈ −1` instead of `+1`, flip the sign of `w` (use `l.cross(nf)`). The constant/tangency tests are sign-agnostic; the sign is pinned by the Task-D round-trip.)*

- [ ] **Step 3: Build + run the test.**
```bash
cmake --build build --target test_face_gradient && ./build/test_face_gradient
```
Expected: all PASS (shape, constants, tangency).

---

### Task C2: `lapFace` builder + C++ tests

**Files:** Modify `src/face_gradient.cpp`, `test/test_face_gradient.cpp`.

- [ ] **Step 1: Add failing tests** to `testGradFace` (or a new `testLapFace`):
```cpp
    Eigen::SparseMatrix<double> K = ops::facegrad::laplacian(m);
    EXPECT(K.rows()==Fn && K.cols()==Fn, "lapFace is [F x F] = [20, 20]");
    EXPECT((K - Eigen::SparseMatrix<double>(K.transpose())).norm() < 1e-10, "lapFace symmetric");
    // == gradFace' W_F gradFace  (single source of truth)
    Eigen::VectorXd A3(3*Fn); { int fi=0; for(auto f:m.mesh().faces()){ double a=m.geometry().faceAreas[f]; A3[3*fi]=a;A3[3*fi+1]=a;A3[3*fi+2]=a; ++fi; } }
    Eigen::SparseMatrix<double> WF(3*Fn,3*Fn); { std::vector<Eigen::Triplet<double>> tw; for(int i=0;i<3*Fn;++i) tw.emplace_back(i,i,A3[i]); WF.setFromTriplets(tw.begin(),tw.end()); }
    Eigen::SparseMatrix<double> GtWG = (Eigen::SparseMatrix<double>(G.transpose())*WF*G).pruned();
    EXPECT((K-GtWG).norm() < 1e-9*K.norm(), "lapFace == gradFace' W_F gradFace");
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Eigen::MatrixXd(K));
    EXPECT(es.eigenvalues().minCoeff() > -1e-9, "lapFace PSD");
    EXPECT(std::abs(es.eigenvalues()(0)) < 1e-9 && es.eigenvalues()(1) > 1e-9, "lapFace kernel is 1-dim (constants)");
```

- [ ] **Step 2: Implement `ops::facegrad::laplacian`** in `src/face_gradient.cpp`, built from `gradient(m)` (single source of truth, mirroring `extrinsicBlockFace`):
```cpp
Eigen::SparseMatrix<double> laplacian(Manifold& m) {
    Eigen::SparseMatrix<double> G = gradient(m);    // [3F x F]
    auto& mesh = m.mesh(); auto& geom = m.geometry(); geom.requireFaceAreas();
    const int Fn = m.nF();
    std::vector<Eigen::Triplet<double>> tw; tw.reserve(3*Fn);
    for (Face f : mesh.faces()){ const int fi=(int)f.getIndex(); const double a=geom.faceAreas[f];
        for (int d=0; d<3; ++d) tw.emplace_back(3*fi+d, 3*fi+d, a); }
    Eigen::SparseMatrix<double> WF(3*Fn,3*Fn); WF.setFromTriplets(tw.begin(),tw.end());
    Eigen::SparseMatrix<double> K = (Eigen::SparseMatrix<double>(G.transpose())*WF*G).pruned();
    K.makeCompressed(); return K;
}
```

- [ ] **Step 3: Build + run.** `cmake --build build --target test_face_gradient && ./build/test_face_gradient` — all PASS.

---

### Task C3: Register operators (OperatorId + cache + facet accessors)

**Files:** `include/nxr/compute.h`, `src/facets.cpp`, `include/nxr/facets.h`. Mirror `DiracFaceIntrinsicD` exactly (traced: enum @ compute.h:86; cache ptr @ ~188; `...Cached_()` @ facets.cpp:267; facet accessor @ facets.cpp:462 + facets.h:185; `has`/`reset` switch arms @ facets.cpp:114/134).

- [ ] **Step 1:** `include/nxr/compute.h` — extend the enum and add cache members + cached-builder declarations:
```cpp
enum class OperatorId { Laplacian, Mass, Hodge, Gradient3D, Dirac, DiracFace,
    DiracD, DiracFaceD, DiracIntrinsicD, DiracFaceIntrinsicD, GradFace, LapFace };
// ...
std::unique_ptr<Eigen::SparseMatrix<double>> cacheGradFace_;  // dual-mesh face gradient [3F×F]
std::unique_ptr<Eigen::SparseMatrix<double>> cacheLapFace_;   // face Laplacian K̃ [F×F]
// ... in the cached-builder declaration block:
const Eigen::SparseMatrix<double>& gradFaceCached_();
const Eigen::SparseMatrix<double>& lapFaceCached_();
```

- [ ] **Step 2:** `src/facets.cpp` — add `has`/`reset` switch arms (near lines 114/134), the two cached builders (near 267), and the two facet accessors (near 462):
```cpp
// has():    case OperatorId::GradFace: return (bool)cacheGradFace_;
//           case OperatorId::LapFace:  return (bool)cacheLapFace_;
// reset():  case OperatorId::GradFace: cacheGradFace_.reset(); break;
//           case OperatorId::LapFace:  cacheLapFace_.reset();  break;
const Eigen::SparseMatrix<double>& Manifold::gradFaceCached_() {
    if (!cacheGradFace_) cacheGradFace_ = std::make_unique<Eigen::SparseMatrix<double>>(ops::facegrad::gradient(*this));
    return *cacheGradFace_; }
const Eigen::SparseMatrix<double>& Manifold::lapFaceCached_() {
    if (!cacheLapFace_) cacheLapFace_ = std::make_unique<Eigen::SparseMatrix<double>>(ops::facegrad::laplacian(*this));
    return *cacheLapFace_; }
// accessors:
const Eigen::SparseMatrix<double>& OperatorsFacet::gradFace() const { return m_.gradFaceCached_(); }
const Eigen::SparseMatrix<double>& OperatorsFacet::lapFace()  const { return m_.lapFaceCached_(); }
```

- [ ] **Step 3:** `include/nxr/facets.h` — declare the two accessors next to `diracFaceIntrinsicD()` (after line 185):
```cpp
    // gradFace(): barycentric dual-mesh gradient of a per-face scalar [3F×F], cached. Green–Gauss.
    const Eigen::SparseMatrix<double>& gradFace() const;
    // lapFace(): face Laplacian K̃ = gradFace' W_F gradFace [F×F], cached. Kernel = constants (genus-0).
    const Eigen::SparseMatrix<double>& lapFace()  const;
```

- [ ] **Step 4: Build the library.** `cmake --build build` — compiles clean.

---

### Task C4: MEX + WASM dispatch + MATLAB binding test

**Files:** `bindings/mex/src/nxr_compute_mex.cpp`, `bindings/wasm/src/nxr_compute_wasm.cpp`, `bindings/mex/test/test_dirac_first_order.m` (or a new `test_face_gradient.m`).

- [ ] **Step 1:** `nxr_compute_mex.cpp` — add two dispatch arms after `diracFaceIntrinsicD` (line 1877) and extend the error string (line 1882):
```cpp
    } else if (family == "gradFace") {
        plhs[0] = eigenSparseToMx(m.operators().gradFace());   // [3F×F]
    } else if (family == "lapFace") {
        plhs[0] = eigenSparseToMx(m.operators().lapFace());    // [F×F]
```

- [ ] **Step 2:** `nxr_compute_wasm.cpp` — add the matching `family == "gradFace" / "lapFace"` arms (mirror the `diracFaceIntrinsicD` arm @ ~351).

- [ ] **Step 3: Write a MATLAB binding test** (append to `bindings/mex/test/test_dirac_first_order.m`): on the bundled closed test mesh, pull `G=operators(h,'gradFace')`, `K=operators(h,'lapFace')`; assert sizes `[3F F]`/`[F F]`, `norm(G*ones(F,1))<1e-10`, `norm(K-K')<1e-10`, and `K ≈ G'*kron(diag(area),I3)*G`.

- [ ] **Step 4: Build the mex** (per the repo's mex build target/script) and run the binding test in MATLAB. All PASS.

---

### Task C5: Install + commit

- [ ] **Step 1: Back up + install the rebuilt mex** into the managed-plugin folder (stale-binary trap, per memory `nxr-bundle-surface-fields`):
```bash
cp ~/.brainstorm/plugins/nxr-compute/*-r2023b/nxr_compute.mexmaca64 \
   ~/.brainstorm/plugins/.../nxr_compute.mexmaca64.bak-$(date +%Y%m%d)-pregradface
# then copy build/.../nxr_compute.mexmaca64 into that -r2023b plugin folder
```
- [ ] **Step 2: Run the full C++ suite** `ctest --test-dir build` — green.
- [ ] **Step 3: Commit on the feature branch** (do NOT cut a `v*` tag):
```bash
git add -A && git commit -m "feat(facegrad): barycentric dual face gradient gradFace + face Laplacian lapFace"
```

---

## PHASE D — Brainstorm face-native Helmholtz rewrite (MATLAB)

### Task D1: `Prepare` — pull the new operators, build the cached pieces

**Files:** Modify `toolbox/math/bst_dirac_helmholtz_face.m`.

- [ ] **Step 1: Rewrite `Prepare`** to cache, per hemisphere: `G̃ = nxr_compute('operators',h,'gradFace')` `[3F×F]`, `K̃ = nxr_compute('operators',h,'lapFace')` `[F×F]`, per-face normals `Nf`, areas `Af`, `SkewG` (apply `n_f×` to each per-face 3-block of `G̃`), and `cholK = decomposition(K̃(free,free),'chol')` with `free = 2:nFh` (pin face 1). Store dual face-adjacency `NbF{hh}` (faces sharing an edge) for core detection. Drop `Dt`/`Gx/Gy/Gz`/vertex `cholK`. Build `SkewG` as: reshape each column's `[3F]` into `[F×3]`, cross with `Nf`, restack — or precompute a sparse `[3F×3F]` block-rotation `Rn` (per face `v ↦ n_f×v`) and set `SkewG = Rn*G̃`.

- [ ] **Step 2: Smoke-check in MATLAB (MCP):** `Op = bst_dirac_helmholtz_face('Prepare', Dirac, LBO, Surf)` returns `numel(Op.G)==2`, `size(Op.G{1},1)==3*numel(Op.fH{1})`, `size(Op.K{1},1)==numel(Op.fH{1})`.

---

### Task D2: `Frame` — face-native decomposition

**Files:** Modify `toolbox/math/bst_dirac_helmholtz_face.m`.

- [ ] **Step 1: Rewrite `Frame(Op, Jf, withCores)`** per hemisphere (`Jf` `[nF×3]`; local faces `fH`):
```matlab
    Jcol = reshape(Jl', [], 1);                 % [3F x 1] stacked (x,y,z per face)
    WF   = repelem(Af, 3);                       % [3F] area weights
    divS = Op.G{hh}'    * (WF .* Jcol);          % G̃' W_F J   -> [F] divergence source
    curlS= Op.SkewG{hh}'* (WF .* Jcol);          % SkewG' W_F J-> [F] curl source
    phi  = i_poisson(Op.cholK{hh}, Op.K{hh}, divS,  Op.free{hh});   % K̃ phi = divS, mean-zero
    psi  = i_poisson(Op.cholK{hh}, Op.K{hh}, curlS, Op.free{hh});
    Virr = reshape(Op.G{hh}    *phi, 3, [])';    % [F x 3]
    Vsol = reshape(Op.SkewG{hh}*psi, 3, [])';
    Vharm= Jl - Virr - Vsol;
```
HarmFrac = `Σ Af·|Vharm|² / Σ Af·|Jl|²`. Scalars now `[nF×1]`: `Ht.Curl(fH)=curlS; Ht.Div(fH)=divS; Ht.Psi(fH)=psi; Ht.Phi(fH)=phi;`. Component fields `[nF×3]` as before. (`i_poisson` pins/mean-zeros on faces now — the existing helper works once `M`→`speye`/area mass on faces and `free` is face-indexed.)

- [ ] **Step 2: Smoke-check** a random face field decomposes and `Virr+Vsol+Vharm == Jf` to 1e-9.

---

### Task D3: Face-domain core detection

**Files:** Modify `toolbox/math/bst_dirac_helmholtz_face.m`; verify `toolbox/math/bst_vortex_persistence.m` accepts a face-neighbor graph via its `'Neighbors'` option.

- [ ] **Step 1:** Rewrite `i_find_cores` to run persistence on the per-FACE potential (`Ht.Psi`/`Ht.Phi`) over the dual adjacency `Op.NbF{hh}`; `pos` = face barycenter. Cores/Sources keep the same struct fields (`iVertex`→reuse as face index or add `iFace`; charge = sign of `curlS`/`divS` at that face).

- [ ] **Step 2: Smoke-check** cores come back as a persistence-sorted struct array with `pos` on face centroids.

---

### Task D4: Tighten the test to the strict round-trip

**Files:** Modify `dev/tests/test_dirac_helmholtz_face.m`.

- [ ] **Step 1:** Replace the planted-skew-gradient block to plant with the NEW operator and gate strictly:
```matlab
    hh=1; fH=Op.fH{hh}; Af=Op.Af{hh};
    % smooth per-face scalar psi0 from face-centroid geometry
    Cf = (Surf.Vertices(Surf.Faces(fH,1),:)+Surf.Vertices(Surf.Faces(fH,2),:)+Surf.Vertices(Surf.Faces(fH,3),:))/3;
    c = round(numel(fH)/2); d2 = sum((Cf-Cf(c,:)).^2,2); psi0 = exp(-d2/(2*0.012^2)); psi0 = psi0-mean(psi0);
    Vsol0 = reshape(Op.SkewG{hh}*psi0, 3, [])';            % planted pure solenoidal on faces
    Jsk = zeros(size(Surf.Faces,1),3); Jsk(fH,:) = Vsol0;
    Hsk = bst_dirac_helmholtz_face('Frame', Op, Jsk);
    nFail = nFail + chk('planted skew-gradient: HarmFrac < 0.02', Hsk.HarmFrac < 0.02);
    nFail = nFail + chk('planted skew-gradient: corr(psi,psi0) > 0.99', abs(corr(Hsk.Psi(fH), psi0)) > 0.99);
```
Update the earlier shape/size checks: `Curl/Div/Psi/Phi` are now `[nF×1]`; `Virr/Vsol` `[nF×3]`.

- [ ] **Step 2: Run** `dev/tests/test_dirac_helmholtz_face.m` via the MCP. Expected: `0 failed`, with the planted field at HarmFrac<0.02, corr>0.99. *(If corr≈−0.99, flip the `w` sign in `gradFace` per the Task-C1 note, rebuild+reinstall the mex, re-run.)*

---

### Task D5: Comparison benchmark + commit

**Files:** `dev/benchmarks/bench_dirac_face_helmholtz.m` (already exists).

- [ ] **Step 1: Run** `R = bench_dirac_face_helmholtz(22.6)` via the MCP. Now apples-to-apples: expect face HarmFrac LOW (genus-0, no spurious harmonic), curl/div finite, cores of comparable count to the vertex pipeline. Inspect the PNG.
- [ ] **Step 2: Run all suites** — `test_dirac_helmholtz_face`, `test_dirac_helmholtz`, `test_vortex_persistence`, `test_vortex_track`, `test_helmholtz_track`, `test_time_derivative`: all `0 failed`.
- [ ] **Step 3: Record the finding** in chat (face vs vertex HarmFrac/cores) and update memory `face-domain-dirac`.
- [ ] **Step 4: Commit on `development`:**
```bash
git add toolbox/math/bst_dirac_helmholtz_face.m dev/tests/test_dirac_helmholtz_face.m dev/benchmarks/bench_dirac_face_helmholtz.m docs/superpowers/
git commit -m "feat(face-helmholtz): face-native Helmholtz via gradFace/lapFace (HarmFrac->0)"
```

---

## Self-review

- **Spec coverage:** `gradFace` (C1) ✓, `lapFace` (C2) ✓, registration (C3) ✓, bindings (C4) ✓, install/commit (C5) ✓; Brainstorm `Prepare`/`Frame`/cores/test/bench (D1–D5) ✓. The §3 contract maps to C1/C2 tests + the D4 round-trip.
- **Placeholder scan:** the only derived constant is `gradFace`'s formula — given in full C++; its one ambiguity (the `n_f×l` sign) is called out with an explicit fix and pinned by the D4 round-trip test. No "TODO"/"handle edge cases".
- **Type/name consistency:** `ops::facegrad::gradient`/`laplacian`; `OperatorId::GradFace`/`LapFace`; facet `gradFace()`/`lapFace()`; MEX/WASM family `"gradFace"`/`"lapFace"`; MATLAB `Op.G`/`Op.K`/`Op.SkewG` used consistently across D1/D2.
- **Risk:** if barycentric conditioning ever fails on real cortex (it shouldn't on ico meshes), the spec's §6 fallback (circumcentric/least-squares) is an operator-only swap behind the same API.
