# Face-domain Helmholtz Phase B (intrinsic) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax. Task 1 is a DECISIVE EXPERIMENT whose result selects the Task-2 path.

**Goal:** Make `bst_dirac_helmholtz_face` a *consistent* intrinsic face-domain Helmholtz using the new `diracFaceIntrinsicD`, gated by a planted pure-skew-gradient field recovering with HarmFrac≈0.

**Architecture:** First swap the extrinsic `diracFaceD` → intrinsic `diracFaceIntrinsicD` in the existing pipeline and measure the planted-field recovery (Task 1). If the intrinsic operator with the vertex cotan Poisson already recovers (HarmFrac→0, corr→1), Phase B is the swap. If not, build the consistent dual: Poisson on the intrinsic face Laplacian `K̃ᵢₙₜ` (scalar block of `D̃ᵢₙₜᵀ W_V D̃ᵢₙₜ`) with a `D̃ᵢₙₜ`-consistent dual reconstruction, derived against the same gate (Task 2).

**Tech Stack:** MATLAB (Brainstorm); nxr `diracFaceIntrinsicD` (installed). Tests via the MATLAB MCP.

**Preconditions:** the Phase-A mex with `diracFaceIntrinsicD` is installed (done); `TutorialAuditory` protocol; toolbox + dev/tests on path.

---

### Task 1: DECISIVE EXPERIMENT — swap to the intrinsic operator + planted-field gate

**Files:** Modify `toolbox/math/bst_dirac_helmholtz_face.m`; Modify `dev/tests/test_dirac_helmholtz_face.m`.

- [ ] **Step 1: Add the planted-field gate to the test.** Append to `test_dirac_helmholtz_face.m` before the final fprintf — a pure skew-gradient face field must recover with HarmFrac≈0 and recovered potential ≈ the planted one:
```matlab
    % --- planted pure-skew-gradient (solenoidal) face field must recover (HarmFrac->0) ---
    hh=1; vH=Op.vH{hh}; fH=Op.fH{hh};
    c = vH(round(numel(vH)/2)); d2 = sum((Surf.Vertices(vH,:)-Surf.Vertices(c,:)).^2,2);
    psi0 = exp(-d2/(2*0.012^2)); psi0 = psi0 - mean(psi0);
    gp = [Op.Gx{hh}*psi0, Op.Gy{hh}*psi0, Op.Gz{hh}*psi0];
    Vsk = cross(Op.Nf{hh}, gp, 2);                          % n x grad(psi0): pure solenoidal
    Jsk = zeros(nF,3); Jsk(fH,:) = Vsk;
    Hsk = bst_dirac_helmholtz_face('Frame', Op, Jsk);
    nFail = nFail + chk('planted skew-gradient: HarmFrac < 0.05', Hsk.HarmFrac < 0.05);
    nFail = nFail + chk('planted skew-gradient: recovered psi corr > 0.95', ...
        abs(corr(Hsk.Psi(vH), psi0)) > 0.95);
```

- [ ] **Step 2: Run, capture the BASELINE (extrinsic) failure** — via the MCP run `dev/tests/test_dirac_helmholtz_face.m`. Expected: the two new checks FAIL (HarmFrac ≫ 0.05, corr ≈ 0.36) — this is the documented extrinsic inconsistency, now an explicit gate.

- [ ] **Step 3: Swap to the intrinsic operator.** In `bst_dirac_helmholtz_face.m` `Prepare`, change the operator pull:
```matlab
        Dt = nxr_compute('operators', h, 'diracFaceIntrinsicD');   % was 'diracFaceD'
```

- [ ] **Step 4: Re-run the gate** — `dev/tests/test_dirac_helmholtz_face.m`. **DECISION:**
  - If both planted-field checks now PASS → the intrinsic operator with the existing vertex cotan Poisson is consistent; Phase B is essentially done. Skip to Task 3.
  - If they still FAIL (expected possibility: `D̃ᵢₙₜ`'s natural Laplacian is the face one, not the vertex cotan) → proceed to Task 2 (full consistent dual). Record the corr/HarmFrac numbers in chat as the decision evidence.

- [ ] **Step 5: Commit (user-gated)** `git add -A && git commit -m "wip(face-helmholtz): intrinsic diracFaceIntrinsicD + planted-field gate"`

---

### Task 2 (CONDITIONAL — only if Task 1 Step 4 still fails): consistent dual pipeline

**Files:** Modify `toolbox/math/bst_dirac_helmholtz_face.m`.

This task is a **derive-validate-iterate loop** against the Task-1 planted-field gate (HarmFrac<0.05, corr>0.95). Use superpowers:systematic-debugging. The candidates, in order:

- [ ] **Candidate A — Poisson on the intrinsic face Laplacian.** Build `K̃ᵢₙₜ` once per hemi in `Prepare`: `Et = Dt' * WV * Dt` (WV = `kron(diag(vertexDualArea), I4)`; vertex dual area = `(1/3)·Σ incident face areas`), then `Kw = Et(1:4:end,1:4:end)` (scalar block, [F×F]). In `Frame`, average the per-vertex `ωV/dvV` to faces (`ωF = Wvf·ωV`, `Wvf` = area-weighted vertex→face mean, the transpose-normalized dual of the vertex pipeline's `Wfv`), then solve `Kw ψ̃ = M_F ωF` on FACES (`M_F` = `diag(faceArea)`; pin one face, mean-zero). Run the gate.

- [ ] **Candidate B (if A's reconstruction is off) — `D̃ᵢₙₜ`-adjoint reconstruction.** Reconstruct the solenoidal field directly from the face stream function via the operator's own adjoint structure: `Vsol = ` the imaginary part of `Dt' * qPsi` where `qPsi` embeds `ψ̃` as the w-part per face — i.e. use `D̃ᵢₙₜᵀ` (not an independent FEM gradient) so the reconstruction is adjoint-consistent with the curl extraction by construction. Validate the round-trip on the planted field.

- [ ] **Candidate C (fallback) — circumcentric immersion.** If barycentric centroids give a poorly-conditioned `K̃ᵢₙₜ`, the operator-level fix is circumcentric duals — out of scope for MATLAB; would loop back to a Phase-A nxr change. Record the finding and stop rather than hack.

- [ ] **Gate (run after each candidate):** `dev/tests/test_dirac_helmholtz_face.m` — the planted skew-gradient must hit HarmFrac<0.05, corr>0.95, AND a planted pure-gradient field must be divergence-dominated. Stop at the first candidate that passes; document which in the Frame comment.

- [ ] **Commit (user-gated)** once a candidate passes.

---

### Task 3: vertex-vs-face comparison on the real frame

**Files:** Modify `dev/benchmarks/bench_dirac_face_helmholtz.m` (already exists from the earlier prototype).

- [ ] **Step 1: Run the comparison** — `R = bench_dirac_face_helmholtz(22.6)` via the MCP. With the consistent intrinsic pipeline, expect: face HarmFrac low (genus-0), curl/div finite, cores of comparable count to the vertex pipeline, and the side-by-side ψ figure broadly matching (now apples-to-apples). Inspect the PNG.

- [ ] **Step 2: Run all suites** — `test_dirac_helmholtz_face`, `test_dirac_helmholtz`, `test_vortex_persistence`, `test_vortex_track`, `test_helmholtz_track`, `test_time_derivative`: all `0 failed`.

- [ ] **Step 3: Record the comparison finding** in chat — the empirical input to whether a full face-based LEADFIELD + inverse (the next phase) is worth building.

- [ ] **Step 4: Commit (user-gated)** `git add -A && git commit -m "feat(face-helmholtz): intrinsic face-domain Helmholtz + comparison"`

---

## Self-review notes

- **Spec coverage:** intrinsic operator wired into the face pipeline (Task 1 Step 3) ✓; consistent dual Poisson + reconstruction as a gated derivation (Task 2) ✓; planted-field gate (Task 1 Step 1) ✓; comparison (Task 3) ✓.
- **Honesty:** Task 2's reconstruction is genuinely exploratory (the spec's flagged open piece); the plan structures it as named candidates A→B→C against an explicit gate rather than pretending it's fixed code.
- **Placeholder scan:** none of the forbidden kind; Candidate steps are concrete with explicit fallbacks, gated by a runnable test.
- **Risk:** if no MATLAB-side candidate hits the gate, the fix is an operator-level (nxr) change — Candidate C says stop-and-report rather than hack.
```
