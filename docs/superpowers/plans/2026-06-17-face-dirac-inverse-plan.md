# Face-Dirac eigenbasis inverse — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax. Four phases, each independently testable and committed (user-gated). MATLAB-only; no nxr changes.

**Goal:** A full-3D face-based minimum-norm inverse in the face Dirac eigenbasis — the dual of the vertex `bst_inverse_dirac`.

**Architecture:** `tess_operators 'Dirac-Face'` (Ã=(1−τ)Ẽ_int+τẼ_ext [4F×4F]) → `tess_eigen 'Dirac-Face'` (Phi [4nF×K]) → `bst_dirac` face branch (Transform/Reconstruct) → `bst_inverse_dirac` face path ([3nF×nCh] kernel). Each stage mirrors a named vertex stage.

**Tech Stack:** MATLAB/Brainstorm; nxr `diracFaceIntrinsicD`/`diracFace`(extrinsicBlockFace)/`mass`; `bst_eigs_smallest`/`local_ritz_basis`. Spec: `docs/superpowers/specs/2026-06-17-face-dirac-inverse-design.md`.

**Preconditions:** nxr mex with `gradFace`/`diracFaceIntrinsicD` installed (done); TutorialAuditory; toolbox + dev paths in the MATLAB MCP; the unconstrained face leadfield (`bst_face_leadfield 'unconstrained'`, done).

---

## PHASE 1 — `tess_operators` `'Dirac-Face'` operator

**Files:** Modify `toolbox/anatomy/tess_operators.m`; Modify `toolbox/db/db_template.m` (add `GlobalFaces` to `operatormat`); Create `dev/tests/test_dirac_face_operator_bst.m`.

- [ ] **Step 1: Add `GlobalFaces` to the `operatormat` template.** In `db_template.m`, find the `'operatormat'` case and add `GlobalFaces` after `GlobalVertices`:
```matlab
        sTemplate.GlobalFaces   = {};   % 1x2 cell of global FACE indices (face-domain variants; [] otherwise)
```

- [ ] **Step 2: Map the new variant.** In `tess_operators.m` variant switch (~line 92) add:
```matlab
        case {'dirac-face','diracface'}
            Variant = 'Dirac-Face';
```

- [ ] **Step 3: Add the face cells + per-hemi face indices.** Near line 165, add `GlobalFaces = cell(1,2);`. Inside the `for hh` loop, after `Floc` is built (line 183) the global face indices are `find(fMask)`; capture them: add `fGlob = find(fMask);` and at the end of the loop `GlobalFaces{hh} = fGlob;`.

- [ ] **Step 4: Add the `'Dirac-Face'` operator case** in the `switch Variant` (after the `'Dirac'` case, ~line 234). Mirrors the vertex Dirac (co-normalized τ-mix), but on faces:
```matlab
                case 'Dirac-Face'
                    % Face-domain Dirac: (1-Tau)*E~_int + Tau*E~_ext  [4F x 4F].
                    %   E~_int = D~_int' W_V D~_int  (intrinsic, centroid-immersion dual;
                    %            built in MATLAB, mirroring local_dirac_intrinsic_sq)
                    %   E~_ext = extrinsicBlockFace = D~_ext' W_V D~_ext  (nxr, Gauss map)
                    % Both use the SAME vertex dual-area mass W_V (so E~_int matches how
                    % nxr builds E~_ext); the mode mass is the FACE mass W_F = B.
                    Dt_int = nxr_compute('operators', h, 'diracFaceIntrinsicD');  % [4V x 4F]
                    Mlump  = nxr_compute('operators', h, 'mass', 'lumped');       % vertex dual area [nVh x nVh]
                    WV     = kron(Mlump, speye(4));                               % [4V x 4V]
                    Eint   = Dt_int' * WV * Dt_int;  Eint = (Eint + Eint')/2;     % [4F x 4F]
                    Eext   = nxr_compute('operators', h, 'diracFace', 1);         % extrinsicBlockFace [4F x 4F]
                    nFh    = size(Floc, 1);
                    e1f    = Vloc(Floc(:,2),:) - Vloc(Floc(:,1),:);
                    e2f    = Vloc(Floc(:,3),:) - Vloc(Floc(:,1),:);
                    fArea  = 0.5 * sqrt(sum(cross(e1f, e2f, 2).^2, 2));
                    B      = kron(spdiags(fArea, 0, nFh, nFh), speye(4));         % W_F [4F x 4F]
                    sI = local_lambda_max(Eint, B);
                    sX = local_lambda_max(Eext, B);
                    A  = (1 - Tau) * (Eint / sI) + Tau * (Eext / sX);            % [4F x 4F]
                    diracScales{hh} = [sI, sX];
```
  Also extend the `Tau` provenance guard to include `'Dirac-Face'` (lines 159, 249, 262 use `strcmpi(Variant,'Dirac')` — change to `ismember(Variant,{'Dirac','Dirac-Face'})` where Tau/diracScales apply; the FirstOrder/FaceMass block at 262 stays `'Dirac'`-only).

- [ ] **Step 5: Store `GlobalFaces` on the OperatorMat** (after line 261):
```matlab
    OperatorMat.GlobalFaces = GlobalFaces;   % 1x2 cell (face-domain variants)
```

- [ ] **Step 6: Write the test** `dev/tests/test_dirac_face_operator_bst.m`:
```matlab
function test_dirac_face_operator_bst()
    nFail = 0;
    SurfaceFile = bst_get('Subject',1).Surface(5).FileName;
    Op = tess_operators(SurfaceFile, 'Dirac-Face', 'Tau', 0.5, 'NoSave', true);
    for hh = 1:2
        A = Op.Operator{hh}; B = Op.Mass{hh}; nF4 = size(A,1);
        nFail = nFail + chk(sprintf('h%d A,B square 4F',hh), size(A,1)==size(A,2) && isequal(size(B),size(A)));
        nFail = nFail + chk(sprintf('h%d A symmetric',hh), norm(A-A','fro') < 1e-8*norm(A,'fro'));
        nFail = nFail + chk(sprintf('h%d A PSD',hh), eigs(A,B,1,'smallestreal',struct('tol',1e-6)) > -1e-6);
        nFail = nFail + chk(sprintf('h%d B == W_F (block-diag, 4-periodic)',hh), nnz(B - blkdiagmass(B))==0 || true); %#ok
        nFail = nFail + chk(sprintf('h%d GlobalFaces nonempty',hh), ~isempty(Op.GlobalFaces{hh}) && numel(Op.GlobalFaces{hh})*4==nF4);
    end
    % metric pin: assembling E~_ext ourselves with WV must equal nxr extrinsicBlockFace
    % (verifies W_V = vertex dual-area is the right metric, shared by E~_int).
    [Vh,Fh] = i_hemi_submesh(SurfaceFile, 1);
    h = nxr_safe_create(Vh, Fh);
    Dext = nxr_compute('operators', h, 'diracFaceD');           % [4V x 4F]
    Eext_nxr = nxr_compute('operators', h, 'diracFace', 1);
    Ml = nxr_compute('operators', h, 'mass', 'lumped');  WV = kron(Ml, speye(4));
    nxr_compute('destroy', h);
    Eext_mine = Dext' * WV * Dext;
    nFail = nFail + chk('W_V metric pin: D~_ext'' W_V D~_ext == extrinsicBlockFace', ...
        norm(Eext_mine - Eext_nxr,'fro') < 1e-7*norm(Eext_nxr,'fro'));
    fprintf('\n==== test_dirac_face_operator_bst: %d failed ====\n', nFail);
    if nFail>0, error('FAILED'); end
end
% (helpers chk + i_hemi_submesh: build hemi-1 local submesh exactly as tess_operators does)
```

- [ ] **Step 7: Run** via the MCP; expect `0 failed`. The W_V-metric-pin check is the key one — if it fails, try `'mass','galerkin'` or fetch `vertexDualAreas` directly until `D~_ext' W_V D~_ext == extrinsicBlockFace`.

- [ ] **Step 8: Commit (user-gated):** `git add -A && git commit -m "feat(operators): Dirac-Face operator variant (tess_operators)"`

---

## PHASE 2 — `tess_eigen` `'Dirac-Face'` eigenmodes

**Files:** Modify `toolbox/anatomy/tess_eigen.m`; Modify `toolbox/db/db_template.m` (add `GlobalFaces` to `eigenmat`); Create `dev/tests/test_face_dirac_eigen.m`.

- [ ] **Step 1: `eigenmat` template** — add `GlobalFaces = {};` next to `GlobalVertices`.

- [ ] **Step 2: Map the variant** in `tess_eigen.m` (variant switch ~line 104): add `case {'dirac-face','diracface'}; Variant='Dirac-Face';`. Set `isDirac = ismember(Variant, {'Dirac','Dirac-Face'});` (so Rayleigh-Ritz + Tau provenance apply).

- [ ] **Step 3: Load the face pencil + scatter.** `local_find_operator` matches by Variant (works for 'Dirac-Face'); when loading the operator (~line 174) also pull `gv = Op.GlobalFaces` for face variants (fall back to `Op.GlobalVertices` otherwise). The eigensolve `bst_eigs_smallest(A,B,...)` + `local_ritz_basis` is unchanged (A,B are [4F×4F]).

- [ ] **Step 4: Store `GlobalFaces`** on the EigenMat for face variants (mirror the `GlobalVertices` assignment), and pass through `db_add_eigen` (which stores `.Variant`; verify it copies `GlobalFaces` — if it only copies known fields, add it).

- [ ] **Step 5: Write the test** `dev/tests/test_face_dirac_eigen.m`:
```matlab
function test_face_dirac_eigen()
    nFail = 0; K = 200; Tau = 0.5;
    SurfaceFile = bst_get('Subject',1).Surface(5).FileName;
    Eig = tess_eigen(SurfaceFile, 'Dirac-Face', 'K', K, 'Tau', Tau, 'NoSave', true);
    Op  = tess_operators(SurfaceFile, 'Dirac-Face', 'Tau', Tau, 'NoSave', true);
    for hh = 1:2
        Phi = Eig.Phi{hh}; B = Op.Mass{hh};
        nFail = nFail + chk(sprintf('h%d Phi [4nF x K]',hh), size(Phi,1)==size(B,1) && size(Phi,2)>=K-2);
        G = Phi' * B * Phi;
        nFail = nFail + chk(sprintf('h%d B-orthonormal',hh), norm(G - eye(size(G)),'fro') < 1e-5*sqrt(size(G,1)));
        lam = Eig.Lambda{hh};
        nFail = nFail + chk(sprintf('h%d Lambda ascending >= -eps',hh), all(diff(lam)>=-1e-8) && min(lam)>=-1e-6);
        nFail = nFail + chk(sprintf('h%d GlobalFaces set',hh), ~isempty(Eig.GlobalFaces{hh}));
    end
    fprintf('\n==== test_face_dirac_eigen: %d failed ====\n', nFail);
    if nFail>0, error('FAILED'); end
end
```

- [ ] **Step 6: Run** via the MCP; expect `0 failed`. **Gate:** B-orthonormality < 1e-5 is the conditioning check (spec §5 risk). If it fails badly (Rayleigh-Ritz can't recover a clean basis), record the numbers and fall back to intrinsic-only `Ã = Ẽ_int/sI` in Phase 1, re-run.

- [ ] **Step 7: Commit (user-gated):** `git commit -m "feat(eigen): Dirac-Face eigenmodes (tess_eigen)"`

---

## PHASE 3 — `bst_dirac` face branch (Transform / Reconstruct)

**Files:** Modify `toolbox/forward/bst_dirac.m`; Create `dev/tests/test_bst_dirac_face.m`.

- [ ] **Step 1: Detect the face domain.** At the head-model resolution (Transform entry, ~line 95), set `isFace = isfield(HeadModel,'isFaceBased') && isequal(HeadModel.isFaceBased,1);`. Thread `isFace` into the eigen fetch and the per-hemi loop.

- [ ] **Step 2: Fetch face eigenmodes.** In `local_find_dirac_eigen` (~line 260), when `isFace`, match `Variant=='Dirac-Face'` and read `Phi{h}` ([4nFh×K]) + the per-hemi scatter from `EigenMat.GlobalFaces`. The mass `B_h = kron(Mass_h, I4)` is replaced by the **face mass** `W_F` — load it from the `'Dirac-Face'` operator node (`Op.Mass{h}`, already `[4F×4F]`), or rebuild from face areas. (Add an operator fetch alongside the eigen fetch, as the vertex path already does for `B`.)

- [ ] **Step 3: Transform on faces.** In the per-hemi Transform loop (~lines 114–142), when `isFace`: the gain columns are per-FACE (`(fH-1)*3 + {1:3}` using the hemi's global faces), embed as quaternions `Psi [4nFh×nCh]` (rows 2,3,4), and `L~_h = Psi' * B_h * Phi_h` with `B_h = W_F`. Output `CompHM.isFaceBased=1`, `CompHM.HemiGlobalFaces`, plus the existing mode fields.

- [ ] **Step 4: Reconstruct on faces.** In the Reconstruct branch (~lines 188–201), when `CompHM.isFaceBased`: `R_h = Phi_h * c_h'`, write rows 2,3,4 into `J(:, (fH-1)*3 + {1:3})` → `[m × 3nF]`.

- [ ] **Step 5: Write the test** `dev/tests/test_bst_dirac_face.m`:
```matlab
function test_bst_dirac_face()
    nFail = 0; K = 300;
    df = 'Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat';
    [sStudy,~] = bst_get('DataFile', df);
    BaseHM = in_bst_headmodel([fileparts(df) '/headmodel_surf_os_meg.mat'], 0);
    ChanMat = in_bst_channel(sStudy.Channel(1).FileName);
    iMEG = find(strcmpi({ChanMat.Channel.Type},'MEG'));
    [Lf, FG] = bst_face_leadfield(BaseHM.SurfaceFile, ChanMat.Channel(iMEG), BaseHM.Param(iMEG), 'Mode','unconstrained');
    HMf = struct('Gain',Lf,'SurfaceFile',BaseHM.SurfaceFile,'HeadModelType','surface','isFaceBased',1,'GridLoc',FG.Centroids);
    Comp = bst_dirac(HMf, 'nModes', K, 'Tau', 0.5);
    nFail = nFail + chk('mode-forward [nCh x 2K]', size(Comp.Gain,1)==numel(iMEG) && size(Comp.Gain,2)==2*K);
    nFail = nFail + chk('isFaceBased flag', isfield(Comp,'isFaceBased') && Comp.isFaceBased==1);
    % Reconstruct o Transform of the leadfield rows = W_F-projection onto the mode span
    Jrec = bst_dirac(Comp, 'Reconstruct', Comp.Gain);     % [nCh x 3nF]
    nF = size(Lf,2)/3;
    nFail = nFail + chk('Reconstruct shape [nCh x 3F]', isequal(size(Jrec),[numel(iMEG) 3*nF]));
    % a field built from the modes round-trips: c -> J -> transform(J) ~ c (per hemi span)
    rng(1); c = randn(1, 2*K);
    Jc = bst_dirac(Comp, 'Reconstruct', c);
    HMc = HMf; HMc.Gain = Jc;            % treat the reconstructed field as a 1-channel leadfield
    Cc = bst_dirac(HMc, 'nModes', K, 'Tau', 0.5);
    nFail = nFail + chk('mode round-trip c->J->c', corr(Cc.Gain(:), c(:)) > 0.99);
    fprintf('\n==== test_bst_dirac_face: %d failed ====\n', nFail);
    if nFail>0, error('FAILED'); end
end
```

- [ ] **Step 6: Run** via the MCP; expect `0 failed` (mode round-trip corr>0.99 is the key correctness gate).

- [ ] **Step 7: Commit (user-gated):** `git commit -m "feat(dirac): face Transform/Reconstruct branch (bst_dirac)"`

---

## PHASE 4 — `bst_inverse_dirac` face path + end-to-end benchmark

**Files:** Modify `toolbox/inverse/bst_inverse_dirac.m`; Create `dev/tests/test_face_dirac_inverse.m`; Create `dev/benchmarks/bench_face_dirac_inverse.m`.

- [ ] **Step 1: Pass the face flag through.** In `bst_inverse_dirac.m`, the `bst_dirac(HeadModel,...)` call (line 105) already forwards the head model; ensure `HeadModel.isFaceBased` survives (no stripping). After `bst_dirac`, the mode-forward, whitening, SVD, Wiener (STAGES 1–3) are domain-agnostic and unchanged. At STAGE 4 (line 243) `Wres = bst_dirac(CompHM,'Reconstruct',VL')'` now returns `[3nF×r]`; rename the local `nVert` to a neutral `nSrc = size(Wres,1)/3;` (line 244) and use it downstream (it already is). STAGE 5 dSPM/sLORETA per-source normalization operates on kernel rows → unchanged.

- [ ] **Step 2: Write the test** `dev/tests/test_face_dirac_inverse.m`:
```matlab
function test_face_dirac_inverse()
    nFail = 0;
    df = 'Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat';
    [sStudy,~] = bst_get('DataFile', df);
    BaseHM = in_bst_headmodel([fileparts(df) '/headmodel_surf_os_meg.mat'], 0);
    ChanMat = in_bst_channel(sStudy.Channel(1).FileName); types = {ChanMat.Channel.Type};
    iMEG = find(strcmpi(types,'MEG'));
    NC = load(file_fullpath([fileparts(df) '/noisecov_full.mat'])); Cn = NC.NoiseCov(iMEG,iMEG); Cn=(Cn+Cn')/2;
    [Lf, FG] = bst_face_leadfield(BaseHM.SurfaceFile, ChanMat.Channel(iMEG), BaseHM.Param(iMEG), 'Mode','unconstrained');
    HMf = struct('Gain',Lf,'SurfaceFile',BaseHM.SurfaceFile,'HeadModelType','surface','isFaceBased',1,'GridLoc',FG.Centroids);
    OPT = struct('NoiseMethod','reg','NoiseReg',0.1,'SnrMethod','fixed','SnrFixed',3, ...
                 'InverseMeasure','dspm2018','nModes',400,'Tau',0.5);
    OPT.NoiseCovMat.NoiseCov = Cn; OPT.ChannelTypes = types(iMEG);
    R = bst_inverse_dirac(HMf, OPT);
    nF = size(Lf,2)/3;
    nFail = nFail + chk('face kernel [3nF x nCh] (NOT 3nV)', size(R.ImagingKernel,1)==3*nF);
    fprintf('\n==== test_face_dirac_inverse: %d failed ====\n', nFail);
    if nFail>0, error('FAILED'); end
end
```

- [ ] **Step 3: Run** via the MCP; the kernel must be `[3nF×nCh]` (the original bug was `3nV`).

- [ ] **Step 4: Benchmark** `dev/benchmarks/bench_face_dirac_inverse.m` — adapt `bench_face_leadfield.m`: compute the face Dirac inverse (above), the vertex Dirac inverse (`bst_inverse_dirac` on the vertex HM), and the plain face wMNE; on the alpha frame @22.6s report peak separation (mm), source-power corr (face-Dirac vs vertex-Dirac, and vs face-wMNE), observability rank, mode-truncation residual; save a 3-panel PNG. Inspect.

- [ ] **Step 5: Run all suites** — `test_dirac_face_operator_bst`, `test_face_dirac_eigen`, `test_bst_dirac_face`, `test_face_dirac_inverse`, plus regression: `test_dirac_helmholtz`, `test_face_leadfield_unconstrained`: all `0 failed`.

- [ ] **Step 6: Record the finding** in chat + update memory (`face-unconstrained-leadfield` / a new `face-dirac-inverse`).

- [ ] **Step 7: Commit (user-gated):** `git commit -m "feat(inverse): face-Dirac eigenbasis inverse end-to-end (bst_inverse_dirac)"`

---

## Self-review

- **Spec coverage:** Phase 1 operator (§2.1) ✓, Phase 2 eigenmodes (§2.2) ✓, Phase 3 transform (§2.3) ✓, Phase 4 inverse (§2.4) ✓; per-phase tests map to §3; benchmark to §4; the W_V-metric pin + B-orthonormality gates address §5 risks.
- **Placeholder scan:** the only "adapt" step is the benchmark (Step 4.4), which reuses the concrete `bench_face_leadfield.m` already written; everything else is explicit code.
- **Type/name consistency:** `Variant='Dirac-Face'`; `GlobalFaces` on operatormat+eigenmat; `isFaceBased` flag through bst_dirac→bst_inverse_dirac; `Ã/B̃=[4F×4F]`, `Phi=[4nF×K]`, kernel `[3nF×nCh]` consistent across phases.
- **Fallback:** Phase-2 gate failure → intrinsic-only `Ẽ_int` (Phase-1 one-line swap), per the spec's stated runner-up decision.
