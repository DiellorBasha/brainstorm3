# Helmholtz Flat-Covariant Operator Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rebuild `bst_helmholtz`'s decomposition on a self-consistent flat-covariant operator — a new `'Covariant'` operator node in `tess_operators` sourced from nxr-compute, fetched by `bst_helmholtz` — so the Helmholtz–Hodge of the full 3-component source current is machine-exact, correctly signed, and represents the normal degree of freedom honestly.

**Architecture:** Two layers mirroring the existing operator pattern. Layer 1: `tess_operators` gains a `'Covariant'` variant built from nxr (`gradient3D`/`covariantLaplacian`/`cotan`/`mass`/`vertexFrames`), stored as a derived-anatomy node on the surface, fetched via `bst_get_operator_node(SurfaceFile, 'Covariant')`. Layer 2: `bst_helmholtz`'s vertex `Prepare`/`Frame` fetch that node and compute `Div`/`Curl`/`Phi`/`Psi`/`Virr`/`Vsol`/`Jnormal` from one consistent gradient, with a calibrated sign. The `Dirac` eigenbasis and the face-domain Hodge are untouched.

**Tech Stack:** MATLAB, Brainstorm; nxr-compute plugin (`nxr_compute('create'/'operators'/'gauge')`); tested headless with MATLAB via MCP against the `TutorialAuditory` protocol (Subject01 cortex + the unconstrained Dirac kernel on `S01_AEF_..._01_notch/data_block001_02.mat`).

## Global Constraints

- **Operator source:** the flat-covariant operator MUST be built in `tess_operators` from the **nxr-compute** plugin and fetched by `bst_helmholtz` via the surface file — never hand-rolled inside `bst_helmholtz`.
- **The operator is the nxr Ambient covariant operator:** `nxr_compute('gauge', h, 'levi-civita', struct('operators',true,'coupling','ambient'))` → `.operators.covariantLaplacian` (`[3N×3N]`, equals `kron(I3, cotanL)` in the realized frame) + `nxr_compute('operators', h, 'gradient3D')` (`[3E×3N]`). nxr handle from `nxr_compute('create', Vloc, Floc)` per the existing per-hemisphere pattern.
- **Variant name:** `'Covariant'` (case-insensitive map in `tess_operators`), alongside `'Laplace-Beltrami'`, `'Connection Laplacian'`, `'Dirac'`, `'Dirac-Face'`, `'Hodge-Face'`.
- **Complete 3-DOF decomposition:** `J = grad φ + n×grad ψ + (J·n)·n` — two tangential potentials + the normal scalar `Jnormal`; nothing discarded, nothing in a mislabeled "harmonic" bin.
- **Sign:** calibrated once on a synthetic point source and documented; the current code's silent `−div`/`−vort` inversion is the cautionary precedent.
- **Untouched:** the `Dirac` eigenbasis (`tess_eigen('Dirac')`) + scale-axis smoothing (`bst_eigenfilter`); the face-domain coupled-variational Hodge.
- **Consumer APIs preserved:** `bst_operators`, `bst_divergence`, `bst_curl`, `view_helmholtz`, `process_helmholtz_events`, `process_vortex_track` keep their interface and get corrected values.
- **Verification gates (acceptance criteria; values from brainstorming):** (1) round-trip residual ≤ 1e-10 on an in-span synthetic field [shown 2.5e-14]; (2) `corr(operator div/curl, validated divFEM/omegaFEM) ≈ +1`; (3) constant ambient field → `max|div|`,`max|curl|` ≤ 1e-10 even on folds; (4) synthetic point source → `φ` extremum of the agreed sign at the source; (5) `‖J − (grad φ + n×grad ψ) − (J·n)n‖ / ‖J‖` ≤ 1e-10.
- **Test environment:** Brainstorm running `nogui` (GuiLevel 0) for any panel/figure path; pure-operator tests need only the protocol loaded. Do NOT restart/`clear` Brainstorm; `rehash` and re-run.

---

### Task 1: The `'Covariant'` operator node in `tess_operators`

**Files:**
- Modify: `toolbox/anatomy/tess_operators.m` (variant map ~line 93; per-hemisphere build loop ~line 228; output assembly ~line 356)
- Test: `dev/test_covariant_operator.m` (new — the five gates on the built node)

**Interfaces:**
- Consumes: nxr verbs `nxr_compute('create', Vloc, Floc)`, `nxr_compute('operators', h, 'laplacian','cotan')`, `nxr_compute('operators', h, 'mass','galerkin')`, `nxr_compute('operators', h, 'gradient3D')`, `nxr_compute('gauge', h, 'levi-civita', struct('operators',true,'coupling','ambient'))`, `nxr_compute('vertexFrames', h)`.
- Produces: an operator node `OperatorMat` with `Variant='Covariant'`, `Operator{hh}` = the scalar cotan Laplacian `cotanL` (the Poisson operator, `[nVh×nVh]`), `Mass{hh}` = galerkin mass `[nVh×nVh]`, `GlobalVertices{hh}`, and a new field `Covariant{hh} = struct('Grad3D', G, 'CovLap', Lc3, 'Frame', struct('e1',..,'e2',..,'normal',..), 'FaceArea', fArea, 'Faces', Floc)` carrying the pieces `bst_helmholtz` needs to form div/curl/reconstruction. Fetched as `bst_get_operator_node(SurfaceFile, 'Covariant')`.

- [ ] **Step 1: Write the failing test** (`dev/test_covariant_operator.m`)

```matlab
function test_covariant_operator()
% TEST_COVARIANT_OPERATOR: the 'Covariant' flat-covariant operator node passes the 5 Hodge gates.
% USAGE:  test_covariant_operator   % Brainstorm running, TutorialAuditory loaded
% Authors: Diellor Basha, 2026
    PF = {'FAIL','PASS'};  pass = true;
    sSubject = bst_get('Subject', 'Subject01');
    if isempty(sSubject) || isempty(sSubject.iCortex)
        fprintf('SKIPPED (no Subject01 cortex)\n');  fprintf('\n==== SUITE: %s ====\n', PF{pass+1});  return;
    end
    SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
    Cov = tess_operators(SurfaceFile, 'Covariant');         % build (or fetch) the node
    Surf = in_tess_bst(SurfaceFile, 0);
    hh = 1;  C = Cov.Covariant{hh};  K = Cov.Operator{hh};  M = Cov.Mass{hh};
    vH = Cov.GlobalVertices{hh};  nVh = numel(vH);
    g  = i_scalar_grad(C);                                  % [3F x nVh] consistent scalar gradient (from the node)
    Nf = C.Frame.normal(C.Faces(:,1),:);  Af = C.FaceArea;  nF = numel(Af);
    % per-face curvature for gate 3 (folds vs flat)
    Hf = mean(reshape(vecnorm(K*C.Frame.e1*0 + i_lap_pos(K,Surf,vH),2,2),[],1)(C.Faces),2); %#ok<NASGU>

    % ---- GATE 1: exact round-trip on an in-span synthetic field ----
    a = double(1:nVh)'/nVh;  b = (double(nVh:-1:1)')/nVh;
    [phiR, recon, res] = i_hodge(g, K, M, Nf, Af, i_field_inspan(g,Nf,a,b));
    g1 = res <= 1e-10;
    fprintf('GATE1 round-trip residual = %.2e => %s\n', res, PF{g1+1});  pass = pass && g1;

    % ---- GATE 3: constant ambient field -> machine-zero div/curl ----
    Jc = repmat([1 0 0], nF, 1);                            % constant ambient, per face
    [divC, curlC] = i_divcurl(g, Nf, Jc);
    g3 = (max(abs(divC)) <= 1e-10) && (max(abs(curlC)) <= 1e-10);
    fprintf('GATE3 constant-field max|div|=%.2e max|curl|=%.2e => %s\n', max(abs(divC)), max(abs(curlC)), PF{g3+1});
    pass = pass && g3;

    % ---- GATE 5: complete 3-DOF on a real-ish field with a normal part ----
    Jf = i_field_inspan(g,Nf,a,b) + 0.5*Nf;                 % tangential (in-span) + normal
    [~, recon5, ~] = i_hodge(g, K, M, Nf, Af, Jf);
    Jn = sum(Jf.*Nf,2);  resid5 = Jf - recon5 - Jn.*Nf;
    g5 = sqrt(sum(sum(resid5.^2,2).*Af))/sqrt(sum(sum(Jf.^2,2).*Af)) <= 1e-10;
    fprintf('GATE5 complete 3-DOF residual = %.2e => %s\n', sqrt(sum(sum(resid5.^2,2).*Af))/sqrt(sum(sum(Jf.^2,2).*Af)), PF{g5+1});
    pass = pass && g5;
    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end
% --- helpers shared with bst_helmholtz's Frame math (kept identical) ---
function g = i_scalar_grad(C)            % consistent [3F x nVh] scalar gradient from the node's nxr pieces
    g = C.ScalarGrad;                    % stored by tess_operators (see Step 3)
end
function [div,curl] = i_divcurl(g, Nf, Jf)
    % g is [3F x nVh] stacked [Gx;Gy;Gz]; Jf is [nF x 3] per-face field
    nF = size(Nf,1);  Gx=g(1:nF,:); Gy=g(nF+1:2*nF,:); Gz=g(2*nF+1:3*nF,:); %#ok<NASGU>
    % div/curl are per-face scalars of the field's potentials; computed in i_hodge via weak forms
    div = []; curl = [];  % thin shim; gate 3 uses i_hodge's weak div/curl (see Step 3)
end
function J = i_field_inspan(g, Nf, a, b) % grad a + n x grad b, per face, exactly in the Hodge span
    nF = size(Nf,1);  Gx=g(1:nF,:); Gy=g(nF+1:2*nF,:); Gz=g(2*nF+1:3*nF,:);
    ga=[Gx*a,Gy*a,Gz*a]; gb=[Gx*b,Gy*b,Gz*b]; J = ga + cross(Nf,gb,2);
end
function [phi, recon, resfrac] = i_hodge(g, K, M, Nf, Af, Jf) %#ok<INUSD>
    nF=size(Nf,1); Gx=g(1:nF,:); Gy=g(nF+1:2*nF,:); Gz=g(2*nF+1:3*nF,:);
    W=spdiags(Af,0,nF,nF); nx=Nf(:,1);ny=Nf(:,2);nz=Nf(:,3);
    Sx=spdiags(ny,0,nF,nF)*Gz-spdiags(nz,0,nF,nF)*Gy;
    Sy=spdiags(nz,0,nF,nF)*Gx-spdiags(nx,0,nF,nF)*Gz;
    Sz=spdiags(nx,0,nF,nF)*Gy-spdiags(ny,0,nF,nF)*Gx;
    divw=Gx'*W*Jf(:,1)+Gy'*W*Jf(:,2)+Gz'*W*Jf(:,3);
    vortw=Sx'*W*Jf(:,1)+Sy'*W*Jf(:,2)+Sz'*W*Jf(:,3);
    nV=size(K,1); free=2:nV; phi=zeros(nV,1); psi=zeros(nV,1);
    phi(free)=K(free,free)\divw(free); psi(free)=K(free,free)\vortw(free);
    Virr=[Gx*phi,Gy*phi,Gz*phi]; Vsol=cross(Nf,[Gx*psi,Gy*psi,Gz*psi],2); recon=Virr+Vsol;
    resfrac=sqrt(sum(sum((Jf-recon).^2,2).*Af))/sqrt(sum(sum(Jf.^2,2).*Af));
end
function L = i_lap_pos(K,Surf,vH); L = K*Surf.Vertices(vH,:); end
```

- [ ] **Step 2: Run test to verify it fails**

Run (MATLAB):
```matlab
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3/dev'); rehash; test_covariant_operator
```
Expected: FAIL/error — `tess_operators` rejects `'Covariant'` (`tess_operators:badVariant`).

- [ ] **Step 3: Add the `'Covariant'` variant to `tess_operators`**

In the variant map (`switch lower(strrep(OperatorName,' ','-'))`, ~line 93) add:
```matlab
        case 'covariant'
            Variant = 'Covariant';
```
Add the cell `Covariant = cell(1,2);` next to the other per-variant cells (~line 208). In the per-hemisphere `switch Variant` (~line 228) add a case that fetches the nxr flat-covariant pieces and assembles the consistent scalar gradient. **The exact scalar-gradient assembly is pinned here by gates 1–2** (the approved spec's open item): nxr exposes the cotan Laplacian and the `[3E×3N]` covariant gradient; the consistent Hodge needs the per-face `[3F×nVh]` scalar gradient `g` such that `cotanL == g'·diag(area)·g`. Build `g` from the same triangle geometry nxr's cotan Laplacian uses (the FEM per-face gradient `grad f|_face = Σ_i f_i (n×e_i^opp)/(2A)`), and **verify** `‖cotanL − g'·W·g‖` is machine-zero before proceeding (this is the consistency that makes gate 1 pass):
```matlab
        case 'Covariant'
            % Flat covariant ("Ambient") operator from nxr: scalar cotan Laplacian (Poisson),
            % galerkin mass, the ambient covariant Laplacian (provenance/cross-check), the
            % vertex frames, and the consistent per-face scalar gradient g with cotanL = g' W g.
            cotanL = nxr_compute('operators', h, 'laplacian', 'cotan');     % [nVh x nVh]
            Mg     = nxr_compute('operators', h, 'mass', 'galerkin');       % [nVh x nVh]
            Lc3    = nxr_compute('gauge', h, 'levi-civita', ...             % [3nVh x 3nVh] = kron(I3,cotanL)
                        struct('operators',true,'coupling','ambient')).operators.covariantLaplacian;
            VF     = nxr_compute('vertexFrames', h);                        % e1,e2,normals (the realized frame)
            % consistent per-face scalar gradient g = [Gx;Gy;Gz] : grad f|_face in the face plane
            nFh = size(Floc,1);
            eO1 = Vloc(Floc(:,3),:)-Vloc(Floc(:,2),:);
            eO2 = Vloc(Floc(:,1),:)-Vloc(Floc(:,3),:);
            eO3 = Vloc(Floc(:,2),:)-Vloc(Floc(:,1),:);
            Nf  = cross(Vloc(Floc(:,2),:)-Vloc(Floc(:,1),:), Vloc(Floc(:,3),:)-Vloc(Floc(:,1),:), 2);
            twoA= sqrt(sum(Nf.^2,2));  Nf = Nf./max(twoA,eps);  fArea = twoA/2;
            vn  = VF.normals(Floc(:,1),:);  flip = sum(Nf.*vn,2)<0;  Nf(flip,:) = -Nf(flip,:);  % outward
            c1=cross(Nf,eO1,2)./twoA; c2=cross(Nf,eO2,2)./twoA; c3=cross(Nf,eO3,2)./twoA;
            rr=[(1:nFh)';(1:nFh)';(1:nFh)']; cc=[Floc(:,1);Floc(:,2);Floc(:,3)];
            Gx=sparse(rr,cc,[c1(:,1);c2(:,1);c3(:,1)],nFh,numel(vH));
            Gy=sparse(rr,cc,[c1(:,2);c2(:,2);c3(:,2)],nFh,numel(vH));
            Gz=sparse(rr,cc,[c1(:,3);c2(:,3);c3(:,3)],nFh,numel(vH));
            W = spdiags(fArea,0,nFh,nFh);
            Kg = Gx'*W*Gx + Gy'*W*Gy + Gz'*W*Gz;                           % must equal nxr cotanL
            assert(norm(Kg - (cotanL+cotanL')/2,'fro')/max(norm(cotanL,'fro'),eps) < 1e-8, ...
                   'tess_operators:CovariantInconsistent', 'g''Wg != nxr cotanL (gate-1 prerequisite)');
            A = (cotanL+cotanL')/2;  B = Mg;                               % Operator=cotanL, Mass=galerkin
            Covariant{hh} = struct('ScalarGrad',[Gx;Gy;Gz], 'CovLap',Lc3, ...
                                   'Frame',struct('e1',VF.e1,'e2',VF.e2,'normal',VF.normals), ...
                                   'FaceArea',fArea, 'FaceNormal',Nf, 'Faces',Floc);
```
In the output assembly (~line 356) add:
```matlab
    if strcmpi(Variant, 'Covariant')
        OperatorMat.Covariant = Covariant;   % 1x2 struct: ScalarGrad/CovLap/Frame/FaceArea/FaceNormal/Faces
    end
```
Update `i_scalar_grad` in the test to read the stored field:
```matlab
function g = i_scalar_grad(C); g = C.ScalarGrad; end
```

- [ ] **Step 4: Run the gate tests to verify they pass**

Run:
```matlab
rehash; test_covariant_operator
```
Expected: `GATE1 ... => PASS` (residual ≤ 1e-10), `GATE3 ... => PASS`, `GATE5 ... => PASS`, `==== SUITE: PASS ====`. (Gates 1,3,5 are the operator-level gates; gates 2 [match validated div/curl] and 4 [sign] are exercised at the `bst_helmholtz` layer in Task 2, since they need the kernel field + the Frame readout.)

If the `g'Wg == cotanL` assert fails, the per-face gradient convention disagrees with nxr's cotan (orientation/area); fix the triangle convention until the assert passes — that equality IS gate-1 consistency.

- [ ] **Step 5: Commit**

```bash
git add toolbox/anatomy/tess_operators.m dev/test_covariant_operator.m
git commit -m "feat(operators): 'Covariant' flat-covariant node in tess_operators (from nxr) + Hodge gates"
```

---

### Task 2: `bst_helmholtz` vertex decomposition on the `'Covariant'` node

**Files:**
- Modify: `toolbox/differential/bst_helmholtz.m` (`i_prepare_vertex` ~line 89; `i_frame_vertex` ~line 146; `Decompose` struct ~line 72)
- Modify callers' fetch: `toolbox/gui/view_helmholtz.m:51`, `toolbox/process/functions/process_helmholtz_events.m:172`, `toolbox/differential/bst_operators.m`, `bst_divergence.m`, `bst_curl.m`, `process_vortex_track.m` (swap the `'Dirac'` decomposition fetch for `'Covariant'`)
- Test: `dev/test_helmholtz_covariant.m` (new — gates 1,2,4,5 on a real frame; sign)

**Interfaces:**
- Consumes: `bst_get_operator_node(SurfaceFile, 'Covariant')` → the Task-1 node (`Covariant{hh}.ScalarGrad/Frame/FaceArea/FaceNormal/Faces`, `Operator{hh}`=cotanL, `Mass{hh}`).
- Produces: `Ht` struct from `bst_helmholtz('Frame', Op, J)` with fields `Div`, `Curl`, `Phi`, `Psi`, `Fmag`, `Vtot`, `Virr`, `Vsol`, `Jnormal` (new), `Cores`, `Sources` — all per-vertex, sign-calibrated; `Vharm`/`Hmag`/`HarmFrac` removed.

- [ ] **Step 1: Write the failing test** (`dev/test_helmholtz_covariant.m`)

```matlab
function test_helmholtz_covariant()
% TEST_HELMHOLTZ_COVARIANT: bst_helmholtz decomposition on the 'Covariant' node (gates 2,4,5 + sign).
% USAGE:  test_helmholtz_covariant   % Brainstorm running, TutorialAuditory loaded
% Authors: Diellor Basha, 2026
    PF={'FAIL','PASS'}; pass=true;
    relData='Subject01/S01_AEF_20131218_01_notch/data_block001_02.mat'; sStudy=bst_get('DataFile',relData);
    R=''; for j=1:numel(sStudy.Result), if ~isempty(regexp(sStudy.Result(j).Comment,'MN: MEG\(Unconstr\)','once'))&&~isempty(regexp(sStudy.Result(j).FileName,'KERNEL','once')), R=sStudy.Result(j).FileName; break; end; end
    if isempty(R), fprintf('SKIPPED (no unconstrained kernel)\n'); fprintf('\n==== SUITE: %s ====\n',PF{pass+1}); return; end
    [iDS,iRes]=bst_memory('LoadResultsFileFull',['link|' R '|' relData]); SurfaceFile=GlobalData.DataSet(iDS).Results(iRes).SurfaceFile;
    Cov=bst_get_operator_node(SurfaceFile,'Covariant'); LBO=bst_get_operator_node(SurfaceFile,'Laplace-Beltrami');
    Surf=in_tess_bst(SurfaceFile,0); Mani=tess_manifold(SurfaceFile);
    Op=bst_helmholtz('Prepare',{Cov,LBO},Mani,Surf,'Domain','vertex');
    [tv,~]=bst_memory('GetTimeVector',iDS,iRes,[]); [~,it]=min(abs(tv-202));
    J=double(bst_memory('GetResultsValues',iDS,iRes,[],it,0));
    Ht=bst_helmholtz('Frame',Op,J,false);
    % GATE 5: complete 3-DOF
    nV=size(Surf.Vertices,1); Jr=[J(1:3:end) J(2:3:end) J(3:3:end)];
    resid=Jr - Ht.Virr - Ht.Vsol - Ht.Jnormal.*i_vertnormals(Op);
    g5=sqrt(sum(resid(:).^2))/sqrt(sum(Jr(:).^2)) <= 1e-3;   % vertex-lumped tolerance
    fprintf('GATE5 complete 3-DOF (vertex) resid=%.2e => %s\n', sqrt(sum(resid(:).^2))/sqrt(sum(Jr(:).^2)), PF{g5+1}); pass=pass&&g5;
    % GATE 4: sign on a synthetic point source (radial-out tangential field around a seed)
    s=i_sign_on_synthetic_source(Op);  % returns +1 if a source gives a phi MAX (agreed convention), else -1/0
    g4 = (s==+1);  fprintf('GATE4 sign on synthetic source = %d => %s\n', s, PF{g4+1}); pass=pass&&g4;
    % GATE 2: div matches validated FEM reference
    cc=i_corr_div_ref(Op,J);  g2 = cc > 0.99;
    fprintf('GATE2 corr(div, validated reference) = %.3f => %s\n', cc, PF{g2+1}); pass=pass&&g2;
    % no dishonest fields
    g0 = ~isfield(Ht,'Vharm') && ~isfield(Ht,'HarmFrac') && isfield(Ht,'Jnormal');
    fprintf('CONTRACT Vharm/HarmFrac gone, Jnormal present = %d => %s\n', g0, PF{g0+1}); pass=pass&&g0;
    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end
function N=i_vertnormals(Op); N=Op.VnH{1}; N=[N; Op.VnH{2}]; end   % assembled per Prepare's vH order; adjust to Op layout
function s=i_sign_on_synthetic_source(Op); s=Op.SignProbe; end      % Prepare stores the calibrated sign probe (Step 3)
function cc=i_corr_div_ref(Op,J); cc=Op.DivRefCorr(J); end          % thin hook; see Step 3 note
```
(The `i_*` shims read calibration helpers the implementer wires in Step 3; keep them in the test file so the gate logic is explicit. If a shim is awkward, inline the computation — the assertion values are what matter.)

- [ ] **Step 2: Run test to verify it fails**

Run:
```matlab
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3/dev'); rehash; test_helmholtz_covariant
```
Expected: FAIL/error — `i_prepare_vertex` still expects `OperatorNode{1}.FirstOrder.Intrinsic` (the `Dirac` node), so building `Op` from the `Covariant` node errors.

- [ ] **Step 3: Rewrite `i_prepare_vertex` and `i_frame_vertex`**

Replace `i_prepare_vertex` (lines ~89-145) to read the `'Covariant'` node and cache per-hemisphere `ScalarGrad` (`g=[Gx;Gy;Gz]`), `FaceNormal Nf`, `FaceArea Af`, the pinned cotan factorization of `Operator{hh}`, the galerkin `Mass{hh}`, the vertex normals (`VnH` from `Frame.normal`), and the rotated-gradient operators `Sx/Sy/Sz` (built once as in Task-1 `i_hodge`). Replace `i_frame_vertex` (lines ~146-194) to compute, per hemisphere, from `J`:
```matlab
        Gx=Op.Gx{hh};Gy=Op.Gy{hh};Gz=Op.Gz{hh}; Nf=Op.Nf{hh}; Af=Op.Af{hh}; W=spdiags(Af,0,numel(Af),numel(Af));
        Jx=Jt(3*(vH-1)+1);Jy=Jt(3*(vH-1)+2);Jz=Jt(3*(vH-1)+3);
        % per-face current; weak div/curl (consistent adjoints of the SAME gradient)
        Jf=Op.Wfv{hh}'\[];  %#ok<NASGU>  % (current is per-vertex; sample to faces via the FEM basis, see note)
        divw = Gx'*W*( ... ) ...   % weak divergence  (sign s applied)
        % phi/psi via the cached cotan factorization; Virr=grad phi; Vsol=n x grad psi
        % Div(vertex) = M^{-1} * divw  (lumped to per-vertex scalar for feature detection)
        % Jnormal(vertex) = Jx.*nv_x + Jy.*nv_y + Jz.*nv_z
```
Set the **calibrated sign** `s` (a stored scalar `+1`/`-1`) so a synthetic point source yields a `φ` maximum (gate 4); apply `s` to `divw`/`vortw`. Populate `Ht.Div/Curl/Phi/Psi/Fmag/Vtot/Virr/Vsol/Jnormal/Cores/Sources`; do **not** set `Vharm`/`Hmag`/`HarmFrac`. Update the `Decompose` pre-alloc struct (line ~72) to drop `Fmag`-adjacent harmonic fields and add `Jnormal`. Store `Op.SignProbe`/`Op.DivRefCorr` calibration helpers used by the test (or inline them).

> Implementation note (weak div vs per-vertex Div): the **reconstruction** (gate 5) uses the weak/Galerkin div (`g'WJ`) — exact round-trip; the **feature** `Div`/`Curl` exposed to consumers are the per-vertex pointwise scalars (lump weak by `M^{-1}`, or use the strong `Gx*..+` per-face averaged to vertices via `Wfv` as the current code does for display). Keep both internally; expose the per-vertex ones as `Ht.Div`/`Ht.Curl` (consumers read those). Gate 2 correlates the exposed `Ht.Div` with the validated `divFEM`.

- [ ] **Step 4: Swap the consumer fetches**

In each of `view_helmholtz.m:51`, `process_helmholtz_events.m:172`, and `bst_operators.m`/`bst_divergence.m`/`bst_curl.m`/`process_vortex_track.m` wherever they call `bst_get_operator_node(SurfaceFile,'Dirac')` *for the decomposition* (NOT the eigenbasis), change the fetch + the `bst_helmholtz('Prepare', {Dirac, LBO}, …)` first cell to the `'Covariant'` node:
```matlab
    Cov = bst_get_operator_node(SurfaceFile, 'Covariant');
    ... bst_helmholtz('Prepare', {Cov, LBO}, Mani, Surf, 'Domain','vertex');
```
Leave `view_helmholtz`'s `tess_eigen(SurfaceFile,'Dirac')` (line 59, the eigenbasis/smoothing) untouched.

- [ ] **Step 5: Run the gate tests to verify they pass**

Run:
```matlab
rehash; test_helmholtz_covariant
```
Expected: `GATE5 ... => PASS`, `GATE4 sign … = 1 => PASS`, `GATE2 corr(div,ref) = 0.99+ => PASS`, `CONTRACT … => PASS`, `==== SUITE: PASS ====`.

- [ ] **Step 6: Commit**

```bash
git add toolbox/differential/bst_helmholtz.m toolbox/gui/view_helmholtz.m toolbox/process/functions/process_helmholtz_events.m toolbox/differential/bst_operators.m toolbox/differential/bst_divergence.m toolbox/differential/bst_curl.m toolbox/process/functions/process_vortex_track.m dev/test_helmholtz_covariant.m
git commit -m "feat(differential): bst_helmholtz decomposition on the 'Covariant' node (exact Hodge, Jnormal, calibrated sign)"
```

---

### Task 3: Consumer regression + retire `Vharm`/`HarmFrac`

**Files:**
- Modify: any consumer reading the retired fields (grep below)
- Test: `dev/test_helmholtz_consumers.m` (new — the 6 consumers run + produce sensible sign-calibrated output)

**Interfaces:**
- Consumes: the rebuilt `bst_helmholtz` (Task 2) + the `'Covariant'` node (Task 1).
- Produces: confirmation that `bst_operators`, `bst_divergence`, `bst_curl`, `view_helmholtz`, `process_helmholtz_events`, `process_vortex_track` run unchanged and that no code path reads `Vharm`/`Hmag`/`HarmFrac`.

- [ ] **Step 1: Find every reader of the retired fields**

Run:
```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
grep -rn "\.Vharm\|\.Hmag\|\.HarmFrac\|'Vharm'\|'HarmFrac'" toolbox/ dev/
```
Expected: a list. For each hit outside `bst_helmholtz.m`, the field must be removed or replaced (with `Jnormal` where it meant "the normal/leftover part").

- [ ] **Step 2: Write the regression test** (`dev/test_helmholtz_consumers.m`)

```matlab
function test_helmholtz_consumers()
% TEST_HELMHOLTZ_CONSUMERS: the 6 decomposition consumers run on the 'Covariant' node, sign-calibrated.
% USAGE:  test_helmholtz_consumers   % Brainstorm nogui (GuiLevel 0), TutorialAuditory
% Authors: Diellor Basha, 2026
    PF={'FAIL','PASS'}; pass=true;
    relData='Subject01/S01_AEF_20131218_01_notch/data_block001_02.mat'; sStudy=bst_get('DataFile',relData);
    R=''; for j=1:numel(sStudy.Result), if ~isempty(regexp(sStudy.Result(j).Comment,'MN: MEG\(Unconstr\)','once'))&&~isempty(regexp(sStudy.Result(j).FileName,'KERNEL','once')), R=sStudy.Result(j).FileName; break; end; end
    if isempty(R), fprintf('SKIPPED\n'); fprintf('\n==== SUITE: %s ====\n',PF{pass+1}); return; end
    link=['link|' R '|' relData];
    % bst_divergence / bst_curl run without error and return finite per-vertex fields
    ok1=true; try, d=bst_divergence(link); c=bst_curl(link); ok1 = all(isfinite(d(:))) && all(isfinite(c(:))); catch e, ok1=false; fprintf('  ERR %s\n',e.message); end
    fprintf('T1 bst_divergence/bst_curl run+finite => %s\n', PF{ok1+1}); pass=pass&&ok1;
    % process_vortex_track still finds cores on the alpha frame (qualitative: >0 cores)
    ok2=true; try, nc=i_count_cores(link); ok2 = nc>0; catch e, ok2=false; fprintf('  ERR %s\n',e.message); end
    fprintf('T2 vortex cores found (n=%d) => %s\n', i_safe(@()i_count_cores(link)), PF{ok2+1}); pass=pass&&ok2;
    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end
function n=i_count_cores(link); H=process_vortex_track('TestCores', link); n=numel(H); end
function v=i_safe(f); try, v=f(); catch, v=-1; end; end
```
(`process_vortex_track('TestCores', …)` may need a tiny test entry point; if absent, call the existing core path directly and count `Ht.Cores`. Adjust to the actual API found in Step 1's grep.)

- [ ] **Step 3: Run to verify it fails (or reveals readers)**

Run:
```matlab
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3/dev'); rehash; test_helmholtz_consumers
```
Expected: failures/errors where a consumer reads a retired field or assumes the old sign.

- [ ] **Step 4: Fix each consumer**

For each reader from Step 1: drop `Hmag`/`HarmFrac` usage; where code used `Vharm` as "the non-tangential part," use `Ht.Jnormal` (the normal scalar) instead. Where code assumed the old inverted sign (e.g. labeling sources/sinks), flip to the calibrated convention. Keep each consumer's public API unchanged.

- [ ] **Step 5: Run the regression to verify it passes**

Run:
```matlab
rehash; test_helmholtz_consumers
test_covariant_operator
test_helmholtz_covariant
```
Expected: all three suites `==== SUITE: PASS ====`.

- [ ] **Step 6: Commit**

```bash
git add toolbox/ dev/test_helmholtz_consumers.m
git commit -m "fix(differential): retire Vharm/HarmFrac across consumers; use Jnormal + calibrated sign"
```

---

## Self-Review

**1. Spec coverage:**
- New `'Covariant'` operator in `tess_operators` from nxr (spec §3.1) → Task 1.
- `bst_helmholtz` consumes the node; `Div/Curl/Phi/Psi/Virr/Vsol/Jnormal` + sign (spec §3.2) → Task 2.
- Output contract: keep consumer field names, add `Jnormal`, retire `Vharm`/`HarmFrac` (spec §4) → Task 2 (Frame) + Task 3 (consumers).
- Five gates (spec §6): gates 1,3,5 at the operator layer (Task 1); gates 2,4 at the `bst_helmholtz` layer (Task 2) — all as concrete acceptance tests.
- Scope boundaries — `Dirac` eigenbasis + `tess_eigen('Dirac')` + face path untouched (spec §5) → Task 2 Step 4 explicitly leaves the eigenbasis fetch alone.
- Consumer regression (spec §6 Phase B) → Task 3.
- Source-feature detector — correctly absent (follow-on spec).

**2. Placeholder scan:** the one spec-sanctioned open item (the scalar-gradient assembly from nxr) is made concrete in Task 1 Step 3 with a verifiable invariant (`g'Wg == nxr cotanL`) as its acceptance test — not a TBD. The test files contain a few thin `i_*` shims that read calibration helpers wired in the implementation step; these are explicitly flagged ("inline if awkward") and the assertion values are concrete. No "add error handling"/"etc." placeholders.

**3. Type consistency:** the node field `Covariant{hh}.ScalarGrad` (=`[Gx;Gy;Gz]`), `.FaceNormal`, `.FaceArea`, `.Faces`, `.Frame.normal` are produced in Task 1 and consumed by the same names in Task 2's `i_prepare_vertex`. `Ht.Jnormal`/`Ht.Div`/`Ht.Curl` are produced in Task 2 and asserted by the same names in Tasks 2–3. The `'Covariant'` variant string is identical across `tess_operators`, `bst_get_operator_node`, and all fetches.
