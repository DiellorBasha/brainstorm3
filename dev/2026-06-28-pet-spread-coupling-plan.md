# Longitudinal Aβ/tau Spread Coupling Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a synthetic coupled-spread generative model (`pet_spread_simulate`) and a reaction-diffusion coupling inverter (`pet_spread_invert`) that recovers the Aβ→tau coupling κ from a longitudinal surface-PET series, validated by a κ-recovery curve.

**Architecture:** Coupled Fisher-KPP reaction-diffusion on the LH cortical manifold (LBO `M`,`K` from `tess_operators`; implicit-Euler diffusion + explicit logistic reaction). The PDE is linear in its parameters, so inversion is a mass-weighted linear regression of the finite-difference time-derivative onto `[τ(1−τ), a·τ(1−τ), −Kτ]`, restricted to the active (non-saturated) range; κ = β₂/β₁.

**Tech Stack:** MATLAB (R2023b), Brainstorm. Reuses `tess_operators` (Laplace-Beltrami `Mass`/`Operator`/`GlobalVertices`), `in_tess_bst`, `pet_suvr` (pseudo-longitudinal check). Testing = synthetic-recovery validation via the MATLAB MCP tools (`check_matlab_code`, `evaluate_matlab_code`). Prototype reference: `dev/benchmarks/proto_pet_spread.m`.

## Global Constraints
- MATLAB only; base MATLAB + Brainstorm (no Image Processing / Statistics toolbox — use `\`, `corrcoef`).
- Never call `clear`; `rehash path` to reload edited `.m`.
- LBO operators are per-hemisphere: `LBO.Mass{1}`/`LBO.Operator{1}`/`LBO.GlobalVertices{1}` = LH (`nL=10242` on ico5). Work in LH-local indices.
- **Sign/discretization consistency is mandatory:** the simulator's implicit step is `(M+dt·D·K)·x_{n+1} = M·(x_n + dt·reaction_n)`. The inverter MUST mirror it: `M·(x_{n+1}−x_n)/dt = β1·M·[reaction basis]_n + β_D·(−K·x_{n+1})`, so recovered `D = β_D > 0`.
- Functions in `toolbox/anatomy/`; validation script in `dev/benchmarks/`. Commit messages end with the two project trailer lines (Co-Authored-By + Claude-Session).
- Each task: write the validation snippet, run it (fail), implement, `check_matlab_code` (clean), run validation (pass), commit.

---

### Task 1: `pet_spread_simulate` — coupled generative model

**Files:**
- Create: `toolbox/anatomy/pet_spread_simulate.m`

**Interfaces:**
- Produces: `[a, tau, info] = pet_spread_simulate(SurfaceFile, Opts)`.
  - `Opts` fields (defaults): `nT=24, dt=1, Da=8e-5, Dt=8e-5, ra=0.6, rt=0.45, kappa=3, seedAmp=0.6, seedA=[], seedT=[]` (seeds are LH-local indices; if empty, derived from Desikan `precuneus L`/`entorhinal L`).
  - `a`, `tau`: `[nL × nT]` field time-series (LH).
  - `info`: `.gv` (LH global vertex indices), `.M`, `.K` (LH LBO operators), `.seedA`, `.seedT`, `.Opts`.

- [ ] **Step 1: Write the failing validation snippet**

```matlab
rehash path;
[sS,~]=bst_get('Subject','sub-MTL0002');
wf=sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},'cortex_white_low\.mat$','once')),1)).FileName;
[a3,t3,info]=pet_spread_simulate(wf, struct('kappa',3));
[a0,t0]=pet_spread_simulate(wf, struct('kappa',0));
assert(mean(a3(:,end)>0.5)>0.3, 'Abeta failed to spread');
assert(mean(t3(:,end)>0.5) > 3*mean(t0(:,end)>0.5), 'coupling did not accelerate tau');
fprintf('PASS: Abeta cover %.0f%%; tau cover kappa3=%.0f%% vs kappa0=%.0f%%\n', ...
    100*mean(a3(:,end)>0.5), 100*mean(t3(:,end)>0.5), 100*mean(t0(:,end)>0.5));
```

- [ ] **Step 2: Run it — expect FAIL** (`Unrecognized function 'pet_spread_simulate'`).

- [ ] **Step 3: Implement** `toolbox/anatomy/pet_spread_simulate.m`:

```matlab
function [a, tau, info] = pet_spread_simulate(SurfaceFile, Opts)
% PET_SPREAD_SIMULATE: synthetic coupled Abeta/tau spread on the cortical manifold (LH).
% Two Fisher-KPP reaction-diffusion fields, LBO diffusion (implicit Euler) + logistic reaction;
% Abeta seeds at precuneus, tau at entorhinal with growth gated by local Abeta (the cascade).
%
% USAGE: [a, tau, info] = pet_spread_simulate(SurfaceFile, Opts)
%
% Author: Diellor Basha, 2026
    D=struct('nT',24,'dt',1,'Da',8e-5,'Dt',8e-5,'ra',0.6,'rt',0.45,'kappa',3,'seedAmp',0.6,'seedA',[],'seedT',[]);
    if nargin<2, Opts=struct(); end
    fn=fieldnames(D); for i=1:numel(fn), if ~isfield(Opts,fn{i}), Opts.(fn{i})=D.(fn{i}); end; end

    sW=in_tess_bst(SurfaceFile); nV=size(sW.Vertices,1);
    LBO=tess_operators(SurfaceFile,'Laplace-Beltrami');
    gv=double(LBO.GlobalVertices{1}(:)); M=LBO.Mass{1}; K=LBO.Operator{1}; nL=numel(gv);
    g2l=zeros(nV,1); g2l(gv)=1:nL;
    sa=Opts.seedA; if isempty(sa), sa=local_seed(sW,g2l,'precuneus L'); end
    st=Opts.seedT; if isempty(st), st=local_seed(sW,g2l,'entorhinal L'); end

    dA=decomposition(M + Opts.dt*Opts.Da*K);
    dT=decomposition(M + Opts.dt*Opts.Dt*K);
    a=zeros(nL,1); a(sa)=Opts.seedAmp; tau=zeros(nL,1); tau(st)=Opts.seedAmp;
    A=zeros(nL,Opts.nT); T=zeros(nL,Opts.nT);
    for n=1:Opts.nT
        ra = Opts.ra * a   .* (1-a);
        rt = Opts.rt * (1+Opts.kappa*a) .* tau .* (1-tau);
        a   = dA\(M*(a   + Opts.dt*ra));  a  =min(max(a,0),1);
        tau = dT\(M*(tau + Opts.dt*rt));  tau=min(max(tau,0),1);
        A(:,n)=a; T(:,n)=tau;
    end
    a=A; tau=T;
    info=struct('gv',gv,'M',M,'K',K,'seedA',sa,'seedT',st,'Opts',Opts);
end

function loc=local_seed(sW, g2l, label)
    ai=find(~cellfun('isempty',regexp({sW.Atlas.Name},'Desikan','once')),1);
    si=find(strcmp({sW.Atlas(ai).Scouts.Label}, label),1);
    gvs=sW.Atlas(ai).Scouts(si).Vertices(:);
    s=sW.Atlas(ai).Scouts(si).Seed; if isempty(s), s=gvs(round(numel(gvs)/2)); end
    loc=g2l(s); if loc==0, loc=g2l(gvs(find(g2l(gvs)>0,1))); end
end
```

- [ ] **Step 4: Lint** — expect clean.
- [ ] **Step 5: Run the Step-1 snippet — expect PASS** (Abeta spreads; tau κ3 ≫ κ0).
- [ ] **Step 6: Commit**

```bash
git add toolbox/anatomy/pet_spread_simulate.m
git commit -m "feat(pet): pet_spread_simulate - coupled Abeta/tau generative model (Fisher-KPP on LBO)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01BTbrKJ8SkHiQg9hX2gkaPi"
```

---

### Task 2: `pet_spread_invert` — reaction-diffusion coupling inversion

**Files:**
- Create: `toolbox/anatomy/pet_spread_invert.m`

**Interfaces:**
- Consumes: `a`, `tau` `[nL × nT]` and `info.M`,`info.K` (or recomputes from `SurfaceFile`).
- Produces: `[est, dinfo] = pet_spread_invert(SurfaceFile, a, tau, dt, Opts)`.
  - `est`: struct `.ra .Da .rt .Dt .kappa`.
  - `Opts`: `.active=[0.05 0.95]` (fit only where the field is in this range), `.ridge=0`, `.diffAt='np1'` (diffusion term evaluated at step n+1 to match the implicit simulator).
  - `dinfo`: `.condX` (design-matrix condition number), `.resid` (relative residual), `.nSamp`.

- [ ] **Step 1: Write the failing validation snippet** (recover known κ from a simulated series)

```matlab
rehash path;
[sS,~]=bst_get('Subject','sub-MTL0002');
wf=sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},'cortex_white_low\.mat$','once')),1)).FileName;
[a,tau]=pet_spread_simulate(wf, struct('kappa',3,'rt',0.45,'Dt',8e-5));
est=pet_spread_invert(wf, a, tau, 1);
fprintf('recovered: kappa=%.2f (true 3) | rt=%.3f (0.45) | Dt=%.2e (8e-5) | Dt>0:%d\n', est.kappa, est.rt, est.Dt, est.Dt>0);
assert(abs(est.kappa-3)/3 < 0.2 && est.Dt>0, 'kappa/Dt not recovered within tolerance');
fprintf('PASS\n');
```

- [ ] **Step 2: Run it — expect FAIL** (`Unrecognized function 'pet_spread_invert'`).

- [ ] **Step 3: Implement** `toolbox/anatomy/pet_spread_invert.m`:

```matlab
function [est, dinfo] = pet_spread_invert(SurfaceFile, a, tau, dt, Opts)
% PET_SPREAD_INVERT: recover the reaction-diffusion parameters {ra,Da,rt,Dt,kappa} from a
% longitudinal Abeta/tau surface series by mass-weighted linear regression (the coupled Fisher-KPP
% is linear in its parameters). Mirrors the implicit-Euler discretization of pet_spread_simulate:
%   M*(x_{n+1}-x_n)/dt = M*[reaction basis]_n + beta_D*(-K*x_{n+1})
% Abeta:  M*da/dt = ra*M*(a(1-a)) - Da*K*a_{n+1}
% tau:    M*dtau/dt = rt*M*(tau(1-tau)) + (kappa*rt)*M*(a*tau(1-tau)) - Dt*K*tau_{n+1}
% kappa = beta2/beta1.
%
% USAGE: [est, dinfo] = pet_spread_invert(SurfaceFile, a, tau, dt, Opts)
%
% Author: Diellor Basha, 2026
    D=struct('active',[0.05 0.95],'ridge',0,'diffAt','np1');
    if nargin<5, Opts=struct(); end
    fn=fieldnames(D); for i=1:numel(fn), if ~isfield(Opts,fn{i}), Opts.(fn{i})=D.(fn{i}); end; end
    if nargin<4||isempty(dt), dt=1; end
    LBO=tess_operators(SurfaceFile,'Laplace-Beltrami'); M=LBO.Mass{1}; K=LBO.Operator{1};

    % ---- Abeta fit: [ra, Da] ----
    [ya,Xa]=local_assemble(a, [], M, K, dt, Opts, false);
    ba=local_solve(Xa, ya, Opts.ridge);
    % ---- tau fit: [rt, kappa*rt, Dt] (coupling basis uses a at step n) ----
    [yt,Xt]=local_assemble(tau, a, M, K, dt, Opts, true);
    bt=local_solve(Xt, yt, Opts.ridge);

    est=struct('ra',ba(1),'Da',ba(2),'rt',bt(1),'Dt',bt(3),'kappa',bt(2)/bt(1));
    dinfo=struct('condX',cond(Xt'*Xt),'resid',norm(Xt*bt-yt)/max(norm(yt),eps),'nSamp',numel(yt));
end

function [y,X]=local_assemble(f, g, M, K, dt, Opts, coupled)
    nT=size(f,2); y=[]; X=[];
    for n=1:nT-1
        fn=f(:,n); fnp=f(:,n+1);
        act = fn>Opts.active(1) & fn<Opts.active(2);     % fit only in the active (non-saturated) band
        if ~any(act), continue; end
        yi  = M*((fnp-fn)/dt);
        x1  = M*(fn.*(1-fn));
        if strcmpi(Opts.diffAt,'np1'), xd=-K*fnp; else, xd=-K*fn; end
        if coupled
            x2 = M*(g(:,n).*fn.*(1-fn));
            Xi = [x1(act), x2(act), xd(act)];
        else
            Xi = [x1(act), xd(act)];
        end
        y=[y; yi(act)]; X=[X; Xi]; %#ok<AGROW>
    end
end

function b=local_solve(X, y, ridge)
    if ridge>0, b=(X'*X + ridge*eye(size(X,2))) \ (X'*y); else, b=X\y; end
end
```

- [ ] **Step 4: Lint** — expect clean.
- [ ] **Step 5: Run the Step-1 snippet — expect PASS** (κ≈3 within 20%, Dt>0). If κ is off, verify `diffAt='np1'` and the active-band mask (saturated samples where `(1−τ)→0` must be excluded).
- [ ] **Step 6: Commit**

```bash
git add toolbox/anatomy/pet_spread_invert.m
git commit -m "feat(pet): pet_spread_invert - reaction-diffusion coupling inversion (recover kappa)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01BTbrKJ8SkHiQg9hX2gkaPi"
```

---

### Task 3: validation — κ-recovery curve + noise/sparse-time robustness

**Files:**
- Create: `dev/benchmarks/validate_pet_spread_coupling.m`

**Interfaces:**
- Consumes: `pet_spread_simulate`, `pet_spread_invert`.
- Produces: `R = validate_pet_spread_coupling()` — recovery table + a figure `dev/benchmarks/pet_spread_kappa_recovery.png` (κ_est vs κ_true; recovery under noise; recovery under sparse sampling).

- [ ] **Step 1: Write the script**

```matlab
function R = validate_pet_spread_coupling()
% VALIDATE_PET_SPREAD_COUPLING: kappa-recovery curve + noise/sparse-time robustness for the
% reaction-diffusion coupling inverter, on synthetic ground truth.
% Author: Diellor Basha, 2026
    [sS,~]=bst_get('Subject','sub-MTL0002');
    wf=sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},'cortex_white_low\.mat$','once')),1)).FileName;
    here=bst_fileparts(mfilename('fullpath'));
    kts=[0 1 3 5]; kEst=zeros(size(kts)); kNoise=zeros(size(kts)); kSparse=zeros(size(kts));
    for i=1:numel(kts)
        [a,tau]=pet_spread_simulate(wf, struct('kappa',kts(i)));
        kEst(i)   = getfield(pet_spread_invert(wf,a,tau,1),'kappa'); %#ok<GFLD>
        an=min(max(a+0.03*sin((1:numel(a))'*7.1).*reshape(ones(size(a)),[],1),0),1); %#ok<NASGU>
        aN=min(max(a+0.03*randnlike(a),0),1); tN=min(max(tau+0.03*randnlike(tau),0),1);
        kNoise(i) = getfield(pet_spread_invert(wf,aN,tN,1),'kappa'); %#ok<GFLD>
        idx=1:4:size(a,2);                                  % sparse: every 4th timepoint
        kSparse(i)= getfield(pet_spread_invert(wf,a(:,idx),tau(:,idx),4),'kappa'); %#ok<GFLD>
    end
    fprintf('\nkappa recovery:\n  true:   %s\n  full:   %s\n  +noise: %s\n  sparse: %s\n', ...
        mat2str(kts), mat2str(round(kEst,2)), mat2str(round(kNoise,2)), mat2str(round(kSparse,2)));
    f=figure('Visible','off','Position',[40 40 560 480]);
    plot(kts,kts,'k--'); hold on;
    plot(kts,kEst,'o-','LineWidth',1.5); plot(kts,kNoise,'s-'); plot(kts,kSparse,'^-');
    xlabel('true \kappa'); ylabel('recovered \kappa'); grid on; axis equal;
    legend({'identity','full series','+3% noise','sparse (1/4 t)'},'Location','northwest');
    title('A\beta\rightarrowtau coupling recovery'); 
    png=fullfile(here,'pet_spread_kappa_recovery.png'); print(f,png,'-dpng','-r120'); close(f);
    fprintf('  figure -> %s\n', png);
    R=struct('kts',kts,'kEst',kEst,'kNoise',kNoise,'kSparse',kSparse);
end

function r=randnlike(x)
    % deterministic pseudo-noise (no Statistics tb / no global rng dependence): hashed sinusoid
    n=numel(x); idx=(1:n)'; r=reshape(sin(idx*12.9898+1.0)*43758.5453,[],1); r=r-floor(r); r=(r-0.5)*2;
    r=reshape(r,size(x));
end
```

- [ ] **Step 2: Lint** — expect clean (the `randnlike` helper avoids `randn`/Date/Stats dependencies).
- [ ] **Step 3: Run it**

```matlab
rehash path; R=validate_pet_spread_coupling();
```
Expected: `full` recovery tracks `true` closely (κ_est ≈ [0 1 3 5]); `+noise` and `sparse` track with larger error but monotone. Sanity gate: full-series κ_est within ±20% of true for κ∈{1,3,5} and ≈0 for κ=0; a saved recovery-curve figure.

- [ ] **Step 4: Commit** (script + figure)

```bash
git add dev/benchmarks/validate_pet_spread_coupling.m dev/benchmarks/pet_spread_kappa_recovery.png
git commit -m "bench(pet): kappa-recovery curve + noise/sparse-time robustness for spread inversion

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01BTbrKJ8SkHiQg9hX2gkaPi"
```

---

### Task 4: pseudo-longitudinal cross-check on the real cohort

**Files:**
- Modify: `dev/benchmarks/validate_pet_spread_coupling.m` (add a second entry point)

**Interfaces:**
- Consumes: `pet_spread_invert`, `pet_suvr`, the real cohort surface SUVR sampling (same mid-surface sampling as `analyze_pet_epicenters`).
- Produces: `Rp = pseudo_longitudinal_check(subjects)` — orders subjects by global SUVR severity to form a pseudo-time series of `[a, tau]` LH maps, runs the inverter, and reports the recovered dynamics (honestly labeled pseudo-time). Figure `dev/benchmarks/pet_spread_pseudolong.png`.

- [ ] **Step 1: Write the failing snippet** (the new function does not exist yet)

```matlab
rehash path;
Rp=pseudo_longitudinal_check({'sub-MTL0166','sub-MTL0079','sub-MTL0018','sub-MTL0054','sub-MTL0268','sub-MTL0284','sub-MTL0311'});
fprintf('pseudo-long kappa=%.2f (sign/finite only; pseudo-time)\n', Rp.kappa);
assert(isfinite(Rp.kappa), 'pseudo-longitudinal inversion returned non-finite kappa');
fprintf('PASS\n');
```

- [ ] **Step 2: Run it — expect FAIL** (`Unrecognized function 'pseudo_longitudinal_check'`).

- [ ] **Step 3: Implement** — append to `dev/benchmarks/validate_pet_spread_coupling.m`:

```matlab
function Rp = pseudo_longitudinal_check(subjects)
% Order subjects by global cortical SUVR severity -> pseudo-time series of LH Abeta & tau surface
% maps -> run the coupling inverter. PSEUDO-TIME (cross-sectional), a sanity cross-check only.
    sDef=bst_get('Subject',0);
    wfT=sDef.Surface(find(~cellfun('isempty',regexp({sDef.Surface.FileName},'cortex_white_low\.mat$','once')),1)).FileName;
    LBO=tess_operators(wfT,'Laplace-Beltrami'); gv=double(LBO.GlobalVertices{1}(:)); nL=numel(gv);
    A=[]; T=[]; sev=[];
    for s=1:numel(subjects)
        m=local_surf_suvr(subjects{s},'PET 18FNAV4694_mean_pvc', wfT, gv);
        n=local_surf_suvr(subjects{s},'PET 18Fflortaucipir_mean_pvc', wfT, gv);
        if isempty(m)||isempty(n), continue; end
        A=[A, m]; T=[T, n]; sev=[sev, mean(m,'omitnan')]; %#ok<AGROW>
    end
    [~,ord]=sort(sev); A=A(:,ord); T=T(:,ord);             % pseudo-time = increasing severity
    rng=max(A(:)); if rng>0, A=A/rng; end; rng=max(T(:)); if rng>0, T=T/rng; end   % normalize to [0,1]
    A=min(max(A,0),1); T=min(max(T,0),1);
    est=pet_spread_invert(wfT, A, T, 1);
    fprintf('\npseudo-longitudinal (%d subjects, severity-ordered): kappa=%.2f rt=%.3f Dt=%.2e\n', size(A,2), est.kappa, est.rt, est.Dt);
    Rp=struct('kappa',est.kappa,'est',est,'nSub',size(A,2));
end

function m=local_surf_suvr(subj, petC, wfT, gv)
    m=[]; [sS,~]=bst_get('Subject',subj); cmt={sS.Anatomy.Comment};
    if ~any(strcmp(cmt,petC)), return; end
    af=@(c) sS.Anatomy(find(strcmp(cmt,c),1)).FileName;
    sf=@(rx) sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},rx,'once')),1)).FileName;
    sW=in_tess_bst(sf('cortex_white_low\.mat$')); sP=in_tess_bst(sf('cortex_pial_low\.mat$'));
    sSuvr=pet_suvr(in_mri_bst(af(petC)), in_mri_bst(af('ASEG')));
    Vmid=0.5*(sW.Vertices+sP.Vertices); vox=cs_convert(sSuvr,'scs','voxel',Vmid);
    cs=size(sSuvr.Cube); ok=all(vox>=1,2)&vox(:,1)<=cs(1)&vox(:,2)<=cs(2)&vox(:,3)<=cs(3);
    v=nan(size(Vmid,1),1); v(ok)=interpn(double(sSuvr.Cube),vox(ok,1),vox(ok,2),vox(ok,3),'linear',NaN);
    m=v(gv); m(~isfinite(m))=0;   % LH vertices on the subject surface align with the template LH index set
end
```

- [ ] **Step 4: Lint** — expect clean.
- [ ] **Step 5: Run the Step-1 snippet — expect PASS** (finite κ). Report the value; interpret cautiously (pseudo-time, 7 subjects). If subject-vs-template LH vertex alignment is invalid (different meshes), note it and restrict to subjects sharing the template tessellation, or map via `tess_interp_tess2tess` (as `analyze_pet_epicenters` does).
- [ ] **Step 6: Commit**

```bash
git add dev/benchmarks/validate_pet_spread_coupling.m
git commit -m "bench(pet): pseudo-longitudinal coupling cross-check on the real cohort (severity-ordered)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01BTbrKJ8SkHiQg9hX2gkaPi"
```

---

## Self-Review

**Spec coverage:** generative model → Task 1; reaction-diffusion inversion (linear regression, κ=β₂/β₁) → Task 2; κ-recovery curve + noise + sparse-time robustness → Task 3; pseudo-longitudinal cross-check → Task 4. All design components covered.

**Placeholder scan:** all steps carry real MATLAB + assertions; no TBD/TODO. (`randnlike` is a deterministic pseudo-noise helper — chosen because `Math.random`/`randn` reproducibility and Statistics-toolbox limits; documented inline.)

**Type consistency:** `pet_spread_simulate` returns `[a,tau,info]` with `info.M/.K/.gv`; `pet_spread_invert(SurfaceFile,a,tau,dt,Opts)` returns `est{.ra .Da .rt .Dt .kappa}` — consumed unchanged in Tasks 3–4. The simulator's implicit discretization and the inverter's `diffAt='np1'` + active-band mask are matched (Global Constraints) so κ recovers.

**Known caveats to watch during execution:**
1. **Active-band masking is essential** — once a field saturates (τ→1), `τ(1−τ)→0` makes the reaction basis vanish; those samples carry no parameter information and add noise. Fit only `0.05<f<0.95`.
2. **Subject↔template LH alignment (Task 4)** — `local_surf_suvr` assumes subject and template share the LH vertex index set (true for ico5 registered subjects). If not, map via `tess_interp_tess2tess`.
3. **Pseudo-time is not real time** — Task 4 is a sanity cross-check; do not over-interpret the recovered κ.
