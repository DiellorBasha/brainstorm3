# Dynamics panel — GUI workflow, live validation & test targets (2026-07-02)

Hands-on walkthrough of every Dynamics-panel capability on the dev test source, **each step
driven live and validated this session**. Findings (✅ works / 🐞 bug / ⚠️ trap) feed the test list
in §5.

**Test source (same as all Dirac-Connectome dev):** protocol `preventad`, subject `sub-MTL0002`,
result `results_DiracEig_KERNEL_260625_2215.mat` — an **unconstrained (nComponents=3) Dirac-dSPM**
source (carries `ImagingKernelMode` + `DiracEigenFile`) over a 62 s / 74381-sample rest recording.

---

## 0. Open the panel

**GUI:** in the database tree, expand the source study, right-click the **Dirac source result**
(Comment contains "dirac") → **`Dynamics (atoms + differential maps)`**. (Menu is gated to
`HeadModelType=='surface'` results whose Comment contains "dirac" — `tree_callbacks.m:1882`.)

This opens two things: a **3-D cortex figure** (blank "Design" surface) and the docked **Dynamics
panel**. The cortex figure carries a `DynamicsOverlay` appdata with the source link + the surface
**Covariant** operator (`D.Cov`, for Helmholtz); the panel state lives in `getappdata(0,'DynamicsTarget')`.

**Programmatic (used for testing):**
```matlab
srcFile = ['link|sub-MTL0002/@raw.../results_DiracEig_KERNEL_260625_2215.mat|' ...
           'sub-MTL0002/@raw.../data_0raw_...meg_resample_notch_high.mat'];
hFig = view_dynamics('FromResult', srcFile);
```
✅ Opens; `DynamicsOverlay.srcResult`, `.Cov` populated.

---

## 1. Panel anatomy (what maps to what)

**Menu bar:** `File` · `Atoms` (▸ `Set operator` [5 radio], ▸ `Show phases`, ▸ `Sort groups`) · `✕` (close session).

**Center:** the **filterbank atom list** (coloured-dot rows; Delete key removes the selected atom).

**Vertical toolbar (top→bottom):**
| Button (tooltip) | Callback | Purpose |
|---|---|---|
| Create atom | `OnCreateAtom` | add a default (diffusion) atom |
| Save filterbank | `OnSaveFilterbank` | write the atom table to disk (prompts path) |
| Localize (toggle) | `OnLocalize` | click a cortex vertex to re-seed the selected atom |
| Threshold | `OnThresholdMenu` | level-set threshold for Scout+Event export (popup) |
| Apply (toggle) | `OnApply` | ON = Preview (filter the real source); OFF = Design (impulse) |
| Analyze | `OnAnalyzeWindow` | frame → **full-time coeff-space scalogram** |
| Localize bands | `OnLocalizeBands` | localize each frame band → a **separate** dynamics table |
| **Helmholtz (filtered)** | `OnHelmholtzFiltered` | Hodge decomp of the **filtered** field at the cursor |
| Show phases (toggle) | `OnShowAll` | show/hide all atom phases |
| Measure | `OnMeasureMenu` | **raw**-source differential map: None/Div/Curl/Potential/Stream |

**Atom section (SOUTH):** `Filter:` dropdown (18 kernels) · per-kernel param sliders · `Direction:`
(Dirac: quaternion preset combo; tangent: angle spinner; scalar: hidden) · `Measure:` Amplitude/dSPM
radio (Dirac-dSPM only) · info readout.

**Frame section:** `A / B / B:A` frame-bound readout · `N` spinner · **`Design tight frame`** · `Show coverage`.

**Filter dropdown — static vs dynamic** (only **static** kernels drive Apply / Helmholtz / the frame;
dynamic ones are time/joint kernels and are rejected there):
- **Static:** Difference of Gaussians, Flat, **Heat/diffusion**, Ideal, Inverse heat, Logarithmic,
  Matérn/SPDE, **Mexican hat**, Power law, Tikhonov.
- **Dynamic (ts/js):** **Diffusion (the DEFAULT)**, Damped wave, Klein-Gordon, Wave, Gabor, Resonator,
  Spatiotemporal 1/f, Traveling wave.

⚠️ **Trap:** a fresh atom defaults to **Diffusion (dynamic)** → Apply/Helmholtz bail with
"…is dynamic; use a static kernel". Switch the Filter to **Heat** (or another static) first.

---

## 2. Operator gating (validated live on this unconstrained source)

`Atoms ▸ Set operator` offers 5 operators; the compat gate keys off `nComponents`:

| Operator | Gate (unconstrained) | Live result |
|---|---|---|
| Geometric (Laplace-Beltrami) | allowed | ✅ realises `[10242×4800]` |
| Connectomic (LB-Connectome) | allowed | 🐞 **"not realisable"** — `bst_eigenfilter: unsupported eigen variant 'LB-Connectome'` |
| Tangent (Connection Laplacian) | blocked | ⚠️ "not realisable" — needs a not-yet-persisted per-vertex tangent frame |
| Dirac | allowed (default) | ✅ realises |
| Dirac (connectome) | allowed | ✅ realises (fiber-spread) |

"not realisable for this atom" is a **caught** `i_atom_realise_core` throw (`panel_bst_dynamics.m:709`)
— graceful message, no crash. See §5 for the two distinct causes.

---

## 3. Step-by-step example workflows (this recording)

### A. Design an atom's impulse response (Design mode)
1. Open the session (§0). `Apply` toggle **OFF** (Design).
2. `Atoms ▸ Set operator ▸ Dirac (connectome)`.
3. Atom list → **Create atom** (adds a default diffusion atom, selected).
4. `Filter:` → **Heat** (static). Adjust the **Rate** slider.
5. Toolbar **Localize** ON → click a left-hemi cortex vertex (re-seeds the atom). Localize OFF.
6. The cortex shows the atom's **impulse response** spreading along fibers from the seed (peak-normalised,
   `source` colormap) + orientation quivers.

### B. Preview the filtered real source (Apply/Preview)
1. With the Heat Dirac-connectome atom + a seed from A, set the time cursor into signal
   (`panel_time('SetCurrentTime',30)` or scrub the recording).
2. **Apply** toggle **ON**. → cortex repaints the **filtered fiber-spread current magnitude** over a 4 s
   window + mid-frame quivers. Info: `Dirac-Connectome | heat [Preview: fiber-spread cortex + quivers;
   source-space only, no sensor]`. ✅ validated.
   (A **Dirac-dSPM** atom instead also overlays a filtered **sensor** forward on the recording figure.)

### C. Design a tight frame (filterbank)
1. Frame section: set `N` = 6. Click **Design tight frame**. → replaces the bank with an itersine tight
   frame of N members; `A/B/B:A` readout ≈ 1 (tight). ✅ validated.
2. `Show coverage` → opens the frame response view `∑ g_m(λ)²`. ✅ validated (toggles on/off).
   ⚠️ Design-frame **resets the seed to vertex 1** — re-Localize before Preview/Helmholtz.

### D. Analyze → windowed scalogram (Morlet full-time, coeff-space)  ★ new
1. Frame from C (or any atom). Set the cursor anywhere.
2. **Analyze**. → computes per-member energies `[3×nT×M]` {Global,LH,RH} in **coefficient space** over
   the **whole recording**, saves a standalone timefreq, and opens it.
3. ✅ validated: **62 s in 0.16 s**; the scalogram's **Time == the recording's own time base** (0–61.98 s,
   74381 pts) so it **co-displays without** the "Time definition not compatible" popup; tight-frame
   residual ≈ 0 %; LH/RH split real (Global = LH+RH).
   ⚠️ cosmetic: the timefreq viewer labels the √λ spatial-scale center as "Hz"; a 1-member frame sits off
   the auto axis → design a multi-member frame (C) to see filled bands.

### E. Localize bands → marker atoms  🐞
1. Frame + seed. **Localize bands**. → localizes each frame band (peak vertex / time window / level set)
   into a **separate** `dynamics_framebands_*.mat` table and opens it.
2. 🐞 **BUG (validated):** opening the table via `view_dynamics(out)` **replaces the global
   `DynamicsTarget` with a source-less session** (new 3-D figure, `DynamicsOverlay` absent). Every
   source-dependent action afterwards (Apply/Analyze/Helmholtz) fails with **"no real source linked"**
   until you reopen the source session. See §5-B1.

### F. Helmholtz on the filtered field (Hodge decomposition)  ★ new
1. **Dirac** or **Dirac (connectome)** atom + a **static** kernel (Heat) + a seed; cursor in signal.
2. **Helmholtz (filtered)**. → builds the g(λ)-filtered current at the cursor frame, runs
   `process_helmholtz`, paints **Curl** (signed, `stat2`). Info: `Helmholtz (filtered) | heat | Curl |
   |Div|max=3.8e+02 |Curl|max=3.8e+02 HarmFrac=0.089`. ✅ validated.
3. The painted scalar follows the **Measure** selection (§G) if one is set (Div/Curl/Potential/Stream),
   else defaults to Curl. The filtered map **differs from the raw** (§G): live `corr≈0.14`.

### G. Raw-source differential maps (Measure menu)
1. **Measure** → pick **Divergence / Curl / Potential / Stream** (or None to restore).
2. → per-frame Hodge decomposition of the **raw** unconstrained current at the cursor, painted signed
   (`stat2`); updates as you scrub. ✅ validated finite: |Div|≈113, |Curl|≈195, |Φ|≈0.069, |Ψ|≈0.048.
   (This is the pre-existing raw path; **Helmholtz (filtered)** in §F is its filtered counterpart.)

---

## 4. End-to-end "happy path" (copy/paste to reproduce a full validated run)
```matlab
% 0. open on the test source, cursor into signal
hFig = view_dynamics('FromResult', srcFile);  panel_time('SetCurrentTime',30);
% 1. Dirac-connectome atom, static heat kernel, left-hemi seed
panel_bst_dynamics('OnCreateAtom');
panel_bst_dynamics('OnSetOperator','Dirac-Connectome');
c = bst_get('PanelControls','Dynamics'); c.jKernel.setSelectedIndex(3);   % Heat
panel_bst_dynamics('OnKernelChange');
st=getappdata(0,'DynamicsTarget'); S=in_tess_bst(st.T.SurfaceFile,0); [~,il]=tess_hemisplit(S);
sv=il(round(numel(il)/2)); st.T.Groups(st.curAtom).vertices=sv; st.atomSeed=sv; setappdata(0,'DynamicsTarget',st);
% 2. Preview -> 3. tight frame N=6 -> 4. Analyze -> 5. Helmholtz
panel_bst_dynamics('i_atom_apply');
c.jFrameN.setValue(int32(6)); panel_bst_dynamics('OnDesignFrame');
panel_bst_dynamics('OnAnalyzeWindow');
% (re-seed after Design frame if needed), then:
panel_bst_dynamics('OnHelmholtzFiltered');
% NB: run 'Localize bands' LAST — it currently hijacks the session (§5-B1).
```

---

## 5. Test targets (from the live findings)

### A. Unit tests (pure, headless — no GUI, no Brainstorm data)
1. **`bst_eigenwavelet('ScalogramEnergy')` ≡ `Scalogram`** on synthetic ax/C/gCell (scalar **and**
   quaternion Phi): energy/residual/centers match to <1e-8. *(Already written:
   `scratchpad/test_scalenergy.m`; promote into the repo test suite.)* Machine-precision this session (4e-16).
2. **Hemi split:** with an explicit `hemi.isL/isR` partition matching two blocks, `ScalogramEnergy`
   equals the block-based fallback; **Global == LH + RH** to ~1e-10 (single-block whole-brain ax → real
   split, not degenerate).
3. **`ScalogramEnergy` empty/degenerate inputs:** all-empty `C` → valid struct with `energy [3×0×M]`;
   single-block ax with `hemi=[]` → all energy in LH (documented fallback).
4. **`process_helmholtz('Compute', J, Cov)` contract:** `J` is `[3nV×nT]` interleaved x,y,z; returns
   finite `Div/Curl/Phi/Psi/…` and `HarmFrac∈[0,1]`; `Curl/Psi/Vsol` sign-invariant under manifold
   normal orientation, `Div/Phi` flip (see `manifold-face-normals-inward`).
5. **`i_project_fulltime` equivalence:** block-wise scalar projection == one-shot `manifold_ft` on a
   short synthetic field (coeffs identical); vector path == `i_vector_coeffs`. Guards memory (never
   materialises `[nV×N]`).

### B. Integration tests (booted Brainstorm + the test source)
1. 🐞 **`OnLocalizeBands` must NOT clobber the source session.** After Localize bands, assert the
   original source figure's `DynamicsOverlay.srcResult` is still intact and a following
   `OnHelmholtzFiltered`/`OnAnalyzeWindow` does **not** report "no real source linked". *(Currently
   FAILS — `view_dynamics(out)` replaces `DynamicsTarget` with a source-less table. Fix: open the
   marker table in a separate view without hijacking `DynamicsTarget`, or restore the source target
   afterward.)*
2. 🐞 **LB-Connectome Design realiser.** `OnSetOperator('LB-Connectome')` is gate-allowed but
   `i_atom_realise_core` throws `bst_eigenfilter: unsupported eigen variant 'LB-Connectome'`. Decide:
   add LB-Connectome support to `bst_eigenfilter` **or** drop it from the impulse-realiser gate. Test
   should pin whichever contract is chosen (and confirm Apply/Analyze — which use the projection path,
   not `bst_eigenfilter` — still work for LB-Connectome).
3. ✅ **Helmholtz filtered vs raw differ.** Dirac-connectome + heat → `OnHelmholtzFiltered`; assert
   `Div/Curl` finite, `HarmFrac<0.5`, and `corr(filteredCurl, rawCurl) < 0.9` (fiber-spread ≠ raw).
4. ✅ **Full-time scalogram co-display.** `OnAnalyzeWindow` → saved timefreq `Time` equals the
   recording time vector (max diff <1e-9), `TF` finite `[3×N×M]`, `RowNames={Global,LH,RH}`, and
   `view_timefreq` opens without throwing the global-time-conflict error.
5. ✅ **Static-kernel gate.** With the **default diffusion** (dynamic) kernel, `OnHelmholtzFiltered`
   and `i_atom_apply` (Dirac) return the "…is dynamic; use a static kernel" info and paint nothing;
   with Heat they succeed.
6. ⚠️ **Design-frame seed reset.** After `OnDesignFrame`, `st.atomSeed`/`Groups.vertices` revert to
   vertex 1 — assert-and-document (or preserve the prior seed).
7. **Operator round-trip.** For each gate-allowed operator, `OnSetOperator` then read
   `i_atom_op(st)` back == the requested variant; gate-blocked ones surface a message, not a crash.

### C. Harness notes
- **Never drive the panel via the MCP `run_matlab_file` tool** — it resets the MATLAB path and **kills
  the running Brainstorm session** (observed this session). Use inline `evaluate_matlab_code`, or a
  headless `brainstorm server` boot inside the test process.
- The Dynamics session **tears down if an action throws hard** (`DynamicsTarget` empties). Integration
  tests should re-assert the session between steps and reopen on teardown.
- Local subfunctions are reachable for testing via the `eval(macro_method)` dispatch, e.g.
  `panel_bst_dynamics('i_atom_realise_core', ax, kernel, kp, seed, seedDir)`.
