# Dirac–Helmholtz analysis of cortical source-vector vortices

**Date:** 2026-06-14
**Author:** Diellor Basha (with Claude)
**Status:** Methods reference (analysis validated; not yet a packaged process)

## Phenomenon

The unconstrained Dirac dSPM inverse (`bst_inverse_dirac`) recovers a full **3-D
source vector field** `J(x,t)` on the cortex. Rendered as a quiver, it shows
**vortices** — the vectors winding around a core — often flipping chirality and
appearing to advect. Test case: `Subject01/S01_AEF_..._01_notch/data_block001_02_band.mat`
(alpha 9–12 Hz), **t=216.147 s, superiorparietal-R** (Desikan–Killiany).

## Core method: the Dirac is the gauge-correct Helmholtz operator

Embed the per-vertex 3-vector as a pure-imaginary quaternion `ψ_J = (0, Jx, Jy, Jz)`
and apply the **intrinsic first-order Dirac** `D_int` (nxr `diracIntrinsicD`,
`[4F×4V]`, immersion/edge-based). Because `leftMulImag(e)·J = −(e·J) + (e×J)`, the
face quaternion `φ = D_int·ψ_J` is the Helmholtz decomposition in one operator:

```
φ per face:   w-part   = Σ(edge·J)  = DIVERGENCE  (source/sink)
              im-part  = Σ(edge×J)  = CURL 3-vector  (·n̂_face = vortex chirality)
```

This is the spin-connection-correct curl/divergence on the folded cortex — the
principled replacement for ad-hoc tangential-projection + per-face curl (which is
frame-naive). Keep the **full 3-D field**; the only `·n̂` is to extract the
chirality scalar from the 3-D curl.

## Pipeline (validated on the superiorparietal-R alpha vortex)

1. **Helmholtz**: `D_int` → per-face curl·n̂ (vortex) and divergence (source/sink).
2. **Scale-resolve** (measured, not guessed): project `J` onto the Dirac
   eigenmodes, sweep a continuous Mexican-hat `g(λ)=(tλ)e^{−tλ}`, reconstruct each
   band, run the Helmholtz, and take the **area-weighted** `∫curl·n̂ dA` / `∫div dA`
   over a **geodesic disk around the curl extremum**. Result here: vortex (curl)
   peaks fine (~mode 276, broad shoulder), divergence coarse (~mode 96), cross
   ~170 → scale-separable. ⚑ Do NOT use hand-picked mode-index bins.
3. **Divergence diagnostic** (normal vs tangential): split `J = J_n + J_t`,
   `D_int` each. The **net** source/sink flux is ~100 % tangential (real in-surface
   feature); the radial part (≈60 % of the field, MEG's weak direction) makes large
   *local* `(J·n̂)×curvature` divergence that **cancels to ~0 net** — it inflated the
   apparent divergence without a net source/sink.
4. **Time-resolve** over the alpha cycle: track area-weighted `∫curl·n̂ dA` and
   `∫div_tan dA`. Both oscillate at alpha; they are in **quadrature** (corr ≈ 0.1).
   Curl⊥div is the signature of a **rotating** pattern (a standing oscillation has
   curl & div *in phase*).
5. **Rotation rate**: phase of the `(curl, div)` trajectory (z-scored, `atan2`,
   unwrap, linear fit — **no temporal Hilbert**). Here: **10.7 Hz, 3.87 deg/ms,
   steady (R²=1.00, a clean circle)** — the alpha rhythm as a rotating wave.
6. **Core tracking**: the `|curl|`-weighted centroid of the *scale-isolated* field
   is robust (the raw `|curl|`-extremum jitters). Here: net displacement 1.4 mm,
   RMS spread 1.6 mm → the vortex **rotates in place, it does not translate**.
7. **Angular-momentum ledger**: orbital `L(t) = ∫(r×J)·n̂ dA` about the core tracks
   the chirality (corr 0.89). Internal **spin = Im(J*×J) = 0** exactly (real
   instantaneous field). Net translation ~0. So `J_total = L_orbital`: a steady
   in-place rotation, no spin, no linear momentum.

**Net picture:** the superiorparietal-R alpha is a **steady in-place rotating
wave** at 10.7 Hz, cycling vortex(CCW) → source → vortex(CW) → sink each period.
A clear vortex riding on a weak, diffuse, genuine in-surface source/sink.

## Figure recipe (Brainstorm-native; required for vector results)

Render in the real `figure_3d`, never hand-rolled trisurf. To show a derived
scalar (divergence/curl) under the actual or scale-filtered vector field:
1. Save the scalar as a temp 1-component results file → `view_surface_data` (gives
   a correct Brainstorm colormap + colorbar; set `cmap_rbw`, absolute off for signed).
2. `setappdata(hFig,'QuiverVectorOverride', V)` with the vector field `V [nVert×3]`
   (the `figure_3d` quiver now reads this — decouples the quiver vector from the
   colored scalar), then `figure_3d('SetShowSourceVectors',...)`.
3. `delete(findobj(hFig,'-regexp','Tag','[Ss]cout'))`, zoom, `out_figure_image`.
⚑ NEVER `close all force` (corrupts GlobalData + scout registry); recover with
`brainstorm stop; brainstorm nogui` (base workspace survives). Don't delete/re-add
colorbars on a live source figure (breaks the shared colormap's threshold marker).

## Implementation dependencies (all committed)

- nxr-compute `diracIntrinsicD` (`515ad15`): intrinsic first-order Dirac.
- `tess_operators` stores `FirstOrder.Intrinsic/.Extrinsic` + `FaceMass` in the
  operator node (`a12f835c`).
- `figure_3d` `QuiverVectorOverride` appdata (`f6094353`).

## Open questions / next

- Is the in-place rotation **specific** to this region/time, or general across the
  alpha-vortex population (beta@central, other bands)? → needs the packaged process.
- The signed `𝔇 = [[0,D†];[D,0]]` propagation (Riesz) branch is wired but unused —
  relevant if/when a vortex that **does** translate is found.
- Package steps 1–7 into `process_dirac_vortex` (scan any region/time, auto-label
  chirality / rotation rate / translate-vs-in-place) once the method is locked.
