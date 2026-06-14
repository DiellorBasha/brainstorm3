# Exploring the quaternion structure of the Dirac operator for vector-field filterbanks

**Date:** 2026-06-13
**Author:** Diellor Basha (with Claude)
**Status:** Exploration / design precursor (not yet an implementation spec)

## Goal

Generalize the scalar Laplace–Beltrami (LBO) spatial filterbank — `W = Φ·diag(g(λ))·Φ'·M·δ_v`,
a vertex delta reconstructed through a kernel `g(λ)` on the eigenvalue axis (heat,
wave, Mexican-hat, DoG, Matérn…) — to the **Dirac square operator** and its eigenmodes,
so we can build *vector* wavelets with knobs for **vertex delta + direction + chirality
+ filter parameters**.

This note maps the quaternion structure as actually implemented, then proposes a
concrete filter family that uses it. Every structural claim below was verified
numerically on the real `Subject01` cortex Dirac operator (`tess_cortex_pial_low`,
K=400, Tau=0.5).

## What the operator actually is (verified)

- **Quaternion field, not a 3-vector field.** The Dirac operator `A` is `[4nV × 4nV]`,
  acting on quaternion-valued functions `ψ: V → ℍ`, each vertex carrying
  `q = w + x·i + y·j + z·k` in a contiguous `[w,x,y,z]` block. Mass `B = kron(Mass, I4)`.
  Construction: `A = (1−τ)·D_int²/s_L + τ·E/s_E` (`tess_operators.m:195–217`), where
  `D_int²` is the intrinsic (immersion/edge) quaternionic Dirac-squared from Crane et al.
  "Spin Transformations of Discrete Surfaces" (`local_dirac_intrinsic_sq`,
  `tess_operators.m:277–302`) and `E` is the extrinsic (Gauss-map) shape term.
- **The source 3-vector is a pure-imaginary quaternion.** `bst_dirac` embeds each
  per-vertex gain `g=[gx,gy,gz]` as `ψ = (0, gx, gy, gz)` (w-slot = 0), projects with
  `L̃_h = Ψ_h'·B_h·Φ_h`, and on reconstruction keeps the imaginary part and **drops w**
  (`bst_dirac.m:17–29, 136–142, 170–201`). So the physical vector field lives in
  `Im(ℍ) = ℝ³` (ambient), while the operator's dynamics couple all four components.
- **Eigenvalues come in exact 4-fold multiplets.** Verified: within-quartet spread /
  median eigenvalue = **1.05e-15**. Quartets: `{0,0,0,0}, {λ₁×4}, {λ₂×4}, …`. The
  scalar part of `D_int²` equals the cotan Laplacian, so the spectrum tracks LBO scale,
  but each scale is **4× degenerate** — that degeneracy is the internal quaternionic
  structure (`tess_eigen.m:199–244`; multiplets handled by Rayleigh–Ritz).

## The structural hook: right-multiplication is an exact internal symmetry

Because the operator is built from quaternion **left**-multiplication (edge vectors act
on the left), global **right**-multiplication `R_q: ψ ↦ ψ·q` commutes with `A` by
associativity. Verified on the real operator for `R_i` (right-mult by `i`):

| Test | Result | Meaning |
|------|--------|---------|
| `‖[A, R_i]‖ / ‖A‖` | **1.3e-16** | `R_i` commutes with the Dirac operator |
| `‖R_iᵀ B R_i − B‖ / ‖B‖` | **0** | `R_i` is exactly B-orthogonal |
| `‖R_i² + I‖` | **0** | `R_i² = −I` → an exact **complex structure** |
| `R_i` on one quartet | eigenvalues **{+i,+i,−i,−i}** | each quartet splits into two helicities |

`R_i, R_j, R_k` all commute with `A` and satisfy the quaternion algebra
(`R_p R_q = R_{qp}`), so they generate an **su(2) of internal rotations that commutes
with the whole operator**. Pick any unit imaginary axis `n̂`; `R_n̂` is a complex
structure splitting every eigenspace into `±i` eigenspaces — **the two chiralities**.

This is the axis the scalar LBO does not have and the connection-Laplacian only has in
a single complex line. The Dirac operator gives, at *every* spatial scale `λ`, a full
internal SU(2) to steer direction and a clean `±` to select handedness.

## A concrete filter family: F = Φ · g(λ) · P(n̂,±) · Φ' · B · (δ_v ⊗ q₀)

Three independent knobs, evaluated in mode space then reconstructed (drop w):

1. **Scale — `g(λ)`** (reuse the existing `eigfilter` library verbatim): heat
   `e^{−tλ}`, wave `cos(t√λ)`, Mexican-hat `tλ·e^{−tλ}`, DoG, Matérn `(κ²+λ)^{−ν}`.
   Sets the spatial localization of the vector wavelet exactly as in the scalar case.
2. **Direction — the seed quaternion `q₀`.** Inject a *quaternion delta*
   `δ_v ⊗ q₀` with `q₀ = (0, d̂)` a chosen launch direction. Because `D_int²` carries
   the spin-connection, `d̂` is **parallel-transported by the curved cortical frame** as
   the kernel spreads — the wavelet is a frame-following vector field, not a constant
   ambient arrow. (This is the curvature-aware payoff over a naïve ambient quiver, and
   is consistent with [[dirac-ambient-flat]]: the basis is ambient-smooth, the transport
   is what bends the arrows.)
3. **Chirality — the right-multiplication projector `P_±(n̂) = (I ∓ i·R_n̂)/2`.**
   Selects one helicity of the quartet. `P_±` is complex, so the output is a complex
   (analytic) vector field: real and imaginary parts are the two quadratures of a
   **circularly-polarized / rotating vector wavelet**, and `P_+` vs `P_−` flip its
   handedness. Ties directly into the existing complex analytic-inverse machinery
   (`bst_eigenmode_analytic_inverse`).

### Worked variants this unlocks
- **Directional heat wavelet:** `g=e^{−tλ}`, seed `(0,d̂)`, no chirality projector —
  a smooth vector blob pointing along the transported `d̂`. The vector analogue of the
  scalar heat wavelet.
- **Chiral Mexican-hat:** `g=tλ·e^{−tλ}` ⊗ `P_+(n̂)` — a band-limited, single-handed
  rotating vector wavelet. A genuinely new primitive: scale-selective *and* helicity-
  selective.
- **Helical traveling wave:** `g=e^{iω(λ)t}` ⊗ `P_±` — a vector wave that both
  propagates (phase on λ) and rotates (helicity), the vector-field generalization of the
  analytic eigenmode wave.
- **Helmholtz/Hodge emphasis:** weight the w-channel (≈ divergence/gradient sector,
  since `scalar(D_int²) = cotanL`) against the imaginary channel to bias curl-free vs
  divergence-free content — a Hodge knob orthogonal to scale.

## The subtlety to get right (and design around)

`R_n̂` acts on the **internal/spinor index**, mixing `w` and the vector part — it is
**not** the ambient SO(3) rotation of the physical arrow. Ambient rotation is
conjugation `u*(·)u = L_{u*}R_u`, and `L` (left-mult) does **not** commute with `A`
(correct: for a fixed embedding the operator is not ambient-rotation invariant). So:

- The exactly-commuting symmetry we filter with is **right-mult (internal)**, which is
  why chirality/helicity is well-defined and basis-independent across multiplets.
- The *physical* reading of `P_±` is a **rotating (circularly polarized) vector field**,
  not a static rotated arrow. "Direction" (the seed `d̂`) and "chirality" (the handedness
  of rotation) are therefore independent, exactly as the user's intuition framed them.
- The eigensolver returns a B-orthonormal basis of each quartet but does **not**
  canonicalize it to `{φ, φi, φj, φk}`. We do **not** need to: `R_n̂` is a global
  `kron(I_nV, ρ_R(n̂))` operator we apply directly to coefficient vectors (verified to
  preserve every multiplet), so `P_±` is well-defined regardless of the intra-multiplet
  basis. No re-canonicalization required. ⚑ This is the implementation shortcut.

## Implementation path (extends the existing suite, minimal new surface)

1. **`bst_dirac_filter(CompHM/EigenNode, δ_v, g, q₀, n̂, helicity)`** — mode-space core:
   build the quaternion delta, `c = Φ'·B·(δ_v⊗q₀)`, apply `diag(g(λ))` and
   `P_±(n̂)`, reconstruct via the existing `bst_dirac(...,'Reconstruct',...)`. Reuses
   `bst_eigfilter_kernel`/`evaluate` unchanged for `g(λ)`.
2. **`ρ_R(n̂)` helper** — the 4×4 right-multiplication matrix; `R_n̂ = kron(I_nV, ρ_R)`.
   (For `n̂=i`: `ρ_R(i)=[0 −1 0 0; 1 0 0 0; 0 0 0 1; 0 0 −1 0]`, verified.)
3. **Viewer** — reuse the `figure_3d` source-vector quiver overlay (now with the
   true-size toggle) to render the resulting 3-vector wavelet; complex output → show the
   real quadrature, with a key to step the helical phase.
4. **Process/UI** — a `process_dirac_wavelet` mirroring `process_eigenmodes_wavelet`,
   options: vertex pick, direction `d̂`, filter type+params (reuse eigfilter), chirality
   `{none, +, −}`, axis `n̂`.

## Suggested first experiment (cheap, decisive)

On `Subject01` low-res cortex: pick a gyral-crown vertex, `g=heat(t)`, seed `d̂` along a
principal curvature direction. Render (a) no projector, (b) `P_+(i)`, (c) `P_−(i)`.
Expect: (a) a transported directional blob; (b,c) opposite-handed rotating vector
wavelets. Quantify handedness by the sign of the surface curl / phase winding around `v`.
This validates the whole chain before any UI work.

## First-order vs squared Dirac: verified findings (2026-06-13)

We use the **squared** operator `A = (1−τ)D_int²/s_L + τ E/s_E` (PSD, clean
generalized eigensolve). Question: does squaring lose chirality? Answer, verified
numerically on the `Subject01` left hemisphere (intrinsic operator, raw `L4`):

- **Exposed the first-order intrinsic Dirac `D`** (Crane, rebuilt verbatim):
  `[4nF × 4nV] = [81920 × 40968]`, rectangular **vertex→face** (a spin-connection
  twisted `d`), with `D·(constant quaternion) = 0` (1e-13) → the expected 4-D harmonic
  kernel. `L4 = D'·MF·D` reproduces the intrinsic squared operator.
- **There are TWO independent chiralities, and squaring treats them oppositely:**
  - *Internal / polarization chirality* `R_n̂` (right-multiplication): `R_i,R_j,R_k`
    commute with `L4` at **1e-16** → eigenspaces are exact **4-fold quaternionic
    multiplets**. This grading **survives squaring** (it commutes with `D` and `D²`).
    This is the chirality the demo recovered.
  - *Spectral-sign chirality* `sign(𝔇)`: the self-adjoint first-order operator is the
    block `𝔇 = [[0, D†];[D, 0]]` on **vertex⊕face** (the quaternionic Hodge–Dirac).
    Its spectrum is symmetric `±μ` with `μ² =` the squared eigenvalues, so **squaring
    is 2-to-1 and folds the sign away** (`μ²` cannot encode `sign μ`). The sign is the
    *relative vertex↔face phase* `(ψ_v, ±ψ_f)`, invisible to the vertex-only `D²`.
  - **The two are independent.** On one clean quartet (`μ²=105.17`), the doubled 8-D
    space factorizes into all four `(sign 𝔇, R_i) ∈ {±1}×{±i}` sectors, 2 real dims
    each; `Γ=𝔇/μ` and `R_i` commute, `Γ²=+I`, `R_i²=−I`.
- **Practical:** the spectral chirality is recoverable **without** an indefinite
  eigensolve — build the squared eigenbasis `Φ` (what we have) plus face partners
  `Ψ = DΦ/μ`; the per-quartet block `[[0, μI];[μI, 0]]` then carries `sign(𝔇)` on top
  of `R_n̂`. ⚑ But a single first-order `D` exists only for the **pure** operator
  (intrinsic XOR extrinsic) — the `τ`-blend is a sum of squares with no first-order
  root. So a "signed branch" of the filterbank must commit to `τ=0` or `τ=1`.

**Implication for the filterbank:** the current squared basis already delivers
scale × direction × internal-helicity. A *first-order branch* (pure intrinsic,
vertex⊕face) would add the genuinely new axis: `sign(𝔇)` projectors = a surface
**Riesz/Hilbert transform**, one-way traveling vector waves (signed `e^{iμt}`,
no ad-hoc Hilbert), and harmonic-spinor / spin-structure / index content.

## Open questions worth probing
- Which axis `n̂` is physically meaningful — the surface normal `N̂(v)` at the seed
  (helicity about the local normal = the natural "cortical handedness"), or a fixed
  ambient axis? `N̂(v)` is the likely answer and makes chirality intrinsic.
- Does the `w`-channel weighting give a usable Hodge (curl/grad) separation in practice,
  or is it contaminated by the extrinsic term `E`? Test with `τ=0` (pure intrinsic).
- Relationship to [[dirac-quaternion-frame-rotation]]: the right-ℍ gauge here is the same
  quaternion-coefficient gauge flagged for source-vector differential analysis.
