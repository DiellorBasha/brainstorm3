
All projects
cortical-flow
In this project, we will develop the Global Cortical-Flow framework (GCF; Phase 1) to track how brain activity propagates across the cortical surface. GCF estimates velocity fields from neurophysiological activity across frequencies and spatial scales on the folded cortex. By quantifying multiple features of cortical dynamics, GCF provides descriptors of brain function that considerably expand the neuroscience toolkit. A key advance is that these measurements support short-term prediction: given an observed pattern, we can estimate its expected displacement over short time windows along the cortical sheet, enabling testable predictions of propagation. This closes the loop between description and prediction, enabling new tests of expected propagation along the functional hierarchy of the cortex across the healthy lifespan (Phase 2) and sensitive detection of abnormal propagation patterns in disease (Phase 3). We hypothesize that cortical activity can be modelled across spatial scales and neurophysiological frequencies as a velocity field, with recurrent source-sink motifs, propagation speed, and trajectories serving as core descriptors of brain dynamics. These descriptors are expected to capture healthy developmental and aging variants, while deviations may yield clinically meaningful markers of disease, including abnormal propagation signatures of epilepsy, including between seizures. Overarching Goal & Specific Aims: We aim to develop and disseminate a scalable and validated cortical-flow framework for dynamic brain mapping, through three Phases: Phase 1: Advance the Global Cortical-Flow framework. We will build and validate methodology for cortical-flow mapping, including joint scale-frequency representations and event-based propagation summaries, to derive robust descriptors of cortical activity propagation that are comparable across individuals and ready for population and clinical applications. All these advances will be implemented and widely disseminated through our established open-source software (Brainstorm7 ). Phase 2: Characterize normative cortical flow across the lifespan. Using our large Lifespan Cohort (n≈1,700; ages 4-88), we will chart the normative distribution of cortical-flow descriptors across development and healthy aging. Analyses will be sex-stratified to detect sex-specific patterns. We will test whether slow spontaneous activity (δ-α, 2-15 Hz) propagates preferentially bottom-up along the cortical functional hierarchy, whereas faster β activity (15-35 Hz) flows top-down in the reverse direction, with modulations in aging. We will also assess whether propagation trajectories and synchronized source-sink motifs provide proxy measures of brain network connectivity. The outcome will be the first population-scale atlas of cortical-flow, openly shared through Brainstorm and neuromaps8 (Summary of Progress) to benchmark discovery and clinical applications. 2 Phase 3: Demonstrate clinical value. We will test whether cortical-flow source-sink motifs and propagation pathways localize epileptogenic regions and improve presurgical hypotheses compared to current standards. This work will leverage the normative atlas from Phase 2 to generate patient-specific deviation maps and will be benchmarked directly against standard clinical workflows. We will leverage two well-characterized cohorts of epilepsy patients with longitudinal follow-up: one pediatric (n=121) and one adult (n=300), already collected at two sites to enable generalizability testing. This phase will evaluate cortical flow as a clinically meaningful tool for surgical planning, with potential to advance research and inform intervention for other conditions such as Alzheimer’s and Parkinson’s disease, where early brain dynamic changes can predict clinical evolution
Show more



How can I help you today?



Start a task in Cowork
Wang GBF paper benchmarking and seizure datasets
Last message 26 minutes ago
MEG source mapping with eigenmodes and harmonic basis functions
Last message 38 minutes ago
Comparing eigenmode analysis across different brain anatomies
Last message 5 hours ago
Building a lightweight MEG forward solver with Zarr storage
Last message 7 hours ago
Can you tell me what are cross...
Last message 10 hours ago
Wavelet packets vs continuous wavelet transforms
Last message 15 hours ago
Heat distance algorithm in geometry central
Last message 19 hours ago
Mass matrix types in differential geometry
Last message yesterday
Layered architecture for NXR geometric precomputations
Last message 5 days ago
Integrating nxr-compute with Brainstorm for eigenmode source mapping
Last message 6 days ago
Laplace-Beltrami operator and material properties
Last message 6 days ago
Choosing a Zarr implementation for MEG/EEG data storage
Last message last week
Forward mapping and source space decimation in MNE
Last message 3 weeks ago
MNE-CPP modules and MNE Browse architecture
Last message 3 weeks ago
Conductivity in Maxwell's equations and weighted Laplace operators
Last message 3 weeks ago
Integrating C++ libraries with Spectra for eigenmode computation
Last message 4 weeks ago
Curvature and eigenmodes bridging dimensions
Last message 4 weeks ago
Cellular vs simplicial complexes
Last message 4 weeks ago
Learning convolutional neural networks for cortical surface analysis
Last message 4 weeks ago
Wave and heat equations comparison
Last message last month
Eigenmode-based MEG analysis repository
Last message last month
Differential geometry and eigenmode analysis on manifolds
Last message last month
Intrinsic metrics on Riemannian manifolds
Last message last month
Fourier analysis for cortical dynamics
Last message 2 months ago
Integrating cortical flow codebases with Claude Code
Last message 2 months ago
Fourier-based level of detail in video
Last message 2 months ago
Cortical flow analysis using geometry processing
Last message 2 months ago
Spatial-temporal frequencies and Fourier analysis in videos
Last message 2 months ago
Surface parametrization for cortical mapping and visualization
Last message 2 months ago
Memory
Only you
Purpose & context Diellor is a researcher and software engineer working at the intersection of computational neuroscience, differential geometry, and signal processing. Core expertise spans MEG/EEG source imaging, spectral geometry, discrete exterior calculus (DEC), and mathematical physics. Work is both theoretical (developing novel frameworks) and applied (building production C++ tooling with WASM/JS bindings). The overarching research direction involves treating cortical dynamics not as scalar fields on a static substrate but as geometric objects with intrinsic structure — manifolds, fiber bundles, tangent fields, and spectral decompositions — amenable to rigorous mathematical analysis. A recurring theme is grounding neuroscientific data analysis in physically and geometrically principled constructions rather than ad hoc signal processing. Current state Diellor is actively developing two interrelated research frameworks: Doubly-spectral source imaging framework: A unified λ–ω (spatial eigenmode × temporal frequency) domain where MEG/EEG source imaging operations — denoising, regularization, neurophysiological constraints — are explicit, separable, and interpretable as operations on a 2D spectral image. Key innovations include eigenmode-direct forward computation (bypassing vertex-wise leadfield for large computational savings), sensor eigenmodes defining a doubly-spectral transfer matrix that empirically determines transmission bandwidth without free parameters, and the neurophysiological dispersion bound (ω ≤ c√λ) as a hard physical mask. The λ–ω image is itself treated as an analyzable object supporting 2D aperiodic-periodic decomposition, traveling wave detection, and principled compression. Connection Laplacian / complex vector wavelet tensor framework: Extends the above to the connection Laplacian on the cortical tangent bundle, yielding vector-valued eigenmodes encoding amplitude, phase, and propagation direction simultaneously. Differential-geometric quantities (divergence, curl, Helmholtz decomposition, phase singularities) are computed as precomputed eigenmode operations with no noise amplification. The spatially varying gain of active cortex is framed as a measurable material property (the wavelet tensor / constitutive relation), with resting-state as the standardized measurement protocol and task conditions as perturbations around the intrinsic baseline. Diellor is also developing a geometry-aware CNN for cortical surface analysis, built on an established preprocessing pipeline (Geometry Central, C++) with pre-computed LBO eigenmodes, trivial connections, Hodge decomposition, DEC operators, and geodesic/heat computations. The CNN design favors physics-informed, non-separable spatiotemporal kernels parameterized by wave speed and damping, unifying heat kernels, difference-of-Gaussians, and wave filters in a single eigenmode-based family. Input data consists of cortical surfaces with time-varying scalar fields (MEG/EEG source maps, dense temporal sampling). Recent concrete engineering decisions include: representing connection Laplacian output as real 2N×2N block expansion by default (complex COO as optional flag), and using a single ConnectionLaplacianOptions struct with defaulted fields for API configuration. On the horizon Two full research proposals have been drafted (doubly-spectral framework; connection Laplacian / wavelet tensor framework) — next steps likely involve formalization, implementation, and validation Open question: fully fixed vs. learnable vs. hybrid filter banks in the geometry-aware CNN Extension of wave propagation kernels to anisotropic case using tangent frame direction fields Deeper integration of the λ–ω image analysis (2D FOOOF generalization, traveling wave detection) into the tooling Key learnings & principles Geometry-first: Physical and geometric constraints should be baked into model architecture and parameterization, not learned away. Unconstrained data-driven approaches are deprioritized given data regimes and the richness of available geometric priors. Resting state as constitutive measurement: Brain states (sleep, wake, arousal) are better understood as phases of a single intrinsic material description rather than different material property maps. Resting/input-free conditions are the standardized protocol for measuring the cortex's constitutive relation; tasks are perturbations δM around that baseline. Separability and explicitness: A central design principle across both frameworks is making implicit assumptions (e.g., regularization as hidden prior) into explicit, inspectable operations in a well-chosen domain. Eigenmode operations over numerical differentiation: Computing differential-geometric quantities via precomputed eigenmode expansions avoids noise amplification that plagues numerical differentiation on noisy neural data. The brain has true spatial extension: Connectivity is not a graph abstraction — axonal and synaptic delays impose genuine propagation limits that must be respected geometrically. Strong preference for analytically grounded constructions (Sturm-Liouville, de Rham complex, fiber bundle geometry) as unifying organizational frameworks rather than treating each method in isolation. Approach & patterns Engages Claude as a collaborative intellectual peer, not a tutor — pushes back on imprecise framing, corrects conflations, and drives toward conceptual clarity Consistently moves from mathematical foundations toward implementation-ready constructions, bridging theory and engineering in the same conversation Prefers to develop intuition for abstract structures through concrete physical examples (electromagnetism, wave optics, neural data) before generalizing Research conversations tend to evolve: beginning with established concepts, then co-developing novel frameworks as the dialogue deepens Saves key outputs as artifacts for later use (research proposals, narrative structures, API designs) Tools & resources C++ / geometry-central: Core computational geometry library; basis of the nxr-compute module and associated WASM/JS bindings MNE-Python, FreeSurfer, fsaverage: Standard neuroimaging stack for MEG/EEG processing and cortical surface reference GBFs (NCC Lab, SUSTech): Reference implementation of LBO eigenmode-based ESI, studied for methodology and dependency patterns Geometry Central preprocessing pipeline: LBO eigenmodes, trivial connections, vector heat method, DEC operators, Hodge decomposition, geodesic and heat-equation tools Libraries surveyed: neuromaps, nibabel, scikit-learn, PyVista (in context of GBFs codebase analysis)

Last updated 1 day ago

Instructions
Add instructions to tailor Claude’s responses

Files
2% of project capacity used

spatiotemporal_filters_proposal.md
735 lines

md



vector_eigenmode_meg_proposal.md
499 lines

md



inverse-methods-as-spectral-filters.md
356 lines

md



diffusionnet-cortical-applications.md
235 lines

md


inverse-methods-as-spectral-filters.md
21.43 KB •356 lines
Formatting may be inconsistent from source

# Inverse Methods as Spectral Filters in Eigenmode Space

> Reference document for the cortical-flow project.
> Maps standard MEG/EEG inverse methods to spectral transfer functions
> in the LBO eigenmode basis, connecting classical source imaging to
> the learned filterbank framework.

---

## 1. Setup: the eigenmode-direct inverse problem

### 1.1 Forward model

The whitened forward model in eigenmode space:

$$\mathbf{d} = \mathbf{L}_\lambda\, \hat{\mathbf{s}} + \mathbf{n}$$

where:
- **d** ∈ ℝ^N — whitened sensor data (N channels)
- **L_λ** = **W** **L** **Φ** ∈ ℝ^{N×K} — whitened eigenmode leadfield
- **Φ** ∈ ℝ^{V×K} — LBO eigenvectors (columns are M-orthonormal eigenmodes)
- **L** ∈ ℝ^{N×V} — vertex-space leadfield from BEM/FEM head model
- **W** = **C_n**^{−½} — whitening matrix from empty-room noise covariance
- **ŝ** ∈ ℝ^K — eigenmode coefficients to estimate
- **n** ~ 𝒩(0, **I**) — whitened noise (identity covariance after whitening)

### 1.2 Key per-eigenmode quantities

For each eigenmode k, define:

- **λ_k** — the k-th LBO eigenvalue (spatial frequency squared)
- **l_k** = L_λ[:,k] — the sensor-space signature of eigenmode k (column of L_λ)
- **σ_k²** = ‖l_k‖² = l_kᵀ l_k — sensor-space power (observability) of eigenmode k
- **φ_k** — the k-th eigenvector (spatial pattern on cortex)

The quantity σ_k² is crucial: it measures how well the sensor array can observe eigenmode k. Low-order eigenmodes (smooth spatial patterns) typically have high σ_k² because they produce coherent signals across many sensors. High-order eigenmodes (fine spatial patterns) have low σ_k² because their signals cancel across neighboring sensors. This falloff is determined entirely by the sensor geometry and head model — it is the sensor array's spatial transfer function.

### 1.3 General linear estimator

Every linear inverse method can be written as:

$$\hat{s}_k = T_k \cdot [\text{data projection onto eigenmode } k]$$

where T_k is the **spectral filter** — a scalar function of the eigenmode index (or equivalently, of the eigenvalue λ_k) that determines how much each eigenmode is preserved or suppressed in the reconstruction. Different inverse methods correspond to different choices of T_k.

### 1.4 Diagonal approximation

When eigenmodes are approximately orthogonal in sensor space (l_jᵀ l_k ≈ 0 for j ≠ k), the inverse problem decouples mode-by-mode and the spectral filter has a clean closed form. This is a reasonable approximation for cortical MEG/EEG because the sensor array is broadband in eigenmode space — it doesn't strongly couple specific eigenmodes. The eigenmode-direct forward model was designed to exploit this near-diagonality.

Under this approximation, the general Tikhonov-type estimator for eigenmode k is:

$$\hat{s}_k = \frac{\sigma_k^2}{\sigma_k^2 + \alpha\, r_k} \cdot \frac{l_k^\top \mathbf{d}}{\sigma_k^2}$$

where:
- The first factor is the spectral filter T_k = σ_k² / (σ_k² + α r_k)
- The second factor is the matched-filter projection of data onto eigenmode k
- r_k is the eigenmode-space regularization penalty (method-specific)
- α is the regularization strength

---

## 2. Method-by-method derivation

### 2.1 Minimum Norm Estimate (MNE)

**Vertex-space formulation:**

$$\hat{\mathbf{s}} = \arg\min_{\mathbf{s}} \|\mathbf{d} - \mathbf{L}\mathbf{s}\|^2 + \alpha\|\mathbf{s}\|^2$$

Penalty is total source power. No spatial preference.

**Eigenmode-space penalty:** r_k = 1 for all k.

The penalty ‖s‖² = ‖Φŝ‖²_M = ŝᵀ ŝ (by M-orthonormality of eigenvectors). All eigenmode coefficients are penalized equally.

**Spectral filter:**

$$T_k^{\text{MNE}} = \frac{\sigma_k^2}{\sigma_k^2 + \alpha}$$

This is a standard Tikhonov filter. For well-observed eigenmodes (σ_k² ≫ α), T_k ≈ 1 — the coefficient passes through. For poorly observed eigenmodes (σ_k² ≪ α), T_k ≈ σ_k²/α — the coefficient is suppressed proportional to its observability.

**Shape:** Lowpass in observability, flat in eigenvalue. The filter depends on σ_k² (sensor geometry), not on λ_k (cortical geometry). MNE has no cortical spatial prior — the regularization is purely driven by what the sensors can see.

---

### 2.2 Weighted Minimum Norm (wMNE)

**Vertex-space formulation:**

$$\hat{\mathbf{s}} = \arg\min_{\mathbf{s}} \|\mathbf{d} - \mathbf{L}\mathbf{s}\|^2 + \alpha\, \mathbf{s}^\top \mathbf{W}_d\, \mathbf{s}$$

where W_d = diag(‖L[:,v]‖²) compensates for the depth bias of MNE by penalizing superficial sources more (they have larger leadfield column norms).

**Eigenmode-space penalty:** r_k = ‖L Φ[:,k]‖² ∝ σ_k² (before whitening) or a related depth-weighted norm.

In eigenmode space, the depth weighting translates to a penalty proportional to the unwhitened sensor-space power of each eigenmode. Eigenmodes dominated by superficial (gyral) sources get penalized more.

**Spectral filter:**

$$T_k^{\text{wMNE}} = \frac{\sigma_k^2}{\sigma_k^2 + \alpha\, w_k}$$

where w_k is the depth-weight for eigenmode k.

**Shape:** Similar to MNE but the depth correction flattens the bias toward well-observed (superficial) eigenmodes.

---

### 2.3 LORETA (Low Resolution Electromagnetic Tomography)

**Vertex-space formulation:**

$$\hat{\mathbf{s}} = \arg\min_{\mathbf{s}} \|\mathbf{d} - \mathbf{L}\mathbf{s}\|^2 + \alpha\|\Delta\, \mathbf{s}\|^2$$

Penalty is the Laplacian of the source distribution — penalizing spatial roughness.

**Eigenmode-space penalty:** In the eigenbasis, the Laplacian is diagonal: Δφ_k = −λ_k φ_k. Therefore:

$$\|\Delta\, \mathbf{s}\|^2 = \sum_k \lambda_k^2\, \hat{s}_k^2$$

So r_k = λ_k².

**Spectral filter:**

$$T_k^{\text{LORETA}} = \frac{\sigma_k^2}{\sigma_k^2 + \alpha\, \lambda_k^2}$$

**Shape:** This is the cleanest result. The filter has two independent suppression mechanisms: poorly observed eigenmodes (small σ_k²) are suppressed by the Tikhonov structure, AND high spatial frequency eigenmodes (large λ_k) are suppressed by the Laplacian penalty. The filter is a 2D function of both observability and spatial frequency.

This is where eigenmode space shines: the Laplacian penalty, which is a complicated sparse matrix operation in vertex space, becomes a simple diagonal weighting by λ_k². The entire LORETA estimator reduces to one scalar formula per eigenmode.

---

### 2.4 Generalized smoothness prior (Lᵝ regularization)

**Vertex-space formulation:**

$$\hat{\mathbf{s}} = \arg\min_{\mathbf{s}} \|\mathbf{d} - \mathbf{L}\mathbf{s}\|^2 + \alpha\|\Delta^{\beta/2}\, \mathbf{s}\|^2$$

Penalty is a fractional power of the Laplacian. β = 0 gives MNE, β = 2 gives LORETA. Intermediate values interpolate.

**Eigenmode-space penalty:** r_k = λ_k^β

**Spectral filter:**

$$T_k^{L^\beta} = \frac{\sigma_k^2}{\sigma_k^2 + \alpha\, \lambda_k^\beta}$$

**Shape:** β controls the steepness of the spatial frequency rolloff. Small β penalizes all eigenmodes nearly equally (approaching MNE). Large β aggressively suppresses high-frequency eigenmodes (stronger smoothness). β = 1 corresponds to a Green's function prior (penalizing the first-order Sobolev norm), which gives a 1/λ decay in the prior variance — matching the empirical observation that cortical source power falls off approximately as 1/λ.

This family unifies MNE and LORETA as endpoints of a continuous spectrum of smoothness priors, all diagonal in eigenmode space.

---

### 2.5 dSPM (Dynamic Statistical Parametric Mapping)

**Vertex-space formulation:** Compute MNE solution, then normalize each source by the standard deviation it would have under pure noise:

$$z_k = \frac{\hat{s}_k^{\text{MNE}}}{\sqrt{\text{var}(\hat{s}_k \mid \text{noise only})}}$$

**Eigenmode-space:** The noise-only variance of the MNE estimator for eigenmode k is proportional to T_k^MNE · σ_k⁻² (the filter gain times the noise projection). The normalization divides out the spatial sensitivity pattern.

**Spectral filter:**

$$T_k^{\text{dSPM}} = \frac{\sigma_k^2}{\sigma_k^2 + \alpha} \cdot \frac{1}{\sigma_k / \sqrt{\sigma_k^2 + \alpha}}$$

$$= \frac{\sigma_k}{\sqrt{\sigma_k^2 + \alpha}}$$

**Shape:** Square root of the MNE filter. The normalization partially compensates for the depth/observability bias: well-observed and poorly-observed eigenmodes are brought closer together than in raw MNE. The output is dimensionless (a z-score), so absolute amplitude information is lost. The filter is still monotonically increasing in σ_k² — it never boosts poorly-observed modes above well-observed ones — but the bias is reduced.

---

### 2.6 sLORETA (Standardized LORETA)

**Vertex-space formulation:** Compute MNE solution, then normalize each source by the square root of the resolution matrix diagonal:

$$z_v = \frac{\hat{s}_v^{\text{MNE}}}{\sqrt{R_{vv}}}$$

where R = C_s Lᵀ (L C_s Lᵀ + α C_n)⁻¹ L C_s is the resolution (model) matrix.

**Eigenmode-space:** In the eigenmode basis (under the diagonal approximation), the resolution matrix diagonal for eigenmode k is:

$$R_{kk} = \frac{\sigma_k^2}{\sigma_k^2 + \alpha} \cdot \frac{\sigma_k^2}{\sigma_k^2 + \alpha} / \text{(eigenmode mixing terms)}$$

Under the diagonal approximation, R_kk ≈ (T_k^MNE)².

**Spectral filter:**

$$T_k^{\text{sLORETA}} = \frac{T_k^{\text{MNE}}}{\sqrt{R_{kk}}} \approx \frac{T_k^{\text{MNE}}}{T_k^{\text{MNE}}} = 1$$

This is the remarkable property of sLORETA: **under the diagonal approximation, the spectral filter is exactly flat.** Every eigenmode passes through with unit gain. The normalization by the resolution matrix perfectly compensates for the Tikhonov suppression.

In practice, sLORETA doesn't give a perfectly flat filter because the diagonal approximation isn't exact — there's eigenmode mixing through the sensor array. But it's close to flat, which is why sLORETA has zero localization bias for point sources: it doesn't preferentially suppress any spatial frequency.

**Shape:** Approximately flat. The output is dimensionless (like dSPM). What sLORETA achieves is not regularization in the spectral sense but rather debiasing — it removes the spectral distortion introduced by the regularization, at the cost of amplifying noise in poorly-observed eigenmodes.

---

### 2.7 eLORETA (Exact LORETA)

**Vertex-space formulation:** Iteratively computes a source-space weight matrix W such that the weighted minimum-norm solution has exactly zero localization bias for any point source at any location.

**Eigenmode-space:** eLORETA finds weights w_k such that:

$$T_k^{\text{eLORETA}} = \frac{\sigma_k^2}{\sigma_k^2 + \alpha\, w_k}$$

produces zero-bias localization. The weights w_k are determined by a fixed-point iteration that accounts for the full (non-diagonal) structure of L_λᵀ L_λ.

**Spectral filter:**

$$T_k^{\text{eLORETA}} = \frac{\sigma_k^2}{\sigma_k^2 + \alpha\, w_k^*}$$

where w_k* are the converged iterative weights. These are not available in closed form but are computable.

**Shape:** Similar to sLORETA's approximately flat response, but achieved through the prior (choosing w_k) rather than through post-hoc normalization. eLORETA preserves amplitude information (unlike sLORETA/dSPM which are dimensionless), and the output is in physical units. The spectral filter is not exactly flat but is flatter than MNE.

---

### 2.8 LCMV Beamformer

**Vertex-space formulation:** For each source location, construct a spatial filter that passes signal from that location while suppressing everything else:

$$\mathbf{w}_v = \frac{\mathbf{C}_d^{-1}\, \mathbf{l}_v}{\mathbf{l}_v^\top\, \mathbf{C}_d^{-1}\, \mathbf{l}_v}$$

where C_d is the data covariance matrix.

**Eigenmode-space:** For eigenmode k:

$$w_k = \frac{\mathbf{C}_d^{-1}\, \mathbf{l}_k}{\mathbf{l}_k^\top\, \mathbf{C}_d^{-1}\, \mathbf{l}_k}$$

The data covariance is C_d = L_λ C_ŝ L_λᵀ + I (whitened), where C_ŝ is the true source covariance. The beamformer uses the empirical data covariance, making it data-adaptive.

**Spectral filter:**

$$T_k^{\text{LCMV}} = \frac{1}{\mathbf{l}_k^\top\, \hat{\mathbf{C}}_d^{-1}\, \mathbf{l}_k}$$

This filter is **data-dependent**: it adapts to the actual covariance structure of the recording. In eigenmode space, the beamformer suppresses eigenmodes that are correlated with other active eigenmodes (because their signals are captured by the dominant components of C_d), and preserves eigenmodes that are uncorrelated with other active sources.

**Shape:** Not a fixed function of λ_k or σ_k² — it depends on which eigenmodes are active in the data. For a single active eigenmode in noise, the beamformer filter converges to the matched filter. For multiple correlated eigenmodes, the filter shows suppression at the correlated modes. This is fundamentally different from the Bayesian methods above: the beamformer is adaptive, while MNE/LORETA/sLORETA apply fixed priors regardless of the data.

---

### 2.9 Minimum-norm with empirical Bayesian covariance

**Vertex-space formulation:** Learn the source covariance C_s from the data (e.g., by fitting a parametric model or using restricted maximum likelihood).

**Eigenmode-space:** If the source covariance is diagonal in eigenmode space:

$$C_{\hat{s}} = \text{diag}(p_k)$$

where p_k are learned from data (e.g., by maximizing marginal likelihood).

**Spectral filter:**

$$T_k^{\text{EB}} = \frac{p_k\, \sigma_k^2}{p_k\, \sigma_k^2 + 1}$$

**Shape:** The prior variance p_k for each eigenmode is learned from data, giving a fully adaptive spectral filter. This is the closest classical method to the learned filterbank — but it's restricted to a diagonal filter (one value per eigenmode), whereas the learned filterbank can also apply cross-eigenmode coupling through its MLP.

---

## 3. Summary table

| Method | Eigenmode penalty r_k | Spectral filter T(λ_k) | Prior on eigenmodes | Output units |
|---|---|---|---|---|
| **MNE** | 1 | σ_k² / (σ_k² + α) | Flat — all modes equal | Physical (A·m) |
| **wMNE** | w_k (depth weight) | σ_k² / (σ_k² + α w_k) | Depth-compensated flat | Physical |
| **LORETA** | λ_k² | σ_k² / (σ_k² + α λ_k²) | Smooth — 1/λ² decay | Physical |
| **Lᵝ general** | λ_k^β | σ_k² / (σ_k² + α λ_k^β) | Smoothness order β | Physical |
| **dSPM** | 1 (then normalize) | σ_k / √(σ_k² + α) | Flat + noise normalization | z-score |
| **sLORETA** | 1 (then normalize) | ≈ 1 (flat after normalization) | Flat + resolution normalization | Dimensionless |
| **eLORETA** | w_k* (iterated) | σ_k² / (σ_k² + α w_k*) | Zero-bias weights | Physical |
| **LCMV** | (data-adaptive) | 1 / (l_kᵀ Ĉ_d⁻¹ l_k) | Data covariance | Physical |
| **Empirical Bayes** | 1/p_k (learned) | p_k σ_k² / (p_k σ_k² + 1) | Learned from data | Physical |
| **Learned filterbank** | (unconstrained) | T(λ_k, ω) — learned 2D | Learned from data, joint λ–ω | Task-dependent |

---

## 4. Visual intuition: filter profiles in eigenvalue space

### 4.1 The observability curve σ_k²

The sensor-space power σ_k² as a function of eigenvalue λ_k is a monotonically decreasing curve, determined entirely by the sensor array geometry and head model. It represents the sensor array's spatial transfer function: how much signal power each eigenmode contributes at the sensors.

For MEG with a typical 300-channel whole-head array:
- Low-order eigenmodes (λ < 100, ~first 20 modes): σ_k² is large. These are smooth, large-scale spatial patterns that coherently drive many sensors.
- Mid-order eigenmodes (100 < λ < 5000, modes 20–200): σ_k² falls off, roughly as 1/λ or faster. Increasingly fine spatial patterns produce increasingly cancelling signals at the sensors.
- High-order eigenmodes (λ > 5000, modes > 200): σ_k² is negligible. These fine spatial patterns are effectively invisible to the sensor array.

This observability curve sets the fundamental bandwidth of the inverse problem. No regularizer can recover eigenmodes that the sensors cannot see.

### 4.2 How each method shapes the filter

All Tikhonov-type methods have the form T_k = σ_k² / (σ_k² + α r_k). The regularization strength α sets an overall threshold, and the penalty r_k tilts the filter:

**MNE** (r_k = 1): The filter follows the observability curve σ_k². Well-observed modes pass; poorly-observed modes are suppressed. The crossover is at σ_k² = α. No spatial frequency preference beyond what the sensors impose.

**LORETA** (r_k = λ_k²): The filter is steeper than MNE — high-frequency eigenmodes are doubly suppressed (by poor observability AND by the Laplacian penalty). The crossover shifts to lower spatial frequencies. LORETA produces smoother source estimates than MNE at the same regularization strength.

**Lᵝ family**: β = 0 gives MNE, β = 2 gives LORETA, and β can be tuned continuously between them. β = 1 (Green's function prior) gives a 1/λ penalty, which matches the empirical 1/f spatial spectrum of cortical activity and is therefore arguably the most neurophysiologically appropriate fixed prior.

**sLORETA**: The post-normalization flattens the filter to approximately 1 everywhere. This means every eigenmode contributes equally to the output — no spatial frequency preference at all. The price: noise from poorly-observed eigenmodes is amplified, and the output is dimensionless.

**Beamformer**: The filter adapts to the data. Eigenmodes carrying strong signal get preserved; eigenmodes carrying mostly interference from correlated sources get suppressed. No fixed spatial frequency profile — the shape changes with every recording.

### 4.3 The learned filterbank in context

The learned filterbank generalizes all of the above. Instead of one filter function T(λ_k) with 1–2 free parameters, it learns a bank of filters with many parameters, combined through nonlinear operations. After training, the effective filter applied to each eigenmode can be extracted and compared directly to the curves above.

If the learned filter looks like σ_k² / (σ_k² + α λ_k²), it has rediscovered LORETA. If it looks like σ_k² / (σ_k² + α λ_k), it has discovered the Green's function prior. If it has a non-monotonic profile (bandpass), it has discovered something no classical method captures — and the location and width of the passband characterize the spatial frequency structure of cortical activity.

The extension to the λ–ω domain adds a second axis. Classical methods apply the same spatial filter at all temporal frequencies (separable). The learned filterbank can apply different spatial filters at different temporal frequencies (non-separable), capturing the dispersion relation and frequency-specific spatial structure.

---

## 5. Implications for the doubly-spectral framework

### 5.1 The regularization axis of λ–ω space

In the doubly-spectral λ–ω image, the λ axis already encodes the spatial regularization through the eigenmode representation. Each column of the λ–ω image (fixed ω, varying λ) is subject to a spectral filter determined by the inverse method:

- With MNE: each column is filtered by σ_k² / (σ_k² + α), independent of ω
- With LORETA: each column is filtered by σ_k² / (σ_k² + α λ_k²), independent of ω
- With a learned non-separable filter: each column gets a different filter depending on ω

### 5.2 The doubly-spectral transfer matrix

The complete input-output relationship of the inverse in the λ–ω domain can be written as a 2D transfer matrix:

$$\hat{S}(\lambda_k, \omega) = T(\lambda_k, \omega) \cdot S_{\text{true}}(\lambda_k, \omega) + N(\lambda_k, \omega)$$

where T(λ_k, ω) is the effective spatiotemporal transfer function and N is the filtered noise. For classical methods, T(λ_k, ω) = T(λ_k) — separable, ω-independent. For the learned filterbank, T(λ_k, ω) is a full 2D function encoding the joint spatiotemporal regularization.

### 5.3 The dispersion bound as a regularizer

The neurophysiological dispersion bound ω ≤ c√λ defines a hard mask in λ–ω space:

$$T^{\text{disp}}(\lambda_k, \omega) = \Theta(c\sqrt{\lambda_k} - \omega)$$

This is a non-separable regularizer that no classical method can express. It says: "suppress any energy at spatial frequency λ and temporal frequency ω if they violate the wave propagation speed limit." This is a physics-based prior that constrains the source estimate to be consistent with cortical wave dynamics.

The dispersion regularizer can be combined with any classical method by multiplication:

$$T^{\text{LORETA+disp}}(\lambda_k, \omega) = \frac{\sigma_k^2}{\sigma_k^2 + \alpha\, \lambda_k^2} \cdot \Theta(c\sqrt{\lambda_k} - \omega)$$

This gives a LORETA estimate that is additionally constrained to respect the dispersion relation — a physically motivated non-separable regularizer that is trivial to implement in the eigenmode-direct framework.

---

## 6. Key references

- **MNE:** Hämäläinen & Ilmoniemi, "Interpreting magnetic fields of the brain." Med. Biol. Eng. Comput., 1994.
- **wMNE:** Lin et al., "Assessing and improving the spatial accuracy in MEG source localization." NeuroImage, 2006.
- **LORETA:** Pascual-Marqui et al., "Low resolution electromagnetic tomography." Int. J. Psychophysiol., 1994.
- **sLORETA:** Pascual-Marqui, "Standardized low-resolution brain electromagnetic tomography." Meth. Find. Exp. Clin. Pharmacol., 2002.
- **eLORETA:** Pascual-Marqui, "Discrete, 3D distributed, linear imaging methods of electric neuronal activity." arXiv, 2007.
- **dSPM:** Dale et al., "Dynamic statistical parametric mapping." Neuron, 2000.
- **LCMV:** Van Veen et al., "Localization of brain electrical activity via linearly constrained minimum variance spatial filtering." IEEE Trans. Biomed. Eng., 1997.
- **Empirical Bayes / MNE:** Wipf & Nagarajan, "A unified Bayesian framework for MEG/EEG source imaging." NeuroImage, 2009.