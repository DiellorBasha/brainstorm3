
All projects
cortical-flow
In this project, we will develop the Global Cortical-Flow framework (GCF; Phase 1) to track how brain activity propagates across the cortical surface. GCF estimates velocity fields from neurophysiological activity across frequencies and spatial scales on the folded cortex. By quantifying multiple features of cortical dynamics, GCF provides descriptors of brain function that considerably expand the neuroscience toolkit. A key advance is that these measurements support short-term prediction: given an observed pattern, we can estimate its expected displacement over short time windows along the cortical sheet, enabling testable predictions of propagation. This closes the loop between description and prediction, enabling new tests of expected propagation along the functional hierarchy of the cortex across the healthy lifespan (Phase 2) and sensitive detection of abnormal propagation patterns in disease (Phase 3). We hypothesize that cortical activity can be modelled across spatial scales and neurophysiological frequencies as a velocity field, with recurrent source-sink motifs, propagation speed, and trajectories serving as core descriptors of brain dynamics. These descriptors are expected to capture healthy developmental and aging variants, while deviations may yield clinically meaningful markers of disease, including abnormal propagation signatures of epilepsy, including between seizures. Overarching Goal & Specific Aims: We aim to develop and disseminate a scalable and validated cortical-flow framework for dynamic brain mapping, through three Phases: Phase 1: Advance the Global Cortical-Flow framework. We will build and validate methodology for cortical-flow mapping, including joint scale-frequency representations and event-based propagation summaries, to derive robust descriptors of cortical activity propagation that are comparable across individuals and ready for population and clinical applications. All these advances will be implemented and widely disseminated through our established open-source software (Brainstorm7 ). Phase 2: Characterize normative cortical flow across the lifespan. Using our large Lifespan Cohort (n≈1,700; ages 4-88), we will chart the normative distribution of cortical-flow descriptors across development and healthy aging. Analyses will be sex-stratified to detect sex-specific patterns. We will test whether slow spontaneous activity (δ-α, 2-15 Hz) propagates preferentially bottom-up along the cortical functional hierarchy, whereas faster β activity (15-35 Hz) flows top-down in the reverse direction, with modulations in aging. We will also assess whether propagation trajectories and synchronized source-sink motifs provide proxy measures of brain network connectivity. The outcome will be the first population-scale atlas of cortical-flow, openly shared through Brainstorm and neuromaps8 (Summary of Progress) to benchmark discovery and clinical applications. 2 Phase 3: Demonstrate clinical value. We will test whether cortical-flow source-sink motifs and propagation pathways localize epileptogenic regions and improve presurgical hypotheses compared to current standards. This work will leverage the normative atlas from Phase 2 to generate patient-specific deviation maps and will be benchmarked directly against standard clinical workflows. We will leverage two well-characterized cohorts of epilepsy patients with longitudinal follow-up: one pediatric (n=121) and one adult (n=300), already collected at two sites to enable generalizability testing. This phase will evaluate cortical flow as a clinically meaningful tool for surgical planning, with potential to advance research and inform intervention for other conditions such as Alzheimer’s and Parkinson’s disease, where early brain dynamic changes can predict clinical evolution
Show more



How can I help you today?



Start a task in Cowork
Wang GBF paper benchmarking and seizure datasets
Last message 27 minutes ago
MEG source mapping with eigenmodes and harmonic basis functions
Last message 39 minutes ago
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


diffusionnet-cortical-applications.md
15.84 KB •235 lines
Formatting may be inconsistent from source

# DiffusionNet: Architecture and Applications to Cortical Source Imaging

> Reference document for the cortical-flow project.
> Based on Sharp, Attaiki, Crane & Ovsjanikov (2022), *ACM Trans. Graph.*
> [Paper](https://nmwsharp.com/media/papers/diffusion-net/DiffusionNet.pdf) · [Code](https://github.com/nmwsharp/diffusion-net)

---

## 1. Core idea

DiffusionNet replaces the fragile, mesh-dependent convolution operations used in geometric deep learning with a single robust primitive: **learned heat diffusion**. Because diffusion is defined by the Laplace–Beltrami operator (LBO) — an intrinsic geometric object — the resulting network automatically generalizes across different meshes, resolutions, and even representations (triangle meshes, point clouds). You can train on one discretization and evaluate on a completely different one.

This property is called **discretization agnosticism**, and it is precisely what is needed for cross-subject cortical analysis, where each brain has a different mesh.

---

## 2. Architecture

A DiffusionNet is a stack of identical blocks (typically 4). Each block has three components applied in sequence:

### 2.1 Learned diffusion layer

The heat equation on a surface S:

$$\frac{\partial u}{\partial t} = \Delta_S\, u$$

evolves a scalar field u by smoothing it over the surface. The **heat operator** for time t is:

$$H_t(u) = \exp(t\,\Delta_S)\, u$$

In the LBO eigenbasis {φ_k} with eigenvalues {λ_k}, this is diagonal:

$$H_t(u) = \sum_k \exp(-\lambda_k\, t)\;\langle u, \phi_k \rangle\; \phi_k$$

Each spectral coefficient is simply multiplied by exp(−λ_k t). Small t leaves high-frequency modes intact (local support); large t suppresses them (global support); t = 0 is the identity; t → ∞ approaches the surface mean.

**What is learned:** A separate diffusion time t_c ≥ 0 for each feature channel c. The network discovers its own receptive field per feature, per layer — from purely local to fully global — as a continuous, differentiable parameter. No manual neighborhood size selection, no pooling hierarchy.

**Computational cost:** Precompute K eigenpairs of the LBO once per mesh. Then each diffusion operation is O(K·V) — project to eigenbasis, pointwise multiply by exp(−λ_k t_c), project back. K = 128 is typical and sufficient.

### 2.2 Spatial gradient features

Diffusion alone produces **isotropic** (radially symmetric) spatial communication. To enable directional sensitivity, DiffusionNet computes the intrinsic surface gradient ∇u of each feature and forms pairwise inner products:

$$g_{ij}(x) = \langle \nabla u_i(x),\; \nabla u_j(x) \rangle$$

These inner products are invariant to the choice of tangent frame (no gauge ambiguity), but encode directional relationships between features — enabling edge detection, orientation-sensitive filters, and anisotropic pattern recognition without ever defining a canonical coordinate system on the surface.

The gradient is computed via standard DEC or FEM gradient operators on the mesh.

### 2.3 Pointwise MLP

A standard multilayer perceptron applied independently at every vertex, with shared weights. Transforms the concatenation of diffused features and gradient features. Handles all channel mixing and nonlinear function approximation, but provides zero spatial communication — that role belongs entirely to the diffusion layer.

### 2.4 Block assembly

One DiffusionNet block:

```
input features (V × D)
    │
    ├──→ learned diffusion (per-channel t_c) ──→ diffused features
    │
    ├──→ gradient ∇u_i, form ⟨∇u_i, ∇u_j⟩  ──→ gradient features
    │
    └──→ concatenate [diffused; gradient]
              │
              ▼
         pointwise MLP (shared weights across vertices)
              │
              ▼
       output features (V × D)
```

Stack 4 blocks. Add a task-specific head (per-vertex softmax for segmentation, global mean-pool + MLP for classification, feature extraction for correspondence).

### 2.5 Theoretical expressivity

**Lemma (Sharp et al. 2022):** Radially symmetric convolutions are contained in the function space defined by diffusion followed by a pointwise map. Specifically, for a signal u on ℝ², the integral over the r-sphere at any point can be recovered from the diffusion history u_t via a pointwise transform. Therefore, convolution with any radial kernel α(r) can be written as a pointwise operation on diffused values — diffusion + pointwise MLP is at least as expressive as geodesic radial convolution.

With gradient features added, the network goes beyond radial symmetry and can learn arbitrary orientation-sensitive filters.

---

## 3. Why it matters for cortical analysis

### 3.1 The cross-subject mesh problem

Every cortical surface reconstruction (FreeSurfer, CAT12, FastSurfer) produces a different mesh per subject — different vertex count (100k–300k), different triangulation, different vertex placement. Standard group analysis handles this by registering all subjects to a common template (fsaverage) and resampling data onto the template mesh. This works but:

- Introduces interpolation error during resampling
- Destroys subject-specific geometric detail
- Couples analysis quality to registration quality
- Makes eigenmode-based analysis problematic (template eigenmodes ≠ subject eigenmodes)

DiffusionNet sidesteps this entirely. Because its operations are defined by the LBO and gradient operator — intrinsic to each surface — **the same trained network applies to any cortical mesh without registration or resampling.** The learned diffusion times and MLP weights parameterize geometric operations (heat flow at a given spatial scale, nonlinear feature transforms) that have the same meaning on every cortex.

### 3.2 Precomputation requirements

Per subject mesh, compute once and cache:

- LBO eigenpairs (φ_k, λ_k) for k = 1..K (K ≈ 128–200)
- Mass matrix M (diagonal, vertex areas)
- Gradient operator G (sparse, from DEC or FEM)

These are already computed in the cortical-flow preprocessing pipeline (geometry-central C++ / eigensolver.cpp). No additional infrastructure needed.

---

## 4. Applications to cortical MEG/EEG source imaging

### 4.1 Spatial filtering front-end for the doubly-spectral framework

**Connection:** DiffusionNet's diffusion layer implements exp(−λ_k t) — a pure lowpass filter in eigenmode space. This is a special case of the spatial filter families in the doubly-spectral framework, corresponding to the **overdamped / pure-diffusion limit** of the wave-propagation kernels.

**Extension — WaveNet variant:** Replace the diffusion transfer function with the wave-propagation kernel:

$$\text{Diffusion:} \quad T_k = \exp(-\lambda_k\, t)$$

$$\text{Wave kernel:} \quad T_k = \exp(-\gamma\, t) \cos(c\sqrt{\lambda_k}\, t)$$

Instead of learning one parameter t per channel, learn (c, γ, t) per channel — wave speed, damping, and time. This gives access to **bandpass spatial filters** (selecting a specific spatial frequency band) and **oscillatory spatial patterns** that pure heat diffusion cannot capture.

The spectral implementation is identical in structure: multiply each eigenmode coefficient by the transfer function T_k, then reconstruct. Everything else (gradient features, MLPs, block stacking) is unchanged.

This creates a geometry-aware CNN whose spatial filters are physically interpretable as damped wave propagation on the cortical surface — directly aligned with the physics of cortical dynamics.

### 4.2 Source classification and decoding on native geometry

**Task:** Given MEG/EEG source-estimated time series on a subject-specific cortical mesh, classify brain states, decode stimuli, or detect events.

**Input:** Per-vertex scalar features at each time point (source amplitudes, spectral power in frequency bands, eigenmode coefficients). Shape: (V × D) where D = number of features.

**Architecture:** Standard DiffusionNet with per-vertex output (segmentation head for spatial ROI identification) or global pooling (classification head for brain state decoding).

**Key advantage:** Train across subjects on their native meshes. A training set can include subjects with 150k and 300k vertices, high- and low-resolution reconstructions, without any preprocessing to harmonize meshes. The geometric invariance handles cross-subject generalization.

**Comparison to standard approach:** Conventional decoding projects all subjects to fsaverage, extracts features at template vertices, and feeds them to a standard classifier. This loses subject-specific geometric structure and introduces registration-dependent noise. DiffusionNet preserves the native geometry while still supporting multi-subject training.

### 4.3 Cortical parcellation without registration

**Task:** Segment the cortical surface into anatomical or functional regions.

**Input:** Per-vertex geometric and functional features (curvature, sulcal depth, cortical thickness, myelin maps, resting-state connectivity profiles).

**Architecture:** DiffusionNet with per-vertex softmax output over region labels.

**Advantage over atlas-based methods:** Desikan-Killiany, Destrieux, and similar parcellations rely on FreeSurfer's registration to transfer labels from a template. DiffusionNet learns the segmentation directly from geometric and functional features, adapting to each subject's individual anatomy. The learned diffusion times allow the network to integrate context at appropriate spatial scales — local features (curvature) at short diffusion times, global context (which lobe am I in?) at long diffusion times.

### 4.4 Learning the wavelet tensor / constitutive relation

**Connection to the connection Laplacian framework:** The wavelet tensor describes the spatially varying gain of active cortex — the material property that relates input (connectivity-mediated drive) to output (local oscillatory response). The resting state is the standardized measurement protocol for this constitutive relation.

**Approach:** Frame the estimation of the wavelet tensor as a learning problem. Given resting-state source data on native cortical geometry, learn per-vertex parameters that predict the observed spatiotemporal dynamics. DiffusionNet (or the wave-kernel extension) provides the architecture. The learned parameters at each vertex become the estimate of the constitutive relation.

**Cross-subject structure:** Because DiffusionNet generalizes across meshes, you can train across subjects and extract a population-level constitutive model while respecting individual geometry. Subject-specific deviations from the population model become interpretable as individual differences in cortical material properties.

### 4.5 Functional correspondence via learned features

**Connection to functional maps:** DiffusionNet was co-authored with Ovsjanikov precisely because its learned features can drive functional map computation. For cortical surfaces, this enables data-driven correspondence as an alternative to atlas registration.

**Approach:** Train DiffusionNet to produce per-vertex feature descriptors that are consistent across subjects (e.g., via a contrastive loss or a functional map loss). These learned descriptors define a functional map C between any pair of subjects — the C matrix from the cross-subject eigenmode discussion, but computed from data-driven features rather than geometric registration alone.

**Implications for the doubly-spectral framework:** The functional map C derived from DiffusionNet features would align per-subject eigenmode bases, enabling rigorous cross-subject comparison in the λ–ω domain without relying on FreeSurfer registration.

---

## 5. Relationship to existing cortical-flow components

| cortical-flow component | DiffusionNet connection |
|---|---|
| LBO eigensolver (eigensolver.cpp) | Provides the precomputed eigenpairs (φ_k, λ_k) that DiffusionNet uses for spectral acceleration of diffusion |
| DEC operators (dec_operators.cpp) | Provide the gradient operator G used for DiffusionNet's spatial gradient features |
| Helmholtz decomposition (helmholtz.cpp) | Could provide irrotational/solenoidal features as input channels to DiffusionNet |
| Eigenmode-direct forward model (doubly-spectral framework) | DiffusionNet's diffusion layer is the heat-kernel special case of the general spatial filter family; wave-kernel extension generalizes it |
| Connection Laplacian | Extension to vector-valued DiffusionNet using connection Laplacian eigenmodes instead of scalar LBO eigenmodes |
| Neurophysiological dispersion bound (ω ≤ c√λ) | Could be enforced as a hard constraint on learned wave-kernel parameters (c, γ, t) to ensure physical plausibility |

---

## 6. Implementation notes

### 6.1 Reference implementation

The [official PyTorch implementation](https://github.com/nmwsharp/diffusion-net) is clean and well-documented. Core files:

- `diffusion_net/layers.py` — DiffusionBlock, LearnedTimeDiffusion, SpatialGradientFeatures
- `diffusion_net/geometry.py` — LBO eigenpair computation, gradient operator assembly
- `diffusion_net/utils.py` — mesh loading, preprocessing

### 6.2 Integration with cortical-flow

Two integration paths:

**Python path (prototyping):** Use the reference implementation directly. Load cortical meshes, precompute LBO eigenpairs and gradient operators using the existing geometry.py utilities or import from cortical-flow's precomputed Zarr stores. Define task-specific heads and train.

**C++ / WASM path (production):** The eigendecomposition and gradient operator assembly already exist in geometry-central C++ (nxr-compute). The learned network weights (diffusion times, MLP parameters) can be exported from PyTorch and evaluated in C++ — the forward pass is just matrix multiplications and per-element exponentials, no autograd needed.

### 6.3 Wave-kernel extension

Minimal modification to the reference code. In `layers.py`, replace:

```python
# Original: heat kernel
diffusion_coefs = torch.exp(-evals * t)

# Wave kernel extension: learn (c, gamma, t) per channel
diffusion_coefs = torch.exp(-gamma * t) * torch.cos(c * torch.sqrt(evals) * t)
```

The eigenvalue array `evals` and eigenvector matrix `evecs` are unchanged. Only the transfer function applied to spectral coefficients changes.

### 6.4 Connection Laplacian extension

For vector-valued data (tangent vector fields on the cortical surface), replace the scalar LBO eigenbasis with the connection Laplacian eigenbasis. The diffusion layer becomes:

$$H_t(\mathbf{v}) = \sum_k \exp(-\mu_k\, t)\;\langle \mathbf{v}, \psi_k \rangle\; \psi_k$$

where (ψ_k, μ_k) are eigenpairs of the connection Laplacian. This naturally handles the gauge ambiguity of tangent vectors — the connection Laplacian eigenmodes define a globally consistent basis for tangent vector fields, avoiding the sign/phase problems that plague vertex-wise comparison of vector data.

This extension would enable DiffusionNet to process cortical current flow vectors (not just scalar source amplitudes), with learned spatial communication that respects the tangent bundle geometry.

---

## 7. Key references

- **DiffusionNet:** Sharp, Attaiki, Crane & Ovsjanikov. "DiffusionNet: Discretization Agnostic Learning on Surfaces." ACM Trans. Graph. 41(4), 2022.
- **Functional maps:** Ovsjanikov, Ben-Chen, Solomon, Butscher & Guibas. "Functional Maps: A Flexible Representation of Maps Between Shapes." ACM Trans. Graph. 31(4), 2012.
- **Vector Heat Method:** Sharp, Soliman & Crane. "The Vector Heat Method." ACM Trans. Graph. 38(3), 2019.
- **Nonmanifold Laplacian:** Sharp & Crane. "A Laplacian for Nonmanifold Triangle Meshes." SGP, 2020.
- **Intrinsic triangulations:** Sharp, Soliman & Crane. "Navigating Intrinsic Triangulations." ACM Trans. Graph. 38(4), 2019.
- **CEPS:** Gillespie, Springborn & Crane. "Discrete Conformal Equivalence of Polyhedral Surfaces." ACM Trans. Graph. 40(4), 2021.
- **Möbius Registration:** Baden, Crane & Kazhdan. "Möbius Registration." SGP, 2018.
- **Dirac operator:** Liu, Jacobson & Crane. "A Dirac Operator for Extrinsic Shape Analysis." SGP, 2017.