
All projects
cortical-flow
In this project, we will develop the Global Cortical-Flow framework (GCF; Phase 1) to track how brain activity propagates across the cortical surface. GCF estimates velocity fields from neurophysiological activity across frequencies and spatial scales on the folded cortex. By quantifying multiple features of cortical dynamics, GCF provides descriptors of brain function that considerably expand the neuroscience toolkit. A key advance is that these measurements support short-term prediction: given an observed pattern, we can estimate its expected displacement over short time windows along the cortical sheet, enabling testable predictions of propagation. This closes the loop between description and prediction, enabling new tests of expected propagation along the functional hierarchy of the cortex across the healthy lifespan (Phase 2) and sensitive detection of abnormal propagation patterns in disease (Phase 3). We hypothesize that cortical activity can be modelled across spatial scales and neurophysiological frequencies as a velocity field, with recurrent source-sink motifs, propagation speed, and trajectories serving as core descriptors of brain dynamics. These descriptors are expected to capture healthy developmental and aging variants, while deviations may yield clinically meaningful markers of disease, including abnormal propagation signatures of epilepsy, including between seizures. Overarching Goal & Specific Aims: We aim to develop and disseminate a scalable and validated cortical-flow framework for dynamic brain mapping, through three Phases: Phase 1: Advance the Global Cortical-Flow framework. We will build and validate methodology for cortical-flow mapping, including joint scale-frequency representations and event-based propagation summaries, to derive robust descriptors of cortical activity propagation that are comparable across individuals and ready for population and clinical applications. All these advances will be implemented and widely disseminated through our established open-source software (Brainstorm7 ). Phase 2: Characterize normative cortical flow across the lifespan. Using our large Lifespan Cohort (n≈1,700; ages 4-88), we will chart the normative distribution of cortical-flow descriptors across development and healthy aging. Analyses will be sex-stratified to detect sex-specific patterns. We will test whether slow spontaneous activity (δ-α, 2-15 Hz) propagates preferentially bottom-up along the cortical functional hierarchy, whereas faster β activity (15-35 Hz) flows top-down in the reverse direction, with modulations in aging. We will also assess whether propagation trajectories and synchronized source-sink motifs provide proxy measures of brain network connectivity. The outcome will be the first population-scale atlas of cortical-flow, openly shared through Brainstorm and neuromaps8 (Summary of Progress) to benchmark discovery and clinical applications. 2 Phase 3: Demonstrate clinical value. We will test whether cortical-flow source-sink motifs and propagation pathways localize epileptogenic regions and improve presurgical hypotheses compared to current standards. This work will leverage the normative atlas from Phase 2 to generate patient-specific deviation maps and will be benchmarked directly against standard clinical workflows. We will leverage two well-characterized cohorts of epilepsy patients with longitudinal follow-up: one pediatric (n=121) and one adult (n=300), already collected at two sites to enable generalizability testing. This phase will evaluate cortical flow as a clinically meaningful tool for surgical planning, with potential to advance research and inform intervention for other conditions such as Alzheimer’s and Parkinson’s disease, where early brain dynamic changes can predict clinical evolution
Show more



How can I help you today?



Start a task in Cowork
Wang GBF paper benchmarking and seizure datasets
Last message 27 minutes ago
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


vector_eigenmode_meg_proposal.md
38.29 KB •499 lines
Formatting may be inconsistent from source

# Vector Eigenmode Source Mapping: A Differential Geometric Framework for Orientation-Consistent MEG/EEG Source Imaging and Cortical Wave Detection

---

## Summary

Magnetoencephalography (MEG) and electroencephalography (EEG) source imaging remains fundamentally limited by the ill-posed nature of the electromagnetic inverse problem and by geometric inconsistencies introduced by the folded cortical surface. Recent work has demonstrated that expressing cortical sources as linear combinations of Laplace-Beltrami Operator (LBO) eigenmodes — Geometric Basis Functions (GBFs) — substantially improves source localization accuracy and biological interpretability (Wang et al., 2026). However, GBF and all existing eigenmode-based approaches model cortical sources as **scalar fields**, projecting the neural current vector onto local surface normals independently at each vertex. This convention introduces systematic sign inconsistencies at opposite sulcal walls, where antiparallel surface normals cause geometrically identical neural activity to appear with opposite sign in the source estimate. Downstream phase-based wave detection — such as optical flow on instantaneous phase fields — then confounds these geometric sign artifacts with true neural propagation, generating spurious traveling wave signatures across sulci.

We propose a principled extension of the GBF framework grounded in differential geometry: a **Vector Eigenmode Source Mapping** approach that models cortical neural currents as vector fields on the cortical manifold, establishes globally consistent orientation frames using trivial connections and the vector heat method (Crane et al., 2010; Sharp et al., 2019), decomposes the current field via the Hodge decomposition into its MEG-visible (irrotational) and MEG-invisible (solenoidal) components, and expands the result in eigenmodes of the **connection Laplacian** — the natural generalization of the scalar LBO to vector-valued fields. The scalar potential extracted from the Hodge decomposition provides an orientation-invariant amplitude field whose instantaneous phase is free of sulcal sign artifacts. Combined with a joint spatiotemporal spectral analysis in the (λ, ω) plane — where λ indexes spatial frequency via LBO eigenvalues and ω indexes temporal frequency via wavelet scale — this framework enables principled discrimination between wave propagation, diffusive spreading, and reaction-diffusion dynamics on the cortical surface, grounded in both Maxwell's equations and known physiological constraints on neural propagation speed.

---

## 1. Introduction

### 1.1 The need for high-fidelity spatiotemporal brain imaging

Brain waves — spatially organized oscillations propagating across cortical networks — underlie fundamental cognitive functions including attention, memory consolidation, sensory processing, and consciousness (Vidaurre et al., 2018; Raut et al., 2021). Understanding these dynamics requires imaging methods with simultaneously high spatial and temporal resolution across the whole brain. MEG and EEG uniquely provide millisecond temporal resolution with whole-brain coverage, making them the primary tools for studying fast electrophysiological dynamics. However, the spatial resolution of sensor-level MEG/EEG is fundamentally limited by the number and geometry of sensors and by the ill-posed nature of the electromagnetic inverse problem (He et al., 2018).

EEG/MEG source imaging (ESI) addresses this limitation by reconstructing the cortical current distribution from sensor measurements using computational forward and inverse modeling. The forward problem — computing the sensor signal from a known cortical current distribution — is uniquely solved by Maxwell's equations given the head geometry (Baillet, 2017). The inverse problem — recovering the cortical current from sensor measurements — is massively underdetermined, with ~10,000 cortical source locations but typically only 64–300 sensors, requiring the introduction of prior information to constrain the solution (He et al., 2018).

### 1.2 Geometric priors and the GBF framework

The geometry of the cortical surface fundamentally constrains the spatial organization of neural dynamics. Cortical folding defines propagation pathways, and cortico-cortical connectivity follows the intrinsic geometry of the cortical manifold. Recent work has demonstrated that whole-brain fMRI dynamics can be efficiently described by the eigenmodes of the Laplace-Beltrami operator (LBO) on the cortical surface — the Geometric Basis Functions — with as few as 200 modes capturing the majority of spatial variance (Pang et al., 2023).

Building on this insight, Wang et al. (2026) introduced the GBF framework for EEG/MEG source imaging, expressing cortical sources as linear combinations of patient-specific LBO eigenmodes:

$$x(t) = \sum_{i} \theta_i(t) \psi_i$$

where $\psi_i$ are the LBO eigenmodes and $\theta_i(t)$ are time-varying scalar coefficients. A MAP inverse with a logarithmic spectral prior $\Sigma^{-1} = \text{diag}[-\beta / \log(\lambda_i)]$ yields a closed-form solution that outperforms MNE, dSPM, sLORETA, eLORETA, LCMV, and wMNE across synthetic benchmarks, task-evoked EEG, resting-state MEG functional connectivity, intracranial stimulation, and epilepsy localization.

GBF represents a significant advance. However, it treats cortical sources as scalar fields — projecting the neural current vector onto local surface normals — and performs wave detection post-hoc using optical flow on scalar phase fields computed in vertex space. Both of these design choices introduce geometric inconsistencies that limit the framework's fidelity for wave propagation analysis.

### 1.3 The orientation problem: a fundamental geometric inconsistency

MEG and EEG are sensitive to current dipoles whose orientation has a tangential component relative to the sensor array. In standard source modeling, dipole orientations are constrained to the local surface normal $\hat{n}(x)$ at each cortical vertex. This convention is geometrically correct locally — it reflects the columnar organization of cortical neurons perpendicular to the cortical surface — but is globally inconsistent.

On the folded cortical surface, opposite walls of a sulcus have antiparallel surface normals by the geometry of folding. Two neural populations with identical activation patterns — same current polarity, same magnitude, same temporal dynamics — but located on opposite sulcal walls will produce source estimates of opposite sign. This is not a failure of the inverse method; it is an unavoidable consequence of using locally defined surface normals as the orientation convention without enforcing global consistency.

The consequences propagate downstream:

1. **Phase analysis**: The Hilbert transform of a sign-flipped signal produces a $\pi$ phase shift. Two vertices with identical dynamics on opposite sulcal walls appear to be in antiphase — a spurious spatial phase discontinuity.

2. **Optical flow wave detection**: Phase gradients computed across sulcal walls generate spurious velocity vectors pointing from positive to negative sulcal walls — artifacts that appear as traveling waves propagating across sulci but reflect geometry, not neural dynamics.

3. **Eigenmode projection**: Sign flips are spatially abrupt, introducing high spatial frequency content into the eigenmode spectrum and contaminating high-$k$ modes with geometric artifacts.

4. **Functional connectivity**: Correlation analyses between opposite sulcal walls are negatively biased, potentially suppressing true functional connectivity between adjacent cortical areas.

This problem is not addressed by GBF or any existing eigenmode-based ESI framework.

---

## 2. Problem Statement

### 2.1 Formal statement

Let $S$ be the cortical surface, a compact two-dimensional Riemannian manifold embedded in $\mathbb{R}^3$ with metric $g$. At each point $x \in S$, the neural current is a vector $\mathbf{J}(x,t) \in \mathbb{R}^3$. Standard ESI methods project $\mathbf{J}$ onto the local surface normal $\hat{n}(x)$, obtaining a scalar source amplitude:

$$s(x,t) = \mathbf{J}(x,t) \cdot \hat{n}(x)$$

The surface normal $\hat{n}(x)$ is defined only up to a global sign — the choice of inward vs. outward pointing is locally arbitrary and globally inconsistent on a folded surface. On opposite sulcal walls, $\hat{n}(x_1) \approx -\hat{n}(x_2)$ for geometrically adjacent points $x_1, x_2$. Hence even if $\mathbf{J}(x_1,t) \approx \mathbf{J}(x_2,t)$, we have $s(x_1,t) \approx -s(x_2,t)$.

Phase-based wave detection computes:

$$\phi(x,t) = \arg\bigl(s(x,t) + i\,\mathcal{H}[s(x,t)]\bigr)$$

where $\mathcal{H}$ is the Hilbert transform. The phase gradient $\nabla\phi$ used for optical flow then contains contributions from both genuine neural propagation and spurious $\pi$ phase jumps at sulcal walls. These cannot be separated without knowledge of the global orientation consistency of $\hat{n}(x)$.

### 2.2 Requirements for a principled solution

A principled solution must:

1. Establish a **globally consistent orientation frame** on the cortical surface, assigning to each vertex a reference direction that varies smoothly across the entire manifold including across sulcal walls.

2. Express the neural current vector $\mathbf{J}(x,t)$ in this globally consistent frame, eliminating locally arbitrary sign conventions.

3. Decompose $\mathbf{J}$ into its **MEG-visible** (irrotational, dipolar) and **MEG-invisible** (solenoidal, closed-loop) components — identifying which component carries measurable electromagnetic signal.

4. Perform phase and wave analysis on a scalar field derived from the MEG-visible component that is **orientation-invariant** — whose sign does not depend on the local choice of surface normal.

5. Extend the eigenmode basis from scalar LBO eigenmodes to **vector eigenmodes** that naturally encode the directional structure of cortical current flow.

---

## 3. Proposed Method

### 3.1 Overview

The proposed framework has five stages, each grounded in a specific differential geometric construction:

```
Stage 1: Cortical geometry
         Mesh + metric → LBO eigenmodes {ψₖ, λₖ}
                       → Connection Laplacian eigenmodes {Ψₖ, μₖ}
         
Stage 2: Globally consistent frame
         Trivial connections + vector heat method
         → Frame field {e₁(x), e₂(x), n̂_consistent(x)}

Stage 3: Vector inverse problem
         Free-orientation MEG inverse
         → Vector source field J(x,t) in consistent frame

Stage 4: Hodge decomposition
         J = ∇f + ∇×ψ + h
         → Scalar potential f(x,t) [MEG-visible, orientation-invariant]

Stage 5: Spatiotemporal spectral analysis
         Morlet CWT of vector eigenmode coefficients
         → Joint (λ, ω, t) tensor
         → Dispersion relation, wave vs. diffusion discrimination
```

### 3.2 Stage 1: Geometric infrastructure

#### 3.2.1 Scalar LBO eigenmodes

The Laplace-Beltrami operator on the cortical surface $S$ is discretized using cotangent weights on the triangular mesh:

$$(\Delta_S f)_i = \frac{1}{2A_i} \sum_{j \in \mathcal{N}(i)} (\cot\alpha_{ij} + \cot\beta_{ij})(f_j - f_i)$$

where $A_i$ is the vertex area, $\alpha_{ij}$ and $\beta_{ij}$ are the angles opposite edge $ij$ in the two incident triangles, and $\mathcal{N}(i)$ is the one-ring neighborhood of vertex $i$. The generalized eigenvalue problem:

$$L\Psi = \Lambda M\Psi$$

where $L$ is the stiffness matrix, $M$ is the mass matrix, and $\Lambda = \text{diag}(\lambda_0, \lambda_1, \ldots, \lambda_K)$, yields scalar eigenmodes $\psi_k : S \to \mathbb{R}$ ordered by spatial frequency. The first $K \approx 200$–300 modes form the spatial basis, as validated by Wang et al. (2026) and Pang et al. (2023).

#### 3.2.2 The connection Laplacian and vector eigenmodes

The scalar LBO operates on functions $f : S \to \mathbb{R}$. Its natural generalization to vector fields is the **connection Laplacian** $\Delta_\nabla$, which operates on sections of the tangent bundle $TS$:

$$\Delta_\nabla \mathbf{V} = \text{div}(\nabla \mathbf{V})$$

where $\nabla$ is the covariant derivative associated with the Levi-Civita connection on $(S, g)$. The connection Laplacian has its own eigenvalue problem:

$$\Delta_\nabla \mathbf{\Psi}_k = -\mu_k \mathbf{\Psi}_k$$

where $\mathbf{\Psi}_k : S \to TS$ are **vector eigenmodes** — smoothly varying vector fields on the cortical surface — and $\mu_k \geq 0$ are the corresponding eigenvalues. On a flat surface, $\mu_k = \lambda_k$ exactly; on the curved, folded cortex, $\mu_k \neq \lambda_k$, with the deviation encoding curvature-dependent coupling between spatial scale and vector field direction.

The vector eigenmodes $\{\mathbf{\Psi}_k\}$ form a complete orthonormal basis for vector fields on $S$:

$$\mathbf{J}(x,t) = \sum_k c_k(t) \, \mathbf{\Psi}_k(x)$$

where $c_k(t) = \langle \mathbf{J}(\cdot,t), \mathbf{\Psi}_k \rangle_{L^2(TS)}$ are the vector eigenmode coefficients.

#### 3.2.3 Discrete implementation

The connection Laplacian is discretized as a complex-valued matrix using the formulation of Crane et al. (2010) and Sharp et al. (2019). For each directed edge $(i,j)$ in the mesh, the parallel transport angle $\phi_{ij}$ — the rotation needed to align the reference frame at vertex $i$ with the frame at vertex $j$ when transporting along edge $(i,j)$ — is computed from the mesh geometry. The discrete connection Laplacian is then:

$$(\tilde{L})_{ij} = \begin{cases} -\frac{1}{2}(\cot\alpha_{ij} + \cot\beta_{ij}) e^{i\phi_{ij}} & j \in \mathcal{N}(i) \\ \sum_{k \in \mathcal{N}(i)} \frac{1}{2}(\cot\alpha_{ik} + \cot\beta_{ik}) & j = i \end{cases}$$

The complex representation encodes 2D tangent vectors as complex numbers — the real and imaginary parts correspond to the two tangential directions in the local frame. This formulation is implemented in **Geometry Central** (Sharp et al., 2019) and is the computational foundation for both the trivial connections and vector heat method computations.

### 3.3 Stage 2: Globally consistent orientation frame via trivial connections

#### 3.3.1 The trivial connection problem

Given the mesh with its intrinsic geometry, we seek a connection $\nabla$ on the tangent bundle $TS$ that is as flat as possible — minimizing the total curvature (holonomy) accumulated when transporting vectors around loops. A connection on a surface is specified by an angle $\phi_{ij}$ for each directed edge $(i,j)$, representing the rotation of the reference frame when transported from $i$ to $j$.

The trivial connection minimizes:

$$E(\{\phi_{ij}\}) = \sum_{\text{faces } f} \left(\sum_{(i,j) \in \partial f} \phi_{ij}\right)^2$$

subject to the topological constraint that the total index of the connection equals the Euler characteristic $\chi(S)$ (Poincaré-Hopf theorem). For the cortical surface, $\chi(S) \approx 2$ (sphere-like topology). The unavoidable holonomy — imposed by the Gauss-Bonnet theorem — is distributed as uniformly as possible across the surface.

This optimization, solved as a linear system by Crane et al. (2010), yields globally consistent parallel transport angles $\{\phi_{ij}^*\}$ and a corresponding globally consistent frame field $\{e_1(x), e_2(x)\}$ at every vertex.

#### 3.3.2 Vector heat method for parallel transport

Given a source vector $\mathbf{v}_0$ at a source vertex $x_0$, the vector heat method (Sharp et al., 2019) computes its parallel transport to all other vertices via a three-step procedure:

**Step 1 — Vector diffusion**: Solve the vector heat equation with the trivial connection:

$$\frac{\partial \mathbf{V}}{\partial t} = \Delta_\nabla \mathbf{V}, \quad \mathbf{V}(x_0, 0) = \mathbf{v}_0, \quad \mathbf{V}(x \neq x_0, 0) = \mathbf{0}$$

at a short time $t = h^2$ where $h$ is the mean edge length. This diffuses the vector consistently with the connection, respecting the intrinsic geometry of the surface.

**Step 2 — Scalar divergence**: Compute the divergence of the diffused vector field to extract a scalar potential.

**Step 3 — Normalization**: Normalize the result to recover unit vectors representing the parallel-transported direction field.

The output is a globally consistent vector field over the entire cortical surface — the reference frame for expressing current vectors at all vertices relative to the source orientation.

#### 3.3.3 Consistent normal orientation

For the normal bundle, a globally consistent sign for $\hat{n}(x)$ is determined by solving:

$$\Delta_S f = 0 \quad \text{with boundary conditions at gyral crowns}$$

where gyral crowns — the most MEG-visible locations with well-defined outward normals — provide anchoring boundary conditions. The smooth solution $f(x)$ determines a sign field $s(x) = \text{sign}(f(x))$ that corrects the locally arbitrary normal convention to a globally consistent one, minimizing sign inconsistencies across sulcal walls.

### 3.4 Stage 3: Vector source reconstruction

#### 3.4.1 Free-orientation inverse in the consistent frame

With the globally consistent frame $\{e_1(x), e_2(x), \hat{n}_{\text{consistent}}(x)\}$ established, the neural current at each vertex is expressed as:

$$\mathbf{J}(x,t) = J_1(x,t)\, e_1(x) + J_2(x,t)\, e_2(x) + J_n(x,t)\, \hat{n}(x)$$

The MEG forward model for each component is linear:

$$\mathbf{d}(t) = L_1 \mathbf{J}_1(t) + L_2 \mathbf{J}_2(t) + L_n \mathbf{J}_n(t)$$

where $L_1, L_2, L_n \in \mathbb{R}^{n_\text{ch} \times n_\text{vert}}$ are lead field matrices for each frame component, computable from the same BEM forward model used by standard pipelines.

The inverse problem recovers all three components jointly. In vector eigenmode space, the source model becomes:

$$\mathbf{J}(x,t) = \sum_k c_k(t)\, \mathbf{\Psi}_k(x)$$

The effective lead field in vector eigenmode space:

$$\tilde{L}_{\text{vec}} = L \cdot \Psi \in \mathbb{R}^{n_\text{ch} \times K}$$

where $\Psi$ concatenates the vector eigenmode spatial patterns evaluated at each vertex. The inverse in this basis is:

$$\hat{\mathbf{c}}(t) = \left(\tilde{L}_{\text{vec}}^T \tilde{L}_{\text{vec}} + \Sigma_c^{-1}\right)^{-1} \tilde{L}_{\text{vec}}^T \mathbf{d}(t)$$

The prior covariance $\Sigma_c$ for vector eigenmodes is indexed by $\mu_k$ — the connection Laplacian eigenvalues — using the same logarithmic spectral prior validated by Wang et al. (2026):

$$\Sigma_c^{-1} = \text{diag}\left[-\beta / \log(\mu_k)\right]$$

This extends the GBF prior from scalar to vector eigenmode space in a geometrically principled way.

#### 3.4.2 Composed sensor-to-eigenmode transform

As in the scalar case, the composed transform from sensor data directly to vector eigenmode coefficients avoids materializing the full vertex-space source estimate:

$$A_{\text{vec}} = \Psi^T M_{\text{inv}} \in \mathbb{R}^{K \times n_\text{ch}}$$

where $M_{\text{inv}}$ is the imaging kernel from the vector inverse operator. The vector eigenmode time series:

$$\mathbf{c}(t) = A_{\text{vec}} \cdot \mathbf{d}(t)$$

is computed by a single matrix multiply at each time step, with no intermediate vertex-space representation required.

### 3.5 Stage 4: Hodge decomposition and orientation-invariant phase

#### 3.5.1 Discrete Hodge decomposition

Using the Discrete Exterior Calculus (DEC) framework (Desbrun et al., 2005), the recovered vector field $\mathbf{J}(x,t)$ is decomposed at each time point into three orthogonal components:

$$\mathbf{J} = \underbrace{\nabla f}_{\text{exact}} + \underbrace{\star d \psi}_{\text{coexact}} + \underbrace{\mathbf{h}}_{\text{harmonic}}$$

where:
- $\nabla f = \delta \alpha$ is the **gradient (exact) component** — irrotational, sourced from a scalar potential $f$. This is the **MEG-visible dipolar component**.
- $\star d\psi$ is the **curl (coexact) component** — divergence-free, forming closed current loops. This is **MEG-invisible** — it produces zero net magnetic flux outside a spherical conductor.
- $\mathbf{h}$ is the **harmonic component** — determined by the topology of the surface ($\dim \mathcal{H} = 2g$ for genus $g$). For a sphere-like cortical surface, this is small.

The scalar potential $f(x,t)$ is recovered by solving:

$$\Delta_S f = \text{div}(\mathbf{J})$$

via the LBO — using the same eigenmode infrastructure already computed. The solution is:

$$f(x,t) = \sum_k \frac{\langle \text{div}(\mathbf{J}(\cdot,t)), \psi_k \rangle}{\lambda_k} \psi_k(x)$$

#### 3.5.2 Orientation-invariant phase analysis

The scalar potential $f(x,t)$ is **orientation-invariant**: it is defined by the divergence of the current field, which is independent of the sign convention for the surface normal. Opposite sulcal walls with antiparallel normals but identical neural current both contribute positive divergence at a current source — there is no sign ambiguity.

The analytic signal of $f$ via the Hilbert transform:

$$z_f(x,t) = f(x,t) + i\,\mathcal{H}[f(x,t)]$$

gives an instantaneous phase $\phi_f(x,t) = \arg(z_f(x,t))$ that is free of the $\pi$ phase jumps induced by sulcal orientation artifacts in the standard signed amplitude approach.

Phase-based wave detection — optical flow, phase gradient velocity, inter-regional phase coherence — applied to $\phi_f(x,t)$ reflects genuine neural propagation rather than geometric artifacts.

### 3.6 Stage 5: Joint spatiotemporal spectral analysis

#### 3.6.1 The vector eigenmode wavelet tensor

The vector eigenmode time series $c_k(t)$ — one per connection Laplacian eigenmode — are transformed via the Morlet continuous wavelet transform:

$$W_k(s,t) = \int c_k(\tau) \, \psi_s^*\!\left(\frac{\tau - t}{s}\right) d\tau$$

where $\psi_s$ is the Morlet wavelet at scale $s$ (corresponding to frequency $f = f_0/s$). The result is a complex tensor:

$$W_k(s,t) \in \mathbb{C}, \quad k = 0,\ldots,K,\quad s \in \mathcal{S},\quad t \in [0,T]$$

This is the **primary data structure** of the framework — a three-index complex tensor jointly resolving spatial scale (via $\mu_k$), temporal frequency (via $s$), and time.

#### 3.6.2 The joint (μ, ω) spectrum and dispersion analysis

The time-integrated power spectrum:

$$P(k, f) = \frac{1}{T} \int_0^T |W_k(f,t)|^2 \, dt$$

gives the **joint spatiotemporal power spectrum** in the $(\mu_k, f)$ plane. This is the primary observable for discriminating dynamical regimes:

**Wave propagation** — energy concentrates on a dispersion curve $\omega = c\sqrt{\mu_k}$ for some propagation speed $c$. The curve shape identifies the wave type (acoustic, dispersive, etc.).

**Heat diffusion** — energy fills a **wedge** in the $(\mu_k, \omega)$ plane, with high-$\mu_k$ modes decaying rapidly in time (broad $\omega$ spectrum) and low-$\mu_k$ modes persisting (narrow $\omega$ spectrum near DC).

**Reaction-diffusion** — energy fills the wedge but with a spectral gap at low spatial frequencies, reflecting the stabilizing reaction term.

#### 3.6.3 Physiological speed limit as a spectral mask

From known neurophysiology — cortico-cortical axonal conduction velocities of 3–8 m/s for myelinated association fibers and synaptic time constants of 2–200 ms — a **physiologically accessible region** in the $(\mu_k, \omega)$ plane is defined:

$$\omega \leq c_{\max} \cdot \sqrt{\mu_k}, \quad \omega \leq \omega_{\max} = 1/\tau_{\min}$$

Any power in the joint spectrum outside this region is guaranteed to be non-neural — it requires faster-than-possible propagation — and can be zeroed without any dynamical model assumption. This provides a principled, model-agnostic spectral mask that complements and extends the eigenmode-based regularization.

The empirical propagation speed:

$$c_{\text{empirical}} = \max_{(k,f)\,:\,P(k,f) > P_{\text{noise}}} \frac{2\pi f}{\sqrt{\mu_k}}$$

is a measurable output of the framework — the fastest neural propagation present in the data, verifiable against known physiological bounds.

#### 3.6.4 Noise calibration and SNR-resolved eigenmode spectrum

The empty-room noise covariance $\Sigma_n$ — standard in MEG preprocessing pipelines — provides a principled noise floor for the joint spectrum without requiring the inverse operator, avoiding circularity:

$$N(k, f) = \frac{\tilde{L}_k^T \, \Sigma_n(f) \, \tilde{L}_k}{\|\tilde{L}_k\|^2}$$

where $\Sigma_n(f)$ is the frequency-resolved noise cross-spectral density and $\tilde{L}_k = L\Psi_k$ is the sensor topography of vector eigenmode $k$. The SNR-resolved eigenmode power spectrum:

$$\text{SNR}(k, f) = \frac{P(k,f)}{N(k,f)}$$

determines the effective number of reliable eigenmodes $K^*(f)$ at each frequency — a data-adaptive quantity that replaces the fixed truncation used in GBF.

---

## 4. Proposed Validation

### 4.1 Synthetic benchmark: sign-flip artifact quantification

**Objective**: Demonstrate that the proposed framework eliminates sulcal sign artifacts in phase-based wave detection.

**Protocol**:
1. Select vertex pairs on opposite sulcal walls of the central sulcus, Sylvian fissure, and calcarine fissure in the fsaverage5 template.
2. Simulate identical neural activations at both walls — same amplitude, same temporal dynamics, same oscillation frequency (10 Hz).
3. Generate synthetic MEG sensor data via BEM forward model with realistic sensor noise (SNR = 10 dB).
4. Apply GBF (scalar, fixed orientation) and the proposed vector eigenmode framework to reconstruct sources.
5. Compute instantaneous phase at both wall vertices using (a) standard signed amplitude and (b) Hodge scalar potential $f(x,t)$.
6. Measure the inter-wall phase difference $\Delta\phi$: for identical activations the true $\Delta\phi = 0$; the scalar approach should produce $\Delta\phi \approx \pi$; the proposed approach should produce $\Delta\phi \approx 0$.

**Metric**: Phase error $|\Delta\phi - 0|$ across 100 randomly selected sulcal wall pairs. Statistical comparison via Wilcoxon signed-rank test.

### 4.2 Synthetic benchmark: traveling wave detection fidelity

**Objective**: Demonstrate superior wave detection accuracy relative to GBF + optical flow.

**Protocol**:
1. Simulate cortical traveling waves on the fsaverage5 surface using the wave equation $\partial^2 u / \partial t^2 = -c^2 \Delta_S u$ at three velocities: 2, 5, and 8 m/s (physiologically motivated bounds).
2. Generate MEG sensor data via BEM forward model at three SNR levels: 0, 10, and 20 dB.
3. Reconstruct sources using GBF (scalar) and the proposed vector eigenmode framework.
4. Apply wave detection via:
   - **GBF**: Hilbert transform on scalar amplitude, optical flow on phase field.
   - **Proposed**: Hilbert transform on Hodge scalar potential $f(x,t)$, optical flow on $\phi_f(x,t)$.
   - **Proposed (spectral)**: Dispersion relation recovery from joint $(\mu_k, \omega)$ spectrum.
5. Compare recovered wave velocity and propagation direction against ground truth.

**Metrics**: Velocity error (m/s), direction error (degrees), detection rate across noise levels.

### 4.3 Comparison with GBF on Meta-Source Benchmark

**Objective**: Confirm that the vector extension does not degrade source localization accuracy relative to scalar GBF.

**Protocol**: Apply the proposed framework to the 200-component Meta-Source Benchmark introduced by Wang et al. (2026), using the same evaluation metrics (NRMSE, localization error, Pearson correlation, cosine similarity, AUC) and noise conditions (Gaussian and realistic noise across SNR = -20 to +20 dB).

**Expected result**: Comparable or superior performance to GBF, particularly in regions of high sulcal complexity (temporal lobe, insula, cingulate).

### 4.4 Task-evoked data: phase-based propagation analysis

**Objective**: Demonstrate improved wave propagation characterization in empirical MEG data.

**Dataset**: HCP MEG dataset (n = 80 participants, 3-Tesla MEG, 248 channels, 1 kHz sampling rate).

**Protocol**:
1. Select motor task data (right-hand finger tapping, well-characterized beta-band (~20 Hz) propagation from motor cortex).
2. Reconstruct sources using GBF and the proposed framework.
3. Compute phase-based wave detection using both scalar amplitude (GBF pipeline) and Hodge scalar potential (proposed pipeline).
4. Identify beta-band traveling wave trajectories and measure propagation velocity.
5. Validate against known motor cortex propagation patterns from ECoG literature (mesoscale traveling waves at ~0.1–1 m/s within motor cortex, longer-range propagation at ~3–8 m/s across cortico-cortical pathways).
6. Compare joint $(\mu_k, \omega)$ spectrum against the physiological speed limit mask — quantify the fraction of detected wave energy within the physiologically accessible region.

**Metric**: Correspondence of detected propagation pathways with known neuroanatomy, velocity agreement with ECoG literature, fraction of energy within physiological bounds.

### 4.5 Resting-state: joint spectral characterization

**Objective**: Characterize the empirical joint $(\mu_k, \omega)$ spectrum of resting-state MEG and discriminate wave from diffusion dynamics.

**Protocol**:
1. Apply the proposed framework to resting-state HCP MEG (n = 80).
2. Compute the group-averaged joint power spectrum $P(k, f)$ in the $(\mu_k, \omega)$ plane.
3. Test for dispersion curve vs. wedge structure:
   - Fit a power-law dispersion relation $\omega = c\mu_k^\alpha$ to the spectral ridge.
   - $\alpha = 0.5$ indicates wave-like dynamics; $\alpha \to 0$ (flat ridge) indicates diffusion-like dynamics.
4. Compare the empirical spectral decay $P_{\text{neural}}(k)$ against the noise floor $N(k)$ after deconvolution by $G(k) = \|\tilde{L}_k\|^2$.
5. Quantify $K^*(f)$ — the frequency-dependent effective eigenmode number — across canonical frequency bands ($\delta$, $\theta$, $\alpha$, $\beta$, $\gamma$).

**Metric**: Dispersion exponent $\alpha$ per frequency band, $K^*(f)$ profile, empirical $c_{\text{empirical}}$ vs. known conduction velocities.

### 4.6 Epilepsy: ictal onset localization

**Objective**: Validate the clinical utility of orientation-consistent source reconstruction for epileptic zone localization.

**Dataset**: HDEEG-IED-SurgOutcome dataset (n = 24, 257-channel HD-EEG, favorable surgical outcomes; Vorderwülbecke et al., 2025) and Huashan Hospital private dataset (n = 4, 256-channel EEG).

**Protocol**:
1. Apply both GBF and the proposed framework to interictal spike EEG.
2. Reconstruct sources at the mid-GFP frame (50% of GFP peak).
3. Measure minimum Euclidean distance from peak activation to resection mask boundary.
4. Additionally apply Hodge decomposition to separate the irrotational (dipolar) component of the spike — hypothesized to be more focal and closer to the true epileptogenic zone.
5. Compare localization error distributions via paired Wilcoxon signed-rank test.

**Expected result**: Improved localization accuracy particularly for sources near deep sulci (hippocampus, insula, cingulate) where sulcal orientation artifacts are most severe.

### 4.7 Intracranial stimulation: wave propagation trajectory validation

**Objective**: Validate the vector eigenmode wave detection against the iES-CCEP ground truth used by Wang et al. (2026).

**Dataset**: iES-CCEP EEG dataset (n = 35, 156-channel EEG, 318 sessions; Parmigiani et al., 2022).

**Protocol**:
1. Reconstruct sources using GBF and the proposed framework.
2. Apply optical flow to (a) signed scalar amplitude (GBF) and (b) Hodge scalar potential (proposed).
3. Compare streamline trajectories against known stimulation-evoked propagation pathways from intracranial recordings (Veit et al., 2021; Lemaréchal et al., 2022).
4. Quantify streamline consistency — fraction of streamlines that follow anatomically plausible white matter pathways — using a white matter tractography atlas.
5. Assess whether spurious sulcal-crossing streamlines (a predicted artifact of the scalar approach) are reduced in the proposed framework.

**Metric**: Streamline anatomical plausibility, reduction in spurious sulcal-crossing trajectories, propagation velocity agreement with iES literature.

---

## 5. Discussion

### 5.1 Relationship to existing work

The proposed framework extends GBF (Wang et al., 2026) in three key dimensions:

**Geometric**: From scalar LBO eigenmodes to vector connection Laplacian eigenmodes, using trivial connections (Crane et al., 2010) to establish globally consistent frames. This is the principled resolution of the sulcal orientation problem.

**Physical**: From projected scalar amplitude to the Hodge scalar potential, correctly separating MEG-visible (irrotational) from MEG-invisible (solenoidal) components of cortical current. This makes the wave analysis physically grounded — it is the dipolar component of the current, not its arbitrary scalar projection, that generates MEG signals.

**Dynamical**: From vertex-space phase analysis to joint $(\mu_k, \omega)$ spectral analysis, enabling discrimination between wave, diffusion, and reaction-diffusion dynamics without assuming a dynamical model, and grounded in physiological propagation speed constraints.

### 5.2 Limitations and open questions

**Computational cost**: The connection Laplacian eigendecomposition and Hodge decomposition add computational overhead relative to scalar GBF. Both are implementable efficiently in Geometry Central, and the composed transform $A_{\text{vec}}$ is a one-time computation per subject. Runtime estimates suggest a factor of 3–5x increase over scalar GBF, which remains tractable for offline analysis.

**Free orientation ill-conditioning**: Allowing free dipole orientation triples the source unknowns. The vector eigenmode basis mitigates this by compressing to $K$ vector modes rather than $3 \times n_\text{vert}}$ scalar unknowns, but the conditioning of the inverse is slightly worse than fixed-orientation GBF. The logarithmic spectral prior on $\mu_k$ is the primary regularization mechanism.

**Subcortical sources**: The current framework is defined on the cortical surface manifold. Subcortical structures — hippocampus, thalamus, basal ganglia — require volumetric extension of the vector eigenmode framework, as noted as a limitation in GBF. Preliminary feasibility is suggested by Wang et al. (2026)'s Extended Data Figure 5.

**Unknown dynamical model**: The joint $(\mu_k, \omega)$ spectrum identifies the dispersion structure empirically, but interpreting it requires some knowledge of which physical model the brain is operating under. The physiological speed limit provides a necessary but not sufficient constraint. Bayesian model comparison between wave, diffusion, and reaction-diffusion hypotheses within the joint spectral framework remains an open methodological question.

### 5.3 Clinical and scientific impact

Clinically, correcting sulcal orientation artifacts is expected to most benefit EZ localization in temporal lobe epilepsy and insula epilepsy — precisely the cases where sources are deepest, sulcal geometry is most complex, and accurate noninvasive localization is most clinically impactful, reducing the need for invasive iEEG implantation.

Scientifically, the joint $(\mu_k, \omega)$ spectral framework provides, for the first time, a principled vocabulary for characterizing the spatiotemporal complexity of cortical dynamics — distinguishing local oscillations, traveling waves, diffusive spreading, and their interactions — using the natural geometric language of the cortical manifold itself.

---

## 6. Conclusion

We have proposed a vector eigenmode framework for MEG/EEG source imaging that addresses the fundamental geometric inconsistency in existing eigenmode-based approaches: the arbitrary sign convention for dipole orientations on opposite sulcal walls. By establishing globally consistent orientation frames using trivial connections and the vector heat method, decomposing cortical currents via the Hodge decomposition to isolate the MEG-visible scalar potential, expanding in eigenmodes of the connection Laplacian rather than the scalar LBO, and analyzing the resulting vector eigenmode time series in the joint $(\mu_k, \omega)$ spectral domain, the framework provides principled solutions to problems that GBF and all existing ESI methods leave unaddressed.

The proposed approach is grounded in differential geometry — specifically the theory of connections on vector bundles over Riemannian manifolds — and is fully implementable using existing tools in Geometry Central and MNE-Python. It is simultaneously a correction to a known artifact, a principled extension of the GBF framework to vector-valued sources, and a new methodology for characterizing cortical wave dynamics from MEG data.

---

## References

Baillet, S. (2017). Magnetoencephalography for brain electrophysiology and imaging. *Nature Neuroscience*, 20(3), 327–339.

Crane, K., Desbrun, M., & Schröder, P. (2010). Trivial connections on discrete surfaces. *Computer Graphics Forum*, 29(5), 1525–1533.

Desbrun, M., Hirani, A. N., Leok, M., & Marsden, J. E. (2005). Discrete exterior calculus. *arXiv preprint math/0508341*.

He, B., Sohrabpour, A., Brown, E., & Liu, Z. (2018). Electrophysiological source imaging: a noninvasive window to brain dynamics. *Annual Review of Biomedical Engineering*, 20, 171–196.

Lemaréchal, J.-D., et al. (2022). A brain atlas of axonal and synaptic delays based on modelling of cortico-cortical evoked potentials. *Brain*, 145(5), 1653–1667.

Pang, J. C., Aquino, K. M., Oldehinkel, M., Robinson, P. A., Fulcher, B. D., Breakspear, M., & Fornito, A. (2023). Geometric constraints on human brain function. *Nature*, 618, 566–574.

Parmigiani, S., et al. (2022). Simultaneous stereo-EEG and high-density scalp EEG recordings to study the effects of intracerebral stimulation parameters. *Brain Stimulation*, 15(3), 664–675.

Raut, R. V., et al. (2021). Global waves synchronize the brain's functional systems with fluctuating arousal. *Science Advances*, 7(30), eabf2709.

Roberts, J. A., et al. (2019). Metastable brain waves. *Nature Communications*, 10(1), 1056.

Sharp, N., Soliman, Y., & Crane, K. (2019). The vector heat method. *ACM Transactions on Graphics*, 38(3), 1–19.

Veit, M. J., et al. (2021). Temporal order of signal propagation within and across intrinsic brain networks. *Proceedings of the National Academy of Sciences*, 118(48), e2105031118.

Vidaurre, D., et al. (2018). Spontaneous cortical activity transiently organises into frequency specific phase-coupling networks. *Nature Communications*, 9(1), 2987.

Vorderwülbecke, B. J., et al. (2025). High-density EEG source localisation of averaged interictal epileptic discharges validated by surgical outcome. *Scientific Data*, 12(1), 1441.

Wang, S., Lou, K., Wei, C., et al. (2026). A geometry aware framework enhances noninvasive mapping of whole human brain dynamics. *arXiv:2604.25592v1*.

---

*Proposal prepared for review. All computational components target implementation in C++ (Geometry Central) for geometric operations and Python (MNE-Python) for MEG preprocessing and forward modeling.*