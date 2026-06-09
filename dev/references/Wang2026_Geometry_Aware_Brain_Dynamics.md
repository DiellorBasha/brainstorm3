---
title: "A geometry aware framework enhances noninvasive mapping of whole human brain dynamics"
authors:
  - Wang, Song
  - Lou, Kexin
  - Wei, Chen
  - Sheng, Zhiyuan
  - Tang, Jiahao
  - Peng, Kaining
  - Shen, Xinke
  - Mei, Shuhao
  - Chen, Liang
  - Gu, Dongfeng
  - Liu, Quanying
year: 2026
journal: "arXiv preprint"
doi: "10.48550/arXiv.2604.25592"
tags:
  - electrical-engineering-signal-processing
  - quantitative-biology-neurons-cognition
  - eeg-source-imaging
  - geometric-eigenmodes
  - inverse-problem
  - cortical-geometry
  - traveling-waves
date_added: "2026-05-02"
zotero_key: "L5ZTV9A3"
---

# A geometry aware framework enhances noninvasive mapping of whole human brain dynamics

**Authors**: Wang S, Lou K, Wei C, Sheng Z, Tang J, Peng K, Shen X, Mei S, Chen L, Gu D, Liu Q
**Year**: 2026 | **Journal**: arXiv preprint
**DOI**: [10.48550/arXiv.2604.25592](https://doi.org/10.48550/arXiv.2604.25592)

---

## TL;DR

This paper introduces Geometric Basis Functions (GBF), a source imaging method that uses participant-specific cortical surface eigenmodes (derived via Laplace-Beltrami decomposition) as spatial priors for the EEG/MEG inverse problem. Validated across synthetic benchmarks, task-evoked data, resting-state connectivity, intracranial stimulation, and clinical epilepsy data, GBF consistently outperforms conventional methods in localization accuracy and spatiotemporal fidelity.

## Key Contributions

- **Geometry-informed inverse solution**: Embeds individual cortical surface eigenmodes directly into the source imaging framework, replacing biologically implausible or generic priors with participant-specific anatomical constraints
- **Meta-Source Benchmark**: Introduces a novel synthetic validation benchmark derived from 26,273 NeuroVault statistical maps annotated with 1,307 Neurosynth cognitive terms, producing 200 functionally grounded source-EEG pairs with known ground truth
- **Multi-domain validation**: Demonstrates superiority across five independent validation domains — synthetic benchmarks, task-evoked EEG, MEG-iEEG functional connectivity, intracranial electrical stimulation, and epileptogenic zone localization
- **Wave propagation tracking**: Combines GBF with optical flow analysis to reconstruct millisecond-scale cortical wave propagation from scalp EEG, validated against known stimulation sites
- **Clinical translation**: Achieves significantly closer localization of epileptogenic zones to surgical resection boundaries compared to standard methods

## Background & Motivation

Noninvasive electrophysiology (EEG/MEG) offers millisecond temporal resolution and whole-brain coverage, making it uniquely suited to studying fast neural dynamics. However, the inverse problem — reconstructing cortical sources from sensor measurements — is fundamentally ill-posed: multiple source configurations can produce identical scalp signals. Current approaches address this through regularization priors, but most are either biologically implausible (e.g., minimum-norm assumptions) or require extensive hand-tuning (e.g., patch-based methods with many free parameters).

A critical missing ingredient is cortical geometry. The brain's folded surface constrains how neural activity propagates and organizes spatially. Recent fMRI work by Pang et al. (2023) demonstrated that whole-brain activity across ~32,000 vertices can be compactly represented by roughly 200 geometric eigenmodes — the eigenfunctions of the [[Laplace-Beltrami Operator]] on the cortical surface. Since EEG/MEG signals predominantly reflect postsynaptic cortical currents, incorporating these geometry-derived basis functions as spatial priors should yield more accurate and neurophysiologically meaningful source estimates.

This paper formalizes that intuition into the GBF framework: extract participant-specific cortical eigenmodes, use them as a spatial basis for the inverse problem, and apply a logarithmic spectral prior that favors smooth, low-frequency modes while retaining higher-frequency anatomical detail for stability.

## Methods

### Geometric basis decomposition

The framework begins by computing the [[Laplace-Beltrami Operator]] on each participant's cortical surface mesh (extracted via FreeSurfer). Solving the eigenvalue problem yields an orthonormal set of eigenmodes (psi_k) with corresponding eigenvalues (lambda_k) sorted by spatial frequency. Low eigenvalues correspond to slowly varying global patterns; high eigenvalues capture localized, fine-grained anatomical variations. Approximately 300 modes are retained empirically — the first 200 capture ~87% of cortical variance, and additional modes increase ill-conditioning without improving effective spatial resolution.

### Inverse problem formulation

Sensor measurements are modeled as y = KA*theta + epsilon, where K is the lead-field matrix (computed via three-layer BEM forward model), A is the geometric eigenmode matrix, theta is the source coefficient vector, and epsilon is Gaussian noise. The key insight is reformulating this as y = L*theta where L = KA, reducing the problem from estimating ~20,000 cortical dipole amplitudes to estimating ~300 eigenmode coefficients.

A logarithmic spectral prior is placed on the coefficients: the inverse covariance is diagonal with entries proportional to -beta/log(lambda_k). This penalizes high-frequency modes while only moderately constraining higher-order contributions — avoiding the over-smoothing produced by power-law or exponential decay schemes. A single global parameter beta controls regularization strength.

The maximum a posteriori (MAP) solution has a closed form: theta_hat = (L^T L + Sigma^{-1})^{-1} L^T y. This is computationally efficient and avoids iterative optimization entirely.

### Meta-Source Benchmark

Recognizing that ground-truth source locations are unknown in real EEG, the authors constructed a comprehensive synthetic benchmark. They selected 26,273 statistical maps from NeuroVault, annotated them with 1,307 Neurosynth cognitive terms, applied iterative PCA (70 components x 200 iterations = 14,000 components), and clustered the results via KNN into 200 unique spatial maps. These were projected to the cortical surface and then to scalp sensors via BEM forward modeling, with added Gaussian and realistic (empirical covariance) noise at controlled SNR levels. The result is 200 source-EEG pairs with functionally grounded, distributed ground-truth source patterns.

### Validation domains

The framework is validated across five domains: (1) the Meta-Source Benchmark with five complementary metrics (NRMSE, localization error, Pearson correlation, cosine similarity, AUC); (2) task-evoked EEG from visual, auditory, somatosensory, and motor paradigms; (3) MEG-derived virtual iEEG connectomes compared against real iEEG from 110 epilepsy patients using amplitude envelope correlation; (4) intracranial electrical stimulation (318 sessions, 35 participants) with simultaneous HD-EEG providing ground-truth stimulation sites; and (5) epileptogenic zone localization in 24 patients with favorable post-surgical outcomes, validated against resection cavity masks.

### Wave propagation analysis

For the stimulation data, GBF-reconstructed source activity is analyzed via optical flow methods (following Roberts et al. 2019). Instantaneous phase is obtained via Hilbert transform, spatial phase gradients are estimated on the cortical surface, and velocity vector fields are derived. Streamline tracing with forward Euler integration (8 mm steps) maps propagation pathways, and DBSCAN clustering identifies convergence hubs.

## Results

### Meta-Source Benchmark

GBF significantly outperformed all conventional methods (MNE, wMNE, sLORETA, eLORETA, dSPM, LCMV) across all five evaluation metrics. At SNR = 5 dB with realistic noise, GBF achieved NRMSE of ~0.14 compared to ~0.18 for MNE and higher for other methods. All pairwise comparisons were statistically significant (paired t-tests and Wilcoxon tests with FDR correction). Regional analysis showed GBF yielded lower NRMSE in approximately 95% of cortical parcels, with residual errors concentrated in deep/ventral regions consistent with known depth-sensitivity limitations.

### Task-evoked EEG

Across all four paradigms, GBF produced well-localized activations in the expected cortical regions: primary visual cortex (peak at 120 ms), bilateral auditory cortices (210 ms), contralateral somatosensory cortex (280 ms), and contralateral motor cortex (beta-band ERD/ERS). Source maps closely matched Neurosynth meta-analytic references. Conventional methods (wMNE, eLORETA) showed greater spatial variability and less precise temporal dynamics.

### MEG-iEEG functional connectivity

GBF-derived virtual iEEG connectomes from 80 HCP MEG participants were compared against real iEEG from 110 epilepsy patients across six frequency bands. GBF achieved the highest cross-modal correlations in all bands, with the best performance in beta (r = 0.42-0.47) and alpha (r = 0.39-0.45). Statistical testing against 5,000 spatially permuted null models confirmed significance at FDR-adjusted P = 4.0 x 10^-4 across all bands.

### Intracranial stimulation

Using 309 quality-checked stimulation sessions from 35 participants, GBF achieved mean localization error of ~20 mm compared to ~30 mm for eLORETA, ~40 mm for LCMV, ~50 mm for MxNE, and ~60+ mm for wMNE. All pairwise comparisons were highly significant (Wilcoxon signed-rank, FDR corrected; p_max = 2.42 x 10^-9). A linear mixed-effects model showed that while localization error increases with stimulation depth for all methods (beta = 0.428, P = 2.83 x 10^-5), GBF maintained a consistent advantage at all depths. Optical flow analysis of stimulation-evoked responses revealed millisecond-scale propagation from stimulation sites to convergence hubs in somatomotor, prefrontal, and dorsal attention regions.

### Epileptogenic zone localization

In 24 patients with favorable surgical outcomes, GBF achieved the smallest mean distance from peak source activity to resection mask boundary (~10-15 mm vs. ~20-25 mm for MNE). GBF vs. MNE was statistically significant (p = 0.012). Additional validation in patients with clinical SOZ annotations showed GBF peaks closer to the annotated onset zone than dSPM or eLORETA. Frequency-domain Granger causality applied to GBF-reconstructed sources revealed directed causal links from the SOZ to downstream network regions.

## Figures

### Figure 1 — Framework overview and validation domains
Illustrates the GBF pipeline: eigenmode extraction from individual cortical meshes via Laplace-Beltrami decomposition, the inverse problem formulation constraining sources as linear combinations of geometric modes, and the five validation domains (benchmark, task, connectivity, stimulation, epilepsy).

### Figure 2 — Meta-Source Benchmark performance
Shows the benchmark data generation pipeline from NeuroVault/Neurosynth, example reconstructions at SNR = 5 dB, and systematic comparison across six methods on five metrics. GBF achieves the lowest NRMSE and highest correlation/similarity scores across noise conditions.
**Stats**: FDR-corrected paired t-tests/Wilcoxon tests, significant across all metrics and pairwise comparisons (Supplementary Tables 10-11)

### Figure 3 — Task-evoked EEG validation
Displays ERP/GFP waveforms and scalp topographies for four paradigms (visual, auditory, somatosensory, motor), with source maps from GBF, wMNE, and eLORETA compared against Neurosynth meta-analytic references. GBF shows tightest spatial correspondence to expected anatomical activations.

### Figure 4 — MEG-iEEG functional connectivity
Presents the cross-modal validation pipeline (110 iEEG patients, 80 HCP MEG participants), frequency-band-specific correlations between GBF-derived and iEEG connectomes, null model statistical testing, and regional edge-level analysis.
**Stats**: GBF vs. null: FDR-adjusted P = 4.0 x 10^-4 across all frequency bands; beta-band GBF r = 0.42-0.47

### Figure 5 — Intracranial stimulation localization and wave propagation
Shows representative GBF localizations, localization error distributions across methods, depth-stratified error analysis, and optical flow streamlines mapping stimulation-evoked wave propagation in N1 and N2 time windows.
**Stats**: GBF vs. all methods: Wilcoxon signed-rank, FDR corrected, p_max = 2.42 x 10^-9; depth model: beta = 0.428, P = 2.83 x 10^-5; GBF vs. eLORETA at mean depth: beta = 6.926, P = 1.86 x 10^-5

### Figure 6 — Epileptogenic zone localization
Displays resection mask reconstruction pipeline, distance distributions from peak source activity to resection boundary across methods, and representative patient cases with clinical SOZ annotations.
**Stats**: GBF vs. MNE distance to resection: p = 0.012, paired two-sided test, n = 24

## Discussion & Implications

GBF represents a shift from generic spatial priors to participant-specific, anatomically grounded constraints for [[EEG Source Imaging]]. By embedding cortical geometry directly into the inverse solution, it achieves biological plausibility without sacrificing computational efficiency — the closed-form MAP solution avoids iterative optimization entirely. The finding that ~200-300 geometric modes provide a compact yet accurate representation of whole-brain electrophysiological dynamics connects with the fMRI eigenmode literature (Pang et al. 2023) and suggests a general principle: cortical geometry fundamentally shapes the vocabulary of large-scale brain activity patterns.

The multi-domain validation strategy is a particular strength. Rather than relying solely on synthetic benchmarks or single experimental paradigms, the paper demonstrates GBF's advantages across the full chain from simulation to clinical application. The intracranial stimulation results are especially compelling because they provide rare ground-truth source locations and enable direct comparison of localization accuracy with depth analysis.

The wave propagation results demonstrate that GBF's temporal fidelity is sufficient to track millisecond-scale dynamics — a capability that opens new possibilities for studying [[Traveling Waves]] and directed cortical communication noninvasively. This connects directly to earlier work on rotating spindle waves ([[Muller2016_Rotating_Waves_Sleep_Spindles]]) and broadens the scope from high-density invasive recordings to scalp-level reconstruction.

## Limitations & Open Questions

- **Cortical-only**: The framework currently focuses on the cortical surface. Preliminary subcortical extension is shown, but hippocampal, thalamic, and basal ganglia source imaging remains unvalidated — critical for many cognitive and clinical applications.
- **Temporal independence**: Each time point is treated independently. Incorporating temporal priors or state-space formulations could further improve reconstruction of fast dynamics.
- **Structural anomalies**: Performance in the presence of cortical lesions, malformations, or significant atrophy is unknown. Adaptive strategies for abnormal geometry are needed for broader clinical applicability.
- **Sample sizes**: The epilepsy validation (n = 24) and some connectivity analyses use relatively small cohorts. Larger, more diverse validations are needed before clinical deployment.
- **Integration with deep learning**: The authors note that GBFs could serve as geometric priors for neural network-based source imaging — an unexplored direction that could combine the anatomical grounding of eigenmodes with the flexibility of learned representations.
- **What determines the optimal number of modes?** The paper uses ~300 empirically, but the relationship between cortical geometry complexity, recording density, and optimal truncation deserves formal treatment.

## Connections

- **Related papers**: [[Muller2016_Rotating_Waves_Sleep_Spindles]] — GBF enables noninvasive tracking of the kind of cortical wave propagation that Muller et al. characterized invasively
- **Key concepts**: [[EEG Source Imaging]], [[Geometric Eigenmodes]], [[Traveling Waves]]
- **Methods referenced**: [[Laplace-Beltrami Operator]], minimum-norm estimation (MNE), sLORETA, eLORETA, dSPM, LCMV beamformer
- **Datasets**: Human Connectome Project (HCP), MNI Open iEEG Atlas, NeuroVault, Neurosynth
- **Key authors**: [[Quanying Liu]] (senior/corresponding author)

## User Annotations

No annotations yet.

## Raw Metadata

<details>
<summary>Full bibliographic details</summary>

- **Item Type**: preprint
- **Repository**: arXiv
- **Archive ID**: 2604.25592
- **DOI**: 10.48550/arXiv.2604.25592
- **URL**: http://arxiv.org/abs/2604.25592
- **Zotero Key**: L5ZTV9A3
- **Date Added**: 2026-05-02
- **Tags**: Electrical Engineering and Systems Science - Signal Processing, Quantitative Biology - Neurons and Cognition

</details>
