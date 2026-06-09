# Spatiotemporal Filter Banks for MEG/EEG: A Principled Framework for Measuring Neural Current Patterns on the Cortical Manifold

---

## Summary

Spectral filtering is the foundational analytical tool of MEG/EEG research. Temporal filters — bandpass, notch, Hilbert-based analytic filters — decompose neural time series into frequency bands associated with canonical neural rhythms. Yet these filters operate exclusively on the temporal axis, treating the spatial distribution of neural activity as secondary. The complementary spatial dimension — which cortical patterns are active — is addressed separately and incompletely by source imaging pipelines that recover static spatial maps but discard the joint spatiotemporal structure of the signal. No existing MEG/EEG analysis framework treats space and time symmetrically as a unified filtering problem.

We propose that the natural unification is achieved through **spatiotemporal filter banks** defined in the joint $(\lambda, \omega)$ spectral domain — where $\lambda$ indexes spatial frequency via the eigenvalues of the Laplace-Beltrami operator (LBO) on the cortical surface, and $\omega$ indexes temporal frequency. Each filter in this bank is a complex-valued operator that selects neural current patterns at a specific combination of spatial scale and temporal frequency, capturing the dynamics — propagation, diffusion, oscillation — that neither purely spatial nor purely temporal filters can resolve. For nonstationary signals, where the active $(\lambda, \omega)$ structure changes over time, the filter bank is extended to a **complex-valued wavelet tensor** $W_k(s,t) \in \mathbb{C}$ — jointly resolving spatial scale, temporal frequency, and time — providing amplitude and phase for each spatiotemporal scale at each moment.

We further propose the development of an open **spatiotemporal filter library** — a principled, computationally efficient collection of filter banks for cortical MEG/EEG data, indexed by spatial eigenmode, temporal frequency, and dynamical regime, implemented in C++ (Geometry Central) and Python (MNE-Python). This library provides, for the first time, a symmetric treatment of space and time in MEG/EEG analysis, grounded in the differential geometry of the cortical surface and the physics of electromagnetic neural measurement.

---

## 1. Introduction

### 1.1 Filtering as the language of neural signal analysis

Signal filtering is the primary analytical vocabulary of systems neuroscience. When a researcher bandpass filters MEG data at 8–12 Hz to isolate alpha oscillations, or applies a Hilbert transform to obtain the oscillation envelope, or uses a wavelet transform to track time-varying spectral power, they are applying a filter — an operator that selects a subset of the signal's content according to some criterion, discarding the rest.

The criterion implicitly defines what the researcher believes is meaningful in the signal. A temporal bandpass filter asserts that the frequency content within the passband is the signal of interest and everything outside is noise or irrelevant. The filter is a hypothesis about signal structure, expressed in mathematical form.

This framing — filtering as hypothesis — is productive because it makes the assumptions explicit and the outputs interpretable. A bandpass-filtered signal at 40 Hz contains the neural activity at that temporal scale. The filter output is a direct measure of the quantity of interest.

The same logic should apply to the spatial dimension. A spatial filter that selects neural activity at a specific cortical spatial scale — smooth global patterns vs. focal patches vs. fine columnar structure — would provide a direct measure of the spatial content of neural dynamics at that scale. Combined with temporal filtering, a **spatiotemporal filter** would measure neural activity at a specific combination of spatial scale and temporal frequency simultaneously.

This is the conceptual starting point of the present proposal. We argue that the natural framework for spatiotemporal filtering of MEG/EEG data is the joint $(\lambda, \omega)$ spectral domain, that the correct spatial spectral basis is provided by the LBO eigenmodes of the individual cortical surface, and that the filter bank over this domain is the principled generalization of standard temporal filtering methods to the full spatiotemporal signal.

### 1.2 The asymmetry problem in current MEG/EEG analysis

Current MEG/EEG analysis treats time and space with profound asymmetry:

**Temporal axis**: Richly developed filter theory. Temporal Fourier transforms, bandpass filters, Hilbert transforms, wavelet transforms, empirical mode decomposition, filter banks — a mature, theoretically grounded toolkit with decades of development.

**Spatial axis**: Ad hoc and model-dependent. Beamformers, minimum norm estimates, sLORETA, dSPM — these are inverse problem solvers, not filters in the signal processing sense. They do not decompose the spatial content of the signal into a basis ordered by spatial frequency. They do not provide a symmetric spatial analog of the temporal bandpass filter. They solve for spatial amplitudes given a forward model, but the spatial structure of the solution is determined by the regularization prior — an engineering choice — rather than by the signal itself.

The spatial analog of temporal filtering has been missing because there was no obvious spatial spectral basis. Cortical space is not a flat plane where Fourier modes are natural — it is a curved, folded, geometrically complex manifold. The eigenvalues and eigenfunctions of the LBO on this manifold provide exactly the spatial spectral basis needed (Pang et al., 2023; Wang et al., 2026), but their role as **spatial filters** rather than merely as basis functions for source reconstruction has not been developed.

### 1.3 Scope and contributions

This proposal develops the spatiotemporal filter framework from first principles, building from familiar temporal filtering through spatial filtering to joint spatiotemporal filtering, and proposes:

1. A formal definition of **LBO eigenmode filters** as the spatial analogs of temporal bandpass filters, characterizing what they measure and what they capture.

2. A **joint $(\lambda, \omega)$ filter bank** that simultaneously selects spatial and temporal frequency content, showing that filters here capture dynamics — propagation velocity, dispersion relations, diffusive spreading — that neither purely spatial nor purely temporal filters access.

3. A **complex-valued wavelet tensor** extension for nonstationary signals, providing full amplitude and phase resolution in the $(\lambda, \omega, t)$ volume.

4. A proposal for an open **spatiotemporal filter library** implementing these filter banks as reusable, composable computational tools for MEG/EEG research.

---

## 2. From Simple to Complex: Building the Spatiotemporal Filter

### 2.1 Temporal filters: the familiar starting point

#### 2.1.1 The Fourier transform as a filter bank

The temporal Fourier transform of a MEG time series $d(t)$ at a single sensor:

$$\tilde{d}(\omega) = \int_{-\infty}^{\infty} d(t)\, e^{-i\omega t}\, dt$$

decomposes the signal into complex-valued frequency components. The magnitude $|\tilde{d}(\omega)|$ is the spectral amplitude at frequency $\omega$; the argument $\arg(\tilde{d}(\omega))$ is the spectral phase. The power spectral density $|\tilde{d}(\omega)|^2$ is the familiar quantity reported in neural oscillation studies.

In filter bank language, the Fourier transform is an infinite bank of **narrowband complex filters**, one per frequency $\omega$, each extracting the amplitude and phase of the signal at that temporal frequency. Applying a bandpass filter is selecting a subset of this bank — retaining components within the passband $[\omega_1, \omega_2]$ and zeroing the rest.

The key properties that make this useful:

- **Parseval's theorem**: Total power is preserved. $\int |d(t)|^2 dt = \int |\tilde{d}(\omega)|^2 d\omega / (2\pi)$.
- **Completeness**: Any signal is exactly reconstructed from its Fourier components.
- **Orthogonality**: Different frequency components are independent — a filter at $\omega_1$ does not contaminate the output at $\omega_2$.
- **Physical interpretability**: Each filter output corresponds to a specific oscillation rate.

#### 2.1.2 The Hilbert transform and analytic signal

For a narrowband signal at frequency $\omega_0$, the analytic signal:

$$z(t) = d(t) + i\,\mathcal{H}[d(t)] = A(t)\, e^{i\phi(t)}$$

provides the instantaneous amplitude $A(t) = |z(t)|$ — the signal envelope — and instantaneous phase $\phi(t) = \arg(z(t))$. The Hilbert transform is a $\pi/2$ phase shifter: it replaces each Fourier component at frequency $\omega$ with a component shifted by $-\pi/2\, \text{sgn}(\omega)$. Together, $d(t)$ and $\mathcal{H}[d(t)]$ form a complex-valued representation of the signal that carries both quadrature components.

The Hilbert transform works cleanly only on narrowband signals. For broadband signals, instantaneous amplitude and phase lose physical meaning. This motivates the wavelet transform.

#### 2.1.3 The Morlet wavelet transform: time-frequency localization

The continuous wavelet transform with a Morlet wavelet:

$$W(s, t) = \int d(\tau)\, \psi_s^*\!\left(\frac{\tau - t}{s}\right) d\tau$$

where $\psi_s(t) = \pi^{-1/4} e^{i\omega_0 t/s} e^{-t^2/(2s^2)}$ is a Gaussian-windowed complex sinusoid at scale $s$ (corresponding to frequency $f = f_0/s$), produces a complex-valued function of time $t$ and scale $s$.

This is the **time-frequency filter bank** — each $(s, t)$ cell contains the complex amplitude $W(s,t) = A(s,t)\, e^{i\phi(s,t)}$, resolving both instantaneous amplitude and phase at each temporal frequency and each time. For nonstationary signals — where the frequency content changes over time — this is the appropriate tool. The Fourier transform and Hilbert transform are special cases: the former integrates over all time (no temporal localization), the latter applies to a pre-selected frequency band (no multi-scale resolution).

The wavelet transform is the **standard tool for nonstationary temporal analysis**. The present proposal argues that the natural generalization to spatiotemporal analysis requires extending this filter bank to include the spatial dimension — and that the correct spatial basis for this extension is the LBO eigenmode basis of the cortical surface.

### 2.2 Spatial filters: the missing symmetric tool

#### 2.2.1 Why spatial filters are harder

Temporal filtering is straightforward because time is a one-dimensional, flat continuum. The Fourier basis $\{e^{i\omega t}\}$ is natural, universal, and independent of the signal content. The filter response is determined purely by the filter design.

Spatial filtering is harder for three reasons:

**Geometry**: Cortical space is a two-dimensional Riemannian manifold — curved, folded, with complex topology. There is no flat spatial Fourier basis. The appropriate spectral basis depends on the geometry of the specific manifold — it is individual-specific.

**Measurement**: MEG sensors do not measure cortical current directly. They measure the magnetic field at sensor locations outside the head. Spatial filtering of the cortical current requires either solving the inverse problem first, or designing filters that operate on sensor data while targeting cortical spatial content.

**Discretization**: The cortical surface is represented as a triangular mesh with ~10,000–160,000 vertices. The discrete approximation to the continuous spatial spectral basis depends on the mesh resolution and quality.

#### 2.2.2 The LBO eigenmode as a spatial filter

The Laplace-Beltrami operator $\Delta_S$ on the cortical surface $S$, discretized via cotangent weights, has eigenmodes $\psi_k : S \to \mathbb{R}$ satisfying:

$$\Delta_S \psi_k = -\lambda_k \psi_k, \quad k = 0, 1, 2, \ldots$$

with $0 = \lambda_0 \leq \lambda_1 \leq \lambda_2 \leq \ldots$. These eigenmodes are the **cortical analogs of Fourier modes** — the natural oscillatory patterns on the cortical manifold ordered by spatial frequency.

A **spatial eigenmode filter** is the projection of a cortical source field $q(x)$ onto eigenmode $k$:

$$\hat{q}_k = \langle q, \psi_k \rangle = \int_S q(x)\, \psi_k(x)\, dA(x)$$

This is a **spatial inner product** — it measures how much of cortical pattern $\psi_k$ is present in the current source distribution $q(x)$. It is the exact spatial analog of the Fourier coefficient $\tilde{d}(\omega) = \int d(t) e^{-i\omega t} dt$ — a projection of the signal onto a basis function.

The filter output $\hat{q}_k$ is a scalar: the amplitude of spatial mode $k$ in the source field. The full spatial filter bank is the set of all projections $\{\hat{q}_k\}_{k=0}^K$, which exactly decomposes the cortical source field into its spatial frequency components.

**What a spatial eigenmode filter captures**:

- **$k = 0$** ($\lambda_0 = 0$): The DC spatial mode — uniform activity across the entire cortical surface. This is the global mean field amplitude.
- **Small $k$** (small $\lambda_k$): Smooth, slowly varying spatial patterns — large-scale bilateral activations, global network patterns, hemisphere-scale gradients.
- **Intermediate $k$**: Mesoscale spatial patterns — area-scale activations, canonical functional network footprints.
- **Large $k$** (large $\lambda_k$): Fine spatial patterns — focal activations, sharp boundaries between active and inactive regions, cortical column-scale structure.

This is directly analogous to temporal frequency content: low-$k$ modes correspond to slow spatial variation just as low $\omega$ corresponds to slow temporal variation. The eigenvalue $\lambda_k$ plays the role of the squared spatial frequency.

#### 2.2.3 Spatial filter properties — the parallel with temporal filters

The LBO eigenmode filter bank satisfies the same fundamental properties as the temporal Fourier filter bank:

**Parseval's theorem** (spatial): $\int_S |q(x)|^2 dA = \sum_k |\hat{q}_k|^2$. Total spatial power is preserved.

**Completeness**: $q(x) = \sum_k \hat{q}_k \psi_k(x)$ exactly, for any square-integrable $q$.

**Orthogonality**: $\langle \psi_k, \psi_{k'} \rangle = \delta_{kk'}$. Different spatial modes are independent.

**Physical interpretability**: Each filter output corresponds to a specific spatial scale on the cortical surface, quantified by $\lambda_k$.

**Spatial lowpass filter**: Truncating the sum at mode $K$ defines a spatial lowpass filter — retaining only the smooth component of the cortical source field. This is the geometric operation underlying GBF source reconstruction (Wang et al., 2026).

**Spatial bandpass filter**: Retaining only modes $k \in [k_1, k_2]$ selects a specific range of cortical spatial scales — analogous to temporal bandpass filtering.

#### 2.2.4 What spatial eigenmode filters do NOT capture

Crucially, a static spatial eigenmode filter applied to $q(x)$ at a single time point $t$ captures the spatial amplitude at that moment — but nothing about dynamics. It cannot distinguish between:

- A source pattern that is stationary in space
- The same pattern propagating as a traveling wave
- The pattern growing and decaying as a heat-like diffusion
- The pattern oscillating in amplitude at temporal frequency $\omega$

All of these scenarios produce the same set of spatial eigenmode amplitudes $\{\hat{q}_k\}$ at any single time point. **Dynamics are invisible to static spatial filters.** This is the fundamental limitation that motivates the joint spatiotemporal filter.

### 2.3 The GBF framework as a spatial filter bank

The GBF framework (Wang et al., 2026) implements exactly the spatial eigenmode filter bank as a source reconstruction method. The source model:

$$x(t) = \sum_{k=1}^K \theta_k(t)\, \psi_k$$

expresses the cortical source field as a sum of spatial eigenmode filters applied to the coefficient time series $\theta_k(t)$. The MAP inverse recovers $\theta_k(t)$ — the output of the spatial filter bank — from sensor data at each time point independently.

GBF is thus a **temporal sequence of static spatial filter applications** — at each $t$, the spatial content of the source field is decomposed into eigenmode amplitudes. The time series $\theta_k(t)$ for each mode $k$ is the output of spatial filter $k$ over time.

But — and this is the key point — GBF applies the spatial filter bank **independently at each time point**. The temporal structure of $\theta_k(t)$ is not used to inform the spatial filtering at any other time point. The time axis is merely an index. The dynamics — how the spatial pattern evolves — are not captured by the filtering operation itself.

This is the precise analog of applying a static spatial Fourier transform frame-by-frame to a video signal without ever using the temporal frequency content. The spatial content is captured; the spatiotemporal structure — what is moving, at what speed, in what direction — is not.

---

## 3. The Joint ($\lambda$, $\omega$) Filter Bank: Filters that Capture Dynamics

### 3.1 Joining the two spectral axes

The temporal filter bank decomposes $d(t)$ over $\omega$. The spatial filter bank decomposes $q(x)$ over $\lambda_k$. The natural unification is a **joint spatiotemporal filter** that selects signals at a specific combination of spatial frequency $\lambda_k$ and temporal frequency $\omega$.

For a spatiotemporal source field $q(x,t)$, the joint filter output is:

$$\hat{Q}(k, \omega) = \int_S \int_{-\infty}^{\infty} q(x,t)\, \psi_k(x)\, e^{-i\omega t}\, dt\, dA(x)$$

This is the **spatiotemporal Fourier-eigenmode transform** — a projection onto the joint basis $\{\psi_k(x) \cdot e^{i\omega t}\}$ that resolves both spatial scale and temporal frequency simultaneously.

The joint power spectrum:

$$P(k, \omega) = |\hat{Q}(k,\omega)|^2$$

is the **joint $(\lambda, \omega)$ spectrum** — the primary observable of the proposed framework.

### 3.2 Why dynamics appear in the joint spectrum

The critical insight is that dynamic phenomena — propagation, diffusion, oscillation — produce **characteristic signatures in the $(\lambda, \omega)$ plane** that are invisible in either the purely spatial or purely temporal spectrum alone.

#### 3.2.1 Traveling waves: energy on a dispersion curve

A cortical traveling wave governed by the wave equation:

$$\frac{\partial^2 q}{\partial t^2} = -c^2 \Delta_S q$$

has eigenmode solutions $q(x,t) = \psi_k(x) \cos(\omega_k t + \phi_k)$ with dispersion relation:

$$\omega_k = c\sqrt{\lambda_k}$$

In the joint $(\lambda, \omega)$ spectrum, the energy of a traveling wave is concentrated **on the curve** $\omega = c\sqrt{\lambda}$ — a parabola in the $(\sqrt{\lambda}, \omega)$ plane, a straight line in the $(\lambda, \omega^2)$ plane. The slope of this curve directly gives the propagation speed $c$.

A purely temporal spectrum at any fixed sensor location would show a peak at $\omega_k$ — but would not reveal $\lambda_k$ or $c$. A purely spatial eigenmode spectrum at any fixed time would show amplitude at mode $k$ — but would not reveal the oscillation frequency or propagation direction. **Only the joint spectrum reveals the dispersion relation** — and thus the propagation speed.

A filter centered at the point $(k^*, \omega^*)$ on the dispersion curve is a **traveling wave filter** — it selects neural activity propagating at speed $c = \omega^*/\sqrt{\lambda_{k^*}}$ with spatial wavelength $\Lambda = 2\pi/\sqrt{\lambda_{k^*}}$.

#### 3.2.2 Heat diffusion: energy in a wedge

A cortical diffusion process governed by the heat equation:

$$\frac{\partial q}{\partial t} = -\alpha \Delta_S q$$

has eigenmode solutions $q(x,t) = \psi_k(x) e^{-\alpha \lambda_k t}$ — exponential decay in time, with decay rate $\alpha \lambda_k$.

In the temporal Fourier domain, exponential decay corresponds to a Lorentzian spectrum with half-width $\alpha \lambda_k$:

$$|\hat{Q}(k,\omega)|^2 \propto \frac{1}{\omega^2 + (\alpha\lambda_k)^2}$$

High-$\lambda_k$ modes have broad Lorentzian spectra (fast decay → broadband temporal content). Low-$\lambda_k$ modes have narrow spectra (slow decay → narrowband near DC). In the joint $(\lambda, \omega)$ plane, diffusion energy fills a **wedge**: all modes contribute at all frequencies up to a mode-dependent cutoff $\omega^* \sim \alpha\lambda_k$.

This wedge shape is diagnostic: heat diffusion cannot produce a dispersion curve signature, and a traveling wave cannot produce a wedge signature. **A joint $(\lambda, \omega)$ filter that selects a specific wedge region is a diffusion filter** — it measures the diffusive spreading of neural activity at a specific spatial scale.

#### 3.2.3 The physiological speed limit as a filter design constraint

From known neurophysiology — cortico-cortical axonal conduction velocities $c_{\max} \approx 3$–8 m/s and synaptic time constants $\tau_{\min} \approx 2$ ms — the physiologically accessible region in the $(\lambda, \omega)$ plane is:

$$\mathcal{R} = \{(\lambda_k, \omega) : \omega \leq c_{\max}\sqrt{\lambda_k},\; \omega \leq \omega_{\max}\}$$

Any filter placed outside $\mathcal{R}$ is guaranteed to capture noise, artifact, or inverse modeling error — not neural signal. The physiological speed limit defines the **admissible filter domain**: only filters within $\mathcal{R}$ are scientifically valid.

This constraint is model-agnostic — it applies regardless of whether the cortical dynamics are waves, diffusion, or something else. It is derived from known biology, not assumed about the dynamics. It provides a principled exclusion criterion that reduces the false positive rate of any spatiotemporal filter analysis.

#### 3.2.4 The dispersion relation as a filter selection criterion

For a system with a known or empirically estimated dispersion relation $\omega = f(\lambda_k)$, the optimal filter bank for detecting that dynamical regime consists of filters placed along the dispersion curve — a **matched filter bank** tuned to the specific dynamics.

For unknown dynamics, the joint spectrum itself reveals the dispersion relation empirically: the ridge of high power in the $(\lambda, \omega)$ plane is the empirical dispersion curve. Filter placement can then be adapted to the data — a **data-adaptive spatiotemporal filter bank**.

This is a fundamentally new analysis paradigm: rather than assuming a temporal frequency band (delta, theta, alpha, etc.) and filtering there regardless of spatial content, the joint spectrum reveals which combinations of spatial scale and temporal frequency carry neural signal, and filters are placed accordingly.

### 3.3 From sensor data to joint spectrum: the measurement pipeline

The joint $(\lambda, \omega)$ spectrum of the cortical source field cannot be computed directly from sensor data without addressing the electromagnetic inverse problem. The composed transform:

$$A = \Phi^T M \in \mathbb{R}^{K \times n_{\text{ch}}}$$

where $\Phi = [\psi_0, \psi_1, \ldots, \psi_K]$ is the eigenmode matrix and $M$ is the MEG imaging kernel, maps directly from sensor data to eigenmode coefficients:

$$\boldsymbol{\theta}(t) = A \cdot \mathbf{d}(t)$$

The $k$-th row of $A$ is the **sensor-space representation** of spatial filter $k$ — the linear combination of sensor channels that extracts the $k$-th LBO eigenmode amplitude from the raw MEG data. This is the spatial filter implemented in sensor space, without explicitly solving the full inverse problem at every time step.

The joint spectrum is then computed from the eigenmode time series $\theta_k(t)$ via temporal Fourier transform:

$$P(k, \omega) = |\mathcal{F}\{\theta_k\}(\omega)|^2$$

The full spatiotemporal filter output — the joint spectrum — is obtained from raw sensor data via two linear operations: the composed sensor-to-eigenmode transform $A$, and the temporal Fourier transform. The spatial filter and temporal filter are **separable** — they commute and can be applied in either order.

---

## 4. Complex-Valued Filters for Nonstationary Signals: The Wavelet Tensor

### 4.1 The stationarity assumption and its failure

The joint $(\lambda, \omega)$ spectrum via temporal Fourier transform assumes **stationarity** — that the temporal frequency content of each eigenmode is constant over the analysis window. For a 10-second MEG recording during a cognitive task, this assumption fails: neural dynamics are nonstationary, with different $(\lambda, \omega)$ combinations active at different times.

This is the exact same limitation that motivated the wavelet transform in temporal analysis. The temporal Fourier transform gives $|\tilde{d}(\omega)|^2$ — global frequency content — while the Morlet wavelet transform gives $|W(s,t)|^2$ — local frequency content at each time. The spatiotemporal analog requires joint $(\lambda, \omega, t)$ localization.

### 4.2 The complex-valued wavelet tensor

For each eigenmode time series $\theta_k(t)$, apply the Morlet continuous wavelet transform:

$$W_k(s,t) = \int \theta_k(\tau)\, \psi_s^*\!\left(\frac{\tau - t}{s}\right) d\tau \in \mathbb{C}$$

where scale $s$ corresponds to temporal frequency $f = f_0/s$. The collection:

$$\mathbf{W} = \{W_k(s,t)\}_{k=0,\ldots,K;\; s \in \mathcal{S};\; t \in [0,T]} \in \mathbb{C}^{K \times |\mathcal{S}| \times T}$$

is the **complex-valued spatiotemporal wavelet tensor** — the primary data structure of the proposed framework.

Each element $W_k(s,t) \in \mathbb{C}$ is a complex number carrying:

- **Amplitude** $|W_k(s,t)|$: How strongly spatial mode $k$ oscillates at temporal frequency $f \sim 1/s$ at time $t$.
- **Phase** $\arg(W_k(s,t))$: The instantaneous phase of spatial mode $k$ at temporal frequency $f$ at time $t$.

This is the **time-varying joint $(\lambda, \omega)$ spectrum** — the full spatiotemporal spectral decomposition of cortical dynamics, resolved in time.

### 4.3 What the complex values encode

The complex value of $W_k(s,t)$ carries information that the real-valued power $|W_k(s,t)|^2$ discards:

#### 4.3.1 Instantaneous dispersion relation

At each time $t$ and scale $s$, the phase $\phi_k(s,t) = \arg(W_k(s,t))$ across eigenmode index $k$ encodes the **instantaneous dispersion relation**. For a traveling wave with dispersion $\omega = c\sqrt{\lambda_k}$, the phase advances with $k$ at a rate determined by $c$ and the time elapsed since the wave onset. Fitting:

$$\phi_k(s,t) = \phi_0(s,t) + c_{\text{empirical}}(s,t)\sqrt{\lambda_k} \cdot \Delta t$$

gives the instantaneous propagation speed $c_{\text{empirical}}(s,t)$ — a time-varying quantity that characterizes the dynamics at each moment.

#### 4.3.2 Inter-eigenmode phase coherence

The complex cross-spectrum between modes $k$ and $k'$:

$$C_{kk'}(s,t) = \frac{W_k(s,t)\, W_{k'}^*(s,t)}{|W_k(s,t)|\, |W_{k'}(s,t)|}$$

measures the **instantaneous phase coherence** between spatial scales $k$ and $k'$ at temporal frequency $f$ and time $t$. For a coherent traveling wave, $C_{kk'}$ is close to unity with a systematic phase lag. For independent noise, $|C_{kk'}|$ is near zero. For diffusive dynamics, $C_{kk'}$ decays with $|k - k'|$ at a rate determined by the diffusivity.

Phase coherence is more robust to the amplitude distortions introduced by the MEG inverse problem than amplitude itself — the imaging kernel $M$ attenuates amplitudes at high spatial frequencies but does not systematically rotate phases. Phase-based analysis therefore provides a noise-robust window into the dynamical structure of the wavelet tensor.

#### 4.3.3 Cross-frequency coupling in eigenmode space

The wavelet tensor simultaneously resolves multiple temporal scales, enabling cross-frequency coupling analysis within eigenmode space. The phase-amplitude coupling between a slow oscillation at scale $s_{\text{slow}}$ and a fast oscillation at scale $s_{\text{fast}}$:

$$\text{PAC}(k, s_{\text{slow}}, s_{\text{fast}}, t) = |W_k(s_{\text{fast}}, t)|\cdot e^{i\phi_k(s_{\text{slow}}, t)}$$

measures whether the amplitude of fast spatial dynamics in mode $k$ is modulated by the phase of slow dynamics in the same mode. This type of cross-frequency coupling — well-documented in invasive recordings — is measurable in the wavelet tensor without any assumption about its spatial localization.

### 4.4 The wavelet tensor as a filter bank

The wavelet tensor implements a **doubly indexed filter bank** — indexed by spatial eigenmode $k$ and temporal scale $s$ — with each filter being a complex-valued spatiotemporal filter:

$$\mathcal{F}_{k,s}[q](t) = W_k(s,t) = \left\langle q(\cdot, \cdot), \psi_k(\cdot) \otimes \psi_s(\cdot - t) \right\rangle$$

where $\otimes$ denotes the outer product in the joint $(x,t)$ domain. Each filter $\mathcal{F}_{k,s}$ is localized in:
- **Spatial frequency**: around $\lambda_k$, with bandwidth determined by the eigenmode spacing
- **Temporal frequency**: around $f_0/s$, with bandwidth determined by the Morlet wavelet uncertainty ($\sigma_f \sim f_0/(2\pi\sigma_t)$)
- **Time**: around $t$, with localization determined by the Gaussian window width

This is a proper filter bank in the signal processing sense — complete (all signal content is captured across the full $k \times s$ grid), approximately orthogonal (different $(k,s)$ pairs are approximately independent), and interpretable (each filter output corresponds to a specific combination of cortical spatial scale and neural oscillation frequency at a specific time).

### 4.5 Filter bank design: choosing the grid in $(\lambda, \omega)$ space

The filter bank grid — the set of $(k, s)$ pairs for which filters are computed — is a design choice with significant practical implications.

**Uniform $k$ spacing**: Compute wavelet coefficients for all $k = 0, 1, \ldots, K$. Complete but computationally expensive. Appropriate when the active spatial scales are unknown.

**Log-spaced $\lambda$ sampling**: Sample eigenmode index $k$ on a logarithmic scale in $\lambda_k$, analogous to the log-frequency spacing standard in auditory filterbanks. Appropriate for signals with scale-invariant (power-law) spatial spectra.

**Dispersion-adapted grid**: Place filter centers along the empirical dispersion curve $\omega = \hat{c}\sqrt{\lambda_k}$ estimated from the data. A **matched filter bank** for the dominant dynamical regime.

**Physiological constraint grid**: Restrict filter centers to the accessible region $\mathcal{R}$, excluding all $(k, \omega)$ combinations outside the physiological speed limit. An **admissible filter bank** that excludes non-neural signal by construction.

**Hierarchical filter bank**: Coarse-scale filters ($k$ small, $s$ large) for detecting global low-frequency dynamics; fine-scale filters for detecting focal high-frequency events. Analogous to a wavelet multiresolution analysis but in 2D $(\lambda, \omega)$ space.

---

## 5. Proposed Spatiotemporal Filter Library

### 5.1 Motivation and design philosophy

The temporal filter tools of neuroscience — MNE-Python's `filter_data`, `tfr_morlet`, EEGLAB's `pop_eegfiltnew`, FieldTrip's `ft_freqanalysis` — are widely used precisely because they are simple to apply, well-documented, and produce interpretable outputs. A researcher can apply a temporal bandpass filter with a single function call without understanding the underlying signal processing theory.

The proposed **spatiotemporal filter library** aims to achieve the same accessibility for joint $(\lambda, \omega)$ filtering. A researcher should be able to apply a spatial eigenmode filter, a joint spatiotemporal filter, or a full wavelet tensor decomposition with comparably simple function calls, receiving outputs that are directly interpretable in terms of cortical spatial scale, temporal frequency, and dynamical regime.

The library is designed around three principles:

**Symmetry**: Spatial and temporal filtering are treated symmetrically. Every temporal filter has a spatial analog and a joint spatiotemporal generalization.

**Composability**: Filters are composable — a spatial filter and a temporal filter can be combined to produce a joint spatiotemporal filter. Filter banks can be nested and hierarchically organized.

**Interpretability**: Every filter output is labeled with its spatial scale ($\lambda_k$ or approximate spatial wavelength $\Lambda_k \sim 1/\sqrt{\lambda_k}$), temporal frequency ($f$), and time ($t$ for wavelet filters). The physical meaning of each output is always explicit.

### 5.2 Core filter classes

#### 5.2.1 Spatial eigenmode filters

The fundamental spatial filtering unit. Input: cortical source field $q(x)$ or sensor data $\mathbf{d}$ with imaging kernel $M$. Output: scalar eigenmode coefficient $\hat{q}_k$.

```python
# Conceptual API
spatial_filter = SpatialEigenmodeFilter(
    surface=cortical_mesh,          # FreeSurfer surface
    n_modes=200,                    # number of LBO eigenmodes
    mode_index=k,                   # specific mode to extract
    imaging_kernel=M                # MEG inverse operator
)
q_hat_k = spatial_filter.apply(sensor_data)  # scalar time series
```

**Filter bank variant**: Apply all $K$ spatial filters simultaneously via the composed transform $A = \Phi^T M$:

```python
spatial_filterbank = SpatialEigenmodeFilterBank(
    surface=cortical_mesh,
    n_modes=K,
    imaging_kernel=M
)
theta = spatial_filterbank.apply(sensor_data)  # (K, T) array
```

#### 5.2.2 Temporal filters (standard, for reference)

Standard temporal filters implemented on eigenmode time series, for direct comparison with sensor-space temporal filtering:

```python
temporal_filter = TemporalBandpassFilter(
    freq_range=(8, 12),  # Hz
    sfreq=1000           # sampling rate
)
theta_alpha = temporal_filter.apply(theta)  # (K, T) bandpassed eigenmode series
```

#### 5.2.3 Joint spatiotemporal Fourier filters

For stationary signals, select a rectangular region in $(\lambda, \omega)$ space:

```python
joint_filter = JointSpectralFilter(
    lambda_range=(lambda_min, lambda_max),   # spatial frequency range
    omega_range=(omega_min, omega_max),      # temporal frequency range
    surface=cortical_mesh,
    imaging_kernel=M
)
output = joint_filter.apply(sensor_data)    # filtered cortical source field
```

#### 5.2.4 Dispersion-matched filter banks

For wave detection, a filter bank matched to a specific dispersion relation:

```python
wave_filterbank = DispersionMatchedFilterBank(
    dispersion='wave',          # omega = c * sqrt(lambda)
    velocity_range=(2, 8),      # m/s physiological range
    n_filters=20,               # number of velocity-tuned filters
    surface=cortical_mesh,
    imaging_kernel=M
)
wave_amplitude, wave_phase = wave_filterbank.apply(sensor_data)
# Output: (n_filters, T) amplitude and phase arrays
```

#### 5.2.5 The spatiotemporal wavelet tensor

The full nonstationary filter bank — the primary analysis tool for task-evoked and resting-state MEG:

```python
wavelet_tensor = SpatiotemporalWaveletTensor(
    surface=cortical_mesh,
    n_modes=K,
    freqs=np.logspace(np.log10(1), np.log10(100), 40),  # Hz
    n_cycles=freqs / 2.,        # frequency-dependent Morlet width
    imaging_kernel=M,
    noise_cov=Sigma_n           # for SNR calibration
)
W = wavelet_tensor.compute(sensor_data)
# W: complex (K, n_freqs, T) tensor
# W.amplitude: |W|
# W.phase: arg(W)
# W.power: |W|^2
# W.snr: |W|^2 / noise_floor(k, f)
```

#### 5.2.6 Physiological constraint mask

Applied to any filter output to restrict to the admissible region:

```python
physio_mask = PhysiologicalConstraintMask(
    c_max=8.0,          # m/s maximum propagation speed
    tau_min=0.002,      # seconds minimum synaptic time constant
    eigenvalues=lambda_k
)
W_masked = physio_mask.apply(W)  # zeros out supraphysiological (k, f) cells
```

### 5.3 Analysis tools built on the filter library

The filter library enables a set of high-level analysis tools, each implemented as operations on the wavelet tensor:

#### 5.3.1 Empirical dispersion relation estimation

```python
dispersion = EmpiricalDispersionEstimator(
    method='spectral_ridge'  # or 'phase_slope', 'group_velocity'
)
omega_k, c_empirical = dispersion.fit(W)
# Returns: dispersion curve omega(k) and estimated propagation speed
```

#### 5.3.2 Wave vs. diffusion discrimination

```python
dynamics_classifier = DynamicsClassifier(
    null_models=['wave', 'diffusion', 'reaction_diffusion'],
    method='spectral_shape'   # fits joint spectrum to each model
)
regime, likelihood_ratio = dynamics_classifier.classify(W)
```

#### 5.3.3 Inter-eigenmode phase coherence

```python
coherence = InterEigenmodeCohrence(
    reference_mode=k_ref,
    method='phase_locking_value'
)
PLV_k = coherence.compute(W)   # (K, n_freqs, T) phase locking values
```

#### 5.3.4 Source recovery at onset

```python
source_recovery = OnsetSourceRecovery(
    onset_detector='wavelet_amplitude_threshold'
)
u_onset = source_recovery.recover(W, surface=cortical_mesh)
# Returns vertex-space source map at wave onset
```

### 5.4 Implementation architecture

The library is implemented in two complementary layers:

**C++ core (Geometry Central-based)**:
- LBO eigendecomposition on cortical mesh
- Connection Laplacian for vector eigenmode extension
- Trivial connections and vector heat method for orientation consistency
- Discrete exterior calculus for Hodge decomposition
- Composed transform $A = \Phi^T M$ construction

**Python interface (MNE-Python-compatible)**:
- MEG preprocessing and forward modeling via MNE-Python
- Filter application and wavelet tensor computation via NumPy/SciPy
- Visualization via MNE-Python surface plotting
- Analysis tools (dispersion estimation, coherence, classification)
- Compatibility with MNE `SourceEstimate` and `EpochsTFR` data structures

The interface between the two layers is the composed transform matrix $A$ — a NumPy array exported from the C++ geometric computation and used directly in the Python analysis pipeline.

### 5.5 Comparison with existing tools

| Feature | MNE-Python temporal filters | GBF (Wang et al., 2026) | Proposed library |
|---|---|---|---|
| Temporal filtering | ✓ Full | ✗ Post-hoc only | ✓ Full |
| Spatial filtering | ✗ None | ✓ Eigenmode basis | ✓ Full eigenmode bank |
| Joint (λ,ω) filtering | ✗ | ✗ | ✓ Core feature |
| Nonstationary analysis | ✓ Wavelet (temporal) | ✗ | ✓ Wavelet tensor |
| Complex-valued output | ✓ Hilbert/wavelet | ✗ | ✓ Throughout |
| Phase information | ✓ Temporal phase | ✗ | ✓ Spatiotemporal phase |
| Dispersion analysis | ✗ | ✗ | ✓ |
| Wave/diffusion discrimination | ✗ | ✗ | ✓ |
| Physiological constraint mask | ✗ | ✗ | ✓ |
| Orientation consistency | N/A | ✗ Scalar only | ✓ Trivial connections |
| Vector source support | N/A | ✗ | ✓ Connection Laplacian |

---

## 6. Proposed Validation

### 6.1 Synthetic validation: filter selectivity and orthogonality

**Objective**: Verify that the spatial eigenmode filter bank has the expected selectivity, orthogonality, and spatial frequency resolution on realistic cortical meshes.

**Protocol**:
1. Generate synthetic cortical source fields as linear combinations of known eigenmodes: $q(x) = \sum_{k \in S} a_k \psi_k(x)$ for a sparse support set $S$.
2. Apply the spatial filter bank and verify: (a) $\hat{q}_k \approx a_k$ for $k \in S$; (b) $\hat{q}_k \approx 0$ for $k \notin S$; (c) Parseval's theorem holds numerically.
3. Generate synthetic traveling waves with known dispersion $\omega = c\sqrt{\lambda_k}$ at three velocities (2, 5, 8 m/s) and three frequencies (10, 20, 40 Hz).
4. Apply the wavelet tensor and verify that the joint power $|W_k(s,t)|^2$ concentrates on the dispersion curve at the correct $(k, \omega)$ location.
5. Generate synthetic heat diffusion at three diffusivities and verify that the joint power fills the predicted wedge shape.

**Metrics**: Reconstruction error $\|\hat{\mathbf{q}} - \mathbf{a}\|_2$, filter cross-talk matrix $|(\Phi^T\Phi - I)|_{\max}$, dispersion curve recovery error.

### 6.2 Comparison with standard temporal filtering on canonical neural rhythms

**Objective**: Demonstrate that joint $(\lambda, \omega)$ filtering provides strictly more information than temporal-only filtering at the same computational cost.

**Protocol**:
1. Apply both standard temporal bandpass filtering (alpha: 8–12 Hz, beta: 13–30 Hz, gamma: 30–80 Hz) and joint $(\lambda, \omega)$ filtering to resting-state MEG (HCP dataset, n = 80).
2. For temporal-only filtering: compute sensor-space power spectral density per band.
3. For joint filtering: compute the wavelet tensor and integrate over the corresponding temporal frequency band to obtain the spatial eigenmode power $P(k) = \int_{f_1}^{f_2} |W_k(f,t)|^2 df\, dt$.
4. Show that the spatial eigenmode power profile $P(k)$ reveals which spatial scales carry each band's power — information unavailable from temporal filtering alone.
5. Test whether the dispersion structure within each band is consistent with known propagation speeds (alpha: ~1–3 m/s for traveling waves; beta: ~3–8 m/s; gamma: ~0.1–1 m/s for local oscillations).

**Metrics**: Spatial spectral power profile per band, dispersion exponent per band, comparison of dispersion velocities with ECoG literature.

### 6.3 Task-evoked filter response: motor beta waves

**Objective**: Demonstrate that the wavelet tensor captures the known spatiotemporal dynamics of motor beta oscillations, providing richer characterization than temporal filtering alone.

**Dataset**: HCP MEG motor task (right-hand finger tapping, n = 80).

**Protocol**:
1. Compute the wavelet tensor for each participant during movement preparation and execution.
2. Extract the beta-band ($13$–30 Hz) slice of the tensor: $W_k(s_\beta, t)$.
3. Track the temporal evolution of $|W_k(s_\beta,t)|$ across eigenmode index $k$ as a function of trial time — show the beta power redistribution from pre-movement (high amplitude, low-$k$ dominant) to movement (beta desynchronization) to post-movement rebound.
4. Compute the instantaneous phase $\phi_k(s_\beta, t)$ across $k$ during the rebound period — test whether the phase advances with $k$ according to a dispersion relation consistent with known motor cortex propagation speeds.
5. Compare with purely temporal beta power analysis (standard event-related desynchronization, ERD) — quantify what additional information the joint analysis provides.

**Metrics**: Correlation of joint spectral dynamics with known ERD/ERS profiles, dispersion relation recovery, propagation speed estimate.

### 6.4 Resting-state: empirical $(\lambda, \omega)$ spectrum characterization

**Objective**: Provide the first systematic characterization of the joint $(\lambda, \omega)$ power spectrum of human resting-state MEG, across participants and frequency bands.

**Protocol**:
1. Compute the time-averaged joint power spectrum $P(k, f) = \langle |W_k(f,t)|^2 \rangle_t$ for each of n = 80 HCP participants.
2. Compute the group-average $\bar{P}(k,f)$ and its variance across participants.
3. For each canonical frequency band, extract the spatial spectral profile $P(k)$ and fit a power law $P(k) \propto \lambda_k^{-\gamma}$ — the spatial spectral exponent $\gamma$ characterizes the cortical spatial smoothness of resting dynamics.
4. Test for dispersion structure: fit $\omega = c\lambda_k^\alpha$ to the spectral ridge in each band. $\alpha = 0.5$ indicates wave-like dynamics; $\alpha < 0.5$ indicates sub-diffusive dynamics.
5. Apply the physiological constraint mask and quantify the fraction of total power within the admissible region $\mathcal{R}$ — a measure of the physical plausibility of the recovered signal.
6. Compare spatial spectral exponents with analogous measures from ECoG and fMRI (Pang et al., 2023) to establish cross-modal consistency.

**Metrics**: Spatial spectral exponent $\gamma$, dispersion exponent $\alpha$, admissible power fraction, inter-subject variability of all measures.

### 6.5 Clinical application: ictal and interictal filter signatures

**Objective**: Demonstrate that different epileptic dynamical regimes (interictal spikes vs. ictal spread) produce distinguishable signatures in the joint $(\lambda, \omega)$ spectrum.

**Dataset**: HDEEG-IED-SurgOutcome dataset (n = 24, 257-channel HD-EEG).

**Protocol**:
1. Compute the wavelet tensor for interictal spike events (time-locked to GFP peak).
2. Compute the wavelet tensor for ictal periods (seizure onset).
3. Compare the joint $(\lambda, \omega)$ signature:
   - Interictal spikes: predicted to be broadband in time (transient), concentrated at specific spatial scales near the epileptogenic zone (focal in $k$-space).
   - Ictal spread: predicted to show a time-evolving dispersion signature as the seizure wavefront propagates — energy sweeping from high-$k$ to low-$k$ as the spatial pattern of activation broadens.
4. Extract the propagation speed from the wavelet tensor during ictal spread: $c(t) = \Delta\omega / \Delta\sqrt{\lambda_k}$ at the spectral ridge.
5. Compare recovered propagation speeds with known ictal spread velocities from ECoG (~1–10 mm/s for slow cortical spread; ~1 m/s for fast propagation along white matter pathways).

**Metrics**: Joint spectral discriminability between interictal and ictal regimes, propagation speed estimate, correspondence with surgical outcome.

### 6.6 Filter library benchmarking

**Objective**: Characterize computational performance and numerical accuracy of the filter library across mesh resolutions and eigenmode numbers.

**Protocol**:
1. Benchmark LBO eigendecomposition time vs. mesh resolution ($n_{\text{vert}} = 2562$, $5124$, $10242$, $20484$ for fsaverage3–6) and mode number $K$ (50, 100, 200, 300, 500).
2. Benchmark wavelet tensor computation time vs. $K$, number of frequency scales $|\mathcal{S}|$, and time samples $T$.
3. Verify numerical accuracy of the composed transform $A = \Phi^T M$ against direct eigenmode projection on reconstructed vertex-space sources.
4. Profile memory usage of the complex wavelet tensor $W \in \mathbb{C}^{K \times |\mathcal{S}| \times T}$ and propose sparse storage strategies for large $K$ and $T$.

**Target performance**: Full wavelet tensor computation for $K = 200$, $|\mathcal{S}| = 40$, $T = 10{,}000$ samples (10 s at 1 kHz) in under 60 seconds on a standard workstation (16-core CPU, 64 GB RAM).

---

## 7. Discussion

### 7.1 Relation to GBF and spatial filtering literature

GBF (Wang et al., 2026) demonstrates that LBO eigenmodes are an effective spatial prior for MEG/EEG source imaging, validating the eigenmode spatial filter concept empirically. The proposed library builds on this foundation by:

1. Making the filter interpretation of eigenmodes explicit and primary — treating GBF not as an inverse method but as a spatial filter bank whose outputs are the starting point for analysis.

2. Adding temporal filter structure via the wavelet tensor — extending GBF from static spatial decomposition to dynamic spatiotemporal decomposition.

3. Making the filter outputs complex-valued — carrying phase as well as amplitude, enabling coherence and phase-based wave analysis that GBF's real-valued source estimates cannot support.

4. Providing the complete filter library infrastructure — reusable, composable, documented tools rather than a single monolithic analysis pipeline.

The relationship to Pang et al. (2023) is also direct: that work demonstrated the representational efficiency of LBO eigenmodes for fMRI spatial data. The present proposal operationalizes this as a filter bank for MEG temporal dynamics — moving from static spatial representation to dynamic spatiotemporal filtering.

### 7.2 Relation to graph signal processing

The spatial eigenmode filter bank on the cortical LBO is a special case of **graph signal processing** (GSP) on the cortical mesh graph (Shuman et al., 2013; Ortega et al., 2018). The LBO cotangent discretization defines a weighted graph Laplacian; its eigenvectors are the graph Fourier basis; the filter bank is the graph Fourier filter bank.

GSP provides a rich theoretical framework for filter design on graphs — including graph wavelets, spectral graph convolutions, and localized graph filters — that is directly applicable to the proposed library. In particular:

- **Spectral graph convolution**: Filters defined as $h(\Lambda)$ — functions of the eigenvalue matrix — can be applied efficiently in eigenmode space without explicitly computing the full spatial filter. This is the graph signal processing generalization of frequency-domain filtering.
- **Localized graph filters**: Filters that are spatially localized on the cortical graph — concentrating the filter's spatial response within a geodesic ball of radius $r$ around a seed vertex — provide a spatially adaptive alternative to globally defined eigenmode filters.
- **Graph wavelets**: Multi-scale spatial analysis tools that are simultaneously localized in cortical space and spatial frequency, analogous to continuous wavelets in time.

The proposed library incorporates these GSP tools as additional filter classes alongside the eigenmode filter bank.

### 7.3 Connection to the vector eigenmode framework

The present proposal describes spatiotemporal filters for **scalar** cortical source fields — the output of fixed-orientation MEG inverse methods. As described in the companion proposal (Vector Eigenmode Source Mapping, this volume), fixed-orientation source estimates suffer from sulcal orientation inconsistencies that introduce geometric artifacts into phase-based analysis.

The spatiotemporal filter library is fully compatible with the vector eigenmode extension: replacing scalar LBO eigenmodes with connection Laplacian eigenmodes and Hodge scalar potentials as the spatial basis provides orientation-consistent filter outputs. The wavelet tensor $W_k(s,t)$ is identical in structure; only the spatial basis functions change. The library architecture is designed to support both scalar and vector eigenmode variants through a common interface.

### 7.4 Limitations

**Stationarity of the spatial basis**: The LBO eigenmodes are computed once from the structural MRI and treated as fixed throughout the MEG recording. This assumes the cortical geometry is the relevant spatial basis for all dynamical states — a reasonable but untested assumption. State-dependent spatial bases (e.g., functional eigenmodes computed from data covariance) are a natural extension.

**Truncation at $K$ modes**: The filter bank is finite — it represents spatial frequencies up to the $K$-th eigenmode. The truncation is justified by MEG spatial resolution limits and empirical spectral decay, but it means spatial structures finer than $\sim 1/\sqrt{\lambda_K}$ are invisible to the filter bank.

**Linear filter assumptions**: All proposed filters are linear — they compute inner products with fixed basis functions. Nonlinear spatiotemporal interactions — such as phase-amplitude coupling between different spatial scales — require nonlinear extensions of the filter framework.

---

## 8. Conclusion

We have proposed a framework for spatiotemporal filtering of MEG/EEG data grounded in the geometry of the cortical surface and the physics of electromagnetic neural measurement. The core argument, built from simple to complex, is as follows:

Temporal filters decompose neural time series by oscillation frequency, providing interpretable, physically meaningful signal components. Spatial eigenmode filters — projections onto LBO eigenmodes of the individual cortical surface — provide the symmetric spatial analog, decomposing cortical source fields by spatial scale. Applied statically, spatial filters characterize the instantaneous spatial content of neural activity but cannot distinguish between static, propagating, oscillating, or diffusing patterns.

Joining the two filter axes into the joint $(\lambda, \omega)$ domain produces spatiotemporal filters whose outputs are sensitive to **dynamics** — the relationship between spatial scale and temporal frequency that characterizes traveling waves (dispersion curve), heat diffusion (wedge), and other dynamical regimes. The physiological speed limit — derived from axonal conduction velocities and synaptic time constants — defines the admissible filter domain and provides a model-agnostic exclusion criterion for non-neural content.

For nonstationary neural signals — the general case in cognitive neuroscience and clinical MEG — the joint filter bank is extended to the **complex-valued spatiotemporal wavelet tensor** $W_k(s,t) \in \mathbb{C}$, jointly resolving spatial scale, temporal frequency, and time with full amplitude and phase information. This tensor is the primary observable: it contains, in compact form, everything that can be measured about the spatiotemporal dynamics of cortical neural activity from MEG sensor data, within the resolution limits of the sensor array and the geometric constraints of the cortical manifold.

The proposed **spatiotemporal filter library** makes these tools accessible — through simple, composable, well-documented function calls compatible with standard MEG/EEG analysis software — to the broad neuroscience community, enabling a symmetric treatment of space and time in MEG/EEG analysis for the first time.

---

## References

Baillet, S. (2017). Magnetoencephalography for brain electrophysiology and imaging. *Nature Neuroscience*, 20(3), 327–339.

Crane, K., Desbrun, M., & Schröder, P. (2010). Trivial connections on discrete surfaces. *Computer Graphics Forum*, 29(5), 1525–1533.

Desbrun, M., Hirani, A. N., Leok, M., & Marsden, J. E. (2005). Discrete exterior calculus. *arXiv preprint math/0508341*.

Gramfort, A., et al. (2013). MEG and EEG data analysis with MNE-Python. *Frontiers in Neuroinformatics*, 7, 267.

He, B., Sohrabpour, A., Brown, E., & Liu, Z. (2018). Electrophysiological source imaging: a noninvasive window to brain dynamics. *Annual Review of Biomedical Engineering*, 20, 171–196.

Morlet, J., Arens, G., Fourgeau, E., & Giard, D. (1982). Wave propagation and sampling theory — Part I: Complex signal and scattering in multilayered media. *Geophysics*, 47(2), 203–221.

Ortega, A., Frossard, P., Kovačević, J., Moura, J. M. F., & Vandergheynst, P. (2018). Graph signal processing: Overview, challenges, and applications. *Proceedings of the IEEE*, 106(5), 808–828.

Pang, J. C., Aquino, K. M., Oldehinkel, M., Robinson, P. A., Fulcher, B. D., Breakspear, M., & Fornito, A. (2023). Geometric constraints on human brain function. *Nature*, 618, 566–574.

Roberts, J. A., Gollo, L. L., Abeysuriya, R. G., Roberts, G., Mitchell, P. B., Woolrich, M. W., & Breakspear, M. (2019). Metastable brain waves. *Nature Communications*, 10(1), 1056.

Sharp, N., Soliman, Y., & Crane, K. (2019). The vector heat method. *ACM Transactions on Graphics*, 38(3), 1–19.

Shuman, D. I., Narang, S. K., Frossard, P., Ortega, A., & Vandergheynst, P. (2013). The emerging field of signal processing on graphs. *IEEE Signal Processing Magazine*, 30(3), 83–98.

Vidaurre, D., et al. (2018). Spontaneous cortical activity transiently organises into frequency specific phase-coupling networks. *Nature Communications*, 9(1), 2987.

Wang, S., Lou, K., Wei, C., et al. (2026). A geometry aware framework enhances noninvasive mapping of whole human brain dynamics. *arXiv:2604.25592v1*.

---

*Proposal prepared for review. Library implementation targets Geometry Central (C++) for geometric computation and MNE-Python (Python) for MEG analysis, with a NumPy/SciPy-based filter computation layer and a public API compatible with standard MNE data structures.*
