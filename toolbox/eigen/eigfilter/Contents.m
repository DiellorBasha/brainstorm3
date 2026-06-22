% EIGFILTER  Spatial-spectral filter kernels for the eigen module (bst_eigfilter)
%
% The spatial "wavelet" kernel library of toolbox/eigen/ (the morlet_wavelet analogue):
% gain shapes h(lambda) along the operator eigenvalue (spatial-frequency) axis, the
% spatial counterpart of bst_freqfilter (temporal frequency). Operator-agnostic -- the
% same kernels apply to any eigenvalue axis (Laplace-Beltrami / Connection Laplacian /
% Dirac). Each kernel is an analytic factory returning a handle g = @(lambda) ...; a
% vector-valued scale parameter returns a cell-array filterbank.
%
% Registry / apply:
%   bst_eigfilter_kernel    - name string -> factory; 'list' / 'info'
%   bst_eigfilter_evaluate  - evaluate a handle (or bank) on eigenvalues
%   bst_eigfilter_compose   - serial composition (pointwise product of gains)
%
% Kernel factories (bst_eigfilter_design_<name>):
%   flat          - g = 1
%   power         - g = lambda^-alpha
%   log           - g = -log(lambda)            (GBF 2026; prior rescales)
%   heat          - g = exp(-t*lambda)          (low-pass / diffusion)
%   inverse_heat  - g = min(exp(+t*lambda), maxgain)   (sharpening)
%   tikhonov      - g = 1/(1+beta*lambda)       (low-pass / membrane)
%   ideal         - g = 1[lo <= lambda <= hi]   (brick-wall band)
%   matern        - g = (kappa^2 + lambda)^-nu  (SPDE / Gaussian field)
%   mexhat        - g = (t*lambda).*exp(-t*lambda)     (band-pass; bank-capable)
%   dog           - g = exp(-t1*lambda) - exp(-t2*lambda)  (band-pass)
%
% Consumed by the file-based eigen filter path: bst_eigen 'filter' -> bst_eigenfilter.
% NOTE: the legacy scalar-LBO eigenmode cluster (bst_eigenmodes_* + process_eigenmodes_*,
% incl. bst_eigenmodes_filter_gain) and the Dirac-wavelet cluster were retired to
% dev/experimental; bst_eigenmode_prior remains only for the deprecated bst_inverse_eigenmodes
% inverses. The live source mapping is bst_dirac / bst_inverse_dirac, which do not use this library.
%
% Authors: Diellor Basha, 2026
