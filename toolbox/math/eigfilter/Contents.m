% EIGFILTER  Spectral-filter kernels for Laplace-Beltrami eigenmodes (bst_eigfilter)
%
% Filters along the eigenvalue (spatial-frequency) axis, the spatial-spectral
% counterpart of bst_freqfilter (temporal frequency). Each kernel is an analytic
% factory returning a function handle g = @(lambda) ...; a vector-valued scale
% parameter returns a cell-array filterbank.
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
% Consumed by: bst_eigenmode_prior (prior R), bst_eigenmodes_filter_gain (analysis).
%
% Authors: Diellor Basha, 2026
