function Gjs = bst_eigfilter_jtv_evaluate(g, domain, Lambda, axis, NFFT)
% BST_EIGFILTER_JTV_EVALUATE: Evaluate a dynamic (time-vertex) eigenfilter kernel into the JOINT-SPECTRAL
% domain [K x NFFT], ready to multiply the output of MANIFOLD_JFT. Handles the eigen-time <-> eigen-
% frequency conversion (the GSPBox 'ts'/'js' distinction): an eigen-time kernel g(lambda,t) is FFT'd along
% the time axis into the joint spectrum; an eigen-frequency kernel g(lambda,omega) is used directly.
%
% USAGE:  Gjs = bst_eigfilter_jtv_evaluate(g, 'ts', Lambda, tLag, NFFT)   % tLag = time-lag axis (s)
%         Gjs = bst_eigfilter_jtv_evaluate(g, 'js', Lambda, omega, NFFT)  % omega = temporal freq (Hz)
% INPUTS:
%   g      : kernel handle g(lambda, x) -> [K x numel(x)]  (from bst_eigfilter_design_<name>)
%   domain : 'ts' (eigen-time) | 'js' (eigen-frequency)
%   Lambda : [K x 1] eigenvalues
%   axis   : [1 x N] 'ts' -> time-lag axis (seconds); 'js' -> temporal-frequency axis (Hz)
%   NFFT   : number of temporal bins (default numel(axis)); the 1/sqrt(NFFT) matches manifold_jft
% OUTPUT:
%   Gjs    : [K x NFFT] joint-spectral gains (multiply manifold_jft(...) elementwise)
%
% SEE ALSO: manifold_jft, manifold_ijft, bst_eigfilter_design_diffusion/wave/dampedwave/kleingordon
%
% Authors: Diellor Basha, 2026
if (nargin < 5) || isempty(NFFT); NFFT = numel(axis); end
G = g(Lambda(:), axis(:)');                          % [K x numel(axis)]
switch lower(domain)
    case 'ts'                                        % eigen-time -> joint-spectral (FFT along time)
        Gjs = fft(G, NFFT, 2) / sqrt(NFFT);
    case 'js'                                        % already eigen-frequency (g defined on omega)
        Gjs = G;
    otherwise
        error('bst_eigfilter_jtv_evaluate: unknown domain ''%s'' (use ''ts'' or ''js'').', domain);
end
end
