function [W, Freqs] = bst_eigenmodes_wavelet(Coeffs, sfreq, Freqs, varargin)
% BST_EIGENMODES_WAVELET: Complex Morlet wavelet tensor of eigenmode coefficients.
%
% USAGE:  [W, Freqs] = bst_eigenmodes_wavelet(Coeffs, sfreq, Freqs)
%         [W, Freqs] = bst_eigenmodes_wavelet(Coeffs, sfreq, [], 'MorletFc',1, 'MorletFwhmTc',3)
%
% DESCRIPTION:
%     Applies the complex Morlet continuous wavelet transform to each eigenmode
%     coefficient time series, producing the time-resolved (lambda, omega, t)
%     tensor W_k(s,t). The result is COMPLEX: |W| is amplitude, arg(W) is phase.
%     This wraps Brainstorm's validated morlet_transform (squared='n', the
%     un-squared complex coefficients), which returns [K x nTime x nFreq]
%     (permuted internally); the conjugate is applied so phase follows the
%     standard positive-rotation analytic-signal convention.
%
% INPUTS:
%     Coeffs : [K x nTime] eigenmode coefficient time series.
%     sfreq  : sampling rate (Hz).
%     Freqs  : [1 x nFreq] frequencies (Hz). If empty, a default log-spaced grid
%              logspace(log10(2), log10(min(100, 0.4*sfreq)), 40) is used.
%
% OPTIONS (name-value):
%     'MorletFc'     : Morlet central frequency (default 1).
%     'MorletFwhmTc' : Morlet time resolution / FWHM in seconds (default 3).
%
% OUTPUTS:
%     W     : [K x nTime x nFreq] complex wavelet tensor.
%     Freqs : [1 x nFreq] frequencies actually used.
%
% SEE ALSO: morlet_transform, bst_eigenmodes_transform, process_eigenmodes_wavelet

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Diellor Basha, 2026

%% ===== PARSE OPTIONS =====
MorletFc = 1; MorletFwhmTc = 3;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'morletfc',     MorletFc     = varargin{i+1};
        case 'morletfwhmtc', MorletFwhmTc = varargin{i+1};
    end
end

Coeffs = double(Coeffs);
[~, nTime] = size(Coeffs);

%% ===== FREQUENCY GRID =====
if nargin < 3 || isempty(Freqs)
    fhiDefault = min(100, 0.4*sfreq);
    if fhiDefault <= 2
        error('bst_eigenmodes_wavelet: sampling rate too low for the default 2 Hz frequency floor; pass an explicit Freqs vector.');
    end
    Freqs = logspace(log10(2), log10(fhiDefault), 40);
end
Freqs = Freqs(:)';

%% ===== MORLET CWT (complex) =====
t = (0:nTime-1) / sfreq;
% morlet_transform computes in [K x nFreq x nTime] and permutes internally,
% so it returns [K x nTime x nFreq]. 'n' keeps the complex coefficients; the
% conjugate aligns phase with the standard positive-rotation analytic-signal convention.
W = conj(morlet_transform(Coeffs, t, Freqs, MorletFc, MorletFwhmTc, 'n'));
end
