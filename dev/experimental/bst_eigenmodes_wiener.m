function Filtered = bst_eigenmodes_wiener(Coeffs, sfreq, Gain, GainFreqs, varargin)
% BST_EIGENMODES_WIENER: Apply a per-mode magnitude gain G(k,f) to coefficients via the FFT.
%
% USAGE:  Filtered = bst_eigenmodes_wiener(Coeffs, sfreq, Gain, GainFreqs)
%         Filtered = bst_eigenmodes_wiener(Coeffs, sfreq, Gain, GainFreqs, 'Mirror', true)
%
% DESCRIPTION:
%     Applies a frequency-dependent, per-eigenmode magnitude response (e.g. a
%     Wiener gain) to a coefficient time series, in the frequency domain. For
%     each mode k, the gain Gain(k,:) defined on the grid GainFreqs is
%     interpolated onto the (folded) FFT frequency grid and applied as a
%     zero-phase magnitude response: Filtered = real(ifft(fft(Coeffs).*H)).
%     Because the gain is interpolated on the folded (absolute) frequency, H is
%     Hermitian-symmetric, so the output is real and phase is preserved. The
%     signal is mirrored at both ends before the FFT (and trimmed after) to
%     suppress circular-convolution edge artifacts, as in bst_bandpass_fft.
%
% INPUTS:
%     Coeffs    : [K x nTime] real coefficient time series (rows = eigenmodes).
%     sfreq     : sampling frequency (Hz).
%     Gain      : [K x nGainFreq] magnitude gain per mode per frequency, in [0,1].
%     GainFreqs : [1 x nGainFreq] frequencies (Hz) for Gain's columns, ascending.
%
% OPTIONS (name-value):
%     'Mirror'  : reflect the signal at both ends before the FFT, default true.
%
% OUTPUT:
%     Filtered  : [K x nTime] real, same size as Coeffs.
%
% SEE ALSO: bst_eigenmodes_noisefloor, process_eigenmodes_wiener, bst_bandpass_fft

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
Mirror = true;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'mirror', Mirror = logical(varargin{i+1});
    end
end

Coeffs    = double(Coeffs);
Gain      = double(Gain);
GainFreqs = double(GainFreqs(:)');     % row

[K, nTime] = size(Coeffs);
if size(Gain, 1) ~= K
    error('bst_eigenmodes_wiener: Gain must have one row per eigenmode (got %d rows for %d modes).', size(Gain,1), K);
end
if numel(GainFreqs) ~= size(Gain, 2)
    error('bst_eigenmodes_wiener: GainFreqs length (%d) must match the number of Gain columns (%d).', numel(GainFreqs), size(Gain,2));
end
if numel(GainFreqs) < 2 || any(diff(GainFreqs) <= 0)
    error('bst_eigenmodes_wiener: GainFreqs must be strictly ascending with at least 2 points.');
end

%% ===== MIRROR (edge handling) =====
domirror = Mirror && (nTime >= 2);
if domirror
    x = [fliplr(Coeffs), Coeffs, fliplr(Coeffs)];   % [K x 3*nTime]
else
    x = Coeffs;
end
N = size(x, 2);

%% ===== BUILD ZERO-PHASE MAGNITUDE RESPONSE H(k,f) =====
% Folded (absolute) frequency of each FFT bin, in [0, Nyquist].
fvec = (0:N-1) * (sfreq / N);
hi   = fvec > (sfreq / 2);
fvec(hi) = sfreq - fvec(hi);                          % fold negative-frequency bins
% Interpolate each mode's gain onto the folded grid (columns = modes).
Hcols = interp1(GainFreqs, Gain.', fvec(:), 'linear');   % [N x K], NaN outside GainFreqs
below = fvec(:) < GainFreqs(1);
above = fvec(:) > GainFreqs(end);
if any(below), Hcols(below, :) = repmat(Gain(:,1).',   sum(below), 1); end
if any(above), Hcols(above, :) = repmat(Gain(:,end).', sum(above), 1); end
Hcols = max(0, min(1, Hcols));                        % clamp to [0,1]
H = Hcols.';                                          % [K x N], Hermitian-symmetric

%% ===== APPLY =====
Y = real(ifft(fft(x, [], 2) .* H, [], 2));

%% ===== TRIM MIRROR =====
if domirror
    Filtered = Y(:, nTime + (1:nTime));
else
    Filtered = Y;
end
end
