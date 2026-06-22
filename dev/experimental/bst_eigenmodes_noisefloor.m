function Out = bst_eigenmodes_noisefloor(Pdata, Nnoise, varargin)
% BST_EIGENMODES_NOISEFLOOR: Joint (lambda,omega) SNR + spectral-subtraction denoising.
%
% USAGE:  Out = bst_eigenmodes_noisefloor(Pdata, Nnoise)
%         Out = bst_eigenmodes_noisefloor(Pdata, Nnoise, 'Alpha',1, 'Floor',0, 'SnrThresh',1)
%
% DESCRIPTION:
%     Combines a data power spectrum and an empty-room (noise) power spectrum,
%     both in the eigenmode x frequency plane, into denoising products. Works on
%     POWER (averaged PSD), never on complex coefficients: the data and noise
%     recordings are different noise realizations, so complex subtraction would
%     add variance, whereas E[|data|^2] - E[|noise|^2] = |signal|^2.
%
%     Products:
%       SNR(k,f)      = Pdata / Nnoise
%       CleanPSD(k,f) = max(Pdata - Alpha*Nnoise, Floor*Nnoise)   (spectral subtraction)
%       Gain(k,f)     = max(max(Pdata - Alpha*Nnoise, 0)/Pdata, GainFloor)  (Wiener gain in [0,1])
%       Kstar(f)      = largest mode index k with SNR(k,f) >= SnrThresh (0 if none)
%
% INPUTS:
%     Pdata  : [K x nFreq] data PSD (power/Hz) per eigenmode per frequency.
%     Nnoise : [K x nFreq] empty-room PSD, same size and frequency grid.
%
% OPTIONS (name-value):
%     'Alpha'     : over-subtraction factor (>=0; values >1 over-subtract), default 1.
%     'Floor'     : spectral floor as a fraction of Nnoise, default 0.
%     'GainFloor' : lower bound for the Wiener gain (in [0,1]), default 0.
%     'SnrThresh' : linear SNR threshold for the reliable-mode cutoff, default 1.
%
% OUTPUTS:
%     Out.SNR, Out.CleanPSD, Out.Gain : [K x nFreq]
%     Out.Kstar                       : [1 x nFreq]
%
% SEE ALSO: bst_eigenmodes_transform, process_eigenmodes_denoise

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
Alpha = 1; Floor = 0; SnrThresh = 1; GainFloor = 0;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'alpha',     Alpha     = varargin{i+1};
        case 'floor',     Floor     = varargin{i+1};
        case 'snrthresh', SnrThresh = varargin{i+1};
        case 'gainfloor', GainFloor = varargin{i+1};
    end
end
if GainFloor < 0 || GainFloor > 1
    error('bst_eigenmodes_noisefloor: GainFloor must be in [0,1].');
end
if Alpha < 0
    error('bst_eigenmodes_noisefloor: Alpha must be >= 0.');
end

Pdata  = double(Pdata);
Nnoise = double(Nnoise);
if ~isequal(size(Pdata), size(Nnoise))
    error('bst_eigenmodes_noisefloor: Pdata and Nnoise must have the same size.');
end

%% ===== COMBINE =====
Nsafe = max(Nnoise, eps);
Out = struct();
Out.SNR      = Pdata ./ Nsafe;
Out.CleanPSD = max(Pdata - Alpha .* Nnoise, Floor .* Nnoise);
Out.Gain     = max( max(Pdata - Alpha .* Nnoise, 0) ./ max(Pdata, eps), GainFloor );

%% ===== RELIABLE-MODE CUTOFF =====
[~, nFreq] = size(Pdata);
Out.Kstar = zeros(1, nFreq);
for f = 1:nFreq
    idx = find(Out.SNR(:, f) >= SnrThresh, 1, 'last');
    if ~isempty(idx)
        Out.Kstar(f) = idx;
    end
end
end
