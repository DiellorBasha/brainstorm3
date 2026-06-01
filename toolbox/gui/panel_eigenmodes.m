function varargout = panel_eigenmodes(varargin)
% PANEL_EIGENMODES: Eigenmode scale lever — a stateful, surface-scoped,
% live-broadcast selection in eigenmode (spatial-frequency) space. Changing the
% selection live-coarsens/smooths a displayed source map via a non-destructive
% display-time band-limited reconstruction.
%
% USAGE:  bstPanel = panel_eigenmodes('CreatePanel')
%         W  = panel_eigenmodes('BuildWeights', shape, kLo, kHi, iCenter, K)
%         panel_eigenmodes('SetBand', kLo, kHi)
%         panel_eigenmodes('SetCurrentMode', k)
%         panel_eigenmodes('SetWindowShape', shape)
%         panel_eigenmodes('SetActive', isActive)
%         W  = panel_eigenmodes('GetWeights')
%         uF = panel_eigenmodes('ApplyToColumn', SurfaceFile, u)
%         panel_eigenmodes('UpdatePanel', hFig)
%
% SEE ALSO: bst_eigenmodes_filter, in_tess_eigenmodes, panel_freq, panel_surface

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

eval(macro_method);
end


%% ===== PURE: build the [1 x K] weight vector from a window shape =====
function W = BuildWeights(shape, kLo, kHi, iCenter, K) %#ok<DEFNU>
    % Clamp all indices into [1,K]
    kLo     = min(max(round(kLo),     1), K);
    kHi     = min(max(round(kHi),     1), K);
    iCenter = min(max(round(iCenter), 1), K);
    if (kHi < kLo)
        [kLo, kHi] = deal(kHi, kLo);
    end
    W = zeros(1, K);
    switch lower(shape)
        case 'single'
            W(iCenter) = 1;
        case 'box'
            W(kLo:kHi) = 1;
        case 'tapered'
            % Tukey window over [kLo,kHi]: cosine shoulders taper to 0 at kLo and kHi;
            % flat interior (=1) over the central portion of the band.
            n = kHi - kLo + 1;
            if (n <= 2)
                W(kLo:kHi) = 1;
            else
                r = 0.5;    % total taper fraction (each shoulder spans r/2 = 25% of the window)
                t = linspace(0, 1, n);
                wt = ones(1, n);
                edge = (r/2);
                iL = t < edge;
                iR = t > (1 - edge);
                wt(iL) = 0.5 * (1 + cos(pi * (2*t(iL)/r - 1)));
                wt(iR) = 0.5 * (1 + cos(pi * (2*t(iR)/r - 2/r + 1)));
                W(kLo:kHi) = wt;
            end
        case 'gain'
            % Gaussian bell centered at iCenter; sigma = half the band width.
            sigma = max((kHi - kLo) / 2, 1);
            k = 1:K;
            W = exp(-0.5 * ((k - iCenter) / sigma).^2);
        otherwise
            error('panel_eigenmodes:BuildWeights: unknown shape "%s".', shape);
    end
end
