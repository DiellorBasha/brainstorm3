function varargout = view_eigenmodes(varargin)
% VIEW_EIGENMODES: Interactively browse Laplace-Beltrami eigenmodes on a surface.
%
% USAGE:  hFig  = view_eigenmodes(SurfaceFile)              % open viewer at mode 1
%         hFig  = view_eigenmodes(SurfaceFile, iMode)       % open viewer at mode iMode
%         disp  = view_eigenmodes('GetModeDisplay', Eig, iMode)     % pure helper (testable)
%         iNext = view_eigenmodes('StepMode', iMode, delta, nModes) % pure helper (testable)
%
% Modeled on view_leadfield_sensitivity.m: the Left/Right arrows step the
% displayed mode, the cortex colormap updates live, and a legend reports the
% mode index, eigenvalue, and approximate spatial wavelength. The figure is
% transient (no database node is created).
%
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

% ===== METHOD DISPATCH (headless-testable pure helpers) =====
if (nargin >= 1) && ischar(varargin{1}) && any(strcmp(varargin{1}, {'GetModeDisplay', 'StepMode'}))
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
% ===== MAIN ENTRY: open the viewer =====
[varargout{1:nargout}] = ViewFigure(varargin{:});
end


%% ===== PURE: per-mode display package =====
function d = GetModeDisplay(Eig, iMode)
    nModes = size(Eig.Vectors, 2);
    iMode  = min(max(round(iMode), 1), max(nModes, 1));   % clamp to [1, nModes]
    data   = double(Eig.Vectors(:, iMode));
    m      = max(abs(data));
    if ~(m > 0)        % degenerate (all-zero) mode: avoid a zero-width colormap range
        m = 1e-30;
    end
    lambda = Eig.Values(iMode);
    if lambda > 0
        wavelength = 2*pi / sqrt(lambda);
        wlStr      = sprintf('%.3g', wavelength);
    else
        wavelength = NaN;
        wlStr      = 'n/a';
    end
    d = struct();
    d.iMode      = iMode;
    d.nModes     = nModes;
    d.Data       = data;
    d.CLim       = [-m, m];          % symmetric: eigenmodes are signed (+/- lobes)
    d.Lambda     = lambda;
    d.Wavelength = wavelength;
    d.Label      = sprintf('Mode %d / %d     lambda = %.4g     wavelength ~ %s', ...
                           iMode, nModes, lambda, wlStr);
end


%% ===== PURE: clamped mode stepping =====
function iNext = StepMode(iMode, delta, nModes)
    iNext = min(max(round(iMode) + round(delta), 1), max(nModes, 1));
end
