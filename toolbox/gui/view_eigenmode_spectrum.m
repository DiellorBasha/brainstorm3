function varargout = view_eigenmode_spectrum(varargin)
% VIEW_EIGENMODE_SPECTRUM: Modal power spectrum of a source map's activations.
%
% USAGE:  hFig = view_eigenmode_spectrum(ResultsFile)
%         pw   = view_eigenmode_spectrum('ComputeModalPower', ThetaCol, Component)
%         ax   = view_eigenmode_spectrum('GetSpectrumAxis', Values, mode)
%         avg  = view_eigenmode_spectrum('GetWindowAverage', Theta, iWin)
%         Th   = view_eigenmode_spectrum('CollapseProject', Eig, ImageGridAmp, nComp, M)
%
% Projects the realized vertex source map onto the surface LBO eigenmodes and
% displays power per mode (Left/Right hemisphere curves) vs eigenvalue (or
% spatial wavelength). The figure is registered in the source map's dataset and
% is driven by Brainstorm's global time cursor (see bst_figures FireCurrentTimeChanged).

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

methodNames = {'ComputeModalPower', 'GetSpectrumAxis', 'GetWindowAverage', ...
               'CollapseProject', 'CreateFigure', 'UpdateFigurePlot', ...
               'CurrentTimeChangedCallback'};
if (nargin >= 1) && ischar(varargin{1}) && ismember(varargin{1}, methodNames)
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
[varargout{1:nargout}] = ViewFigure(varargin{:});
end


%% ===== PURE: |theta|^2 split by hemisphere component =====
function pw = ComputeModalPower(ThetaCol, Component)
    p = abs(ThetaCol(:)) .^ 2;
    Component = Component(:);
    pw.left  = p(Component == 1);
    pw.right = p(Component == 2);
end


%% ===== PURE: spectrum x-axis (eigenvalue or spatial wavelength) =====
function ax = GetSpectrumAxis(Values, mode)
    Values = Values(:);
    switch lower(mode)
        case 'eigenvalue'
            ax.x     = Values;
            ax.label = 'Eigenvalue \lambda';
        case 'wavelength'
            x = nan(size(Values));
            pos = (Values > 0);
            x(pos) = 2 * pi ./ sqrt(Values(pos));
            ax.x     = x;
            ax.label = 'Spatial wavelength \approx 2\pi/\surd\lambda';
        otherwise
            error('Unknown spectrum axis mode: %s', mode);
    end
end


%% ===== PURE: mean modal power over a sample window =====
function avg = GetWindowAverage(Theta, iWin)
    if isempty(iWin)
        iWin = 1:size(Theta, 2);
    end
    avg = mean(abs(Theta(:, iWin)) .^ 2, 2);
end
