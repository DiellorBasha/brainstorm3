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


%% ===== GUI: open the surface figure and wire up mode browsing =====
function hFig = ViewFigure(SurfaceFile, iMode)
    hFig = [];
    if (nargin < 2) || isempty(iMode)
        iMode = 1;
    end
    % Load eigenmodes embedded in the surface file
    [Eig, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if ~isComputed || isempty(Eig) || ~isfield(Eig, 'Vectors') || isempty(Eig.Vectors) || ~isfield(Eig, 'Values')
        bst_error(['No eigenmodes found on this surface.' 10 ...
                   'Right-click the cortex and run "Compute eigenmodes" first.'], ...
                   'View eigenmodes', 0);
        return;
    end
    nModes = size(Eig.Vectors, 2);
    iMode  = StepMode(iMode, 0, nModes);   % clamp into range

    % Open a transient surface figure (no database node)
    [hFig, iDS, iFig] = view_surface(SurfaceFile, 0, [], 'NewFigure', 0); %#ok<ASGLU>
    if isempty(hFig)
        bst_error('Could not open the surface figure.', 'View eigenmodes', 0);
        return;
    end
    % Register the source colormap so the overlay is colormapped (+ colorbar)
    bst_colormaps('AddColormapToFigure', hFig, 'source');
    % Figure name + default view
    set(hFig, 'Name', ['Eigenmodes: ' SurfaceFile]);
    figure_3d('SetStandardView', hFig, 'left');
    % Legend (bottom-left), mirrors view_leadfield_sensitivity
    hLabel = uicontrol('Style', 'text', 'String', '...', 'Units', 'Pixels', ...
        'Position', [6 0 520 20], 'HorizontalAlignment', 'left', ...
        'FontUnits', 'points', 'FontSize', bst_get('FigFont'), ...
        'ForegroundColor', [.9 .9 .9], 'BackgroundColor', [0 0 0], 'Parent', hFig);
    % Hook the keyboard callback (preserve the original for unhandled keys)
    KeyPressFcn_bak = get(hFig, 'KeyPressFcn');
    set(hFig, 'KeyPressFcn', @KeyPress_Callback);
    % Initial render
    UpdateMode();

    % ===== NESTED: render the current mode onto the surface =====
    function UpdateMode()
        d = GetModeDisplay(Eig, iMode);
        TessInfo = getappdata(hFig, 'Surface');
        TessInfo(1).ColormapType        = 'source';   % render eigenmode with the source palette (not anatomy)
        TessInfo(1).DataSource.Type     = 'Source';
        TessInfo(1).DataSource.FileName = SurfaceFile;   % real DB file -> satisfies CLim guard
        TessInfo(1).Data                = d.Data;
        TessInfo(1).DataMinMax          = d.CLim;
        TessInfo(1).DataLimitValue      = d.CLim;
        setappdata(hFig, 'Surface', TessInfo);
        panel_surface('UpdateSurfaceColormap', hFig);
        set(hLabel, 'String', d.Label);
    end

    % ===== NESTED: keyboard navigation =====
    function KeyPress_Callback(h, keyEvent)
        switch (keyEvent.Key)
            case 'leftarrow'
                iMode = StepMode(iMode, -1, nModes);   UpdateMode();
            case 'rightarrow'
                iMode = StepMode(iMode, +1, nModes);   UpdateMode();
            case 'pageup'
                iMode = StepMode(iMode, +10, nModes);  UpdateMode();
            case 'pagedown'
                iMode = StepMode(iMode, -10, nModes);  UpdateMode();
            case 'h'
                java_dialog('msgbox', ...
                    ['Eigenmode viewer shortcuts:' 10 10 ...
                     '   Left / Right arrow  :  previous / next mode' 10 ...
                     '   Page Up / Page Down :  +/- 10 modes' 10 ...
                     '   H                   :  this help'], 'Eigenmode viewer');
            otherwise
                if ~isempty(KeyPressFcn_bak)
                    KeyPressFcn_bak(h, keyEvent);
                end
        end
    end
end
