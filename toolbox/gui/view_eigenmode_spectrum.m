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


%% ===== CORE: realized vertex field -> eigenmode coefficients =====
% ImageGridAmp: [nComp*nVert x nTime] source matrix (kernel already applied).
% nComp: 1 (constrained, signed) or 2/3 (unconstrained -> rms magnitude per vertex).
% M: [nVert x nVert] sparse mass matrix for the eigenbasis Eig.Vectors.
function Theta = CollapseProject(Eig, ImageGridAmp, nComp, M)
    if (nComp == 3) || (nComp == 2)
        S = bst_source_orient([], nComp, [], ImageGridAmp, 'rms');   % [nVert x nTime]
    else
        S = ImageGridAmp;                                            % [nVert x nTime]
    end
    Theta = bst_eigenmodes_project(Eig, S, M);                       % [nModes x nTime]
end


%% ===== IO: load realized source map + eigenmodes -> coefficients =====
% Returns Theta [nModes x nTime] and Info (Values/Component/CompRank/Time/SurfaceFile).
% Returns Theta=[] (after a bst_error) when the surface has no eigenmodes.
function [Theta, Info] = GetActivationCoeffs(ResultsFile)
    Theta = [];
    Info  = [];
    % Realized vertex source map (kernel applied on the fly).
    ResultsMat = in_bst_results(ResultsFile, 1);
    if ~isfield(ResultsMat, 'SurfaceFile') || isempty(ResultsMat.SurfaceFile)
        bst_error('This source file has no associated surface.', 'Eigenspectrum', 0);
        return;
    end
    % Eigenmodes on that surface.
    [Eig, isComputed] = in_tess_eigenmodes(ResultsMat.SurfaceFile);
    if ~isComputed || isempty(Eig) || ~isfield(Eig, 'Vectors') || isempty(Eig.Vectors)
        bst_error(['No eigenmodes found on this surface.' 10 ...
                   'Right-click the cortex and run "Compute eigenmodes" first.'], 'Eigenspectrum', 0);
        return;
    end
    % Mass matrix consistent with the stored mass type.
    % Older eigenmode files may not store MassType; default to barycentric.
    if ~isfield(Eig, 'MassType') || isempty(Eig.MassType)
        Eig.MassType = 'barycentric';
    end
    sSurf = in_tess_bst(ResultsMat.SurfaceFile);
    [~, M] = tess_laplacian(sSurf.Vertices, sSurf.Faces, 'MassType', Eig.MassType);
    % Project. Only surface models with 1 (constrained) or 2/3 (unconstrained)
    % orientations per vertex are supported; reject anything else (e.g. mixed
    % head models, nComponents==0) before the projection would error.
    nComp = ResultsMat.nComponents;
    nVert = size(Eig.Vectors, 1);
    if ~ismember(nComp, [1 2 3]) || (size(ResultsMat.ImageGridAmp, 1) ~= nComp * nVert)
        bst_error(['Eigenspectrum supports surface source models with 1 or 3 ' 10 ...
                   'orientations per vertex. This result is not compatible ' 10 ...
                   '(for example, a mixed head model).'], 'Eigenspectrum', 0);
        return;
    end
    Theta = CollapseProject(Eig, ResultsMat.ImageGridAmp, nComp, M);
    Info = struct('Values',      Eig.Values(:), ...
                  'Component',   Eig.Component(:), ...
                  'CompRank',    Eig.CompRank(:), ...
                  'Time',        ResultsMat.Time, ...
                  'SurfaceFile', ResultsMat.SurfaceFile);
end


%% ===== GUI: open the spectrum figure registered in the source map's dataset =====
function hFig = ViewFigure(ResultsFile)
    hFig = [];
    bst_progress('start', 'Eigenspectrum', 'Loading source activations...');
    try
        % Coefficients first (fails fast with a friendly error if no eigenmodes).
        [Theta, Info] = GetActivationCoeffs(ResultsFile);
        if isempty(Theta)
            bst_progress('stop');
            return;
        end
        % Load the results into a dataset so the global time cursor spans this file.
        [iDS, ~] = bst_memory('LoadResultsFile', ResultsFile);
        if isempty(iDS)
            bst_progress('stop');
            return;
        end
        % Create a managed figure of our type, always a fresh one.
        FigureId = db_template('FigureId');
        FigureId.Type     = 'EigenSpectrum';
        FigureId.SubType  = '';
        FigureId.Modality = '';
        [hFig, ~] = bst_figures('CreateFigure', iDS, FigureId, 'AlwaysCreate', ResultsFile);
        if isempty(hFig)
            bst_progress('stop');
            return;
        end
        % Stash everything the redraw needs.
        setappdata(hFig, 'Theta',       Theta);
        setappdata(hFig, 'SpecInfo',    Info);
        setappdata(hFig, 'ResultsFile', ResultsFile);
        setappdata(hFig, 'AxisMode',    'eigenvalue');
        setappdata(hFig, 'isStatic',    size(Theta, 2) <= 2);
        % First plot + show + select.
        UpdateFigurePlot(hFig);
        set(hFig, 'Visible', 'on');
        bst_figures('SetCurrentFigure', hFig, '2D');
    catch ME
        bst_progress('stop');
        bst_error(['Could not open the eigenspectrum:' 10 ME.message], 'Eigenspectrum', 0);
        return;
    end
    bst_progress('stop');
end


%% ===== GUI: build the bare figure (called by bst_figures CreateFigure) =====
function hFig = CreateFigure(FigureId)
    hFig = figure( ...
        'Visible',       'off', ...
        'NumberTitle',   'off', ...
        'IntegerHandle', 'off', ...
        'MenuBar',       'none', ...
        'Toolbar',       'figure', ...
        'DockControls',  'on', ...
        'Units',         'pixels', ...
        'Color',         [.9 .9 .9], ...
        'Pointer',       'arrow', ...
        'Tag',           'EigenSpectrum', ...
        'Name',          'Eigenspectrum', ...
        'CloseRequestFcn', @(h,ev)bst_figures('DeleteFigure', h, ev), ...
        'KeyPressFcn',     @FigureKeyPressedCallback);
    setappdata(hFig, 'FigureId', FigureId);
    setappdata(hFig, 'isStatic', 0);
    % Single axes for the spectrum curves.
    axes('Parent', hFig, 'Tag', 'AxesEigenSpectrum', 'Units', 'normalized', ...
         'Position', [0.12 0.13 0.82 0.78]);
end


%% ===== GUI: redraw the spectrum at the current global time =====
function UpdateFigurePlot(hFig)
    global GlobalData;
    Theta    = getappdata(hFig, 'Theta');
    Info     = getappdata(hFig, 'SpecInfo');
    AxisMode = getappdata(hFig, 'AxisMode');
    if isempty(Theta)
        return;
    end
    % Current time -> column index.
    Time = Info.Time;
    if isempty(GlobalData) || isempty(GlobalData.UserTimeWindow.CurrentTime) || (numel(Time) < 2)
        iTime = 1;
    else
        iTime = bst_closest(GlobalData.UserTimeWindow.CurrentTime, Time);
    end
    iTime = max(1, min(iTime, size(Theta, 2)));
    % Power split + axis + full-window average overlay.
    pw  = ComputeModalPower(Theta(:, iTime), Info.Component);
    ax  = GetSpectrumAxis(Info.Values, AxisMode);
    avg = GetWindowAverage(Theta, []);
    Comp = Info.Component(:);
    xL = ax.x(Comp == 1);  xR = ax.x(Comp == 2);
    aL = avg(Comp == 1);   aR = avg(Comp == 2);
    % Draw (collect handles so the legend only labels curves actually present).
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'AxesEigenSpectrum');
    cla(hAxes);
    hold(hAxes, 'on');
    hLines = []; labels = {};
    if ~isempty(xL)
        hLines(end+1) = plot(hAxes, xL, pw.left, '-',  'Color', [0.85 0.2 0.2], 'LineWidth', 1.5); labels{end+1} = 'Left (t)';
        hLines(end+1) = plot(hAxes, xL, aL,      '--', 'Color', [0.85 0.2 0.2], 'LineWidth', 0.75); labels{end+1} = 'Left (avg)';
    end
    if ~isempty(xR)
        hLines(end+1) = plot(hAxes, xR, pw.right, '-',  'Color', [0.2 0.3 0.85], 'LineWidth', 1.5); labels{end+1} = 'Right (t)';
        hLines(end+1) = plot(hAxes, xR, aR,       '--', 'Color', [0.2 0.3 0.85], 'LineWidth', 0.75); labels{end+1} = 'Right (avg)';
    end
    hold(hAxes, 'off');
    xlabel(hAxes, ax.label);
    ylabel(hAxes, 'Modal power |\theta_k|^2');
    if ~isempty(hLines)
        legend(hAxes, hLines, labels, 'Location', 'northeast');
    end
    if (numel(Time) >= 1)
        tStr = sprintf('%1.3f s', Time(iTime));
    else
        tStr = 'n/a';
    end
    title(hAxes, sprintf('Eigenspectrum  |  t = %s  |  L=%d R=%d modes  |  [e/w axis, arrows step time]', ...
          tStr, numel(xL), numel(xR)), 'Interpreter', 'tex');
end


%% ===== GUI: global time cursor moved (called by bst_figures FireCurrentTimeChanged) =====
function CurrentTimeChangedCallback(hFig)
    if getappdata(hFig, 'isStatic')
        return;
    end
    UpdateFigurePlot(hFig);
end


%% ===== GUI: keys — 'e'/'w' toggle axis; arrows step the global time cursor =====
function FigureKeyPressedCallback(hFig, ev)
    global GlobalData;
    switch (ev.Key)
        case 'e'
            setappdata(hFig, 'AxisMode', 'eigenvalue');
            UpdateFigurePlot(hFig);
        case 'w'
            setappdata(hFig, 'AxisMode', 'wavelength');
            UpdateFigurePlot(hFig);
        case {'leftarrow', 'rightarrow', 'uparrow', 'downarrow', 'pageup', 'pagedown'}
            if isempty(GlobalData) || isempty(GlobalData.UserTimeWindow.SamplingRate)
                return;
            end
            sr = GlobalData.UserTimeWindow.SamplingRate;
            t  = GlobalData.UserTimeWindow.CurrentTime;
            switch (ev.Key)
                case {'leftarrow', 'downarrow'},  t = t - sr;
                case {'rightarrow', 'uparrow'},   t = t + sr;
                case 'pageup',                    t = t + 10 * sr;
                case 'pagedown',                  t = t - 10 * sr;
            end
            panel_time('SetCurrentTime', t);   % moves the one global clock -> redraws all figures
    end
end
