function varargout = view_eigenmode_spectrum(varargin)
% VIEW_EIGENMODE_SPECTRUM: Modal power spectrum of a source map's activations.
%
% USAGE:  hFig = view_eigenmode_spectrum(ResultsFile)
%         pw   = view_eigenmode_spectrum('ComputeModalPower', ThetaCol, Component)
%         ax   = view_eigenmode_spectrum('GetSpectrumAxis', Values, mode)
%         avg  = view_eigenmode_spectrum('GetWindowAverage', Theta, iWin)
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
               'CreateFigure', 'UpdateFigurePlot', 'CurrentTimeChangedCallback'};
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


%% ===== GUI: open the spectrum figure registered in the source map's dataset =====
function hFig = ViewFigure(ResultsFile)
    global GlobalData;
    hFig = [];
    bst_progress('start', 'Eigenspectrum', 'Loading source activations...');
    try
        % Load the source result into a dataset (kernel link OK, including raw).
        [iDS, iResult] = bst_memory('LoadResultsFile', ResultsFile);
        if isempty(iDS)
            bst_progress('stop');
            return;
        end
        % Materialize the kernel/grid matrices: LoadResultsFile loads only metadata,
        % so ImagingKernel/ImageGridAmp and nComponents are empty until this runs
        % (mirrors panel_surface's UpdateSurfaceData path for the cortex display).
        if isempty(GlobalData.DataSet(iDS).Results(iResult).ImageGridAmp) && ...
           isempty(GlobalData.DataSet(iDS).Results(iResult).ImagingKernel)
            bst_memory('LoadResultsMatrix', iDS, iResult);
        end
        % Surface associated with the source model.
        SurfaceFile = GlobalData.DataSet(iDS).Results(iResult).SurfaceFile;
        if isempty(SurfaceFile)
            bst_progress('stop');
            bst_error('This source file has no associated surface.', 'Eigenspectrum', 0);
            return;
        end
        % Eigenmodes on that surface.
        [Eig, isComputed] = in_tess_eigenmodes(SurfaceFile);
        if ~isComputed || isempty(Eig) || ~isfield(Eig, 'Vectors') || isempty(Eig.Vectors)
            bst_progress('stop');
            bst_error(['No eigenmodes found on this surface.' 10 ...
                       'Right-click the cortex and run "Compute eigenmodes" first.'], 'Eigenspectrum', 0);
            return;
        end
        % Older eigenmode files may not store MassType; default to barycentric.
        if ~isfield(Eig, 'MassType') || isempty(Eig.MassType)
            Eig.MassType = 'barycentric';
        end
        % Reject head models we cannot project (e.g. mixed, nComponents==0).
        nComp = GlobalData.DataSet(iDS).Results(iResult).nComponents;
        if ~ismember(nComp, [1 2 3])
            bst_progress('stop');
            bst_error(['Eigenspectrum supports surface source models with 1 or 3' 10 ...
                       'orientations per vertex (not mixed head models).'], 'Eigenspectrum', 0);
            return;
        end
        % Mass matrix consistent with the stored mass type.
        sSurf = in_tess_bst(SurfaceFile);
        [~, M] = tess_laplacian(sSurf.Vertices, sSurf.Faces, 'MassType', Eig.MassType);
        nVert = size(Eig.Vectors, 1);
        % Probe the current time: confirms the source grid matches the eigenmode surface
        % AND that the realized (lazy) field is available for this result.
        Sprobe = bst_memory('GetResultsValues', iDS, iResult, [], 'CurrentTimeIndex');
        if isempty(Sprobe) || (size(Sprobe, 1) ~= nVert)
            bst_progress('stop');
            bst_error(['The source grid does not match the eigenmode surface' 10 ...
                       '(expected ' num2str(nVert) ' vertices).'], 'Eigenspectrum', 0);
            return;
        end
        % Static if the result has <= 2 time samples.
        TimeVector = bst_memory('GetTimeVector', iDS, iResult, 'UserTimeWindow');
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
        % Stash what the lazy redraw needs (no precomputed all-time matrix).
        setappdata(hFig, 'iDS',         iDS);
        setappdata(hFig, 'iResult',     iResult);
        setappdata(hFig, 'Eig',         Eig);
        setappdata(hFig, 'MassMatrix',  M);
        setappdata(hFig, 'SpecInfo',    struct('Values', Eig.Values(:), ...
                                               'Component', Eig.Component(:), ...
                                               'CompRank', Eig.CompRank(:), ...
                                               'SurfaceFile', SurfaceFile));
        setappdata(hFig, 'ResultsFile', ResultsFile);
        setappdata(hFig, 'AxisMode',    'eigenvalue');
        setappdata(hFig, 'isStatic',    numel(TimeVector) <= 2);
        setappdata(hFig, 'AvgCache',    []);
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


%% ===== GUI: redraw the spectrum at the current global time (lazy) =====
function UpdateFigurePlot(hFig)
    global GlobalData;
    iDS     = getappdata(hFig, 'iDS');
    iResult = getappdata(hFig, 'iResult');
    Eig     = getappdata(hFig, 'Eig');
    M       = getappdata(hFig, 'MassMatrix');
    Info    = getappdata(hFig, 'SpecInfo');
    AxisMode = getappdata(hFig, 'AxisMode');
    if isempty(iDS) || isempty(Eig)
        return;
    end
    % Realized vertex scalar field at the current time (lazy; valid for raw).
    S = bst_memory('GetResultsValues', iDS, iResult, [], 'CurrentTimeIndex');   % [nVert x 1]
    if isempty(S)
        return;
    end
    ThetaCol = bst_eigenmodes_project(Eig, S, M);          % [nModes x 1]
    pw  = ComputeModalPower(ThetaCol, Info.Component);
    ax  = GetSpectrumAxis(Info.Values, AxisMode);
    avg = GetWindowAverageLazy(hFig, iDS, iResult, Eig, M); % [nModes x 1] or [] if unavailable
    Comp = Info.Component(:);
    xL = ax.x(Comp == 1);  xR = ax.x(Comp == 2);
    % Draw (collect handles so the legend only labels curves actually present).
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'AxesEigenSpectrum');
    cla(hAxes);
    hold(hAxes, 'on');
    hLines = []; labels = {};
    if ~isempty(xL)
        hLines(end+1) = plot(hAxes, xL, pw.left, '-', 'Color', [0.85 0.2 0.2], 'LineWidth', 1.5); labels{end+1} = 'Left (t)';
        if ~isempty(avg)
            hLines(end+1) = plot(hAxes, xL, avg(Comp == 1), '--', 'Color', [0.85 0.2 0.2], 'LineWidth', 0.75); labels{end+1} = 'Left (avg)';
        end
    end
    if ~isempty(xR)
        hLines(end+1) = plot(hAxes, xR, pw.right, '-', 'Color', [0.2 0.3 0.85], 'LineWidth', 1.5); labels{end+1} = 'Right (t)';
        if ~isempty(avg)
            hLines(end+1) = plot(hAxes, xR, avg(Comp == 2), '--', 'Color', [0.2 0.3 0.85], 'LineWidth', 0.75); labels{end+1} = 'Right (avg)';
        end
    end
    hold(hAxes, 'off');
    xlabel(hAxes, ax.label);
    ylabel(hAxes, 'Modal power |\theta_k|^2');
    if ~isempty(hLines)
        legend(hAxes, hLines, labels, 'Location', 'northeast');
    end
    if ~isempty(GlobalData) && ~isempty(GlobalData.UserTimeWindow.CurrentTime)
        tStr = sprintf('%1.3f s', GlobalData.UserTimeWindow.CurrentTime);
    else
        tStr = 'n/a';
    end
    title(hAxes, sprintf('Eigenspectrum  |  t = %s  |  L=%d R=%d modes  |  [e/w axis, arrows step time]', ...
          tStr, numel(xL), numel(xR)), 'Interpreter', 'tex');
end


%% ===== GUI: window-averaged modal power over the loaded page (cached) =====
function avg = GetWindowAverageLazy(hFig, iDS, iResult, Eig, M)
    global GlobalData;
    avg = [];
    tw = [];
    if ~isempty(GlobalData) && ~isempty(GlobalData.UserTimeWindow.Time)
        tw = GlobalData.UserTimeWindow.Time;
    end
    % Reuse the cached average while the loaded time window is unchanged.
    cache = getappdata(hFig, 'AvgCache');
    if ~isempty(cache) && isequal(cache.tw, tw)
        avg = cache.avg;
        return;
    end
    try
        Swin = bst_memory('GetResultsValues', iDS, iResult, [], 'UserTimeWindow');   % [nVert x nWin]
        if ~isempty(Swin)
            Theta = bst_eigenmodes_project(Eig, Swin, M);
            avg = GetWindowAverage(Theta, []);
        end
    catch
        avg = [];   % overlay omitted if the window cannot be read
    end
    setappdata(hFig, 'AvgCache', struct('tw', tw, 'avg', avg));
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
