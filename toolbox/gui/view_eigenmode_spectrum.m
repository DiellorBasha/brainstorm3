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

methodNames = {'ComputeModalPower', 'GetSpectrumAxis', 'GetWindowAverage', 'BuildModeSpectrum', ...
               'CreateFigure', 'UpdateFigurePlot', 'CurrentTimeChangedCallback', ...
               'SetAxisMode', 'SetPowerMode'};
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


%% ===== PURE: build the [2 x K] L/R power matrix + shared x-vector for an axis mode =====
% Maps the per-mode coefficients to a 2-signal (Left/Right) x K-rank layout that the
% figure_timeseries engine can plot. The x-vector is shared by both hemispheres:
%   'index'      -> within-hemisphere rank k (exact)
%   'eigenvalue' -> per-rank eigenvalue (averaged across hemispheres present)
%   'wavelength' -> 2*pi/sqrt(eigenvalue)
% Missing (hemisphere lacks a given rank) entries are NaN so the line skips them.
% powerMode selects the amplitude shown: 'power' (|theta|^2), 'magnitude' (|theta|),
% or 'log' (10*log10 power, dB). Default 'power'.
function spec = BuildModeSpectrum(Info, ThetaCol, xmode, powerMode)
    if (nargin < 4) || isempty(powerMode)
        powerMode = 'power';
    end
    Comp = Info.Component(:);
    Rank = Info.CompRank(:);
    Val  = Info.Values(:);
    a    = abs(ThetaCol(:));
    switch lower(powerMode)
        case 'power'
            p = a .^ 2;            spec.ylabel = 'Modal power |\theta_k|^2';
        case 'magnitude'
            p = a;                 spec.ylabel = 'Modal magnitude |\theta_k|';
        case 'log'
            pw = a .^ 2; pw(pw <= 0) = eps;
            p = 10 * log10(pw);    spec.ylabel = 'Modal power (dB)';
        otherwise
            error('Unknown power mode: %s', powerMode);
    end
    K    = max(Rank);
    F    = nan(2, K);
    lam  = nan(1, K);
    for k = 1:K
        iL = (Comp == 1) & (Rank == k);
        iR = (Comp == 2) & (Rank == k);
        if any(iL), F(1, k) = p(find(iL, 1)); end
        if any(iR), F(2, k) = p(find(iR, 1)); end
        lk = Val(Rank == k);
        if ~isempty(lk), lam(k) = mean(lk); end
    end
    switch lower(xmode)
        case 'index'
            spec.x     = 1:K;
            spec.label = 'Mode index k';
        case 'eigenvalue'
            spec.x     = lam;
            spec.label = 'Eigenvalue \lambda';
        case 'wavelength'
            w = nan(1, K);
            pos = (lam > 0);
            w(pos) = 2 * pi ./ sqrt(lam(pos));
            spec.x     = w;
            spec.label = 'Spatial wavelength \approx 2\pi/\surd\lambda';
        otherwise
            error('Unknown x-axis mode: %s', xmode);
    end
    % figure_timeseries requires a finite, ASCENDING x-vector (it sets XLim from the
    % ends). Drop unusable modes (e.g. wavelength of lambda<=0 -> NaN) and sort
    % ascending, reordering the power columns to match (wavelength is descending in rank).
    xv    = spec.x(:)';
    valid = ~isnan(xv) & ~isinf(xv);
    xv    = xv(valid);
    F     = F(:, valid);
    [xv, iSort] = sort(xv, 'ascend');
    spec.x = xv;
    spec.F = F(:, iSort);
    spec.rowLabels = {'Left', 'Right'};
end


%% ===== GUI: open the spectrum figure registered in the source map's dataset =====
function hFig = ViewFigure(ResultsFile)
    global GlobalData;
    hFig = [];
    bst_progress('start', 'Eigenspectrum', 'Loading source activations...');
    try
        % --- Surface + eigenmodes first (metadata only; does NOT load recordings) ---
        ResHdr = in_bst_results(ResultsFile, 0, 'SurfaceFile');
        if ~isfield(ResHdr, 'SurfaceFile') || isempty(ResHdr.SurfaceFile)
            bst_progress('stop');
            bst_error('This source file has no associated surface.', 'Eigenspectrum', 0);
            return;
        end
        SurfaceFile = ResHdr.SurfaceFile;
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
        % --- Run the standard "cortical activations" machinery ---
        % Opening the cortex figure sets up the dataset (raw paging/navigation, kernel,
        % time-stepping), so the eigenspectrum reads the SAME per-page source map the cortex
        % colormap uses, on the fly -- no full-recording load, correct page advance.
        hFig3d = view_surface_data(SurfaceFile, ResultsFile);
        if isempty(hFig3d)
            bst_progress('stop');
            return;
        end
        % Dataset/result indices the cortex figure just set up (returns existing, no reload).
        [iDS, iResult] = bst_memory('LoadResultsFile', ResultsFile);
        if isempty(iDS)
            bst_progress('stop');
            return;
        end
        % Reject head models we cannot project (e.g. mixed, nComponents==0).
        nComp = GlobalData.DataSet(iDS).Results(iResult).nComponents;
        if ~ismember(nComp, [1 2 3])
            bst_progress('stop');
            bst_error(['Eigenspectrum supports surface source models with 1 or 3' 10 ...
                       'orientations per vertex (not mixed head models).'], 'Eigenspectrum', 0);
            return;
        end
        % Mass matrix: reuse the one stored with the eigenmodes (computed once at
        % "Compute eigenmodes"); fall back to building it for older files without it.
        nVert = size(Eig.Vectors, 1);
        if isfield(Eig, 'MassMatrix') && ~isempty(Eig.MassMatrix) && isequal(size(Eig.MassMatrix), [nVert, nVert])
            M = Eig.MassMatrix;
        else
            sSurf = in_tess_bst(SurfaceFile);
            [~, M] = tess_laplacian(sSurf.Vertices, sSurf.Faces, 'MassType', Eig.MassType);
        end
        % Probe the current time: confirm the source grid matches the eigenmode surface.
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
        setappdata(hFig, 'PowerMode',   'power');
        setappdata(hFig, 'isStatic',    numel(TimeVector) <= 2);
        % figure_timeseries display state: its engine provides butterfly/column,
        % log/linear scales, gain and zoom buttons; we drive the plotting with our
        % per-time modal data over a mode x-axis.
        TsInfo = db_template('TsInfo');
        TsInfo.DisplayMode   = 'butterfly';
        TsInfo.LinesLabels   = {'Left', 'Right'};
        TsInfo.LinesColor    = {[0.85 0.2 0.2; 0.2 0.3 0.85]};   % 1 x nAxes cell of [nLines x 3]
        TsInfo.ShowEvents    = 0;
        TsInfo.AutoScaleY    = 1;
        TsInfo.DefaultFactor = 1;
        TsInfo.XScale        = 'linear';
        TsInfo.YScale        = 'linear';
        setappdata(hFig, 'TsInfo', TsInfo);
        % Let bst_figures('ReloadFigures') redraw us (e.g. butterfly<->column toggle).
        setappdata(hFig, 'ReloadCall', {'view_eigenmode_spectrum', 'UpdateFigurePlot', hFig, 0});
        % First plot + show + select.
        UpdateFigurePlot(hFig, 0);
        set(hFig, 'Visible', 'on');
        bst_figures('SetCurrentFigure', hFig, '2D');
    catch ME
        bst_progress('stop');
        bst_error(['Could not open the eigenspectrum:' 10 ME.message], 'Eigenspectrum', 0);
        return;
    end
    bst_progress('stop');
end


%% ===== GUI: build the figure (called by bst_figures CreateFigure) =====
% Build it THROUGH figure_timeseries so it inherits all the appdata
% (isPlotEditToolbar, GraphSelection, Colormap, ...) and callbacks (resize,
% mouse, scroll, display-config menu) the reused engine relies on. Then we just
% rename it and point keys at our handler (x-axis modes), delegating the rest.
function hFig = CreateFigure(FigureId)
    hFig = figure_timeseries('CreateFigure', FigureId);
    set(hFig, 'Name', 'Eigenspectrum', 'KeyPressFcn', @FigureKeyPressedCallback);
    % Wrap the engine resize: it lays the axes out full-height (its ShowEvents=0
    % path leaves no room for our x-tick labels/title); we then carve margins.
    set(hFig, bst_get('ResizeFunction'), @(h,ev)ResizeCallback(h, ev));
end


%% ===== GUI: resize = engine layout + margins for the mode x-axis label/title =====
function ResizeCallback(hFig, ev)
    figure_timeseries('ResizeCallback', hFig, ev);   % buttons + axes (full height)
    AdjustAxesMargins(hFig);
end

function AdjustAxesMargins(hFig)
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'AxesGraph');
    % Only the single-axes (butterfly) layout needs fixing; the multi-axes (split)
    % layout is already stacked with proper margins by the engine.
    if numel(hAxes) ~= 1 || ~ishandle(hAxes)
        return;
    end
    Scaling = bst_get('InterfaceScaling') / 100;
    marginBottom = 45 * Scaling;   % room for x-tick labels + axis label
    marginTop    = 30 * Scaling;   % room for the title
    figPos = get(hFig, 'Position');
    set(hAxes, 'Units', 'pixels');
    p = get(hAxes, 'Position');    % keep engine's left/width (legend + button column)
    p(2) = marginBottom;
    p(4) = max(figPos(4) - marginBottom - marginTop, 1);
    set(hAxes, 'Position', p);
end


%% ===== GUI: redraw the eigenspectrum at the current global time =====
% Mirrors "Display on cortex": project the SINGLE current-time per-vertex scalar
% field onto the eigenmodes, then render through the figure_timeseries engine so
% the standard buttons (butterfly/column, log/linear scales, gain, zoom) all apply.
% isFastUpdate=1 just refreshes the line data (fast time-stepping); =0 full replot.
function UpdateFigurePlot(hFig, isFastUpdate)
    global GlobalData;
    if (nargin < 2) || isempty(isFastUpdate)
        isFastUpdate = 0;
    end
    iDS      = getappdata(hFig, 'iDS');
    iResult  = getappdata(hFig, 'iResult');
    Eig      = getappdata(hFig, 'Eig');
    M        = getappdata(hFig, 'MassMatrix');
    Info     = getappdata(hFig, 'SpecInfo');
    AxisMode = getappdata(hFig, 'AxisMode');
    if isempty(iDS) || isempty(Eig)
        return;
    end
    [~, iFig] = bst_figures('GetFigure', hFig);
    if isempty(iFig)
        return;
    end
    % Current cortical activation (single time index) -> eigenmode coefficients.
    S = bst_memory('GetResultsValues', iDS, iResult, [], 'CurrentTimeIndex');   % [nVert x 1]
    if isempty(S)
        return;
    end
    PowerMode = getappdata(hFig, 'PowerMode');
    if isempty(PowerMode), PowerMode = 'power'; end
    ThetaCol = bst_eigenmodes_project(Eig, S, M);              % [nModes x 1]
    spec = BuildModeSpectrum(Info, ThetaCol, AxisMode, PowerMode);   % .F [2 x K], .x [1 x K], .label, .ylabel
    % The shared x-vector is the mode axis (NOT seconds); store it for figure_timeseries callbacks.
    setappdata(hFig, 'TimeVector', spec.x);
    % Layout from the menu's display mode:
    %   'butterfly' -> both hemispheres overlaid in one axes
    %   'column'    -> split into two stacked axes, each autoscaled to fill
    TsInfo  = getappdata(hFig, 'TsInfo');
    isSplit = strcmpi(TsInfo.DisplayMode, 'column');
    if isSplit
        F = {spec.F(1, :), spec.F(2, :)};                 % 2 axes (Left, Right)
        TsInfo.LinesColor  = {[0.85 0.2 0.2], [0.2 0.3 0.85]};
        TsInfo.LinesLabels = {'Left', 'Right'};
        TsInfo.AxesLabels  = {'Left', 'Right'};
    else
        F = {spec.F};                                     % 1 axes, 2 overlaid lines
        TsInfo.LinesColor  = {[0.85 0.2 0.2; 0.2 0.3 0.85]};
        TsInfo.LinesLabels = {'Left', 'Right'};
        TsInfo.AxesLabels  = {};
    end
    % Render each axes as butterfly (autoscale-fill); restore the menu-facing display
    % mode afterwards so the butterfly/column radio stays in sync.
    menuMode = TsInfo.DisplayMode;
    TsInfo.DisplayMode = 'butterfly';
    setappdata(hFig, 'TsInfo', TsInfo);
    figure_timeseries('PlotFigure', iDS, iFig, F, spec.x, isFastUpdate, []);
    TsInfo = getappdata(hFig, 'TsInfo');
    TsInfo.DisplayMode = menuMode;
    setappdata(hFig, 'TsInfo', TsInfo);
    % Adapt the time-series defaults to a modal spectrum: x-axis label + Y label.
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'AxesGraph');
    for k = 1:numel(hAxes)
        xlabel(hAxes(k), spec.label, 'Interpreter', 'tex');
        ylabel(hAxes(k), spec.ylabel, 'Interpreter', 'tex');
    end
    % Hide the vertical "current time" cursor: it is positioned in seconds and is
    % meaningless on a mode/eigenvalue axis.
    % Hide the seconds time-cursor on every axes (Handles is a struct array in split mode).
    Handles = GlobalData.DataSet(iDS).Figure(iFig).Handles;
    if isfield(Handles, 'hCursor')
        hc = [Handles.hCursor];
        hc = hc(ishandle(hc));
        if ~isempty(hc)
            set(hc, 'Visible', 'off');
        end
    end
    % Keep only the Left/Right data lines in the legend (exclude cursor/baseline).
    for h = findobj(hAxes, 'Type', 'line')'
        if ~strcmp(get(h, 'Tag'), 'DataLine')
            set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
        end
    end
    % Title with the current time + key hints.
    if ~isempty(GlobalData.UserTimeWindow.CurrentTime)
        tStr = sprintf('%1.3f s', GlobalData.UserTimeWindow.CurrentTime);
    else
        tStr = 'n/a';
    end
    if ~isempty(hAxes)
        % Title on the topmost axes (split mode has two stacked axes).
        ytop = zeros(1, numel(hAxes));
        for k = 1:numel(hAxes)
            pp = get(hAxes(k), 'Position');
            ytop(k) = pp(2);
        end
        [~, iTop] = max(ytop);
        title(hAxes(iTop), sprintf('Eigenspectrum  |  t = %s  |  x-axis: e=\\lambda  k=index  w=wavelength', tStr), ...
              'Interpreter', 'tex');
    end
    % Single-axes (butterfly): the engine lays it out full-height -> carve margins.
    % Multi-axes (split): the engine already stacks them with proper margins.
    AdjustAxesMargins(hFig);
end


%% ===== GUI: x-axis mode (eigenvalue/index/wavelength) — from menu or e/k/w keys =====
function SetAxisMode(hFig, mode)
    setappdata(hFig, 'AxisMode', mode);
    UpdateFigurePlot(hFig, 0);
end


%% ===== GUI: amplitude mode (power/magnitude/log) — from the Amplitude menu =====
function SetPowerMode(hFig, mode)
    setappdata(hFig, 'PowerMode', mode);
    UpdateFigurePlot(hFig, 0);
end


%% ===== GUI: global time cursor moved (called by bst_figures FireCurrentTimeChanged) =====
function CurrentTimeChangedCallback(hFig)
    if getappdata(hFig, 'isStatic')
        return;
    end
    UpdateFigurePlot(hFig, 1);   % fast update: refresh line data only
end


%% ===== GUI: keys — e/k/w switch x-axis; everything else -> figure_timeseries =====
function FigureKeyPressedCallback(hFig, ev)
    switch (ev.Key)
        case 'e'
            SetAxisMode(hFig, 'eigenvalue');
        case 'k'
            SetAxisMode(hFig, 'index');
        case 'w'
            SetAxisMode(hFig, 'wavelength');
        otherwise
            % Standard shortcuts (arrows = step time, gain, etc.).
            figure_timeseries('FigureKeyPressedCallback', hFig, ev);
    end
end
