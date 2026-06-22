function varargout = view_eigenmode_spectrum(varargin)
% VIEW_EIGENMODE_SPECTRUM: Modal spectrum of a Dirac source map at the current time.
%
% USAGE:  hFig = view_eigenmode_spectrum(ResultsFile)
%         spec = view_eigenmode_spectrum('BuildModeSpectrum', Info, ThetaCol, xmode, powerMode)
%
% A modular eigenspectrum viewer for the canonical Dirac eigenbasis. It shows the
% SAME mode-coefficient object as view_eigen_timeseries -- c(t) = ImagingKernelMode
% * M(GoodChannel,t) -- but as a SINGLE time point: modal power |c_k|^2 per Dirac
% eigenvalue lambda, split into Left/Right hemisphere curves, advancing with the
% global time cursor. It opens NO 3D figure (it is not conflated with the cortex
% map; the cursor is shared globally).
%
% SOURCE TOGGLE (key 'a' / 'd'):
%   'mode'   (amplitude) -- c = ImagingKernelMode * M : the kernel-native mode
%            coefficients, identical to view_eigen_timeseries. Purely spectral.
%   'vertex' (dSPM/etc)  -- J = ImagingKernel * M (the per-vertex reconstruction
%            for the result's measure, e.g. dSPM), then projected back onto the
%            Dirac eigenbasis ( c_k = <J, phi_k>_B ). Because dSPM normalizes per
%            vertex, this spectrum differs from the amplitude one.
%
% X-AXIS (key 'e' / 'k' / 'w'): eigenvalue lambda (default) / mode index / spatial
% wavelength 2*pi/sqrt(lambda) -- all derived from the canonical Dirac lambda.
%
% SEE ALSO: view_eigen_timeseries, bst_inverse_dirac, bst_dirac, tess_eigen

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

methodNames = {'BuildModeSpectrum', 'GetSpectrumAxis', ...
               'CreateFigure', 'UpdateFigurePlot', 'CurrentTimeChangedCallback', ...
               'SetAxisMode', 'SetSpectrumSource'};
if (nargin >= 1) && ischar(varargin{1}) && ismember(varargin{1}, methodNames)
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
[varargout{1:nargout}] = ViewFigure(varargin{:});
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


%% ===== PURE: build the [2 x K] L/R power matrix + shared x-vector for an axis mode =====
% Maps per-mode coefficients to a 2-signal (Left/Right) x K-rank layout for the
% figure_timeseries engine. The x-vector is shared by both hemispheres:
%   'index'      -> within-hemisphere rank k (exact)
%   'eigenvalue' -> per-rank eigenvalue (averaged across hemispheres present)
%   'wavelength' -> 2*pi/sqrt(eigenvalue)
% Missing entries are NaN so the line skips them. powerMode: 'power' (|c|^2),
% 'magnitude' (|c|), or 'log' (10*log10 power, dB). Default 'power'.
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
            p = a .^ 2;            spec.ylabel = 'Modal power |c_k|^2';
        case 'magnitude'
            p = a;                 spec.ylabel = 'Modal magnitude |c_k|';
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
    % figure_timeseries requires a finite, ASCENDING x-vector. Drop unusable modes
    % (e.g. wavelength of lambda<=0 -> NaN) and sort ascending, reordering columns.
    xv    = spec.x(:)';
    valid = ~isnan(xv) & ~isinf(xv);
    xv    = xv(valid);
    F     = F(:, valid);
    [xv, iSort] = sort(xv, 'ascend');
    spec.x = xv;
    spec.F = F(:, iSort);
    spec.rowLabels = {'Left', 'Right'};
end


%% ===== GUI: open the spectrum figure (no 3D map) =====
function hFig = ViewFigure(ResultsFile)
    hFig = [];
    bst_progress('start', 'Eigenspectrum', 'Computing eigenmode coefficients...');
    try
        % --- resolve the result + its recording (same path as view_eigen_timeseries) ---
        [~, DataFile] = file_resolve_link(ResultsFile);
        if isempty(DataFile)
            bst_progress('stop');
            bst_error(['No recording is associated with this result.' 10 ...
                'Open the Dirac source link on a data block (a recordings node).'], 'Eigenspectrum', 0);
            return;
        end
        R = in_bst_results(ResultsFile, 0, 'ImagingKernelMode', 'Eigenvalues', 'ModeHemisphere', ...
                           'GoodChannel', 'ImagingKernel', 'DiracEigenFile', 'SurfaceFile', 'Function');
        if ~isfield(R, 'ImagingKernelMode') || isempty(R.ImagingKernelMode)
            bst_progress('stop');
            bst_error(['This source file has no persisted Dirac eigenmode kernel.' 10 ...
                'Recompute with "Compute sources: Dirac eigenmodes".'], 'Eigenspectrum', 0);
            return;
        end
        % --- mode-coefficient time series c(t) = Kmode * M(GoodChannel,:) (== view_eigen_timeseries) ---
        DataMat = in_bst_data(DataFile, 'F', 'Time');
        gc = R.GoodChannel; if isempty(gc), gc = 1:size(DataMat.F,1); end
        M    = double(DataMat.F(gc, :));                 % [nCh x nTime]
        Camp = double(R.ImagingKernelMode) * M;          % [nMode x nTime] amplitude mode coeffs
        Time = DataMat.Time;
        % --- dataset binding for the GLOBAL time cursor (loads no figure) ---
        [iDS, iResult] = bst_memory('LoadResultsFile', ResultsFile);
        if isempty(iDS)
            bst_progress('stop');
            return;
        end
        % --- within-hemisphere rank (kernel order is ascending lambda per hemi) ---
        hemi = double(R.ModeHemisphere(:));
        CompRank = zeros(numel(hemi), 1);
        for hh = 1:2
            idx = find(hemi == hh);
            CompRank(idx) = 1:numel(idx);
        end
        Info = struct('Values', double(R.Eigenvalues(:)), 'Component', hemi, 'CompRank', CompRank);

        % --- create a managed figure of our type (no 3D map) ---
        FigureId = db_template('FigureId');
        FigureId.Type     = 'EigenSpectrum';
        FigureId.SubType  = '';
        FigureId.Modality = '';
        [hFig, ~] = bst_figures('CreateFigure', iDS, FigureId, 'AlwaysCreate', ResultsFile);
        if isempty(hFig)
            bst_progress('stop');
            return;
        end
        % --- stash what the per-time redraw needs ---
        setappdata(hFig, 'iDS',           iDS);
        setappdata(hFig, 'iResult',       iResult);
        setappdata(hFig, 'Camp',          Camp);                 % [nMode x nTime] amplitude coeffs
        setappdata(hFig, 'SensorData',    M);                    % [nCh   x nTime] (for the dSPM projection)
        setappdata(hFig, 'Time',          Time(:)');
        setappdata(hFig, 'Info',          Info);
        setappdata(hFig, 'Kernel',        struct('ImagingKernel', double(R.ImagingKernel), ...
                                                  'DiracEigenFile', R.DiracEigenFile, ...
                                                  'Function', R.Function));
        setappdata(hFig, 'Basis',         []);                   % lazily loaded on first dSPM use
        setappdata(hFig, 'ResultsFile',   ResultsFile);
        setappdata(hFig, 'AxisMode',      'eigenvalue');         % canonical Dirac lambda axis
        setappdata(hFig, 'SpectrumSource','mode');               % 'mode' (amplitude) | 'vertex' (dSPM)
        setappdata(hFig, 'isStatic',      numel(Time) <= 2);
        % figure_timeseries display state (butterfly/column, scales, gain, zoom).
        TsInfo = db_template('TsInfo');
        TsInfo.DisplayMode   = 'butterfly';
        TsInfo.LinesLabels   = {'Left', 'Right'};
        TsInfo.LinesColor    = {[0.85 0.2 0.2; 0.2 0.3 0.85]};
        TsInfo.ShowEvents    = 0;
        TsInfo.AutoScaleY    = 1;     % spectra span orders of magnitude -> autoscale by default
        TsInfo.DefaultFactor = 1;
        TsInfo.XScale        = 'linear';
        TsInfo.YScale        = 'linear';
        setappdata(hFig, 'TsInfo', TsInfo);
        % Let bst_figures('ReloadFigures') redraw us (butterfly<->column toggle).
        setappdata(hFig, 'ReloadCall', {'view_eigenmode_spectrum', 'UpdateFigurePlot', hFig, 0});
        % First plot + show + select. (Controls are (re)created inside UpdateFigurePlot,
        % since figure_timeseries('PlotFigure') wipes figure-level uicontrols each redraw.)
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
% Build THROUGH figure_timeseries so it inherits the engine's appdata and callbacks
% (resize, mouse, scroll, display config); then rename + point keys at our handler.
function hFig = CreateFigure(FigureId)
    hFig = figure_timeseries('CreateFigure', FigureId);
    set(hFig, 'Name', 'Eigenspectrum', 'KeyPressFcn', @FigureKeyPressedCallback);
    set(hFig, bst_get('ResizeFunction'), @(h,ev)ResizeCallback(h, ev));
end


%% ===== GUI: resize = engine layout + margins for the mode x-axis label/title =====
function ResizeCallback(hFig, ev)
    figure_timeseries('ResizeCallback', hFig, ev);
    AdjustAxesMargins(hFig);
    PositionControls(hFig);
end


%% ===== GUI: on-figure controls (visible source + axis toggles, time readout) =====
function AddControls(hFig)
    if ~isempty(getappdata(hFig, 'hCtrlSource')) && ishandle(getappdata(hFig, 'hCtrlSource'))
        return;
    end
    src = getappdata(hFig, 'SpectrumSource'); if isempty(src), src = 'mode'; end
    ax  = getappdata(hFig, 'AxisMode');       if isempty(ax),  ax  = 'eigenvalue'; end
    hS = uicontrol(hFig, 'Style', 'popupmenu', 'String', {'Source: Amplitude', 'Source: dSPM'}, ...
        'Value', SrcToVal(src), 'Units', 'pixels', 'FontSize', 9, 'BackgroundColor', [1 1 1], ...
        'Tag', 'EigenSrcCtrl', ...
        'TooltipString', 'Spectrum source: amplitude mode coefficients, or the dSPM vertex map projected onto the Dirac modes', ...
        'Callback', @(h,e)SetSpectrumSource(hFig, ValToSrc(get(h, 'Value'))));
    hA = uicontrol(hFig, 'Style', 'popupmenu', 'String', {'Axis: Eigenvalue', 'Axis: Mode index', 'Axis: Wavelength'}, ...
        'Value', AxisToVal(ax), 'Units', 'pixels', 'FontSize', 9, 'BackgroundColor', [1 1 1], ...
        'Tag', 'EigenAxisCtrl', ...
        'TooltipString', 'Spectrum x-axis (all derived from the canonical Dirac eigenvalue)', ...
        'Callback', @(h,e)SetAxisMode(hFig, ValToAxis(get(h, 'Value'))));
    hT = uicontrol(hFig, 'Style', 'text', 'String', 't = n/a', 'Units', 'pixels', 'FontSize', 9, ...
        'BackgroundColor', get(hFig, 'Color'), 'ForegroundColor', [0 0 0], ...
        'HorizontalAlignment', 'right', 'Tag', 'EigenTimeCtrl');
    setappdata(hFig, 'hCtrlSource', hS);
    setappdata(hFig, 'hCtrlAxis',   hA);
    setappdata(hFig, 'hCtrlTime',   hT);
    PositionControls(hFig);
end

function PositionControls(hFig)
    hS = getappdata(hFig, 'hCtrlSource');
    hA = getappdata(hFig, 'hCtrlAxis');
    hT = getappdata(hFig, 'hCtrlTime');
    if isempty(hS) || ~ishandle(hS), return; end
    p = get(hFig, 'Position');  W = p(3);  H = p(4);
    y = H - 22;
    set(hS, 'Position', [6,   y, 150, 20]);
    set(hA, 'Position', [160, y, 150, 20]);
    if ~isempty(hT) && ishandle(hT), set(hT, 'Position', [max(W-130,316), y, 122, 16]); end
end

% source/axis <-> popup value
function v = SrcToVal(s),   if strcmpi(s, 'vertex'), v = 2; else, v = 1; end, end
function s = ValToSrc(v),   if v == 2, s = 'vertex'; else, s = 'mode'; end, end
function v = AxisToVal(a),  switch lower(a), case 'index', v = 2; case 'wavelength', v = 3; otherwise, v = 1; end, end
function a = ValToAxis(v),  switch v, case 2, a = 'index'; case 3, a = 'wavelength'; otherwise, a = 'eigenvalue'; end, end

function AdjustAxesMargins(hFig)
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'AxesGraph');
    if numel(hAxes) ~= 1 || ~ishandle(hAxes)
        return;   % multi-axes (split) layout is already stacked with proper margins
    end
    Scaling = bst_get('InterfaceScaling') / 100;
    marginBottom = 45 * Scaling;
    marginTop    = 30 * Scaling;
    figPos = get(hFig, 'Position');
    set(hAxes, 'Units', 'pixels');
    p = get(hAxes, 'Position');
    p(2) = marginBottom;
    p(4) = max(figPos(4) - marginBottom - marginTop, 1);
    set(hAxes, 'Position', p);
end


%% ===== GUI: redraw the eigenspectrum at the current global time =====
% Single-time-point view of the SAME coefficients view_eigen_timeseries shows.
% isFastUpdate=1 keeps the current zoom (fast time-stepping); =0 full replot.
function UpdateFigurePlot(hFig, isFastUpdate)
    global GlobalData;
    if (nargin < 2) || isempty(isFastUpdate)
        isFastUpdate = 0;
    end
    iDS      = getappdata(hFig, 'iDS');
    Info     = getappdata(hFig, 'Info');
    AxisMode = getappdata(hFig, 'AxisMode');
    Source   = getappdata(hFig, 'SpectrumSource');
    Time     = getappdata(hFig, 'Time');
    if isempty(iDS) || isempty(Info) || isempty(Time)
        return;
    end
    [~, iFig] = bst_figures('GetFigure', hFig);
    if isempty(iFig)
        return;
    end
    % --- current time index into our stored coefficient/data matrices ---
    ct = GlobalData.UserTimeWindow.CurrentTime;
    if isempty(ct)
        it = 1;
    else
        [~, it] = min(abs(Time - ct));
    end
    % --- modal coefficients at this time: amplitude (mode) or dSPM (vertex->mode) ---
    if strcmpi(Source, 'vertex')
        K = getappdata(hFig, 'Kernel');
        basis = LoadBasis(hFig);
        Mt = getappdata(hFig, 'SensorData');
        J  = K.ImagingKernel * Mt(:, it);                 % [3*nVert x 1] vertex field (measure = Function)
        ThetaCol = ProjectField(J, basis);                % [nMode x 1] Dirac modal coeffs
        srcTag = sprintf('vertex %s', K.Function);
    else
        Camp = getappdata(hFig, 'Camp');
        ThetaCol = Camp(:, it);                           % [nMode x 1] amplitude mode coeffs
        srcTag = 'mode (amplitude)';
    end
    spec = BuildModeSpectrum(Info, ThetaCol, AxisMode, 'power');
    setappdata(hFig, 'TimeVector', spec.x);   % the x-axis is the mode axis (NOT seconds)

    % --- layout from the display mode (butterfly=overlaid, column=split L/R) ---
    TsInfo  = getappdata(hFig, 'TsInfo');
    isSplit = strcmpi(TsInfo.DisplayMode, 'column');
    if isSplit
        F = {spec.F(1, :), spec.F(2, :)};
        TsInfo.LinesColor  = {[0.85 0.2 0.2], [0.2 0.3 0.85]};
        TsInfo.LinesLabels = {'Left', 'Right'};
        TsInfo.AxesLabels  = {'Left', 'Right'};
    else
        F = {spec.F};
        TsInfo.LinesColor  = {[0.85 0.2 0.2; 0.2 0.3 0.85]};
        TsInfo.LinesLabels = {'Left', 'Right'};
        TsInfo.AxesLabels  = {};
    end
    % Capture zoom on a fast time-step so it isn't reset every frame.
    keepZoom = isFastUpdate;
    xlimOld = {}; ylimOld = {};
    if keepZoom
        hAxOld = findobj(hFig, '-depth', 1, 'Tag', 'AxesGraph');
        if ~isempty(hAxOld)
            xlimOld = get(hAxOld, 'XLim');  if ~iscell(xlimOld), xlimOld = {xlimOld}; end
            ylimOld = get(hAxOld, 'YLim');  if ~iscell(ylimOld), ylimOld = {ylimOld}; end
        end
    end
    % Render each axes as butterfly (autoscale-fill); restore the menu-facing mode after.
    menuMode = TsInfo.DisplayMode;
    TsInfo.DisplayMode = 'butterfly';
    setappdata(hFig, 'TsInfo', TsInfo);
    figure_timeseries('PlotFigure', iDS, iFig, F, spec.x, isFastUpdate, []);
    TsInfo = getappdata(hFig, 'TsInfo');
    TsInfo.DisplayMode = menuMode;
    setappdata(hFig, 'TsInfo', TsInfo);
    % Axis labels: mode x-axis + modal-power y-axis (annotated with the source).
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'AxesGraph');
    yl = [spec.ylabel '  [' srcTag ']'];
    for k = 1:numel(hAxes)
        xlabel(hAxes(k), spec.label, 'Interpreter', 'tex');
        ylabel(hAxes(k), yl, 'Interpreter', 'tex');
    end
    % Hide the seconds time-cursor (meaningless on a mode/eigenvalue axis).
    Handles = GlobalData.DataSet(iDS).Figure(iFig).Handles;
    if isfield(Handles, 'hCursor')
        hc = [Handles.hCursor];
        hc = hc(ishandle(hc));
        if ~isempty(hc), set(hc, 'Visible', 'off'); end
    end
    % Keep only the Left/Right data lines in the legend.
    for h = findobj(hAxes, 'Type', 'line')'
        if ~strcmp(get(h, 'Tag'), 'DataLine')
            set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
        end
    end
    % (Re)create the on-figure toggles -- PlotFigure above wipes figure-level
    % uicontrols, so they must be rebuilt on every redraw (like the engine's buttons).
    AddControls(hFig);
    % Time readout in the top-right control (no axes title -> no overlap with the toggles).
    if ~isempty(GlobalData.UserTimeWindow.CurrentTime)
        tStr = sprintf('%1.3f s', GlobalData.UserTimeWindow.CurrentTime);
    else
        tStr = 'n/a';
    end
    hT = getappdata(hFig, 'hCtrlTime');
    if ~isempty(hT) && ishandle(hT)
        set(hT, 'String', ['t = ' tStr]);
    end
    AdjustAxesMargins(hFig);
    PositionControls(hFig);
    % Restore the user's zoom on a time-step: X always; Y unless autoscale is on.
    if keepZoom && ~isempty(xlimOld)
        hAxNew = findobj(hFig, '-depth', 1, 'Tag', 'AxesGraph');
        ts = getappdata(hFig, 'TsInfo');
        keepY = isempty(ts) || ~ts.AutoScaleY;
        if numel(hAxNew) == numel(xlimOld)
            for k = 1:numel(hAxNew)
                set(hAxNew(k), 'XLim', xlimOld{k});
                if keepY, set(hAxNew(k), 'YLim', ylimOld{k}); end
            end
        end
    end
end


%% ===== GUI: x-axis mode (eigenvalue/index/wavelength) =====
function SetAxisMode(hFig, mode)
    setappdata(hFig, 'AxisMode', mode);
    hA = getappdata(hFig, 'hCtrlAxis');
    if ~isempty(hA) && ishandle(hA), set(hA, 'Value', AxisToVal(mode)); end
    UpdateFigurePlot(hFig, 0);
end


%% ===== GUI: spectrum source (mode=amplitude / vertex=dSPM) =====
function SetSpectrumSource(hFig, src)
    setappdata(hFig, 'SpectrumSource', src);
    hS = getappdata(hFig, 'hCtrlSource');
    if ~isempty(hS) && ishandle(hS), set(hS, 'Value', SrcToVal(src)); end
    UpdateFigurePlot(hFig, 0);
end


%% ===== GUI: global time cursor moved =====
function CurrentTimeChangedCallback(hFig)
    if getappdata(hFig, 'isStatic')
        return;
    end
    UpdateFigurePlot(hFig, 1);   % fast update
end


%% ===== GUI: keys — e/k/w x-axis; a/d source; else -> figure_timeseries =====
function FigureKeyPressedCallback(hFig, ev)
    switch (ev.Key)
        case 'e', SetAxisMode(hFig, 'eigenvalue');
        case 'k', SetAxisMode(hFig, 'index');
        case 'w', SetAxisMode(hFig, 'wavelength');
        case 'a', SetSpectrumSource(hFig, 'mode');     % amplitude mode coefficients
        case 'd', SetSpectrumSource(hFig, 'vertex');   % dSPM (vertex) -> Dirac modes
        otherwise
            figure_timeseries('FigureKeyPressedCallback', hFig, ev);
    end
end


%% ===== load the Dirac eigenbasis (Phi) + operator mass (B), cached =====
function basis = LoadBasis(hFig)
    basis = getappdata(hFig, 'Basis');
    if ~isempty(basis), return; end
    K = getappdata(hFig, 'Kernel');
    if ~isfield(K, 'DiracEigenFile') || isempty(K.DiracEigenFile)
        error(['This result does not reference a Dirac eigenbasis (DiracEigenFile);' 10 ...
               'the dSPM (vertex) spectrum needs it. Recompute the Dirac sources.']);
    end
    E = in_bst_eigen(K.DiracEigenFile);
    if ~isfield(E, 'Phi') || isempty(E.Phi) || ~isfield(E, 'OperatorFile')
        error('Dirac eigen node is missing Phi/OperatorFile.');
    end
    O = in_bst_operator(E.OperatorFile);
    if ~isfield(O, 'Mass') || numel(O.Mass) ~= 2
        error('Dirac operator node is missing the 1x2 Mass (B).');
    end
    basis.Phi = E.Phi;                 % {1x2} [4Vh x K]
    basis.gv  = E.GlobalVertices;      % {1x2} global vertex indices
    basis.B   = O.Mass;                % {1x2} [4Vh x 4Vh] = kron(Mass_h, I4)
    setappdata(hFig, 'Basis', basis);
end


%% ===== project a per-vertex 3-vector field onto the Dirac eigenbasis =====
% c_k = <J, phi_k>_B with J embedded as a pure-imaginary quaternion field
% (psi = [0, Jx, Jy, Jz]); stacked [L block; R block] to match ModeHemisphere.
function col = ProjectField(J, basis)
    nVert = numel(J) / 3;
    Jr = reshape(J, 3, nVert);                 % rows x/y/z, cols vertices
    col = [];
    for h = 1:2
        vH  = basis.gv{h}(:);
        Phi = double(basis.Phi{h});
        B   = basis.B{h};
        nVh = numel(vH);
        psi = zeros(4 * nVh, 1);
        psi(2:4:end) = Jr(1, vH);
        psi(3:4:end) = Jr(2, vH);
        psi(4:4:end) = Jr(3, vH);
        col = [col; Phi' * (B * psi)];         %#ok<AGROW>
    end
    col = real(col);   % the relative-Dirac operator is real-symmetric -> real coeffs
end
