function varargout = view_eigenmodes_timeseries(varargin)
% VIEW_EIGENMODES_TIMESERIES: Plot eigenmode coefficients theta_k(t) over time.
%
% USAGE:  hFig = view_eigenmodes_timeseries(DataFile)
%         [iRows, Labels, Hemi] = view_eigenmodes_timeseries('GetBandTraces', Component, CompRank, kLo, kHi)
%         colors = view_eigenmodes_timeseries('HemiColors', Hemi)
%         view_eigenmodes_timeseries('ModesChangedCallback', hFig)
%
% One trace per eigenmode coefficient (sensor->mode transform), for the paired
% ranks in the EigenModes panel band. Each paired rank yields a left and a right
% trace. The figure tracks the panel band live via bst_figures('FireModesChanged').
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

if (nargin >= 1) && ischar(varargin{1}) && any(strcmp(varargin{1}, {'GetBandTraces','HemiColors','ModesChangedCallback','CurrentTimeChangedCallback','ReloadCallback'}))
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
[varargout{1:nargout}] = ViewFigure(varargin{:});
end


%% ===== GUI: build the eigenmode coefficient time series figure =====
function hFig = ViewFigure(DataFile)
    global GlobalData;
    hFig = [];
    if isempty(DataFile) || ~ischar(DataFile)
        bst_error('No data file provided.', 'Eigenmode time series', 0);
        return;
    end
    % ----- Study + head model + surface -----
    [sStudy, ~] = bst_get('AnyFile', DataFile);
    if isempty(sStudy) || ~isfield(sStudy, 'iHeadModel') || isempty(sStudy.iHeadModel) ...
            || (sStudy.iHeadModel < 1) || (length(sStudy.HeadModel) < sStudy.iHeadModel)
        bst_error('No head model available for this study.', 'Eigenmode time series', 0);
        return;
    end
    HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;
    HeadModelMat  = in_bst_headmodel(HeadModelFile, 0, 'HeadModelType', 'SurfaceFile');
    if ~strcmpi(HeadModelMat.HeadModelType, 'surface')
        bst_error('Eigenmode transform requires a surface head model.', 'Eigenmode time series', 0);
        return;
    end
    SurfaceFile = HeadModelMat.SurfaceFile;

    % ----- Eigenmodes -----
    [Eig, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if ~isComputed
        bst_error(['No eigenmodes on this surface.' 10 'Run "Compute eigenmodes" first.'], 'Eigenmode time series', 0);
        return;
    end

    % ----- Constrained gain (fixed orientation: [nch x nVert]) -----
    HM   = in_bst_headmodel(HeadModelFile, 1);
    Gain = double(HM.Gain);
    if size(Gain, 2) ~= size(Eig.Vectors, 1)
        bst_error(sprintf(['Head model has %d vertices but eigenmodes have %d.' 10 ...
            'Recompute the head model.'], size(Gain,2), size(Eig.Vectors,1)), 'Eigenmode time series', 0);
        return;
    end

    % ----- Load the recordings dataset (raw or imported); read the CURRENT window -----
    % This is the lazy path used by "Display on cortex": GetRecordingsValues loads
    % the current page on demand for raw files and returns the displayed window.
    iDS = bst_memory('LoadDataFile', DataFile);
    if isempty(iDS)
        bst_error('Could not load the recordings.', 'Eigenmode time series', 0);
        return;
    end
    Channels    = GlobalData.DataSet(iDS).Channel;
    ChannelFlag = GlobalData.DataSet(iDS).Measures.ChannelFlag;
    if isempty(Channels)
        bst_error('No channels found for this recording.', 'Eigenmode time series', 0);
        return;
    end
    if isempty(ChannelFlag)
        ChannelFlag = ones(length(Channels), 1);
    end
    iCh = good_channel(Channels, ChannelFlag, 'MEG');
    if isempty(iCh)
        iCh = good_channel(Channels, ChannelFlag, 'EEG');
    end
    if isempty(iCh)
        bst_error('No good MEG or EEG channels found.', 'Eigenmode time series', 0);
        return;
    end

    % ----- Transform kernel over ALL raw modes (complete pairing; rank-safe pinv) -----
    K_raw = double(Eig.nModes);
    Phi   = double(Eig.Vectors(:, 1:K_raw));
    [Kernel, ~] = bst_eigenmodes_transform(Gain(iCh, :), Phi);   % [K_raw x nCh]

    % ----- Current-window recordings -> coefficients (unscaled: isGradMagScale=0) -----
    F = bst_memory('GetRecordingsValues', iDS, iCh, 'UserTimeWindow', 0);   % [nCh x nTime]
    [TimeVector, ~] = bst_memory('GetTimeVector', iDS, [], 'UserTimeWindow');
    Theta = Kernel * F;                                                     % [K_raw x nTime]

    % ----- Cache (kernel + channels let us recompute on window change) -----
    cache = struct( ...
        'SurfaceFile', SurfaceFile, ...
        'DataFile',    DataFile, ...
        'Kernel',      Kernel, ...
        'GoodChannel', iCh, ...
        'Component',   Eig.Component(1:K_raw), ...
        'CompRank',    Eig.CompRank(1:K_raw), ...
        'Theta',       Theta, ...
        'TimeVector',  TimeVector, ...
        'WindowTime',  GlobalData.DataSet(iDS).Measures.Time);  % [t0 t1]; CurrentTimeChangedCallback uses it to detect a window change

    % ----- Ensure the lever is initialized for this surface (paired ranks) -----
    Kp = double(max(cache.CompRank));
    st = panel_eigenmodes('GetState');
    if ~file_compare(st.SurfaceFile, SurfaceFile) || (st.nModes ~= Kp)
        panel_eigenmodes('ResetState', SurfaceFile, Kp);
    end
    band = panel_eigenmodes('GetState');
    band = band.Band;
    cache.Band = band;   % last-rendered band; lets a reload replot without the panel

    % ----- First plot -----
    [Fband, Labels, colors] = BandData(cache, band);
    if isempty(Fband)
        bst_error('No eigenmodes in the selected band.', 'Eigenmode time series', 0);
        return;
    end
    hFig = view_timeseries_matrix(DataFile, {Fband}, cache.TimeVector, '', {'Eigenmode coefficients'}, Labels, colors, []);
    if isempty(hFig)
        return;
    end
    set(hFig, 'Name', ['Eigenmode time series: ' SurfaceFile]);
    setappdata(hFig, 'EigenTimeSeries', cache);

    % ----- Show + sync the panel -----
    gui_brainstorm('ShowToolTab', 'EigenModes');
    try
        panel_eigenmodes('RefreshControls');
    catch
        % Non-fatal: the panel still works, controls just won't pre-sync.
    end
end


%% ===== band -> {F, Labels, colors} from a cache's current Theta =====
function [F, Labels, colors] = BandData(cache, band)
    [iRows, Labels, Hemi] = GetBandTraces(cache.Component, cache.CompRank, band(1), band(2));
    if isempty(iRows)
        F = []; colors = {};
        return;
    end
    F      = cache.Theta(iRows, :);
    colors = HemiColors(Hemi);
end


%% ===== Redraw cache.Band into the existing figure (the new TsInfo.DisplayMode
%%       is honoured because view_timeseries_matrix preserves it on reuse) =====
function PlotBand(hFig, cache)
    [F, Labels, colors] = BandData(cache, cache.Band);
    if isempty(F)
        return;
    end
    view_timeseries_matrix(cache.DataFile, {F}, cache.TimeVector, '', {'Eigenmode coefficients'}, Labels, colors, hFig);
    % Restore our title (view_timeseries_matrix calls UpdateFigureName, which
    % overwrites it). The EigenTimeSeries appdata survives the redraw; re-assert
    % it defensively in case view_timeseries_matrix's behaviour ever changes.
    setappdata(hFig, 'EigenTimeSeries', cache);
    set(hFig, 'Name', ['Eigenmode time series: ' cache.SurfaceFile]);
end


%% ===== Re-slice the panel band and redraw (lever change) =====
function RefreshTraces(hFig)
    cache = getappdata(hFig, 'EigenTimeSeries');
    if isempty(cache)
        return;
    end
    st = panel_eigenmodes('GetState');
    if ~file_compare(st.SurfaceFile, cache.SurfaceFile)
        return;   % panel currently driving a different surface
    end
    cache.Band = st.Band;
    PlotBand(hFig, cache);
end


%% ===== Figure reload (display-mode toggle, montage change, "reload all"):
%%       replot the cached coefficients; never reload the raw file via view_matrix =====
function ReloadCallback(hFig) %#ok<DEFNU>
    cache = getappdata(hFig, 'EigenTimeSeries');
    if isempty(cache) || ~isfield(cache, 'Band') || isempty(cache.Band)
        return;
    end
    PlotBand(hFig, cache);
end


%% ===== PURE: band (paired-rank) -> raw-column traces =====
% For each paired rank k in kLo:kHi, emit its left column(s) then right column(s).
% Labels carry an L/R suffix only when the data has two components.
function [iRows, Labels, Hemi] = GetBandTraces(Component, CompRank, kLo, kHi)
    Component = Component(:);
    CompRank  = CompRank(:);
    isPaired  = any(Component == 2);
    iRows = [];
    Labels = {};
    Hemi = [];
    for k = kLo:kHi
        % Left (component 1) then right (component 2)
        for c = find(CompRank == k & Component == 1)'
            iRows(end+1) = c; %#ok<AGROW>
            Hemi(end+1)  = 1; %#ok<AGROW>
            if isPaired
                Labels{end+1} = sprintf('Mode %d L', k); %#ok<AGROW>
            else
                Labels{end+1} = sprintf('Mode %d', k);   %#ok<AGROW>
                Hemi(end)     = 0;
            end
        end
        for c = find(CompRank == k & Component == 2)'
            iRows(end+1) = c; %#ok<AGROW>
            Hemi(end+1)  = 2; %#ok<AGROW>
            Labels{end+1} = sprintf('Mode %d R', k); %#ok<AGROW>
        end
    end
end


%% ===== Lever changed: re-slice the band (no data re-read) =====
function ModesChangedCallback(hFig) %#ok<DEFNU>
    RefreshTraces(hFig);
end


%% ===== Displayed window changed (e.g. raw page scroll): re-read + recompute =====
function CurrentTimeChangedCallback(hFig, iDS) %#ok<DEFNU>
    global GlobalData;
    cache = getappdata(hFig, 'EigenTimeSeries');
    if isempty(cache)
        return;
    end
    if (nargin < 2) || isempty(iDS) || (iDS < 1) || (iDS > numel(GlobalData.DataSet))
        return;
    end
    % Only re-read when the displayed window bounds actually changed (a raw page
    % scroll). Moving the time cursor within the same page leaves Measures.Time
    % unchanged, so this is a cheap no-op on every cursor move.
    Win = GlobalData.DataSet(iDS).Measures.Time;
    if isequal(Win, cache.WindowTime)
        return;
    end
    F = bst_memory('GetRecordingsValues', iDS, cache.GoodChannel, 'UserTimeWindow', 0);
    [TimeVector, ~] = bst_memory('GetTimeVector', iDS, [], 'UserTimeWindow');
    cache.Theta      = cache.Kernel * F;
    cache.TimeVector = TimeVector;
    cache.WindowTime = Win;
    setappdata(hFig, 'EigenTimeSeries', cache);
    RefreshTraces(hFig);
end


%% ===== PURE: hemisphere id -> RGB cell array (left warm, right cool) =====
function colors = HemiColors(Hemi)
    Hemi = Hemi(:)';
    cL = [0.85 0.33 0.10];   % left  = warm orange
    cR = [0.00 0.45 0.74];   % right = cool blue
    c0 = [0.20 0.20 0.20];   % single-component = neutral
    colors = cell(1, numel(Hemi));
    for i = 1:numel(Hemi)
        switch Hemi(i)
            case 1, colors{i} = cL;
            case 2, colors{i} = cR;
            otherwise, colors{i} = c0;
        end
    end
end
