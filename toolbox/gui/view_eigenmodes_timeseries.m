function varargout = view_eigenmodes_timeseries(varargin)
% VIEW_EIGENMODES_TIMESERIES: Plot eigenmode coefficients theta_k(t) over time.
%
% USAGE:  hFig = view_eigenmodes_timeseries(ResultsFile)   % an "eigenmode_harmonic" node
%         [iRows, Labels, Hemi] = view_eigenmodes_timeseries('GetBandTraces', Component, CompRank, kLo, kHi)
%         colors = view_eigenmodes_timeseries('HemiColors', Hemi)
%         view_eigenmodes_timeseries('ModesChangedCallback', hFig)
%
% One trace per eigenmode coefficient theta_k(t), read from a Harmonic results
% node's whitened kernel (EigenKernel = M̃) so the traces are consistent with the
% displayed Harmonic cortex map. Traces cover the paired ranks in the EigenModes
% panel band (left + right per rank) and track the band live via FireModesChanged.
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


%% ===== GUI: build the eigenmode coefficient time series from a Harmonic node =====
function hFig = ViewFigure(ResultsFile)
    global GlobalData;
    hFig = [];
    if isempty(ResultsFile) || ~ischar(ResultsFile)
        bst_error('Open this from an "Eigenmode HARMONIC" results node.', 'Eigenmode time series', 0);
        return;
    end
    % ----- Load the Harmonic node (kernel only: LoadFull=0 does NOT multiply by data) -----
    ResMat = in_bst_results(ResultsFile, 0, 'Function', 'EigenKernel', 'GoodChannel', 'DataFile', 'SurfaceFile');
    if isempty(ResMat) || ~isfield(ResMat, 'Function') || ~strcmpi(ResMat.Function, 'eigenmode_harmonic') ...
            || ~isfield(ResMat, 'EigenKernel') || isempty(ResMat.EigenKernel)
        bst_error(['This is not an Eigenmode HARMONIC results node.' 10 ...
                   'Compute sources with method "Harmonic (eigenmodes)" first.'], 'Eigenmode time series', 0);
        return;
    end
    Mtilde      = double(ResMat.EigenKernel);     % [K x nGoodCh]
    GoodChannel = ResMat.GoodChannel;             % index vector into the channel file
    DataFile    = ResMat.DataFile;
    SurfaceFile = ResMat.SurfaceFile;
    if isempty(DataFile)
        bst_error('This Harmonic node has no associated recordings to project.', 'Eigenmode time series', 0);
        return;
    end

    % ----- Eigenmodes for paired-rank trace mapping (first K columns) -----
    [Eig, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if ~isComputed
        bst_error('No eigenmodes on this surface.', 'Eigenmode time series', 0);
        return;
    end
    K = size(Mtilde, 1);
    if size(Eig.Vectors, 2) < K
        bst_error('Eigenmode count is smaller than the Harmonic kernel rank.', 'Eigenmode time series', 0);
        return;
    end

    % ----- Load the recordings dataset (raw or imported); read the CURRENT window -----
    iDS = bst_memory('LoadDataFile', DataFile);
    if isempty(iDS)
        bst_error('Could not load the recordings.', 'Eigenmode time series', 0);
        return;
    end
    F = bst_memory('GetRecordingsValues', iDS, GoodChannel, 'UserTimeWindow', 0);   % [nGoodCh x nTime]
    [TimeVector, ~] = bst_memory('GetTimeVector', iDS, [], 'UserTimeWindow');
    Theta = Mtilde * F;                                                             % [K x nTime]

    % ----- Cache (M̃ + GoodChannel let SyncWindow recompute on a window change) -----
    cache = struct( ...
        'SurfaceFile', SurfaceFile, ...
        'DataFile',    DataFile, ...
        'Kernel',      Mtilde, ...
        'GoodChannel', GoodChannel, ...
        'Component',   Eig.Component(1:K), ...
        'CompRank',    Eig.CompRank(1:K), ...
        'Theta',       Theta, ...
        'TimeVector',  TimeVector, ...
        'WindowTime',  GlobalData.DataSet(iDS).Measures.Time);  % [t0 t1]; window-change detection

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


%% ===== Re-read the displayed window into the cache IF its bounds changed (a raw
%%       page change). Returns the (possibly updated) cache and a changed flag. =====
function [cache, changed] = SyncWindow(hFig, iDS, cache)
    global GlobalData;
    changed = false;
    if isempty(iDS) || (iDS < 1) || (iDS > numel(GlobalData.DataSet))
        return;
    end
    if isempty(GlobalData.DataSet(iDS).DataFile) || ~file_compare(GlobalData.DataSet(iDS).DataFile, cache.DataFile)
        return;   % not our recordings dataset
    end
    Win = GlobalData.DataSet(iDS).Measures.Time;
    if isequal(Win, cache.WindowTime)
        return;   % same window: nothing to re-read
    end
    F = bst_memory('GetRecordingsValues', iDS, cache.GoodChannel, 'UserTimeWindow', 0);
    [TimeVector, ~] = bst_memory('GetTimeVector', iDS, [], 'UserTimeWindow');
    cache.Theta      = cache.Kernel * F;
    cache.TimeVector = TimeVector;
    cache.WindowTime = Win;
    setappdata(hFig, 'EigenTimeSeries', cache);
    changed = true;
end


%% ===== Figure reload (display-mode toggle, montage change, raw PAGE change,
%%       "reload all"): re-read the current window if it changed, then replot the
%%       cached coefficients. Never reloads the raw file via view_matrix. =====
function ReloadCallback(hFig, iDS) %#ok<DEFNU>
    if (nargin < 2), iDS = []; end
    cache = getappdata(hFig, 'EigenTimeSeries');
    if isempty(cache) || ~isfield(cache, 'Band') || isempty(cache.Band)
        return;
    end
    cache = SyncWindow(hFig, iDS, cache);
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


%% ===== Time cursor moved: re-read + replot only if it crossed into a new raw
%%       page (window bounds changed). A move within the page is a cheap no-op. =====
function CurrentTimeChangedCallback(hFig, iDS) %#ok<DEFNU>
    if (nargin < 2), iDS = []; end
    cache = getappdata(hFig, 'EigenTimeSeries');
    if isempty(cache) || ~isfield(cache, 'Band') || isempty(cache.Band)
        return;
    end
    [cache, changed] = SyncWindow(hFig, iDS, cache);
    if changed
        PlotBand(hFig, cache);
    end
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
