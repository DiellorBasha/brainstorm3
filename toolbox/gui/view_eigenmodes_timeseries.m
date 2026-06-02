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

if (nargin >= 1) && ischar(varargin{1}) && any(strcmp(varargin{1}, {'GetBandTraces','HemiColors','ModesChangedCallback'}))
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
[varargout{1:nargout}] = ViewFigure(varargin{:});
end


%% ===== GUI: build the eigenmode coefficient time series figure =====
function hFig = ViewFigure(DataFile)
    hFig = [];
    if isempty(DataFile)
        bst_error('No data file provided.', 'Eigenmode time series', 0);
        return;
    end
    % ----- Study + head model + surface -----
    [sStudy, iStudy] = bst_get('AnyFile', DataFile); %#ok<ASGLU>
    if isempty(sStudy) || ~isfield(sStudy, 'iHeadModel') || isempty(sStudy.iHeadModel) || (sStudy.iHeadModel < 1)
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

    % ----- Channels + recordings -----
    ChannelFile = bst_get('ChannelFileForStudy', sStudy.FileName);
    if isempty(ChannelFile)
        bst_error('No channel file found.', 'Eigenmode time series', 0);
        return;
    end
    ChannelMat = in_bst_channel(ChannelFile);
    DataMat    = in_bst_data(DataFile);
    if isstruct(DataMat.F)
        bst_error('Eigenmode time series requires imported (non-raw) recordings.', 'Eigenmode time series', 0);
        return;
    end
    if isfield(DataMat, 'ChannelFlag') && ~isempty(DataMat.ChannelFlag)
        ChannelFlag = DataMat.ChannelFlag;
    else
        ChannelFlag = ones(length(ChannelMat.Channel), 1);
    end
    iCh = good_channel(ChannelMat.Channel, ChannelFlag, 'MEG');
    if isempty(iCh)
        iCh = good_channel(ChannelMat.Channel, ChannelFlag, 'EEG');
    end
    if isempty(iCh)
        bst_error('No good MEG or EEG channels found.', 'Eigenmode time series', 0);
        return;
    end

    % ----- Transform: full coefficient matrix Theta [K_raw x nTime] -----
    nCh   = numel(iCh);
    K_raw = min(nCh, double(Eig.nModes));
    Phi   = double(Eig.Vectors(:, 1:K_raw));
    [Kernel, ~] = bst_eigenmodes_transform(Gain(iCh, :), Phi);   % [K_raw x nCh]
    Theta = Kernel * double(DataMat.F(iCh, :));                  % [K_raw x nTime]

    % ----- Cache everything the live refresh needs -----
    cache = struct( ...
        'SurfaceFile', SurfaceFile, ...
        'DataFile',    DataFile, ...
        'Theta',       Theta, ...
        'Component',   Eig.Component(1:K_raw), ...
        'CompRank',    Eig.CompRank(1:K_raw), ...
        'TimeVector',  DataMat.Time);

    % ----- Ensure the lever is initialized for this surface (paired ranks) -----
    Kp = double(max(cache.CompRank));
    st = panel_eigenmodes('GetState');
    if ~file_compare(st.SurfaceFile, SurfaceFile) || (st.nModes ~= Kp)
        panel_eigenmodes('ResetState', SurfaceFile, Kp);
    end
    band = panel_eigenmodes('GetState');
    band = band.Band;

    % ----- First plot -----
    [iRows, Labels, Hemi] = GetBandTraces(cache.Component, cache.CompRank, band(1), band(2));
    if isempty(iRows)
        bst_error('No eigenmodes in the selected band.', 'Eigenmode time series', 0);
        return;
    end
    F      = cache.Theta(iRows, :);
    colors = HemiColors(Hemi);
    hFig = view_timeseries_matrix(DataFile, {F}, cache.TimeVector, '', {'Eigenmode coefficients'}, Labels, colors, []);
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
