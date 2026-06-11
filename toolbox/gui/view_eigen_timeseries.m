function hFig = view_eigen_timeseries(ResultsFile, hFigIn)
% VIEW_EIGEN_TIMESERIES: Dirac eigenmode-coefficient time series of a source result.
%
% USAGE:  hFig = view_eigen_timeseries(ResultsFile)
%
% DESCRIPTION:
%     Displays the inverse source estimate expressed in the Dirac eigenbasis --
%     the mode-coefficient time series  c(t) = ImagingKernelMode * M(GoodChannel,t)
%     -- as a standard time-series figure (rows = Dirac modes, ordered by
%     eigenvalue lambda). This is the intermediate spectral view that sits between
%     the sensor time series and the cortical reconstruction.
%
%     It is PURELY a time series (no 3D map is opened): the rows replace sensor
%     channels and everything else (time stepping, zoom, the global time cursor)
%     behaves exactly like the sensor viewer, via view_timeseries_matrix. Opening
%     the cortical source map ("View sources") separately, the single global time
%     cursor advances both figures together.
%
%     Requires a Dirac source result with the persisted eigenmode kernel
%     (ImagingKernelMode + Eigenvalues), produced by process_inverse_dirac.
%
% SEE ALSO: bst_inverse_dirac, process_inverse_dirac, view_timeseries_matrix

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

    if (nargin < 2), hFigIn = []; end
    isFirst = isempty(hFigIn);          % force default column layout only on first open
    hFig = [];
    bst_progress('start', 'Dirac modes', 'Computing eigenmode coefficients...');
    try
        % --- resolve the result + its data file (handles shared-kernel links) ---
        [~, DataFile] = file_resolve_link(ResultsFile);
        if isempty(DataFile)
            bst_progress('stop');
            bst_error(['No recording is associated with this result.' 10 ...
                'Open the Dirac source link on a data block (a recordings node).'], 'Dirac modes', 0);
            return;
        end
        R = in_bst_results(ResultsFile, 0, 'ImagingKernelMode', 'Eigenvalues', 'ModeHemisphere', 'GoodChannel');
        if ~isfield(R, 'ImagingKernelMode') || isempty(R.ImagingKernelMode)
            bst_progress('stop');
            bst_error(['This source file has no persisted Dirac eigenmode kernel.' 10 ...
                'Recompute with "Compute sources: Dirac eigenmodes" (it now stores the mode kernel).'], 'Dirac modes', 0);
            return;
        end

        % --- mode-coefficient time series: c(t) = Kmode * M(GoodChannel,:) ---
        DataMat = in_bst_data(DataFile, 'F', 'Time');
        gc = R.GoodChannel; if isempty(gc), gc = 1:size(DataMat.F,1); end
        c  = double(R.ImagingKernelMode) * double(DataMat.F(gc, :));    % [nMode x nTime]

        % --- order modes by eigenvalue (ascending = coarse -> fine spatial scale) ---
        lam = double(R.Eigenvalues(:));
        [lamS, ord] = sort(lam, 'ascend');
        c = c(ord, :);
        hemi = ones(numel(ord),1); if ~isempty(R.ModeHemisphere), hemi = R.ModeHemisphere(ord); end
        labels = cell(numel(lamS), 1);
        for i = 1:numel(lamS)
            tag = 'L'; if hemi(i) == 2, tag = 'R'; end
            labels{i} = sprintf('%s \\lambda=%.2g', tag, lamS(i));
        end

        % --- display via the standard time-series engine (decoupled; global-cursor linked) ---
        % Units: these are amplitude-min-norm mode coefficients in the B-orthonormal
        % (mass-weighted) eigenbasis -- a signed amplitude in arbitrary units, NOT a
        % calibrated source unit (pA.m); labelled 'a.u.'.
        [hFig, iDS, iFig] = view_timeseries_matrix(ResultsFile, {c}, DataMat.Time, [], ...
            {'Dirac mode coefficients (\lambda \uparrow)'}, labels, [], hFigIn, [], 'a.u.');
        if isempty(hFig), bst_progress('stop'); return; end
        set(hFig, 'Name', 'Dirac eigenmode time series');
        % Reload-safe: on a display toggle (butterfly<->column) Brainstorm calls
        % ReloadFigures, which without a ReloadCall reloads the whole recording. Point
        % it back here so the toggle just recomputes c(t) (a cheap matrix product) and
        % re-plots into the same figure.
        setappdata(hFig, 'ReloadCall', {'view_eigen_timeseries', ResultsFile, hFig});
        % Default to COLUMN display (rows ordered by lambda) on first open only;
        % preserve the user's choice on later toggles/reloads.
        if isFirst
            TsInfo = getappdata(hFig, 'TsInfo');
            TsInfo.DisplayMode = 'column';
            setappdata(hFig, 'TsInfo', TsInfo);
            figure_timeseries('PlotFigure', iDS, iFig, {c}, DataMat.Time, 0, []);
        end
        % In column mode the eigenvalue axis is a true CONTINUUM, not categorical
        % channels: replace the 800 per-row "lambda=..." labels with sparse ticks
        % (eigenvalue values) and a single axis title. (Butterfly keeps the a.u.
        % amplitude axis.)
        TsInfo = getappdata(hFig, 'TsInfo');
        if strcmpi(TsInfo.DisplayMode, 'column')
            local_eigen_axis(hFig, iDS, iFig, lamS);
        end
    catch ME
        bst_progress('stop');
        bst_error(['Could not open the Dirac eigenmode time series:' 10 ME.message], 'Dirac modes', 0);
        return;
    end
    bst_progress('stop');
end


%% ===== continuous EIGENVALUE y-axis (sparse ticks + one title) for column mode =====
function local_eigen_axis(hFig, iDS, iFig, lamS)
    global GlobalData;
    try
        H = GlobalData.DataSet(iDS).Figure(iFig).Handles(1);
        if ~isfield(H, 'ChannelOffsets') || isempty(H.ChannelOffsets), return; end
        off = H.ChannelOffsets(:);                 % row k (lambda-sorted) sits at off(k)
        n = numel(off);
        if n < 2, return; end
        ti  = unique(round(linspace(1, n, min(8, n))));   % ~8 ticks across the spectrum
        yt  = off(ti);
        ytl = arrayfun(@(i) sprintf('%.2g', lamS(i)), ti(:), 'UniformOutput', 0);
        [yt, oo] = sort(yt(:));  ytl = ytl(oo);    % YTick must be ascending
        hAx = findobj(hFig, '-depth', 1, 'Tag', 'AxesGraph');
        if isempty(hAx), return; end
        set(hAx(1), 'YTickMode', 'manual', 'YTickLabelMode', 'manual', 'YTick', yt, 'YTickLabel', ytl);
        ylabel(hAx(1), 'Eigenvalue (\lambda)');
    catch
        % non-fatal: leave the default axis if the figure internals differ
    end
end
