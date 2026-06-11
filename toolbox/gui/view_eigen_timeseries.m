function hFig = view_eigen_timeseries(ResultsFile, hFigIn, DisplayKind)
% VIEW_EIGEN_TIMESERIES: Dirac eigenmode-coefficient view of a source result.
%
% USAGE:  hFig = view_eigen_timeseries(ResultsFile)                       % image (default)
%         hFig = view_eigen_timeseries(ResultsFile, [], 'image')          % (lambda x time) image
%         hFig = view_eigen_timeseries(ResultsFile, [], 'trace')          % stacked time series
%         hFig = view_eigen_timeseries(ResultsFile, hFig, DisplayKind)    % re-use a figure
%
% DESCRIPTION:
%     Displays the inverse source estimate expressed in the Dirac eigenbasis --
%     the mode-coefficient time series  c(t) = ImagingKernelMode * M(GoodChannel,t)
%     -- as the intermediate spectral view between the sensor time series and the
%     cortical reconstruction. Two display kinds:
%
%       'image' (default): a (lambda x time) pixel field via view_image_reg /
%               figure_image -- the continuum-honest view. The eigenvalue axis is
%               a true scale continuum (one pixel row per mode), so dynamic
%               signatures are directly visible: traveling waves as tilted stripes,
%               diffusion as a downward drift of the energy centroid.
%       'trace': stacked mode-coefficient traces (column mode), one row per mode.
%
%     Both are PURELY a coefficient view (no 3D map is opened): the global time
%     cursor advances the cortical source map separately.
%
%     PER-HEMISPHERE ORDER. Dirac modes are per-hemisphere objects (the operator
%     and eigensolve are block-diagonal by hemisphere; each mode is supported on
%     one hemisphere only). The kernel stacks them [L block ; R block], each
%     ascending in lambda. This view preserves that block order (it does NOT
%     globally interleave the two hemispheres), so a hemisphere-local pattern
%     stays a contiguous band. Rows 1..nL are left, nL+1..end are right, with a
%     divider drawn between them.
%
%     Requires a Dirac source result with the persisted eigenmode kernel
%     (ImagingKernelMode + Eigenvalues), produced by process_inverse_dirac.
%
% SEE ALSO: bst_inverse_dirac, process_inverse_dirac, view_image_reg, view_timeseries_matrix

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
    if (nargin < 3) || isempty(DisplayKind), DisplayKind = 'image'; end
    DisplayKind = lower(DisplayKind);
    isFirst = isempty(hFigIn);          % force default column layout only on first trace open
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
        Time = DataMat.Time;

        % --- PER-HEMISPHERE order: L block (asc lambda), then R block (asc lambda) ---
        % Mirrors the kernel's native [L;R] stacking; does NOT globally interleave.
        lam  = double(R.Eigenvalues(:));
        hemi = ones(numel(lam),1);
        if ~isempty(R.ModeHemisphere), hemi = double(R.ModeHemisphere(:)); end
        ordL = find(hemi == 1); [~, sL] = sort(lam(ordL), 'ascend'); ordL = ordL(sL);
        ordR = find(hemi == 2); [~, sR] = sort(lam(ordR), 'ascend'); ordR = ordR(sR);
        ord  = [ordL; ordR];
        c    = c(ord, :);
        lamS = lam(ord);
        hemiS= hemi(ord);
        nL   = numel(ordL);

        switch DisplayKind
            case 'trace'
                hFig = local_plot_trace(ResultsFile, c, Time, lamS, hemiS, nL, hFigIn, isFirst);
            otherwise   % 'image'
                hFig = local_plot_image(ResultsFile, c, Time, lamS, nL, hFigIn);
        end
        if isempty(hFig), bst_progress('stop'); return; end
        % Reload-safe: re-render into the SAME figure (recomputes c(t), a cheap
        % matrix product) instead of reloading the whole recording.
        setappdata(hFig, 'ReloadCall', {'view_eigen_timeseries', ResultsFile, hFig, DisplayKind});
    catch ME
        bst_progress('stop');
        bst_error(['Could not open the Dirac eigenmode view:' 10 ME.message], 'Dirac modes', 0);
        return;
    end
    bst_progress('stop');
end


%% ===== (lambda x time) IMAGE via figure_image =====
function hFig = local_plot_image(ResultsFile, c, Time, lamS, nL, hFigIn)
    nMode = size(c, 1);
    % Image volume [N1 x N2 x Ntime x Nfreq]: rows = modes, columns(time) on dim 3.
    F = reshape(c, [nMode, 1, numel(Time), 1]);
    Labels = cell(1, 4);
    Labels{1} = lamS(:);            % numeric -> continuous lambda Y axis (figure_image)
    Labels{2} = [];
    Labels{3} = Time(:)';
    Labels{4} = [];
    lamChar = char(955);            % literal Greek lambda (figure_image uses interpreter 'none')
    % Signed/symmetric diverging colormap ('stat2' = mandrill, isAbsoluteValues=0).
    [hFig, ~, ~] = view_image_reg(F, Labels, [1, 3], ...
        {['Eigenvalue (' lamChar ')'], 'Time (s)'}, ResultsFile, hFigIn, 'stat2', 1, [], 'a.u.');
    if isempty(hFig), return; end
    set(hFig, 'Name', 'Dirac eigenmode image');
    % Keep the user's time/lambda zoom across redraws (selecting a point fires a
    % time-cursor redraw that would otherwise snap back to the full range).
    setappdata(hFig, 'KeepZoomOnUpdate', 1);
    local_hemi_decor_image(hFig, nL, nMode, numel(Time));
end


%% ===== stacked TRACE view via figure_timeseries =====
function hFig = local_plot_trace(ResultsFile, c, Time, lamS, hemiS, nL, hFigIn, isFirst)
    labels = cell(numel(lamS), 1);
    for i = 1:numel(lamS)
        tag = 'L'; if hemiS(i) == 2, tag = 'R'; end
        labels{i} = sprintf('%s \\lambda=%.2g', tag, lamS(i));
    end
    % Units: amplitude-min-norm mode coefficients in the B-orthonormal (mass-weighted)
    % eigenbasis -- a signed amplitude in arbitrary units, NOT a calibrated source unit.
    [hFig, iDS, iFig] = view_timeseries_matrix(ResultsFile, {c}, Time, [], ...
        {'Dirac mode coefficients (\lambda \uparrow, L|R)'}, labels, [], hFigIn, [], 'a.u.');
    if isempty(hFig), return; end
    set(hFig, 'Name', 'Dirac eigenmode time series');
    % Default to COLUMN display on first open only; preserve the user's later choice.
    if isFirst
        TsInfo = getappdata(hFig, 'TsInfo');
        TsInfo.DisplayMode = 'column';
        setappdata(hFig, 'TsInfo', TsInfo);
        figure_timeseries('PlotFigure', iDS, iFig, {c}, Time, 0, []);
    end
    TsInfo = getappdata(hFig, 'TsInfo');
    if strcmpi(TsInfo.DisplayMode, 'column')
        local_eigen_axis_trace(hFig, iDS, iFig, lamS, nL);
    end
end


%% ===== hemisphere divider + L/R labels on the image (persist across redraws) =====
function local_hemi_decor_image(hFig, nL, nMode, nTime)
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'AxesImage');
    if isempty(hAxes), return; end
    hAxes = hAxes(1);
    % figure_image's UpdateFigurePlot deletes only ImageSurf/*Marker tags on redraw,
    % so these (different tags) survive a time-cursor refresh. Clear stale ones first.
    delete(findobj(hAxes, 'Tag', 'EigenHemiDiv'));
    delete(findobj(hAxes, 'Tag', 'EigenHemiLbl'));
    if nL < 1 || nL >= nMode, return; end    % single hemisphere -> no divider
    yB = nL + 1;                              % grid boundary between L (1..nL) and R (nL+1..)
    line([1, nTime+1], [yB yB], [2 2], 'Color', [0 0 0], 'LineWidth', 1.5, ...
        'Parent', hAxes, 'Tag', 'EigenHemiDiv', 'Clipping', 'on');
    xT = 1 + 0.012 * nTime;
    text(xT, 1 + 0.5*nL,        2, 'L', 'Parent', hAxes, 'Color', [0 0 0], 'FontWeight', 'bold', ...
        'VerticalAlignment', 'middle', 'Tag', 'EigenHemiLbl', 'Clipping', 'on');
    text(xT, yB + 0.5*(nMode-nL), 2, 'R', 'Parent', hAxes, 'Color', [0 0 0], 'FontWeight', 'bold', ...
        'VerticalAlignment', 'middle', 'Tag', 'EigenHemiLbl', 'Clipping', 'on');
end


%% ===== continuous EIGENVALUE y-axis (sparse per-hemisphere ticks + divider) for column trace =====
function local_eigen_axis_trace(hFig, iDS, iFig, lamS, nL)
    global GlobalData;
    try
        H = GlobalData.DataSet(iDS).Figure(iFig).Handles(1);
        if ~isfield(H, 'ChannelOffsets') || isempty(H.ChannelOffsets), return; end
        off = H.ChannelOffsets(:);                 % row k (per-hemi order) sits at off(k)
        n = numel(off);
        if n < 2, return; end
        hAx = findobj(hFig, '-depth', 1, 'Tag', 'AxesGraph');
        if isempty(hAx), return; end
        hAx = hAx(1);
        % ~4 sparse lambda ticks PER hemisphere block (so the L|R reset is honest)
        blocks = {1:min(nL, n)};
        if nL < n, blocks{end+1} = (nL+1):n; end
        yt = []; ytl = {};
        for b = 1:numel(blocks)
            idx = blocks{b};
            ti  = unique(round(linspace(idx(1), idx(end), min(4, numel(idx)))));
            yt  = [yt; off(ti)];                                                       %#ok<AGROW>
            ytl = [ytl; arrayfun(@(i) sprintf('%.2g', lamS(i)), ti(:), 'UniformOutput', 0)]; %#ok<AGROW>
        end
        [yt, oo] = sort(yt(:));  ytl = ytl(oo);    % YTick must be ascending
        set(hAx, 'YTickMode', 'manual', 'YTickLabelMode', 'manual', 'YTick', yt, 'YTickLabel', ytl);
        ylabel(hAx, 'Eigenvalue (\lambda)');
        % L|R divider
        delete(findobj(hAx, 'Tag', 'EigenHemiDiv'));
        if (nL >= 1) && (nL < n)
            yBnd = mean(off([nL, nL+1]));
            xl = get(hAx, 'XLim');
            line(xl, [yBnd yBnd], 'Color', [.4 .4 .4], 'LineStyle', '--', 'Parent', hAx, 'Tag', 'EigenHemiDiv');
        end
    catch
        % non-fatal: leave the default axis if the figure internals differ
    end
end
