function varargout = view_eigenmodes(varargin)
% VIEW_EIGENMODES: Browse Laplace-Beltrami eigenmodes on a surface.
%
% USAGE:  hFig = view_eigenmodes(SurfaceFile)
%         [Grid, K, Info] = view_eigenmodes('BuildPairedGrid', Eigenmodes)
%         col  = view_eigenmodes('SynthColumn', PairedGrid, W)
%         view_eigenmodes('ModesChangedCallback', hFig)
%
% The viewer initialises the eigenmode lever (panel_eigenmodes) for the
% surface's K_paired modes and registers a single-frame Source result whose
% ImageGridAmp is re-synthesised live whenever the lever changes via
% bst_figures('FireModesChanged') -> ModesChangedCallback.
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

if (nargin >= 1) && ischar(varargin{1}) && any(strcmp(varargin{1}, {'BuildPairedGrid','SynthColumn','ModesChangedCallback'}))
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
[varargout{1:nargout}] = ViewFigure(varargin{:});
end


%% ===== PURE: paired display grid (column k = each component's CompRank==k mode) =====
function [Grid, K, Info] = BuildPairedGrid(Eig)
    nV = size(Eig.Vectors, 1);
    nK = size(Eig.Vectors, 2);
    if isfield(Eig, 'CompRank') && ~isempty(Eig.CompRank)
        CompRank  = Eig.CompRank(:);
        Component = Eig.Component(:);
    else
        CompRank  = (1:nK)';
        Component = ones(nK, 1);
    end
    K = max(CompRank);
    Grid = zeros(nV, K);
    for k = 1:K
        cols = find(CompRank == k);
        if ~isempty(cols)
            Grid(:, k) = sum(Eig.Vectors(:, cols), 2);   % disjoint support across components
        end
    end
    Info = struct('K', K, 'Component', Component, 'CompRank', CompRank, 'Values', Eig.Values(:));
end


%% ===== PURE: synthesized display column = PairedGrid * W(:) =====
function col = SynthColumn(PairedGrid, W) %#ok<DEFNU>
    W = W(:);
    if (numel(W) ~= size(PairedGrid,2))
        error('view_eigenmodes:SynthColumn: weight length (%d) must equal K_paired (%d).', numel(W), size(PairedGrid,2));
    end
    col = PairedGrid * W;
end


%% ===== GUI: display the paired modes as a transient registered Source result =====
function hFig = ViewFigure(SurfaceFile, ~)
    hFig = [];
    % Load eigenmodes
    [Eig, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if ~isComputed || isempty(Eig) || ~isfield(Eig, 'Vectors') || isempty(Eig.Vectors) || ~isfield(Eig, 'Values')
        bst_error(['No eigenmodes found on this surface.' 10 ...
                   'Right-click the cortex and run "Compute eigenmodes" first.'], 'View eigenmodes', 0);
        return;
    end
    % Resolve subject + intra study
    [sSubject, iSubject] = bst_get('SurfaceFile', SurfaceFile); %#ok<ASGLU>
    [sStudy, iStudy] = bst_get('AnalysisIntraStudy', iSubject);
    if isempty(iStudy)
        bst_error('Could not find the intra-subject study.', 'View eigenmodes', 0);
        return;
    end
    % Build paired display grid (mode k shows every component's rank-k mode)
    [Grid, Kp, Info] = BuildPairedGrid(Eig);

    % Initialise the lever for this surface and select mode 1
    panel_eigenmodes('ResetState', SurfaceFile, Kp);
    panel_eigenmodes('SetWindowShape', 'single');
    panel_eigenmodes('SetCurrentMode', 1);
    W0 = panel_eigenmodes('GetWeights');

    % Build a single-frame Source result (two-sample static — valid for source display)
    col0 = SynthColumn(Grid, W0);
    ResMat = db_template('resultsmat');
    ResMat.ImageGridAmp  = [col0 col0];     % 2 samples => valid static result
    ResMat.ImagingKernel = [];
    ResMat.nComponents   = 1;
    ResMat.Time          = [0 1];
    ResMat.SurfaceFile   = SurfaceFile;
    ResMat.HeadModelType = 'surface';
    ResMat.nAvg          = 1;
    ResMat.Leff          = 1;
    ResMat.ColormapType  = 'stat2';   % diverging, non-absolute (signed +/- lobes) without touching 'source'
    ResMat.Comment       = sprintf('Eigenmode viewer (%d modes/component)', Kp);
    ResMat = bst_history('add', ResMat, 'eigenmodes_view', 'Transient eigenmode viewer result');

    % Save to the intra study and register
    StudyDir   = bst_fileparts(file_fullpath(sStudy.FileName));
    OutputFile = bst_process('GetNewFilename', StudyDir, 'results_eigenview');
    bst_save(OutputFile, ResMat, 'v6');
    db_add_data(iStudy, OutputFile, ResMat);

    % Display via the standard surface-data path (colormap UI works natively)
    hFig = view_surface_data(SurfaceFile, file_short(OutputFile));
    if isempty(hFig)
        % Display failed: remove the transient result we just registered
        try
            file_delete(file_fullpath(OutputFile), 1);
            db_reload_studies(iStudy);
        catch
            % Non-fatal
        end
        bst_error('Could not open the surface figure.', 'View eigenmodes', 0);
        return;
    end
    set(hFig, 'Name', ['Eigenmodes: ' SurfaceFile]);

    % Tag the figure so FireModesChanged can route to ModesChangedCallback
    setappdata(hFig, 'EigenView', struct('SurfaceFile', SurfaceFile, 'PairedGrid', Grid, 'Info', Info, 'ResultsFile', file_short(OutputFile)));

    % Bottom-left legend (handle stored in EigenView so ModesChangedCallback can update it)
    hLabel = uicontrol('Style', 'text', 'String', '...', 'Units', 'Pixels', ...
        'Position', [6 0 560 20], 'HorizontalAlignment', 'left', ...
        'FontUnits', 'points', 'FontSize', bst_get('FigFont'), ...
        'ForegroundColor', [.9 .9 .9], 'BackgroundColor', [0 0 0], 'Parent', hFig);
    ev = getappdata(hFig, 'EigenView'); ev.LabelHandle = hLabel; setappdata(hFig, 'EigenView', ev);
    % Custom keyboard stepping (drives the lever)
    KeyPressFcn_bak = get(hFig, 'KeyPressFcn');
    set(hFig, 'KeyPressFcn', @KeyPress_Callback);
    % Auto-remove the transient result when the figure is destroyed
    set(hFig, 'DeleteFcn', @(h,e) CleanupResult());
    % Show + populate the panel
    gui_brainstorm('ShowToolTab', 'EigenModes');
    panel_eigenmodes('UpdatePanel', hFig);
    % Initial repaint + label (driven by ModesChangedCallback so it stays consistent)
    ModesChangedCallback(hFig);

    % ===== NESTED: keyboard navigation via the lever =====
    function KeyPress_Callback(h, keyEvent)
        cur = 1;
        try
            cur = panel_eigenmodes('GetCurrentMode');
        catch
        end
        switch (keyEvent.Key)
            case 'leftarrow',  panel_eigenmodes('SetCurrentMode', cur - 1);
            case 'rightarrow', panel_eigenmodes('SetCurrentMode', cur + 1);
            case 'pageup',     panel_eigenmodes('SetCurrentMode', cur + 10);
            case 'pagedown',   panel_eigenmodes('SetCurrentMode', cur - 10);
            case 'h'
                java_dialog('msgbox', ['Eigenmode viewer:' 10 '  Left/Right: step mode' 10 '  PgUp/PgDn: +/-10' 10 '  Use the "Spatial scale (eigenmodes)" panel to superpose a range.'], 'Eigenmode viewer');
            otherwise
                if ~isempty(KeyPressFcn_bak), KeyPressFcn_bak(h, keyEvent); end
        end
    end

    % ===== NESTED: delete the transient result on figure close =====
    function CleanupResult()
        try
            file_delete(file_fullpath(OutputFile), 1);
            db_reload_studies(iStudy);
        catch
            % Non-fatal: leave the node if cleanup fails (user can delete it)
        end
    end
end


%% ===== Re-synthesize the viewer's displayed column on a lever change =====
function ModesChangedCallback(hFig) %#ok<DEFNU>
    global GlobalData;
    ev = getappdata(hFig, 'EigenView');
    if isempty(ev), return; end
    W = panel_eigenmodes('GetWeights');
    if isempty(W) || (numel(W) ~= size(ev.PairedGrid,2)), return; end
    col = SynthColumn(ev.PairedGrid, W);
    [iDS, iResult] = bst_memory('GetDataSetResult', ev.ResultsFile);
    if isempty(iDS), return; end
    GlobalData.DataSet(iDS).Results(iResult).ImageGridAmp = [col col];
    panel_surface('UpdateSurfaceData', hFig);
    % Update bottom-left legend
    if isfield(ev, 'LabelHandle') && ishandle(ev.LabelHandle)
        K = size(ev.PairedGrid, 2);
        nKeep = nnz(W > 1e-6);
        lo = find(W > 1e-6, 1, 'first'); hi = find(W > 1e-6, 1, 'last');
        cur = 1;
        try, cur = panel_eigenmodes('GetCurrentMode'); catch, end
        % Representative eigenvalue for the current paired rank
        lamStr = '';
        if isfield(ev, 'Info') && isfield(ev.Info, 'Values') && isfield(ev.Info, 'CompRank')
            ik = find(ev.Info.CompRank(:) == cur, 1);
            if ~isempty(ik), lamStr = sprintf('    lambda = %.4g', ev.Info.Values(ik)); end
        end
        if (nKeep <= 1)
            str = sprintf('Mode %d / %d%s', cur, K, lamStr);
        else
            str = sprintf('Modes %d-%d / %d   (%d modes)', lo, hi, K, nKeep);
        end
        set(ev.LabelHandle, 'String', str);
    end
end
