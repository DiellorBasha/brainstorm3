function varargout = view_eigenmodes(varargin)
% VIEW_EIGENMODES: Browse Laplace-Beltrami eigenmodes on a surface.
%
% USAGE:  hFig = view_eigenmodes(SurfaceFile)
%         [Grid, K, Info] = view_eigenmodes('BuildPairedGrid', Eigenmodes)
%
% The viewer displays each component's rank-k mode together (mode k shows both
% hemispheres) as a registered Brainstorm Source result, so the standard colormap
% UI applies; Left/Right arrows step modes.
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

if (nargin >= 1) && ischar(varargin{1}) && strcmp(varargin{1}, 'BuildPairedGrid')
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
    [Grid, K, Info] = BuildPairedGrid(Eig);

    % Build a Source result (modes as the "time" axis)
    ResMat = db_template('resultsmat');
    ResMat.ImageGridAmp  = Grid;
    ResMat.ImagingKernel = [];
    ResMat.nComponents   = 1;
    ResMat.Time          = 1:K;
    ResMat.SurfaceFile   = SurfaceFile;
    ResMat.HeadModelType = 'surface';
    ResMat.nAvg          = 1;
    ResMat.Leff          = 1;
    ResMat.ColormapType  = 'stat2';   % diverging, non-absolute (signed +/- lobes) without touching 'source'
    ResMat.Comment       = sprintf('Eigenmode viewer (%d modes/component, %d component(s))', K, Eig.nComponents);
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

    % Current mode (closure state), starting at mode 1
    curMode = 1;
    % Bottom-left legend
    hLabel = uicontrol('Style', 'text', 'String', '...', 'Units', 'Pixels', ...
        'Position', [6 0 560 20], 'HorizontalAlignment', 'left', ...
        'FontUnits', 'points', 'FontSize', bst_get('FigFont'), ...
        'ForegroundColor', [.9 .9 .9], 'BackgroundColor', [0 0 0], 'Parent', hFig);
    % Custom keyboard stepping (drives the global time = mode index) + legend
    KeyPressFcn_bak = get(hFig, 'KeyPressFcn');
    set(hFig, 'KeyPressFcn', @KeyPress_Callback);
    % Auto-remove the transient result when the figure is destroyed
    set(hFig, 'DeleteFcn', @(h,e) CleanupResult());
    % Initial position + legend
    SetMode(1);

    % ===== NESTED: move to mode k and refresh the legend =====
    function SetMode(k)
        curMode = min(max(round(k), 1), K);
        panel_time('SetCurrentTime', curMode);   % Time vector is 1:K, so time==mode index
        % Per-component eigenvalue(s) for this mode (rank == curMode)
        lv = Info.Values(Info.CompRank == curMode);
        if numel(lv) >= 2
            lamStr = sprintf('lambda = [%.4g, %.4g]', lv(1), lv(2));
        elseif ~isempty(lv)
            lamStr = sprintf('lambda = %.4g', lv(1));
        else
            lamStr = 'lambda = n/a';
        end
        set(hLabel, 'String', sprintf('Mode %d / %d     %s', curMode, K, lamStr));
    end

    % ===== NESTED: keyboard navigation =====
    function KeyPress_Callback(h, keyEvent)
        switch (keyEvent.Key)
            case 'leftarrow',  SetMode(curMode - 1);
            case 'rightarrow', SetMode(curMode + 1);
            case 'pageup',     SetMode(curMode + 10);
            case 'pagedown',   SetMode(curMode - 10);
            case 'h'
                java_dialog('msgbox', ...
                    ['Eigenmode viewer shortcuts:' 10 10 ...
                     '   Left / Right arrow  :  previous / next mode' 10 ...
                     '   Page Up / Page Down :  +/- 10 modes' 10 ...
                     '   H                   :  this help'], 'Eigenmode viewer');
            otherwise
                if ~isempty(KeyPressFcn_bak)
                    KeyPressFcn_bak(h, keyEvent);
                end
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
