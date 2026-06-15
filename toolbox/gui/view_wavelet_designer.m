function hFig = view_wavelet_designer(NodeFile)
% VIEW_FILTER_DESIGNER: Open a transient Dirac/eigenmode filterbank design session.
%
% USAGE:  view_wavelet_designer(EigenFile)        % design a new bank on this eigenbasis
%         view_wavelet_designer(WaveletFile)   % re-open a saved bank for editing
%
% Opens a figure_3d cortex preview (a temporary unconstrained results file) and docks
% panel_wavelet_designer. Clicking the cortex seeds a delta; the panel filters the seed
% live and pushes the field back to the preview. On Save/Cancel/figure-close the panel
% undocks and the preview figure + spectrum strip close (a transient session).
%
% SEE ALSO: panel_wavelet_designer, tess_eigen, db_add_wavelet, view_eigfilter_response
%
% Authors: Diellor Basha, 2026

    hFig = [];

    % --- resolve the node: eigen vs filterbank ---
    [~,~,~,iE] = bst_get('EigenFile', NodeFile);
    if ~isempty(iE)
        EigenMat  = load(file_fullpath(NodeFile));
        EigenFile = NodeFile;  loadBank = [];
    else
        [~,~,~,iFb] = bst_get('WaveletFile', NodeFile);
        if isempty(iFb)
            bst_error('Node is neither an eigen nor a filterbank file.', 'Filter designer', 0);
            return;
        end
        loadBank  = load(file_fullpath(NodeFile));
        EigenFile = loadBank.ParentEigen;
        EigenMat  = load(file_fullpath(EigenFile));
    end
    SurfaceFile = EigenMat.ParentSurface;

    % --- singleton: close any existing session ---
    gui_hide('WaveletDesigner');
    view_eigfilter_response('close');

    % --- preview figure: a transient unconstrained results node on the parent cortex,
    %     registered in the global default study (always exists) so view_surface_data
    %     can resolve it; removed on teardown ---
    nV = double(max(cellfun(@(x) max(x(:)), EigenMat.GlobalVertices)));

    % --- local cortical frame (U,V,N) per vertex from the manifold node ---
    bst_progress('text', 'Loading cortical frame (manifold)...');
    ManifoldMat = tess_manifold(SurfaceFile);   % find-or-load-or-create
    Frame = view_manifold('DeriveVertexFrame', ManifoldMat.Embedded, ManifoldMat.Gauge, nV);

    [tmpResultsFile, iStudyPrev] = i_create_preview(SurfaceFile, nV);
    [hFig, iDS, iResult] = i_open_preview(SurfaceFile, tmpResultsFile);
    if isempty(hFig); return; end
    setappdata(hFig, 'WaveletDesignerTemp',  tmpResultsFile);
    setappdata(hFig, 'WaveletDesignerStudy', iStudyPrev);
    setappdata(hFig, 'WaveletDesignerDS',    [iDS iResult]);

    % --- context callbacks handed to the panel ---
    ctxFn = struct();
    ctxFn.PushField = @(J) i_push_field(hFig, J);
    ctxFn.Close     = @() i_teardown(hFig);

    % --- build + dock the panel, make it the active tab ---
    bstPanel = panel_wavelet_designer('CreatePanel', EigenMat, EigenFile, hFig, ctxFn, Frame);
    gui_show(bstPanel, 'BrainstormTab', 'tools');
    try, gui_brainstorm('SetSelectedTab', 'WaveletDesigner', 0); catch, end %#ok<CTCH>

    % --- show the full 3D source vectors by default (Dirac design is vector-native) ---
    try, figure_3d('SetShowSourceVectors', hFig, 1, 1); catch, end %#ok<CTCH>

    % --- link: closing the figure ends the session; clicks seed the delta ---
    set(hFig, 'CloseRequestFcn', @(h,e) i_teardown(h));
    setappdata(hFig, 'WaveletDesignerPick', @(iVertex) panel_wavelet_designer('SetSeedVertex', 'WaveletDesigner', iVertex));

    % --- if re-opening a saved bank, restore its design + seed ---
    if ~isempty(loadBank)
        panel_wavelet_designer('LoadBank', 'WaveletDesigner', loadBank);
    end
end


%% ===== TRANSIENT PREVIEW RESULTS NODE (zero-filled unconstrained field) =====
function [tmpFile, iStudy] = i_create_preview(SurfaceFile, nV)
    R = db_template('resultsmat');
    R.ImageGridAmp  = zeros(3*nV, 2);          % [3nV x nTime] vector field (zeros)
    R.ImagingKernel = [];
    R.nComponents   = 3;                        % unconstrained -> vector display + quiver
    R.Time          = [0 1];
    R.HeadModelType = 'surface';
    R.SurfaceFile   = file_short(SurfaceFile);
    R.GoodChannel   = [];
    R.Comment       = 'Filter designer preview';
    iStudy  = -3;                               % global default study (always exists)
    tmpFile = db_add(iStudy, R);
end


%% ===== OPEN PREVIEW + LOCATE ITS DATASET/RESULT =====
function [hFig, iDS, iResult] = i_open_preview(SurfaceFile, tmpResultsFile)
    hFig = []; iDS = []; iResult = [];
    [hFig, iDS] = view_surface_data(SurfaceFile, tmpResultsFile, [], 'NewFigure');
    if isempty(hFig); return; end
    [iDS2, iResult] = bst_memory('GetDataSetResult', tmpResultsFile);
    if ~isempty(iDS2); iDS = iDS2; end
end


%% ===== PUSH A FILTERED FIELD TO THE PREVIEW =====
function i_push_field(hFig, J)
    global GlobalData;
    if isempty(hFig) || ~ishandle(hFig); return; end
    dsr = getappdata(hFig, 'WaveletDesignerDS');
    if isempty(dsr); return; end
    iDS = dsr(1); iResult = dsr(2);
    % normalize to [3nV x 2] (two identical time frames so the time slider is happy)
    J = real(J(:));
    Amp = [J, J];
    GlobalData.DataSet(iDS).Results(iResult).ImageGridAmp = Amp;
    % reset cached scaling so the colormap rescales to the new field magnitude
    TessInfo = getappdata(hFig, 'Surface');
    for k = 1:numel(TessInfo); TessInfo(k).DataMinMax = []; end
    setappdata(hFig, 'Surface', TessInfo);
    % the quiver shows the true vector field (per-vertex 3-vector)
    setappdata(hFig, 'QuiverVectorOverride', reshape(J, 3, [])');
    % repaint scalar magnitude + colormap
    panel_surface('UpdateSurfaceData', hFig);
    panel_surface('UpdateSurfaceColormap', hFig);
    % re-plot the source-vector quiver from the new field (reads the override)
    try, figure_3d('SetShowSourceVectors', hFig, 1, 1); catch, end %#ok<CTCH>
end


%% ===== TEARDOWN (idempotent) =====
function i_teardown(hFig)
    if isempty(hFig) || ~ishandle(hFig)
        gui_hide('WaveletDesigner'); view_eigfilter_response('close'); return;
    end
    if isappdata(hFig, 'WaveletDesignerClosing'); return; end     % guard re-entry
    setappdata(hFig, 'WaveletDesignerClosing', true);
    gui_hide('WaveletDesigner');
    view_eigfilter_response('close');
    tmpFile = getappdata(hFig, 'WaveletDesignerTemp');
    iStudy  = getappdata(hFig, 'WaveletDesignerStudy');
    % close the figure via Brainstorm's bookkeeping
    try, bst_figures('DeleteFigure', hFig, []); catch, delete(hFig); end %#ok<CTCH>
    % remove the transient preview node + file
    if ~isempty(tmpFile)
        try, file_delete(file_fullpath(tmpFile), 1); catch, end %#ok<CTCH>
        if ~isempty(iStudy)
            try, db_reload_studies(iStudy); catch, end %#ok<CTCH>
        end
    end
end
