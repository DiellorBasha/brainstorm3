function varargout = view_dynamics( varargin )
% VIEW_DYNAMICS: Display a nested atom table on the cortex + open the dynamics panel.
%
% Loads a dynamics table (db_template('dynamicsmat')), opens its cortex surface, draws each
% atom group's spatial occurrences as colored markers (one color per group, offset along the
% vertex normal to clear the surface), and opens panel_bst_dynamics (the Record-style atom
% panel: group tree + occurrence list + File/Atoms menus). Extended/temporal groups (windows
% with no vertex) draw no markers but appear in the tree for navigation.
%
% USAGE:
%   [hFig, T] = view_dynamics(DynamicsFile)            % load a table, open cortex + panel
%   [hFig, T] = view_dynamics('FromResult', ResultsFile) % reuse/create a table for a Dirac result
%               view_dynamics('Redraw', hFig, T)       % redraw markers from an edited table
%
% SEE ALSO: panel_bst_dynamics, bst_dynamics, process_source_atoms
%
% Authors: Diellor Basha, 2026

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

    % --- Redraw verb: view_dynamics('Redraw', hFig, T) ---
    if (nargin >= 1) && ischar(varargin{1}) && strcmp(varargin{1}, 'Redraw')
        Redraw(varargin{2:end});
        return;
    end

    % --- PickScalar verb (pure, for tests): view_dynamics('PickScalar', Ht, Op) ---
    if (nargin >= 1) && ischar(varargin{1}) && strcmp(varargin{1}, 'PickScalar')
        varargout{1} = i_pick_scalar(varargin{2}, varargin{3});
        return;
    end

    % --- per-frame overlay verbs: view_dynamics('Overlay'|'RefreshOverlay', hFig) ---
    if (nargin >= 2) && ischar(varargin{1}) && any(strcmp(varargin{1}, {'Overlay','RefreshOverlay'}))
        i_dynamics_overlay(varargin{2});
        return;
    end

    % --- atom-field preview verbs (paint a realised W[nGv x nT] through the overlay) ---
    % AtomFrameScalar is pure (for tests): scatter frame iT of W onto a full-surface vector.
    if (nargin >= 1) && ischar(varargin{1}) && strcmp(varargin{1}, 'AtomFrameScalar')
        varargout{1} = i_atom_frame_scalar(varargin{2}, varargin{3}, varargin{4}, varargin{5});
        return;
    end
    if (nargin >= 2) && ischar(varargin{1}) && strcmp(varargin{1}, 'SetAtomField')
        SetAtomField(varargin{2:end});
        return;
    end
    if (nargin >= 2) && ischar(varargin{1}) && strcmp(varargin{1}, 'SetFilteredField')
        SetFilteredField(varargin{2:end});
        return;
    end
    if (nargin >= 2) && ischar(varargin{1}) && strcmp(varargin{1}, 'ClearAtomField')
        ClearAtomField(varargin{2});
        return;
    end
    if (nargin >= 2) && ischar(varargin{1}) && strcmp(varargin{1}, 'SetFilteredSensors')
        SetFilteredSensors(varargin{2:end});  return;
    end
    if (nargin >= 2) && ischar(varargin{1}) && strcmp(varargin{1}, 'ClearFilteredSensors')
        ClearFilteredSensors(varargin{2});    return;
    end

    % --- FromResult verb: open (reuse, else create) a dynamics table for a Dirac result ---
    SrcResult = '';
    if (nargin >= 1) && ischar(varargin{1}) && strcmp(varargin{1}, 'FromResult')
        SrcResult = varargin{2};
        DynamicsFile = AtomsFromResult(SrcResult);
        if isempty(DynamicsFile), return; end
    else
        DynamicsFile = varargin{1};
    end
    T = bst_dynamics('Load', DynamicsFile);
    % An empty table is fine: the panel opens and "Detect" creates the first atoms.
    SurfaceFile = T.SurfaceFile;
    if isempty(SurfaceFile) && ~isempty(T.Groups), SurfaceFile = T.Groups(1).SurfaceFile; end
    if isempty(SurfaceFile)
        error('view_dynamics: table has no SurfaceFile.');
    end

    % ===== PROVENANCE: source result + recording for the linked figures =====
    if isempty(SrcResult)
        for g = 1:numel(T.Groups)
            if ~isempty(T.Groups(g).ResultsFile), SrcResult = T.Groups(g).ResultsFile;  break; end
        end
    end
    DataFile = T.DataFile;
    if isempty(DataFile)
        for g = 1:numel(T.Groups)
            if ~isempty(T.Groups(g).DataFile), DataFile = T.Groups(g).DataFile;  break; end
        end
    end

    % ===== OPEN THE LINKED FIGURES =====
    % Spatial view = the source 3D cortex (atom markers are drawn on THIS cortex, and the
    % differential overlay is painted on it); temporal view = the recording time series. Both
    % follow the global time cursor. Fall back to a bare surface when the table has no
    % source/recording provenance.
    if ~isempty(SrcResult)
        hFig = i_open_source_figure(SrcResult);   % clean cortex; no raw recording timeseries (Design)
        if isempty(hFig), return; end
    else
        hFig = view_surface(SurfaceFile);
    end
    setappdata(hFig, 'DynamicsFile', DynamicsFile);
    Redraw(hFig, T);

    % ===== OPEN THE PANEL (docked tools tab, like Record) =====
    % The Dynamics panel owns the differential overlay on the source figure (operator combobox
    % + per-frame CustomOverlayFcn); reopen it fresh for this figure.
    try, gui_hide('Dynamics'); catch, end %#ok<CTCH>
    bstPanel = panel_bst_dynamics('CreatePanel');
    gui_show(bstPanel, 'BrainstormTab', 'tools');
    try, gui_brainstorm('SetSelectedTab', 'Dynamics', 0); catch, end %#ok<CTCH>
    panel_bst_dynamics('SetTarget', hFig, T);

    % Auto-teardown: closing the linked source figure (directly, or via "Close all figures")
    % tears the Dynamics session down too -- so the docked panel never orphans. figure_3d uses
    % CloseRequestFcn for its own cleanup, leaving the figure DeleteFcn free for this hook.
    set(hFig, 'DeleteFcn', @(h,e)panel_bst_dynamics('OnFigureDeleted', h));

    if (nargout >= 1), varargout{1} = hFig; end
    if (nargout >= 2), varargout{2} = T;    end
end


%% ===== GET A DYNAMICS TABLE FOR A DIRAC SOURCE RESULT (reuse newest, else empty) =====
% Returns the newest dynamics_*.mat in the result's study folder, or creates an EMPTY
% table carrying provenance (DataFile / SurfaceFile). The panel's "Detect windows" is the
% first atom-creation step -- no atoms are auto-populated here.
function DynamicsFile = AtomsFromResult(ResultsFile)
    DynamicsFile = '';
    R = in_bst_results(ResultsFile, 0, 'DataFile', 'SurfaceFile');
    if isempty(R.DataFile)
        bst_error('This result has no associated recording (DataFile).', 'Atoms', 0);  return;
    end
    studyDir = bst_fileparts(file_fullpath(R.DataFile));
    d = dir(fullfile(studyDir, 'dynamics_*.mat'));
    if ~isempty(d)
        [~, ix] = max([d.datenum]);
        DynamicsFile = fullfile(studyDir, d(ix).name);  return;
    end
    % None yet: create an empty table with provenance; Detect creates the first atoms.
    T = bst_dynamics('New', 'atoms');
    T.DataFile = R.DataFile;  T.SurfaceFile = R.SurfaceFile;
    DynamicsFile = bst_process('GetNewFilename', studyDir, 'dynamics');
    bst_dynamics('Save', DynamicsFile, T);
end


%% ===== REDRAW MARKERS + REGIONS FROM A (possibly edited) TABLE =====
function Redraw(hFig, T, showPhase)
    if (nargin < 3) || isempty(showPhase), showPhase = [1 1 1 1]; end
    if isempty(hFig) || ~ishandle(hFig), return; end
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    if isempty(hAxes), return; end
    hAxes = hAxes(1);
    set(hAxes, 'NextPlot', 'add');   % low-level line()/patch() never resets the cortex
    % Cache the surface (normals + faces) per figure to avoid reloading on every edit
    Surf = getappdata(hFig, 'DynamicsSurf');
    if isempty(Surf)
        SurfaceFile = T.SurfaceFile;  if isempty(SurfaceFile), SurfaceFile = T.Groups(1).SurfaceFile; end
        Surf = in_tess_bst(SurfaceFile, 0);
        setappdata(hFig, 'DynamicsSurf', Surf);
    end
    hasNorm = isfield(Surf, 'VertNormals') && ~isempty(Surf.VertNormals);
    nVert   = size(Surf.Vertices, 1);
    % Clear old markers, regions, selection
    delete(findobj(hAxes, '-regexp', 'Tag', '^AtomMarker'));
    delete(findobj(hAxes, '-regexp', 'Tag', '^AtomRegion'));
    delete(findobj(hAxes, 'Tag', 'AtomSel'));
    % One marker set + region patches per spatial group
    GroupsPosOff = cell(1, numel(T.Groups));
    for g = 1:numel(T.Groups)
        G = T.Groups(g);
        pk = i_phase_index(i_phase_type(G));
        if (pk >= 1) && ~showPhase(pk)
            GroupsPosOff{g} = zeros(0,3);   % phase filtered out: no marker, no region
            continue;
        end
        if isempty(G.pos)
            GroupsPosOff{g} = zeros(0,3);   % temporal-only group (no localized occurrence)
            continue;
        end
        nOcc = size(G.pos, 1);
        col  = G.color;  if isempty(col), col = [1 0 1]; end
        % per-occurrence offset seed position; unlocalized (NaN) rows stay NaN -> no marker
        po = nan(nOcc, 3);
        for o = 1:nOcc
            p = G.pos(o, :);
            if any(~isfinite(p)), continue; end
            if hasNorm && (o <= numel(G.vertices)) && isfinite(G.vertices(o))
                nrm = Surf.VertNormals(G.vertices(o), :);
            else
                nrm = p ./ max(sqrt(sum(p.^2)), eps);
            end
            po(o, :) = p + 0.002 * nrm;     % ~2 mm out, clears the surface
        end
        GroupsPosOff{g} = po;
        line(po(:,1), po(:,2), po(:,3), 'Parent', hAxes, ...
            'Marker','o', 'MarkerFaceColor',col, 'MarkerEdgeColor',[0 0 0], ...
            'MarkerSize',6, 'LineStyle','none', 'Tag', sprintf('AtomMarker%d', g));
        % captured regions: a translucent patch of faces fully inside region{o}
        if isfield(G, 'region') && ~isempty(G.region)
            for o = 1:min(nOcc, numel(G.region))
                rv = G.region{o};
                if isempty(rv), continue; end
                inReg = false(nVert, 1);  inReg(rv) = true;
                fIn = all(inReg(Surf.Faces), 2);
                if ~any(fIn), continue; end
                patch('Faces', Surf.Faces(fIn,:), 'Vertices', Surf.Vertices, 'Parent', hAxes, ...
                    'FaceColor', col, 'FaceAlpha', 0.35, 'EdgeColor', 'none', ...
                    'Tag', sprintf('AtomRegion%d_%d', g, o));
            end
        end
    end
    setappdata(hFig, 'GroupsPosOff', GroupsPosOff);
    % Selection marker (hidden until an occurrence is picked)
    line(NaN, NaN, NaN, 'Parent', hAxes, ...
        'Marker','o', 'MarkerSize',13, 'MarkerEdgeColor',[1 1 0], ...
        'LineWidth',2, 'LineStyle','none', 'Visible','off', 'Tag','AtomSel');
end


% Phase name -> filter index (peak=1 trough=2 rising=3 falling=4; 0 = not a phase group).
function k = i_phase_index(ph)
    if isempty(ph), k = 0; return; end
    switch lower(char(ph))
        case 'peak',    k = 1;
        case 'trough',  k = 2;
        case 'rising',  k = 3;
        case 'falling', k = 4;
        otherwise,      k = 0;
    end
end

% Readable phase type from a group's label suffix ('<band>_peak' etc.); '' if none.
% phase is now a NUMERIC value (radians), so the type is parsed from the label.
function t = i_phase_type(G)
    t = '';
    if isempty(G.label), return; end
    tok = regexp(G.label, '_(peak|trough|rising|falling)$', 'tokens', 'once');
    if ~isempty(tok), t = tok{1}; end
end


%% ===== operator name -> per-vertex scalar field (one of the Compute outputs) =====
function scal = i_pick_scalar(Ht, Op)
    switch Op
        case 'Divergence', scal = Ht.Div;
        case 'Curl',       scal = Ht.Curl;
        case 'Potential',  scal = Ht.Phi;
        case 'Stream',     scal = Ht.Psi;
        otherwise, error('view_dynamics:badOp', 'Unknown differential operator: %s', Op);
    end
end


%% ===== atom-field preview: paint a realised W[nGv x nT] through the overlay =====
% The panel's live preview. W is the realised atom field on the eigenbasis support gv (global
% vertex indices), already normalized; isSigned picks the colormap (peak/diverging vs density/seq).
function SetAtomField(hFig, W, gv, isSigned, V3)
    if nargin < 5, V3 = []; end
    if isempty(hFig) || ~ishandle(hFig), return; end
    D = getappdata(hFig, 'DynamicsOverlay');  if isempty(D), return; end
    D.AtomField  = W;  D.AtomGV = gv(:);  D.AtomSigned = logical(isSigned);  D.Op = 'atom';  D.AtomWin = [];
    setappdata(hFig, 'DynamicsOverlay', D);
    if isSigned, bst_colormaps('AddColormapToFigure', hFig, 'stat2');
    else,        bst_colormaps('AddColormapToFigure', hFig, 'source'); end
    i_dynamics_overlay(hFig);
    % Append quiver block: if V3 is non-empty and 3-column, set override and enable source vectors
    if (nargin >= 5) && ~isempty(V3) && (size(V3,2) == 3)
        setappdata(hFig, 'QuiverVectorOverride', V3);
        try, figure_3d('SetShowSourceVectors', hFig, D.iTess, 1); catch, end %#ok<CTCH>
    else
        setappdata(hFig, 'QuiverVectorOverride', []);
        try, figure_3d('SetShowSourceVectors', hFig, D.iTess, 0); catch, end %#ok<CTCH>
    end
end

% Preview: paint the filtered real source Ffilt[nV x nWin] over the recording window `win` (sample
% indices). Same scatter-paint as SetAtomField; `win` lets the overlay map the global cursor to the
% right window column so the filtered map scrubs with time.
function SetFilteredField(hFig, W, gv, win, isSigned)
    if isempty(hFig) || ~ishandle(hFig), return; end
    D = getappdata(hFig, 'DynamicsOverlay');  if isempty(D), return; end
    D.AtomField  = W;  D.AtomGV = gv(:);  D.AtomSigned = logical(isSigned);  D.Op = 'atom';  D.AtomWin = win(:)';
    setappdata(hFig, 'DynamicsOverlay', D);
    if isSigned, bst_colormaps('AddColormapToFigure', hFig, 'stat2');
    else,        bst_colormaps('AddColormapToFigure', hFig, 'source'); end
    i_dynamics_overlay(hFig);
end

% Drop the atom field and restore the native source paint (the 'none' state).
function ClearAtomField(hFig)
    if isempty(hFig) || ~ishandle(hFig), return; end
    D = getappdata(hFig, 'DynamicsOverlay');  if isempty(D), return; end
    D.Op = 'none';  D.AtomField = [];  D.AtomGV = [];
    setappdata(hFig, 'DynamicsOverlay', D);
    TessInfo = getappdata(hFig, 'Surface');
    if ~isempty(TessInfo) && (D.iTess <= numel(TessInfo))
        TessInfo(D.iTess).ColormapType = 'source';
        setappdata(hFig, 'Surface', TessInfo);
    end
    bst_colormaps('AddColormapToFigure', hFig, 'source');
    panel_surface('UpdateSurfaceData', hFig, D.iTess);
    panel_surface('UpdateSurfaceColormap', hFig);
    % Clear quiver override and disable source vectors
    setappdata(hFig, 'QuiverVectorOverride', []);
    try, D = getappdata(hFig,'DynamicsOverlay'); if ~isempty(D) && isfield(D,'iTess'), figure_3d('SetShowSourceVectors', hFig, D.iTess, 0); end, catch, end %#ok<CTCH>
end

%% ===== filtered-sensor overlay: stash Dfilt on the recording figure + trigger the draw hook =====
function SetFilteredSensors(hRec, Dfilt, tWin)
    if isempty(hRec) || ~ishandle(hRec), return; end
    setappdata(hRec, 'FilteredSensorsOverlay', struct('Dfilt',Dfilt, 'tWin',tWin));
    figure_timeseries('DrawFilteredOverlay', hRec);
end
function ClearFilteredSensors(hRec)
    if isempty(hRec) || ~ishandle(hRec), return; end
    setappdata(hRec, 'FilteredSensorsOverlay', []);
    figure_timeseries('DrawFilteredOverlay', hRec);
end

% Scatter frame iT of W (rows aligned with gv) onto a full-surface [nV x 1] vector (index-clamped).
function scal = i_atom_frame_scalar(W, gv, iT, nV)
    if isempty(iT) || iT < 1, iT = 1; end
    iT = min(iT, size(W, 2));
    scal = zeros(nV, 1);
    scal(gv) = W(:, iT);
end


%% ===== open a figure_3d SOURCE figure on an unconstrained result + install the overlay =====
function hFig = i_open_source_figure(SrcResult)
    global GlobalData;
    hFig = [];
    [iDS, iResult] = bst_memory('GetDataSetResult', SrcResult);
    if isempty(iResult)
        try, [iDS, iResult] = bst_memory('LoadResultsFileFull', SrcResult); catch, iResult = []; end %#ok<CTCH>
    end
    if isempty(iResult)
        bst_error('Could not load this source map.', 'Dynamics source view', 0);  return;
    end
    R = GlobalData.DataSet(iDS).Results(iResult);
    if isempty(R.nComponents) || (R.nComponents ~= 3)
        bst_error('The dynamics differential view requires an unconstrained (3-component) source.', 'Dynamics source view', 0);  return;
    end
    SurfaceFile = R.SurfaceFile;
    bst_progress('start', 'Dynamics source view', 'Loading the covariant operator...');
    Cov  = tess_operators(SurfaceFile, 'Covariant');           % find-or-create
    Surf = in_tess_bst(SurfaceFile, 0);  nV = size(Surf.Vertices, 1);
    bst_progress('stop');
    % Design view = a CLEAN cortex. We open with view_surface_data (it loads the source result for
    % Apply AND builds the proper colorbar/data scaffolding the overlay paint relies on), then BLANK
    % the displayed source so nothing of the real map shows. The only field ever painted is the
    % atom's impulse response (Design) or the filtered real source (Apply). No raw recording
    % timeseries is opened (Apply later shows the filtered sensor view instead).
    [hFig, iDSf] = view_surface_data(SurfaceFile, SrcResult, [], 'NewFigure');
    if isempty(hFig), return; end
    iTess = i_find_tess(hFig);
    TI = getappdata(hFig,'Surface');
    if iTess <= numel(TI)
        TI(iTess).DataThreshold = 0;
        TI(iTess).Data = zeros(size(TI(iTess).Data,1), 1);   % flat -> the real source is not shown
        setappdata(hFig,'Surface',TI);
    end
    D = struct('Cov',Cov, 'Op','none', 'AtomField',[], 'AtomGV',[], 'AtomSigned',false, 'AtomWin',[], ...
               'Cache',containers.Map('KeyType','double','ValueType','any'), ...
               'srcDS',iDSf, 'srcResult',iResult, 'iTess',iTess, 'nV',nV);
    setappdata(hFig, 'DynamicsOverlay', D);
    setappdata(hFig, 'CustomOverlayFcn', @(h) i_dynamics_overlay(h));   % fires per frame
    try, figure_3d('SetShowSourceVectors', hFig, iTess, 0); catch, end %#ok<CTCH>  % no raw quivers in Design
    panel_surface('UpdateSurfaceColormap', hFig);                                  % repaint the blanked field
end


%% ===== the CustomOverlayFcn: ephemeral per-frame differential scalar =====
function i_dynamics_overlay(hFig)
    if isempty(hFig) || ~ishandle(hFig), return; end
    D = getappdata(hFig, 'DynamicsOverlay');  if isempty(D), return; end
    % atom-field mode: paint the precomputed realised atom W(:,iT) (the panel's live preview).
    % W is already normalized (density or peak) by the caller; D.AtomSigned picks the colormap.
    if strcmpi(D.Op, 'atom')
        if isempty(D.AtomField), return; end
        TessInfo = getappdata(hFig, 'Surface');
        if isempty(TessInfo) || (D.iTess > numel(TessInfo)) || ~ishandle(TessInfo(D.iTess).hPatch), return; end
        [~, iT] = bst_memory('GetTimeVector', D.srcDS, D.srcResult, 'CurrentTimeIndex');
        if isfield(D,'AtomWin') && ~isempty(D.AtomWin)            % Preview: map global cursor -> window column
            iw = find(D.AtomWin == iT, 1);  if isempty(iw), iw = 1; end;  iT = iw;
        end
        TessInfo(D.iTess).Data         = i_atom_frame_scalar(D.AtomField, D.AtomGV, iT, D.nV);
        TessInfo(D.iTess).DataMinMax   = [min(D.AtomField(:)), max(D.AtomField(:))];
        if D.AtomSigned, TessInfo(D.iTess).ColormapType = 'stat2'; else, TessInfo(D.iTess).ColormapType = 'source'; end
        setappdata(hFig, 'Surface', TessInfo);
        panel_surface('UpdateSurfaceColormap', hFig);
        return;
    end
    % 'none' = no differential overlay: leave whatever the native pipeline painted (the RMS-norm
    % source scalar + 'source' colormap + raw source-vector quivers). Nothing to do.
    if strcmpi(D.Op, 'none'), return; end
    TessInfo = getappdata(hFig, 'Surface');
    if isempty(TessInfo) || (D.iTess > numel(TessInfo)) || ~ishandle(TessInfo(D.iTess).hPatch), return; end
    [~, iT] = bst_memory('GetTimeVector', D.srcDS, D.srcResult, 'CurrentTimeIndex');
    if isempty(iT) || iT < 1, iT = 1; end
    if ~isKey(D.Cache, iT)
        Jt = double(bst_memory('GetResultsValues', D.srcDS, D.srcResult, [], iT, 0));   % raw 3-vector
        if size(Jt,1) ~= 3*D.nV, return; end
        D.Cache(iT) = process_helmholtz('Compute', Jt, D.Cov);        % one call, all fields
        setappdata(hFig, 'DynamicsOverlay', D);
    end
    bst_colormaps('AddColormapToFigure', hFig, 'stat2');   % lazy: signed differential -> diverging colorbar
    scal = i_pick_scalar(D.Cache(iT), D.Op);
    TessInfo(D.iTess).Data         = scal;
    TessInfo(D.iTess).DataMinMax   = i_minmax(scal);
    TessInfo(D.iTess).ColormapType = 'stat2';      % the four differential operators are signed
    setappdata(hFig, 'Surface', TessInfo);
    panel_surface('UpdateSurfaceColormap', hFig);
end


%% ===== find the Source tess slot painted by view_surface_data =====
function iTess = i_find_tess(hFig)
    iTess = 1;  TessInfo = getappdata(hFig, 'Surface');
    for i = 1:numel(TessInfo)
        if ~isempty(TessInfo(i).DataSource) && strcmpi(TessInfo(i).DataSource.Type, 'Source'), iTess = i; return; end
    end
end


%% ===== symmetric color limits for a signed field (stat2) =====
function mm = i_minmax(s)
    m = max(abs(s(:)));  if isempty(m) || m <= 0, m = eps; end;  mm = [-m m];
end
