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
    if isempty(T.Groups)
        error('view_dynamics: no atom groups in %s', DynamicsFile);
    end
    SurfaceFile = T.SurfaceFile;
    if isempty(SurfaceFile), SurfaceFile = T.Groups(1).SurfaceFile; end
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
    % Spatial view = the Helmholtz source 3D (atom markers are drawn on THIS cortex); temporal
    % view = the recording time series. Both follow the global time cursor. Fall back to a bare
    % surface when the table has no source/recording provenance.
    if ~isempty(SrcResult)
        hFig = view_helmholtz(SrcResult);
        if ~isempty(DataFile)
            try, view_timeseries(DataFile); catch, end %#ok<CTCH>
        end
    else
        hFig = view_surface(SurfaceFile);
    end
    setappdata(hFig, 'DynamicsFile', DynamicsFile);
    Redraw(hFig, T);

    % ===== OPEN THE PANEL (docked tools tab, like Record/Helmholtz) =====
    try, gui_hide('Dynamics'); catch, end %#ok<CTCH>
    bstPanel = panel_bst_dynamics('CreatePanel');
    gui_show(bstPanel, 'BrainstormTab', 'tools');
    try, gui_brainstorm('SetSelectedTab', 'Dynamics', 0); catch, end %#ok<CTCH>
    panel_bst_dynamics('SetTarget', hFig, T);

    if (nargout >= 1), varargout{1} = hFig; end
    if (nargout >= 2), varargout{2} = T;    end
end


%% ===== GET A DYNAMICS TABLE FOR A DIRAC SOURCE RESULT (reuse newest, else detect) =====
function DynamicsFile = AtomsFromResult(ResultsFile)
    DynamicsFile = '';
    R = in_bst_results(ResultsFile, 0, 'DataFile');
    if isempty(R.DataFile)
        bst_error('This result has no associated recording (DataFile).', 'Source atoms', 0);  return;
    end
    studyDir = bst_fileparts(file_fullpath(R.DataFile));
    % Reuse the newest dynamics_*.mat already in the study folder, if any
    d = dir(fullfile(studyDir, 'dynamics_*.mat'));
    if isempty(d)
        % None yet: run the detector (default alpha) to create one
        bst_progress('start', 'Source atoms', 'Detecting source atoms...');
        bst_process('CallProcess', 'process_source_atoms', ResultsFile, [], 'freqband','alpha', 'npeaks',3);
        bst_progress('stop');
        d = dir(fullfile(studyDir, 'dynamics_*.mat'));
    end
    if isempty(d)
        bst_error(['No source atoms were created. The recording needs band events first -- ' 10 ...
                   'run "Detect bursts (phase polarity)" on it, then retry.'], 'Source atoms', 0);
        return;
    end
    [~, ix] = max([d.datenum]);
    DynamicsFile = fullfile(studyDir, d(ix).name);
end


%% ===== REDRAW MARKERS FROM A (possibly edited) TABLE =====
function Redraw(hFig, T)
    if isempty(hFig) || ~ishandle(hFig), return; end
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    if isempty(hAxes), return; end
    hAxes = hAxes(1);
    set(hAxes, 'NextPlot', 'add');   % low-level line() never resets the cortex
    % Cache the surface (normals) per figure to avoid reloading on every edit
    Surf = getappdata(hFig, 'DynamicsSurf');
    if isempty(Surf)
        SurfaceFile = T.SurfaceFile;  if isempty(SurfaceFile), SurfaceFile = T.Groups(1).SurfaceFile; end
        Surf = in_tess_bst(SurfaceFile, 0);
        setappdata(hFig, 'DynamicsSurf', Surf);
    end
    hasNorm = isfield(Surf, 'VertNormals') && ~isempty(Surf.VertNormals);
    % Clear old markers + selection
    delete(findobj(hAxes, '-regexp', 'Tag', '^AtomMarker'));
    delete(findobj(hAxes, 'Tag', 'AtomSel'));
    % One marker set per spatial group
    GroupsPosOff = cell(1, numel(T.Groups));
    for g = 1:numel(T.Groups)
        G = T.Groups(g);
        if isempty(G.pos)
            GroupsPosOff{g} = zeros(0,3);   % temporal-only group (window): no markers
            continue;
        end
        vtx = double(G.vertices(:));
        if hasNorm, nrm = Surf.VertNormals(vtx, :);
        else,       nrm = G.pos ./ max(sqrt(sum(G.pos.^2,2)), eps); end
        po = G.pos + 0.002 * nrm;           % ~2 mm out, clears the surface
        GroupsPosOff{g} = po;
        col = G.color;  if isempty(col), col = [1 0 1]; end
        line(po(:,1), po(:,2), po(:,3), 'Parent', hAxes, ...
            'Marker','o', 'MarkerFaceColor',col, 'MarkerEdgeColor',[0 0 0], ...
            'MarkerSize',6, 'LineStyle','none', 'Tag', sprintf('AtomMarker%d', g));
    end
    setappdata(hFig, 'GroupsPosOff', GroupsPosOff);
    % Selection marker (hidden until an occurrence is picked)
    line(NaN, NaN, NaN, 'Parent', hAxes, ...
        'Marker','o', 'MarkerSize',13, 'MarkerEdgeColor',[1 1 0], ...
        'LineWidth',2, 'LineStyle','none', 'Visible','off', 'Tag','AtomSel');
end
