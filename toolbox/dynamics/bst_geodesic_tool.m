function varargout = bst_geodesic_tool( varargin )
% BST_GEODESIC_TOOL: interactive heat-distance disk on the cortex (the dynamics source-axis
% point+extent primitive). A seed vertex (center) + a geodesic radius (extent) define a disk
% via the heat-distance field; the tool draws it as a TRANSIENT overlay (never a scout) and the
% atom panel snapshots it into an atom. Reuses the tess_scout_area engine (toolbox/anatomy).
%
% USAGE:
%   bst_geodesic_tool('Toggle', onoff)        % enter/exit the cortical-pick mode
%   tf  = bst_geodesic_tool('IsActive')
%   bst_geodesic_tool('OnClick', hFig)        % figure_3d: clicked-vertex -> Seed + Draw
%   ok  = bst_geodesic_tool('OnScroll', n)    % figure_3d: grow/shrink + Draw; returns consumed
%   bst_geodesic_tool('Seed', SurfaceFile, vi)% headless: compute the disk around vi
%   bst_geodesic_tool('Grow', scrollCount)    % headless: re-threshold at a new radius
%   st  = bst_geodesic_tool('GetState')       % struct(seed,phi,radius,vertices,pos,SurfaceFile[,hFig]) | []
%   bst_geodesic_tool('Draw', hFig) / ('Clear', hFig)
%
% SEE ALSO: tess_scout_area, panel_bst_dynamics
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

eval(macro_method);
end


%% ===== persistent state cache (one active disk) =====
function out = Cache(in)
    persistent C;
    if (nargin >= 1), C = in; end
    out = C;
end


%% ===== TOGGLE: enter/exit the dynamics cortical-pick mode =====
function Toggle(onoff) %#ok<DEFNU>
    hFigures = bst_figures('GetFiguresForScouts');
    if isempty(hFigures), hFigures = bst_figures('GetFigureWithSurfaces'); end
    if onoff
        if isempty(hFigures)
            java_dialog('warning', 'Open a 3D cortex figure first.', 'Region tool');  return;
        end
        for hFig = hFigures(:)'
            setappdata(hFig, 'isDynamicsGeodesicPick', 1);
            setappdata(hFig, 'isSelectingCorticalSpot', 0);   % mutual exclusion with scout pick
            set(hFig, 'Pointer', 'cross');
        end
        Cache([]);                                            % fresh session
        SurfaceFile = i_tool_surface(hFigures);
        if ~isempty(SurfaceFile)
            bst_progress('start', 'Region tool', 'Pre-factorizing the heat operator...');
            try, tess_scout_area('prewarm', SurfaceFile); catch, end %#ok<CTCH>
            bst_progress('stop');
        end
    else
        for hFig = hFigures(:)'
            setappdata(hFig, 'isDynamicsGeodesicPick', 0);
            set(hFig, 'Pointer', 'arrow');
            Clear(hFig);
        end
    end
end


%% ===== IS ACTIVE =====
function tf = IsActive() %#ok<DEFNU>
    tf = false;
    hFigures = bst_figures('GetFiguresForScouts');
    for hFig = hFigures(:)'
        if isappdata(hFig, 'isDynamicsGeodesicPick') && getappdata(hFig, 'isDynamicsGeodesicPick')
            tf = true;  return;
        end
    end
end


%% ===== SEED: compute the geodesic disk around vi (headless core) =====
function Seed(SurfaceFile, vi) %#ok<DEFNU>
    if isempty(SurfaceFile) || isempty(vi), return; end
    R0 = 0.003;     % initial geodesic radius [m] = 3 mm
    isProg = ~bst_progress('isVisible');
    if isProg, bst_progress('start', 'Region tool', 'Computing geodesic distance...'); end
    [verts, phi] = tess_scout_area(SurfaceFile, vi, R0);
    if isProg, bst_progress('stop'); end
    Surf = in_tess_bst(SurfaceFile, 0);
    Cache(struct('seed',double(vi), 'phi',phi, 'radius',R0, 'vertices',verts, ...
                 'pos',Surf.Vertices(vi,:), 'SurfaceFile',SurfaceFile, ...
                 'Vertices',Surf.Vertices, 'Faces',Surf.Faces));
end


%% ===== GROW: re-threshold cached phi at a new radius (headless core) =====
function Grow(scrollCount) %#ok<DEFNU>
    c = Cache();
    if isempty(c), return; end
    STEP = 0.003;   % 3 mm per scroll tick
    R = max(STEP, c.radius - double(scrollCount) * STEP);   % scroll up (<0) grows
    c.vertices = tess_scout_area(c.SurfaceFile, c.seed, R, c.phi);   % reuse cached distance
    c.radius = R;
    Cache(c);
end


%% ===== GET STATE =====
function st = GetState() %#ok<DEFNU>
    st = Cache();
end


%% ===== ON CLICK (figure_3d): resolve clicked vertex -> seed + draw =====
function OnClick(hFig) %#ok<DEFNU>
    TessInfo = getappdata(hFig, 'Surface');  iTess = getappdata(hFig, 'iSurface');
    if isempty(iTess) || isempty(TessInfo) || isempty(TessInfo(iTess).hPatch), return; end
    if strcmpi(TessInfo(iTess).Name, 'Anatomy'), return; end
    [~, vout, vi] = select3d(TessInfo(iTess).hPatch);
    if isempty(vout) || isempty(vi), return; end
    Seed(TessInfo(iTess).SurfaceFile, vi);
    c = Cache();  c.hFig = hFig;  Cache(c);    % remember the figure for OnScroll redraws
    Draw(hFig);
end


%% ===== ON SCROLL (figure_3d): grow/shrink + redraw; returns consumed flag =====
function handled = OnScroll(scrollCount) %#ok<DEFNU>
    handled = false;
    if ~IsActive(), return; end
    c = Cache();
    if isempty(c) || ~isfield(c,'hFig') || ~ishandle(c.hFig), return; end
    Grow(scrollCount);
    Draw(c.hFig);
    handled = true;
end


%% ===== DRAW the transient disk overlay =====
function Draw(hFig) %#ok<DEFNU>
    c = Cache();
    if isempty(c) || isempty(hFig) || ~ishandle(hFig), return; end
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');  if isempty(hAxes), return; end
    hAxes = hAxes(1);  set(hAxes, 'NextPlot', 'add');     % low-level: avoid the axes-reset trap
    Clear(hFig);
    inR = false(size(c.Vertices,1), 1);  inR(c.vertices) = true;
    fIn = all(inR(c.Faces), 2);
    if ~any(fIn), return; end
    patch('Faces', c.Faces(fIn,:), 'Vertices', c.Vertices, 'Parent', hAxes, ...
        'FaceColor', [0.2 0.7 1.0], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'Tag', 'GeodesicToolDisk');
    % keep the Dynamics Source block fields in sync with the live disk (if the panel is open)
    if ~isempty(bst_get('PanelControls', 'Dynamics'))
        try, panel_bst_dynamics('SyncSource'); catch, end %#ok<CTCH>
    end
end


%% ===== CLEAR overlay =====
function Clear(hFig) %#ok<DEFNU>
    if isempty(hFig) || ~ishandle(hFig), return; end
    delete(findobj(hFig, 'Tag', 'GeodesicToolDisk'));
end


%% ===== cortex surface backing the figures (for prewarm) =====
function SurfaceFile = i_tool_surface(hFigures)
    SurfaceFile = '';
    for hFig = hFigures(:)'
        TessInfo = getappdata(hFig, 'Surface');
        if isempty(TessInfo), continue; end
        for i = 1:numel(TessInfo)
            sf = TessInfo(i).SurfaceFile;
            if ~isempty(sf) && strcmpi(file_gettype(sf), 'cortex'), SurfaceFile = sf;  return; end
        end
    end
end
