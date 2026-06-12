function hFig = view_connection_phase(SurfaceFile, varargin)
% VIEW_CONNECTION_PHASE: Connection-Laplacian eigenmode viewer (phase + field).
%
% USAGE:  hFig = view_connection_phase(SurfaceFile)
%         hFig = view_connection_phase(SurfaceFile, 'nModes', 50, 'MaxArrows', 2000)
%
% DESCRIPTION:
%     Opens an opaque cortex figure showing a connection-Laplacian eigenmode:
%       - the FreeSurfer-gauge phase as a cyclic surface colormap (the between-
%         subject location coordinate);
%       - the gauge-independent 3D tangent field as headless glyphs (the
%         reference-free intrinsic winding);
%       - the field singularities as markers.
%     Keyboard: m / Shift+m  previous / next mode rank (live recompute);
%               g  toggle glyphs;  p  toggle singularities;
%               Left / Right  fewer / more glyphs;  other keys drive the
%               standard Brainstorm view controls.
%     Data pipeline: bst_conn_eigenmodes_ensure -> nxr vertexFrame ->
%     manifold node gauge frame -> bst_conn_phase.
%
%     Requires the nxr-compute plugin and a FreeSurfer-registered cortex.
%
% SEE ALSO: bst_conn_phase, bst_conn_eigenmodes_ensure, tess_manifold, view_eigenmodes

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

%% ===== PARSE INPUTS =====
nModes    = 50;
MaxArrows = 2000;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'nmodes',    nModes = varargin{i+1};
        case 'maxarrows', MaxArrows = varargin{i+1};
    end
end

bst_progress('start', 'Connection phase viewer', 'Preparing eigenmodes and frames...');

try

%% ===== DATA PIPELINE =====
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
if ~isOk
    error('view_connection_phase:nxr', 'nxr-compute required: %s', errMsg);
end
bst_plugin('Load', 'nxr-compute');

% SurfaceFile must be a DB-registered surface (view_surface_data requires the
% surface to live in the protocol's anatomy folder). This matches the
% view_eigenmodes contract.
ConnEig = bst_conn_eigenmodes_ensure(SurfaceFile, nModes);
TessMat = in_tess_bst(SurfaceFile);
Vtx = TessMat.Vertices;
Fcs = double(TessMat.Faces);
mctx   = nxr.manifold.context(Vtx, Fcs);
vFrame = nxr.manifold.measure.vertexFrame(mctx);
[Uv, Vv] = GetManifoldFrame(SurfaceFile);    % manifold gauge frame (vertex domain)
FsFrame  = struct('e1', Uv, 'e2', Vv);

R = bst_conn_phase(ConnEig, vFrame, 'Rank', 1, 'FsFrame', FsFrame, 'nSing', 2);

%% ===== SURFACE FIGURE + PHASE COLORMAP =====
% Register a transient Source result whose ImageGridAmp is the phase (mirrors
% view_eigenmodes). Off-support / undefined-phase vertices set to 0.
phase = R.Phase;
phase(~isfinite(phase)) = 0;
[sStudy, iStudy] = GetHostStudy(SurfaceFile);
if isempty(iStudy) || isempty(sStudy)
    error('view_connection_phase:study', 'Could not find a study to host the transient phase result.');
end
ResMat = db_template('resultsmat');
ResMat.ImageGridAmp  = [phase, phase];
ResMat.ImagingKernel = [];
ResMat.nComponents   = 1;
ResMat.Time          = [0 1];
ResMat.SurfaceFile   = SurfaceFile;
ResMat.HeadModelType = 'surface';
ResMat.nAvg          = 1;
ResMat.Leff          = 1;
ResMat.ColormapType  = 'connphase';
ResMat.Comment       = 'Connection phase (Fiedler, FS gauge)';
ResMat = bst_history('add', ResMat, 'connphase_view', 'Transient connection-phase viewer result');
StudyDir   = bst_fileparts(file_fullpath(sStudy.FileName));
OutputFile = bst_process('GetNewFilename', StudyDir, 'results_connphase');
bst_save(OutputFile, ResMat, 'v6');
db_add_data(iStudy, OutputFile, ResMat);

EnsureCyclicColormap();

hFig = view_surface_data(SurfaceFile, file_short(OutputFile));
if isempty(hFig)
    % Display failed: remove the transient result we just registered
    try
        file_delete(file_fullpath(OutputFile), 1);
        db_reload_studies(iStudy);
    catch
        % Non-fatal
    end
    error('view_connection_phase:figure', 'Could not open the surface figure.');
end
hAxes = GetAxes3D(hFig);
panel_surface('SetSurfaceSmooth', hFig, 1, 0, 0);
panel_surface('SetSurfaceEdges',  hFig, 1, 1);
panel_surface('SetDataThreshold', hFig, 1, 0);   % phase is angular: show all vertices (no amplitude threshold)
set(hFig, 'Name', ['Connection phase: ' SurfaceFile]);

% Auto-remove the transient result on close
set(hFig, 'DeleteFcn', @(h,e) CleanupResult(OutputFile, iStudy));

%% ===== STORE STATE =====
% Trim the eigenmode operators out of the cached state: the complex sparse
% ConnLaplacian + MassMatrix bloat figure memory and no consumer needs them
% (bst_conn_phase Recompute uses only Vectors/Component/CompRank).
ConnEigLite = ConnEig;
for f = {'ConnLaplacian','MassMatrix'}
    if isfield(ConnEigLite, f{1}), ConnEigLite = rmfield(ConnEigLite, f{1}); end
end
st = struct();
st.SurfaceFile = SurfaceFile;
st.ResultsFile = file_short(OutputFile);
st.ConnEig     = ConnEigLite;
st.vFrame      = vFrame;
st.FsFrame     = FsFrame;
st.R           = R;
st.Rank        = 1;
st.nModes      = max(ConnEig.CompRank(:));
st.MaxArrows   = MaxArrows;
st.Vtx         = Vtx;
st.Fcs         = Fcs;
st.showGlyphs  = true;
st.showSing    = true;
% Capture Brainstorm's key handler so ours can chain unhandled keys to it
% (preserving the figure_3d view shortcuts / rotation controls).
st.KeyPressFcn_bak = get(hFig, 'KeyPressFcn');
% Cache the dataset/result indices for live phase updates on mode change.
[st.iDS, st.iResult] = bst_memory('GetDataSetResult', file_short(OutputFile));
setappdata(hFig, 'ConnPhase', st);

DrawOverlays(hFig, hAxes, st);
set(hFig, 'KeyPressFcn', @(h,ev) KeyPress_Callback(h, ev));

bst_progress('stop');

catch ME
    bst_progress('stop');
    rethrow(ME);
end
end


%% ===== ENSURE A CYCLIC COLORMAP TYPE EXISTS =====
% Registers a DEDICATED 'connphase' colormap type so other source figures are not
% affected. bst_colormaps('GetColormap'/'SetColormap') error on unknown types
% (they index GlobalData.Colormaps.(type)), so the field is created directly here
% the same way bst_colormaps('Initialize') seeds the built-in types.
function EnsureCyclicColormap()
    global GlobalData;
    N = 256;
    CMap = hsv2rgb([linspace(0, 1, N)', ones(N, 1), ones(N, 1)]);
    % Seed the dedicated type from the source defaults if it does not exist yet.
    if isempty(GlobalData) || isempty(GlobalData.Colormaps) || ~isfield(GlobalData.Colormaps, 'connphase')
        sColormap = bst_colormaps('GetColormap', 'source');   % template struct (source always exists)
        GlobalData.Colormaps.connphase = sColormap;
    else
        sColormap = GlobalData.Colormaps.connphase;
    end
    sColormap.Name             = 'Connection phase (cyclic)';
    sColormap.CMap             = CMap;
    sColormap.isAbsoluteValues = 0;
    sColormap.DisplayColorbar  = 1;
    sColormap.MaxMode          = 'custom';
    sColormap.MinValue         = -pi;
    sColormap.MaxValue         = pi;
    bst_colormaps('SetColormap', 'connphase', sColormap);
end


%% ===== DRAW OVERLAYS (field glyphs + singularity markers) =====
function DrawOverlays(hFig, hAxes, st)
    if isempty(hAxes)
        hAxes = GetAxes3D(hFig);
    end
    delete(findobj(hAxes, '-depth', 1, '-regexp', 'Tag', '^connPhase(Field|Sing)'));
    Vtx = st.Vtx;  Fcs = st.Fcs;  R = st.R;

    % Hold the axes so successive plot calls (quiver3, plot3) do not replace
    % each other (Brainstorm surface axes default to NextPlot='replace').
    prevHold = get(hAxes, 'NextPlot');
    hold(hAxes, 'on');

    % --- Field glyphs (per-vertex, subsampled to MaxArrows on the support) ---
    if st.showGlyphs
        supp = find(any(R.Field ~= 0, 2));
        if ~isempty(supp)
            step = max(1, ceil(numel(supp) / st.MaxArrows));
            idx  = supp(1:step:end);
            meanEdge = MeanEdgeLength(Vtx, Fcs);
            len  = 0.8 * meanEdge;
            n    = st.vFrame.normals(idx, :);
            unitN = n ./ max(sqrt(sum(n.^2,2)), eps);
            w    = R.Field(idx, :);
            dirw = w ./ max(sqrt(sum(w.^2, 2)), eps);
            B    = Vtx(idx, :) + (0.1 * meanEdge) .* unitN;   % lift off the surface
            Vec  = dirw .* len;
            quiver3(B(:,1), B(:,2), B(:,3), Vec(:,1), Vec(:,2), Vec(:,3), 0, ...
                'Parent', hAxes, 'Color', [0.95 0.85 0.1], 'LineWidth', 1.0, ...
                'ShowArrowHead', 'off', 'Tag', 'connPhaseField');
        end
    end

    % --- Singularity markers ---
    if st.showSing && ~isempty(R.Singularities)
        P = Vtx(R.Singularities, :);
        plot3(P(:,1), P(:,2), P(:,3), 'o', 'Parent', hAxes, ...
            'MarkerFaceColor', [0.9 0.1 0.1], 'MarkerEdgeColor', [.2 .2 .2], ...
            'MarkerSize', 10, 'LineStyle', 'none', 'Tag', 'connPhaseSing');
    end

    % Restore original hold state
    set(hAxes, 'NextPlot', prevHold);
end


%% ===== KEYBOARD =====
% m / Shift+m : previous / next mode rank (recomputes phase + field)
% g           : toggle field glyphs
% p           : toggle singularity markers
% Left/Right  : fewer / more glyphs
% Any other key chains to Brainstorm's figure_3d handler (view shortcuts, etc.).
function KeyPress_Callback(hFig, ev)
    st = getappdata(hFig, 'ConnPhase');
    if isempty(st)
        return;
    end
    switch ev.Key
        case 'm'
            if ismember('shift', ev.Modifier)
                st.Rank = min(st.nModes, st.Rank + 1);
            else
                st.Rank = max(1, st.Rank - 1);
            end
            st = Recompute(hFig, st);
        case 'g'
            st.showGlyphs = ~st.showGlyphs;
        case 'p'
            st.showSing = ~st.showSing;
        case 'rightarrow'
            st.MaxArrows = ceil(st.MaxArrows * 1.5);
        case 'leftarrow'
            st.MaxArrows = max(50, floor(st.MaxArrows / 1.5));
        otherwise
            % Chain unhandled keys to Brainstorm's handler (preserves view shortcuts).
            if isfield(st, 'KeyPressFcn_bak') && ~isempty(st.KeyPressFcn_bak)
                st.KeyPressFcn_bak(hFig, ev);
            end
            return;
    end
    setappdata(hFig, 'ConnPhase', st);
    DrawOverlays(hFig, GetAxes3D(hFig), st);
end


%% ===== RECOMPUTE PHASE/FIELD FOR A NEW MODE RANK =====
function st = Recompute(hFig, st)
    global GlobalData;
    st.R = bst_conn_phase(st.ConnEig, st.vFrame, 'Rank', st.Rank, 'FsFrame', st.FsFrame, 'nSing', 2);
    phase = st.R.Phase;
    phase(~isfinite(phase)) = 0;
    GlobalData.DataSet(st.iDS).Results(st.iResult).ImageGridAmp = [phase, phase];
    panel_surface('UpdateSurfaceData', hFig);
end


%% ===== MEAN EDGE LENGTH =====
function L = MeanEdgeLength(Vtx, Fcs)
    e1 = Vtx(Fcs(:,2),:) - Vtx(Fcs(:,1),:);
    e2 = Vtx(Fcs(:,3),:) - Vtx(Fcs(:,2),:);
    e3 = Vtx(Fcs(:,1),:) - Vtx(Fcs(:,3),:);
    L = mean([sqrt(sum(e1.^2,2)); sqrt(sum(e2.^2,2)); sqrt(sum(e3.^2,2))]);
end


%% ===== MANIFOLD GAUGE FRAME (per vertex) =====
% NOTE: view_connection_phase is slated for rework/deprecation. It derives the
% per-vertex gauge frame ad hoc here; the canonical handler for manifold data is
% tess_manifold (find-or-load-or-create), which we go through for the node.
% Per-vertex trivial-gauge frame: U=real(grid.*rot), V=imag(grid.*rot),
% grid=Embedded.vertex.grid, rot=Gauge.vertex.rotation.
function [U, V] = GetManifoldFrame(SurfaceFile)
    M  = tess_manifold(SurfaceFile);    % canonical: returns the complete ManifoldMat
    nV = size(in_tess_bst(SurfaceFile).Vertices, 1);
    U  = zeros(nV, 3);  V = zeros(nV, 3);
    for hh = 1:2
        vH = double(M.Embedded(hh).GlobalVertices(:));
        cR = M.Embedded(hh).vertex.grid .* M.Gauge(hh).vertex.rotation(:);
        U(vH,:) = real(cR);
        V(vH,:) = imag(cR);
    end
end


%% ===== GET 3D AXES =====
% Returns the main 3D axes of a Brainstorm surface figure.
% Prefers the 'Axes3D'-tagged axes (present when the figure was created via
% figure_3d in a GUI session); falls back to the first non-Colorbar axes so
% that the viewer also works correctly in nogui / headless MATLAB sessions.
function hAxes = GetAxes3D(hFig)
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    if isempty(hAxes)
        allAx = findobj(hFig, 'Type', 'axes');
        for ia = 1:numel(allAx)
            if ~strcmpi(get(allAx(ia), 'Tag'), 'Colorbar')
                hAxes = allAx(ia);
                return;
            end
        end
    end
end


%% ===== FIND A STUDY TO HOST THE TRANSIENT RESULT =====
% Mirrors view_eigenmodes: host the transient result in the subject's
% intra-subject (analysis) study.
function [sStudy, iStudy] = GetHostStudy(SurfaceFile)
    sStudy = [];
    iStudy = [];
    [~, iSubject] = bst_get('SurfaceFile', SurfaceFile);
    if ~isempty(iSubject)
        [sStudy, iStudy] = bst_get('AnalysisIntraStudy', iSubject);
    end
    if isempty(iStudy) || isempty(sStudy)
        % Fall back to the current study
        [sStudy, iStudy] = bst_get('Study');
    end
end


%% ===== DELETE THE TRANSIENT RESULT ON FIGURE CLOSE =====
function CleanupResult(OutputFile, iStudy)
    try
        file_delete(file_fullpath(OutputFile), 1);
        db_reload_studies(iStudy);
    catch
        % Non-fatal: leave the node if cleanup fails (user can delete it)
    end
end
