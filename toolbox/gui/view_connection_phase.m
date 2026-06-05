function hFig = view_connection_phase(SurfaceFile, varargin)
% VIEW_CONNECTION_PHASE: Connection-Laplacian eigenmode viewer (phase + field).
%
% USAGE:  hFig = view_connection_phase(SurfaceFile)
%         hFig = view_connection_phase(SurfaceFile, 'nModes', 50, 'MaxArrows', 2000)
%
% DESCRIPTION:
%     Opens an opaque cortex figure showing a connection-Laplacian eigenmode's
%     FreeSurfer-gauge phase as a cyclic surface colormap (the between-subject
%     location coordinate). Field glyphs and singularity markers are added by
%     later tasks. Data pipeline: bst_conn_eigenmodes_ensure -> nxr vertexFrame ->
%     bst_tangent_face2vertex (FS frame) -> bst_conn_phase.
%
%     Requires the nxr-compute plugin and a FreeSurfer-registered cortex.
%
% SEE ALSO: bst_conn_phase, bst_conn_eigenmodes_ensure, tess_tangents, view_eigenmodes

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
Nv  = TessMat.VertNormals;

mctx   = nxr.manifold.context(Vtx, Fcs);
vFrame = nxr.manifold.measure.vertexFrame(mctx);
[Uf, ~]  = tess_tangents(SurfaceFile, 'NoSave', 1);
[Uv, Vv] = bst_tangent_face2vertex(Fcs, Uf, Nv);
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
hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D'); %#ok<NASGU>
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
setappdata(hFig, 'ConnPhase', st);

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
