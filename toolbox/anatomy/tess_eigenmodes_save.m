function TessMat = tess_eigenmodes_save(SurfaceFile, Eigenmodes, Vertices, Faces, isInteractive)
% TESS_EIGENMODES_SAVE: Save precomputed eigenmodes into a Brainstorm surface file.
%
% USAGE:  TessMat = tess_eigenmodes_save(SurfaceFile, Eigenmodes, Vertices, Faces)
%         TessMat = tess_eigenmodes_save(SurfaceFile, Eigenmodes, Vertices, Faces, isInteractive)
%
% DESCRIPTION:
%     Saves precomputed eigenmodes into a Brainstorm surface file (TessMat).
%     The eigenmodes are stored as a field in the .mat file alongside
%     Vertices, Faces, Atlas, etc. This avoids recomputation and allows
%     downstream processes to access eigenmodes directly from the surface.
%
%     If the mesh was repaired during eigenmode computation (non-manifold
%     defects fixed), the repaired Vertices and Faces replace the original
%     ones in the surface file. VertConn, VertNormals, and Curvature are
%     recomputed to match.
%
% INPUTS:
%     SurfaceFile   : Relative path to the surface file in the Brainstorm
%                     database (e.g., 'sub-0002/tess_cortex_pial_low.mat')
%     Eigenmodes    : Structure from tess_eigenmodes with fields:
%                     .Vectors, .Values, .nModes, .MassType, etc.
%     Vertices      : [nV x 3] vertex positions (from tess_eigenmodes output,
%                     potentially repaired)
%     Faces         : [nF x 3] face indices (from tess_eigenmodes output,
%                     potentially repaired)
%     isInteractive : (optional) If true, print messages. Default: true.
%
% OUTPUTS:
%     TessMat : Updated surface structure with Eigenmodes field.
%
% SEE ALSO: tess_eigenmodes, tess_eigenmodes_load

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

if nargin < 5
    isInteractive = true;
end

%% ===== VALIDATE INPUTS =====
if ~isstruct(Eigenmodes) || ~isfield(Eigenmodes, 'Vectors') || ~isfield(Eigenmodes, 'Values')
    error('Eigenmodes must be a structure with Vectors and Values fields.');
end
if size(Eigenmodes.Vectors, 1) ~= size(Vertices, 1)
    error('Eigenmodes.Vectors (%d rows) must match Vertices (%d rows).', ...
        size(Eigenmodes.Vectors, 1), size(Vertices, 1));
end

%% ===== LOAD SURFACE FILE =====
% Resolve full path
if ~isempty(strfind(SurfaceFile, filesep)) && exist(SurfaceFile, 'file') %#ok<STREMP>
    SurfaceFileFull = SurfaceFile;
else
    SurfaceFileFull = file_fullpath(SurfaceFile);
end

if ~exist(SurfaceFileFull, 'file')
    error('Surface file not found: %s', SurfaceFileFull);
end

% Load the surface
TessMat = load(SurfaceFileFull);

%% ===== UPDATE MESH IF REPAIRED =====
nV_orig = size(TessMat.Vertices, 1);
nV_new = size(Vertices, 1);

if nV_new ~= nV_orig || size(Faces, 1) ~= size(TessMat.Faces, 1)
    if isInteractive
        fprintf('BST> Mesh was repaired (%d->%d V, %d->%d F). Updating surface.\n', ...
            nV_orig, nV_new, size(TessMat.Faces, 1), size(Faces, 1));
    end
    TessMat.Vertices = Vertices;
    TessMat.Faces = Faces;
    % Recompute derived fields
    TessMat.VertConn = tess_vertconn(Vertices, Faces);
    TessMat.VertNormals = tess_normals(Vertices, Faces, TessMat.VertConn);
    % Clear fields that depend on vertex count (will be recomputed on demand)
    TessMat.Curvature = [];
    TessMat.SulciMap = [];
    TessMat.tess2mri_interp = [];
    % Note: Atlas scout vertices may reference removed vertices.
    % For now, clear atlas. A more sophisticated approach would remap.
    if isfield(TessMat, 'Atlas') && ~isempty(TessMat.Atlas)
        warning('BST:AtlasCleared', ...
            'Atlas cleared because mesh vertices changed during repair.');
        TessMat.Atlas = [];
    end
end

%% ===== STORE EIGENMODES =====
% Store as single precision to save space (300 modes x 15k vertices = 18 MB
% in single vs 36 MB in double). The eigenvalues stay double.
EigenmodesStore = struct();
EigenmodesStore.Vectors     = single(Eigenmodes.Vectors);
EigenmodesStore.Values      = Eigenmodes.Values(:);  % column vector, double
EigenmodesStore.nModes      = Eigenmodes.nModes;
EigenmodesStore.MassType    = Eigenmodes.MassType;
EigenmodesStore.Sigma       = Eigenmodes.Sigma;
EigenmodesStore.Tolerance   = Eigenmodes.Tolerance;
EigenmodesStore.nRemoved    = Eigenmodes.nRemoved;
EigenmodesStore.ComputeTime = Eigenmodes.ComputeTime;
EigenmodesStore.ComputeDate = datestr(now, 'yyyy-mm-dd HH:MM:SS');

TessMat.Eigenmodes = EigenmodesStore;

%% ===== SAVE =====
bst_save(SurfaceFileFull, TessMat, 'v7');

% Add history entry
TessMat = bst_history('add', TessMat, 'eigenmodes', ...
    sprintf('Computed %d Laplace-Beltrami eigenmodes (%s mass, sigma=%.1e)', ...
    Eigenmodes.nModes, Eigenmodes.MassType, Eigenmodes.Sigma));
bst_save(SurfaceFileFull, TessMat, 'v7');

if isInteractive
    fprintf('BST> Saved %d eigenmodes to: %s\n', Eigenmodes.nModes, SurfaceFile);
end

end
