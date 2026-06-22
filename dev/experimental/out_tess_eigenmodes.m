function TessMat = out_tess_eigenmodes(SurfaceFile, Eigenmodes, Vertices, Faces, isInteractive)
% OUT_TESS_EIGENMODES: Save precomputed eigenmodes into a Brainstorm surface file.
%
% USAGE:  TessMat = out_tess_eigenmodes(SurfaceFile, Eigenmodes, Vertices, Faces)
%         TessMat = out_tess_eigenmodes(SurfaceFile, Eigenmodes, Vertices, Faces, isInteractive)
%
% SEE ALSO: in_tess_eigenmodes, tess_eigenmodes

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
SurfaceFileFull = file_fullpath(SurfaceFile);
if ~file_exist(SurfaceFileFull)
    error('Surface file not found: %s', SurfaceFileFull);
end
TessMat = load(SurfaceFileFull);

%% ===== UPDATE MESH IF REPAIRED =====
if (size(Vertices, 1) ~= size(TessMat.Vertices, 1)) || (size(Faces, 1) ~= size(TessMat.Faces, 1))
    if isInteractive
        fprintf('BST> Mesh was repaired (%d->%d V, %d->%d F). Updating surface.\n', ...
            size(TessMat.Vertices, 1), size(Vertices, 1), size(TessMat.Faces, 1), size(Faces, 1));
    end
    TessMat.Vertices    = Vertices;
    TessMat.Faces       = Faces;
    TessMat.VertConn    = tess_vertconn(Vertices, Faces);
    TessMat.VertNormals = tess_normals(Vertices, Faces, TessMat.VertConn);
    TessMat.Curvature   = [];
    TessMat.SulciMap    = [];
    TessMat.tess2mri_interp = [];
    if isfield(TessMat, 'Atlas') && ~isempty(TessMat.Atlas)
        warning('BST:AtlasCleared', 'Atlas cleared because mesh vertices changed during repair.');
        TessMat.Atlas = [];
    end
end

%% ===== STORE EIGENMODES (single precision vectors) =====
EigenmodesStore = struct();
EigenmodesStore.Vectors     = single(Eigenmodes.Vectors);
EigenmodesStore.Values      = Eigenmodes.Values(:);
EigenmodesStore.nModes      = Eigenmodes.nModes;
if isfield(Eigenmodes, 'Component'),   EigenmodesStore.Component   = Eigenmodes.Component(:);   end
if isfield(Eigenmodes, 'CompRank'),    EigenmodesStore.CompRank    = Eigenmodes.CompRank(:);    end
if isfield(Eigenmodes, 'nComponents'), EigenmodesStore.nComponents = Eigenmodes.nComponents;    end
EigenmodesStore.MassType    = Eigenmodes.MassType;
% Persist the whole-mesh operators (and their types) for reuse (projection, filtering,
% heat kernels, ...). Sparse matrices stay double (no sparse-single); negligible vs Vectors.
if isfield(Eigenmodes, 'MassMatrix')    && ~isempty(Eigenmodes.MassMatrix),    EigenmodesStore.MassMatrix    = Eigenmodes.MassMatrix;    end
if isfield(Eigenmodes, 'Laplacian')     && ~isempty(Eigenmodes.Laplacian),     EigenmodesStore.Laplacian     = Eigenmodes.Laplacian;     end
if isfield(Eigenmodes, 'LaplacianType') && ~isempty(Eigenmodes.LaplacianType), EigenmodesStore.LaplacianType = Eigenmodes.LaplacianType; end
EigenmodesStore.Sigma       = Eigenmodes.Sigma;
EigenmodesStore.Tolerance   = Eigenmodes.Tolerance;
EigenmodesStore.nRemoved    = Eigenmodes.nRemoved;
EigenmodesStore.ComputeTime = Eigenmodes.ComputeTime;
EigenmodesStore.ComputeDate = datestr(now, 'yyyy-mm-dd HH:MM:SS');
TessMat.Eigenmodes = EigenmodesStore;

%% ===== HISTORY + SINGLE SAVE =====
TessMat = bst_history('add', TessMat, 'eigenmodes', ...
    sprintf('Computed %d Laplace-Beltrami eigenmodes (%s mass, sigma=%.1e)', ...
        Eigenmodes.nModes, Eigenmodes.MassType, Eigenmodes.Sigma));
bst_save(SurfaceFileFull, TessMat, 'v7');

if isInteractive
    fprintf('BST> Saved %d eigenmodes to: %s\n', Eigenmodes.nModes, SurfaceFile);
end
end
