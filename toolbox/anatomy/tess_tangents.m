function [U, V] = tess_tangents(SurfaceFile, varargin)
% TESS_TANGENTS: Globally consistent per-face tangent frame field.
%
% DEPRECATED: superseded by the manifold frame (view_manifold / tess_frame).
%     Retained only for the per-face callers (bst_wavefront_track,
%     tess_nxr_populate) until per-face manifold frames are available
%     (nxr Gauge.face.rotation is currently deferred for the trivial gauge).
%
% USAGE:  [U, V] = tess_tangents(SurfaceFile)              % compute, store, return
%         [U, V] = tess_tangents(SurfaceFile, 'NoSave', 1) % compute + return only
%
% DESCRIPTION:
%     Computes a globally consistent tangent frame field on a cortical surface
%     using nxr-compute's trivial connection. Two singularities are placed at the
%     north and south poles of each hemisphere's FreeSurfer registration sphere
%     (Reg.Sphere). Because the two hemispheres are disconnected genus-0 spheres,
%     the field is solved PER HEMISPHERE (each: Euler characteristic 2, two +1
%     singularities). The result is a per-FACE orthonormal tangent frame (U, V).
%
%     The frame is stored on the surface file as TessMat.TangentFrame (Domain
%     'face'). Transferring it to a per-vertex frame (via the Hodge star / DEC)
%     is future work; the storage format reserves a 'Domain' field for it.
%
% INPUT:
%     - SurfaceFile : Brainstorm cortex surface with a FreeSurfer registration
%                     sphere (TessMat.Reg.Sphere.Vertices).
% OPTIONS:
%     - NoSave : (logical) If true, do not write TangentFrame back to the file.
%                Default: false (store).
% OUTPUT:
%     - U : [nF x 3] per-face first tangent direction (e1).
%     - V : [nF x 3] per-face orthogonal tangent (e2); cross(U,V) is aligned to
%           the face normal (right-handed frame).
%
% Requires the nxr-compute plugin (no MATLAB fallback for the trivial connection).
%
% SEE ALSO: tess_normals, tess_hemisplit, tess_addsphere

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
NoSave = false;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'nosave', NoSave = varargin{i+1};
    end
end

%% ===== LOAD SURFACE =====
TessFile = file_fullpath(SurfaceFile);
TessMat  = in_tess_bst(SurfaceFile);
Vtx = TessMat.Vertices;
Fcs = double(TessMat.Faces);
nF  = size(Fcs, 1);

%% ===== REQUIRE FREESURFER REGISTRATION SPHERE =====
if ~isfield(TessMat, 'Reg') || ~isstruct(TessMat.Reg) || ...
   ~isfield(TessMat.Reg, 'Sphere') || ~isfield(TessMat.Reg.Sphere, 'Vertices') || ...
   isempty(TessMat.Reg.Sphere.Vertices)
    error('tess_tangents:noRegSphere', ...
        ['Surface has no FreeSurfer registration sphere (Reg.Sphere.Vertices). ' ...
         'Import with surface registration or run tess_addsphere first.']);
end
Sphere = TessMat.Reg.Sphere.Vertices;

%% ===== REQUIRE IMPORT-TIME HEMISPHERE LABELS (no geometric re-split) =====
% Use the left/right labels recorded at import (the 'Structures' atlas). Do NOT
% re-split the mesh via connected components or coordinates: the cortex is
% already split and labeled at import time, so we only READ the labels here.
hasLabels = false;
if isfield(TessMat, 'Atlas') && ~isempty(TessMat.Atlas)
    iStruct = find(strcmpi({TessMat.Atlas.Name}, 'Structures'), 1);
    if ~isempty(iStruct) && ~isempty(TessMat.Atlas(iStruct).Scouts)
        scouts  = TessMat.Atlas(iStruct).Scouts;
        labels  = {scouts.Label};
        regions = {scouts.Region};
        reg1    = cellfun(@(c) c(1), regions(~cellfun(@isempty, regions)), 'UniformOutput', false);
        hasL = any(strcmpi(labels, 'lh')) || any(strcmpi(reg1, 'L'));
        hasR = any(strcmpi(labels, 'rh')) || any(strcmpi(reg1, 'R'));
        hasLabels = hasL && hasR;
    end
end
if ~hasLabels
    error('tess_tangents:noHemisphereLabels', ...
        ['Surface has no ''Structures'' atlas with left/right hemisphere labels. ' ...
         'tess_tangents uses the import-time hemisphere labels and never re-splits the ' ...
         'mesh geometrically (no connected-components / coordinate split). Re-import the ' ...
         'FreeSurfer surface so the ''Structures'' atlas is present.']);
end

%% ===== ENSURE NXR-COMPUTE PLUGIN =====
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
if ~isOk
    error('tess_tangents:nxrUnavailable', ...
        'tess_tangents requires the nxr-compute plugin: %s', errMsg);
end

%% ===== HEMISPHERE SPLIT (from import labels, not geometry) =====
% Labels are guaranteed present (checked above), so tess_hemisplit reads the
% 'Structures' atlas and does not use its geometric region-grow / y-split fallback.
[rH, lH, isConnected] = tess_hemisplit(TessMat);
if isConnected
    error('tess_tangents:connectedHemispheres', ...
        ['The two hemispheres are connected (or could not be cleanly separated). ' ...
         'tess_tangents solves each hemisphere as an independent genus-0 sphere and ' ...
         'requires disconnected hemispheres. Use a surface with separated hemispheres.']);
end
hemis    = {lH(:)', rH(:)'};
hemiTags = {'L', 'R'};

%% ===== FACE NORMALS (for frame orientation) =====
fn = cross(Vtx(Fcs(:,2),:) - Vtx(Fcs(:,1),:), Vtx(Fcs(:,3),:) - Vtx(Fcs(:,1),:));
fn = fn ./ max(sqrt(sum(fn.^2, 2)), eps);

%% ===== PER-HEMISPHERE TRIVIAL CONNECTION =====
U = zeros(nF, 3);
V = zeros(nF, 3);
assigned  = false(nF, 1);
singVerts = [];
singIdx   = [];
hemiOf    = {};
for h = 1:numel(hemis)
    vH = hemis{h};
    if isempty(vH)
        continue;
    end
    [Uh, Vh, fMask, polesGlobal] = solve_hemisphere(Vtx, Fcs, Sphere, vH);
    if any(assigned & fMask)
        error('tess_tangents:overlap', 'A face was assigned to both hemispheres.');
    end
    U(fMask, :)     = Uh;
    V(fMask, :)     = Vh;
    assigned(fMask) = true;
    singVerts = [singVerts; polesGlobal(:)];          %#ok<AGROW>
    singIdx   = [singIdx;   [1; 1]];                  %#ok<AGROW>
    hemiOf    = [hemiOf, hemiTags(h), hemiTags(h)];   %#ok<AGROW>
end

%% ===== COVERAGE CHECK =====
if ~all(assigned)
    error('tess_tangents:unassignedFaces', ...
        '%d of %d faces were not assigned to a hemisphere.', nnz(~assigned), nF);
end

%% ===== ORIENT FRAME (right-handed w.r.t. face normal) =====
flip = sum(cross(U, V, 2) .* fn, 2) < 0;
V(flip, :) = -V(flip, :);

%% ===== STORE =====
if ~NoSave
    Sing = struct();
    Sing.Vertices   = singVerts;
    Sing.Indices    = singIdx;
    Sing.Hemisphere = hemiOf;
    TF = struct();
    TF.Domain        = 'face';
    TF.U             = single(U);
    TF.V             = single(V);
    TF.Singularities = Sing;
    TF.Method        = 'nxr trivial-connection (FreeSurfer poles)';
    TessMat.TangentFrame = TF;
    TessMat = bst_history('add', TessMat, 'tangents', ...
        'Computed per-face tangent frame field (nxr trivial connection, FreeSurfer poles).');
    bst_save(TessFile, TessMat, 'v7');
end
end


%% ========================================================================
function [Uh, Vh, fMask, polesGlobal] = solve_hemisphere(Vtx, Fcs, Sphere, vH)
% Solve the trivial-connection frame on one hemisphere submesh.
%   Returns the per-face direction (Uh) and orthogonal (Vh) vectors for the
%   faces selected by fMask (into the full face array), plus the two pole
%   vertex indices (global) used as singularities.
nVtot = size(Vtx, 1);
isVH  = false(nVtot, 1);
isVH(vH) = true;
% Faces wholly inside this hemisphere (valid: hemispheres are disconnected).
fMask = all(isVH(Fcs), 2);
% Local re-indexing (global vertex -> 1-based local index).
map = zeros(nVtot, 1);
map(vH) = 1:numel(vH);
Floc = map(Fcs(fMask, :));
Vloc = Vtx(vH, :);
% Poles from the registration sphere restricted to this hemisphere.
sph = Sphere(vH, :);
[~, iN] = max(sph(:, 3));   % north pole (local index)
[~, iS] = min(sph(:, 3));   % south pole (local index)
% Trivial connection: two +1 singularities, sum = chi = 2.
mctx = nxr.manifold.context(Vloc, Floc);
r = nxr.manifold.interpolate.trivial(mctx, [iN; iS], [1; 1]);
if ~r.gaussBonnetSatisfied
    error('tess_tangents:gaussBonnet', ...
        'Gauss-Bonnet not satisfied for hemisphere (chi=%.3f, expected 2).', ...
        r.eulerCharacteristic);
end
Uh = r.directionVectors;
Vh = r.orthogonalVectors;
polesGlobal = vH([iN, iS]);
end
