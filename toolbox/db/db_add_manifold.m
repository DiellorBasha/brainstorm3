function iManifold = db_add_manifold(iSubject, ParentSurfaceFile, ManifoldMat, Comment)
% DB_ADD_MANIFOLD: Save a manifold_*.mat under a parent surface and register it as a child node.
%
% USAGE:  iManifold = db_add_manifold(iSubject, ParentSurfaceFile, ManifoldMat, Comment)
%
% INPUT:
%    - iSubject          : Index of the subject (0 = default subject)
%    - ParentSurfaceFile : Relative or full path to the parent surface file
%    - ManifoldMat       : Structure to save (db_template('manifoldmat') schema)
%    - Comment           : Description string for this manifold node (default: 'Manifold')
% OUTPUT:
%    - iManifold : Index of the new manifold entry in sSubject.Surface(iSurface).Manifold

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

% Default comment
if (nargin < 4) || isempty(Comment)
    Comment = 'Manifold';
end

% Get protocol subjects database
ProtocolSubjects = bst_get('ProtocolSubjects');

% Resolve subject (handle default subject iSubject==0)
if (iSubject == 0)
    sSubject = ProtocolSubjects.DefaultSubject;
else
    sSubject = ProtocolSubjects.Subject(iSubject);
end

% Find the parent surface in the subject's surface list
iSurface = find(strcmpi({sSubject.Surface.FileName}, file_short(ParentSurfaceFile)), 1);
if isempty(iSurface)
    error('db_add_manifold:noSurface', 'Parent surface not found in subject: %s', ParentSurfaceFile);
end

% Build a unique manifold filename in the same anatomy folder as the parent surface
ProtocolInfo = bst_get('ProtocolInfo');
c = clock;
strTime = sprintf('%02.0f%02.0f%02.0f_%02.0f%02.0f', c(1)-2000, c(2:5));
OutputFile = ['manifold_' strTime '.mat'];
OutputFileFull = file_unique(bst_fullfile(ProtocolInfo.SUBJECTS, bst_fileparts(sSubject.FileName), OutputFile));

% Set required fields on the manifold structure
ManifoldMat.ParentSurface = file_short(ParentSurfaceFile);
ManifoldMat.Comment = Comment;

% Replace (not duplicate) any existing manifold of the same gauge. The Comment
% encodes the gauge (e.g. 'Manifold (gauge=trivial)'), so a recompute of the same
% gauge supersedes the old node: delete its file and drop its DB entry before
% appending the new one. Distinct-gauge manifolds (different Comment) coexist.
% Without this, every ForceRecompute leaves a stale node as Manifold(1), which the
% tess_manifold cache path would then load in preference to the fresh result.
if ~isempty(sSubject.Surface(iSurface).Manifold)
    oldComments = {sSubject.Surface(iSurface).Manifold.Comment};
    iSupersede  = find(strcmpi(oldComments, Comment));
    for k = iSupersede(:)'
        oldFile = file_fullpath(sSubject.Surface(iSurface).Manifold(k).FileName);
        if exist(oldFile, 'file')
            file_delete(oldFile, 1);
        end
    end
    sSubject.Surface(iSurface).Manifold(iSupersede) = [];
end

% Save to disk
bst_save(OutputFileFull, ManifoldMat, 'v7');

% Normalize all surfaces against db_template('Surface') so the struct array stays
% homogeneous (required by db_surface_sort and MATLAB's struct-array assignment
% rules). db_template is the single source of truth for the Manifold/Eigen/Operator
% empty schemas, so they cannot drift. Existing protocols are normally backfilled at
% load time by db_update (DB v5.04); this remains correct (and a near no-op) even if a
% non-normalized surface is passed. A fresh array is rebuilt rather than assigned
% element-by-element, because widening the field set of one element of an existing
% struct array in place throws "Subscripted assignment between dissimilar structures".
templateSurface = db_template('Surface');
tFields = fieldnames(templateSurface);
if ~isempty(sSubject.Surface) && ~all(isfield(sSubject.Surface, tFields))
    sNew = cell(1, numel(sSubject.Surface));
    for k = 1:numel(sSubject.Surface)
        s = struct_copy_fields(sSubject.Surface(k), templateSurface, 0);
        extraFields = setdiff(fieldnames(s), tFields, 'stable');
        sNew{k} = orderfields(s, [tFields; extraFields]);
    end
    sSubject.Surface = reshape([sNew{:}], size(sSubject.Surface));
end

% Append a child entry to the parent surface's Manifold list.
% Cache the discriminating spec (Gauge) so bst_get can find this node from the cache alone.
newEntry.FileName = file_short(OutputFileFull);
newEntry.Comment  = Comment;
newEntry.Gauge = '';
if isfield(ManifoldMat, 'Provenance') && isstruct(ManifoldMat.Provenance) && isfield(ManifoldMat.Provenance, 'Gauge')
    newEntry.Gauge = ManifoldMat.Provenance.Gauge;
end
sSubject.Surface(iSurface).Manifold(end+1) = newEntry;
iManifold = numel(sSubject.Surface(iSurface).Manifold);

% Write the modified subject back via the canonical round-trip
if (iSubject == 0)
    ProtocolSubjects.DefaultSubject = sSubject;
else
    ProtocolSubjects.Subject(iSubject) = sSubject;
end
bst_set('ProtocolSubjects', ProtocolSubjects);

% Refresh tree and save database
panel_protocols('UpdateNode', 'Subject', iSubject);
db_save();

end
