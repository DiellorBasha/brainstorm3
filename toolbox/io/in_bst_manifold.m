function ManifoldMat = in_bst_manifold(ManifoldFile, varargin)
% IN_BST_MANIFOLD: Read a manifold_*.mat node (nxr geometry backbone) in Brainstorm format.
%
% USAGE:  ManifoldMat = in_bst_manifold(ManifoldFile, FieldsList) : Read the specified fields
%         ManifoldMat = in_bst_manifold(ManifoldFile)             : Read all the fields
%
% INPUT:
%    - ManifoldFile : Absolute or relative path to the manifold_*.mat file to read
%    - FieldsList   : List of field names to read from the file (default: all)
%
% DESCRIPTION:
%     Loader (the "in_bst" side) for the manifold_ derived-anatomy node produced by
%     tess_manifold / db_add_manifold. It only reads the file off disk; it neither
%     resolves the DB cache (that is bst_get('ManifoldFile', ...)) nor registers
%     anything (that is db_add_manifold). The geometry backbone is stored alongside its
%     parent surface in the anatomy tree; relative paths are resolved against the
%     protocol SUBJECTS folder.
%
%     NOTE: tess_manifold is the canonical find-or-load-or-create entry point and should
%     be preferred when a manifold may need to be computed. Use in_bst_manifold when the
%     node is known to exist and only a plain read is wanted.
%
% OUTPUT:
%    - ManifoldMat : struct matching db_template('manifoldmat'): Comment, ParentSurface,
%                    Topology(1x2), Embedded(1x2), Intrinsic(1x2), Extrinsic(1x2),
%                    Gauge(1x2), Provenance.
%
% SEE ALSO: tess_manifold, db_add_manifold, in_bst_eigen, in_bst_operator, db_template

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

%% ===== RESOLVE FILENAME (relative -> absolute) =====
ManifoldFile = file_fullpath(ManifoldFile);
if ~file_exist(ManifoldFile)
    error(['Manifold file not found:' 10 file_short(ManifoldFile) 10 ...
           'You should reload this protocol (right-click > reload).']);
end

%% ===== LOAD =====
if (nargin < 2)
    % Read all fields
    ManifoldMat = load(ManifoldFile);
    % Type-signature guard (only meaningful on a full read)
    if ~isfield(ManifoldMat, 'Topology')
        error('Not a valid manifold node (no Topology field): %s', file_short(ManifoldFile));
    end
else
    % Read the requested fields only
    FieldsToRead = varargin;
    warning off MATLAB:load:variableNotFound
    ManifoldMat = load(ManifoldFile, FieldsToRead{:});
    warning on MATLAB:load:variableNotFound
    % Ensure every requested field exists (empty if absent on disk)
    for i = 1:numel(FieldsToRead)
        if ~isfield(ManifoldMat, FieldsToRead{i})
            ManifoldMat.(FieldsToRead{i}) = [];
        end
    end
end
