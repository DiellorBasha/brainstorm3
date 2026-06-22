function OperatorMat = in_bst_operator(OperatorFile, varargin)
% IN_BST_OPERATOR: Read an operator_*.mat node (per-hemisphere discrete operators) in Brainstorm format.
%
% USAGE:  OperatorMat = in_bst_operator(OperatorFile, FieldsList) : Read the specified fields
%         OperatorMat = in_bst_operator(OperatorFile)             : Read all the fields
%
% INPUT:
%    - OperatorFile : Absolute or relative path to the operator_*.mat file to read
%    - FieldsList   : List of field names to read from the file (default: all)
%
% DESCRIPTION:
%     Loader (the "in_bst" side) for the operator_ derived-anatomy node produced by
%     tess_operators / db_add_operator. It only reads the file off disk; it neither
%     resolves the DB cache (that is bst_get('OperatorFile', ...)) nor registers
%     anything (that is db_add_operator). The operator pencil (A, B) is stored
%     alongside its parent surface in the anatomy tree; relative paths are resolved
%     against the protocol SUBJECTS folder.
%
% OUTPUT:
%    - OperatorMat : struct matching db_template('operatormat'): Comment, ParentSurface,
%                    Variant, Operator(1x2), Mass(1x2), GlobalVertices(1x2),
%                    GlobalFaces(1x2), and (Dirac/Hodge-Face variants) FirstOrder /
%                    FaceMass / FaceAux, Provenance.
%
% SEE ALSO: tess_operators, db_add_operator, in_bst_eigen, in_bst_manifold, db_template

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
OperatorFile = file_fullpath(OperatorFile);
if ~file_exist(OperatorFile)
    error(['Operator file not found:' 10 file_short(OperatorFile) 10 ...
           'You should reload this protocol (right-click > reload).']);
end

%% ===== LOAD =====
if (nargin < 2)
    % Read all fields
    OperatorMat = load(OperatorFile);
    % Type-signature guard (only meaningful on a full read)
    if ~isfield(OperatorMat, 'Operator') || ~isfield(OperatorMat, 'Mass')
        error('Not a valid operator node (no Operator/Mass field): %s', file_short(OperatorFile));
    end
else
    % Read the requested fields only
    FieldsToRead = varargin;
    warning off MATLAB:load:variableNotFound
    OperatorMat = load(OperatorFile, FieldsToRead{:});
    warning on MATLAB:load:variableNotFound
    % Ensure every requested field exists (empty if absent on disk)
    for i = 1:numel(FieldsToRead)
        if ~isfield(OperatorMat, FieldsToRead{i})
            OperatorMat.(FieldsToRead{i}) = [];
        end
    end
end
