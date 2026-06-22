function EigenMat = in_bst_eigen(EigenFile, varargin)
% IN_BST_EIGEN: Read an eigen_*.mat node (per-hemisphere eigenbasis) in Brainstorm format.
%
% USAGE:  EigenMat = in_bst_eigen(EigenFile, FieldsList) : Read the specified fields
%         EigenMat = in_bst_eigen(EigenFile)             : Read all the fields
%
% INPUT:
%    - EigenFile  : Absolute or relative path to the eigen_*.mat file to read
%    - FieldsList : List of field names to read from the file (default: all)
%
% DESCRIPTION:
%     Loader (the "in_bst" side) for the eigen_ derived-anatomy node produced by
%     tess_eigen / db_add_eigen. It only reads the file off disk; it neither resolves
%     the DB cache (that is bst_get('EigenFile', ...)) nor registers anything (that is
%     db_add_eigen). The eigenbasis is stored alongside its parent surface in the
%     anatomy tree; relative paths are resolved against the protocol SUBJECTS folder.
%
%     The referenced operator node (EigenMat.OperatorFile), which carries the mass B,
%     is NOT auto-loaded: call in_bst_operator(EigenMat.OperatorFile) when the mass is
%     needed (single-responsibility, one loader per node type).
%
% OUTPUT:
%    - EigenMat : struct matching db_template('eigenmat'): Comment, ParentSurface,
%                 OperatorFile, Variant, Phi(1x2), Lambda(1x2), K, GlobalVertices(1x2),
%                 GlobalFaces(1x2), Provenance.
%
% SEE ALSO: tess_eigen, db_add_eigen, in_bst_operator, in_bst_manifold, db_template

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
EigenFile = file_fullpath(EigenFile);
if ~file_exist(EigenFile)
    error(['Eigen file not found:' 10 file_short(EigenFile) 10 ...
           'You should reload this protocol (right-click > reload).']);
end

%% ===== LOAD =====
if (nargin < 2)
    % Read all fields
    EigenMat = load(EigenFile);
    % Type-signature guard (only meaningful on a full read)
    if ~isfield(EigenMat, 'Phi')
        error('Not a valid eigen node (no Phi field): %s', file_short(EigenFile));
    end
else
    % Read the requested fields only. If nModes is requested, also read the
    % legacy field K so the shim below can map it (mirrors in_bst_results, which
    % reads companion fields when one is requested).
    FieldsToRead = varargin;
    LoadFields   = FieldsToRead;
    if ismember('nModes', FieldsToRead) && ~ismember('K', LoadFields)
        LoadFields{end+1} = 'K';
    end
    warning off MATLAB:load:variableNotFound
    EigenMat = load(EigenFile, LoadFields{:});
    warning on MATLAB:load:variableNotFound
    % Ensure every requested field exists (empty if absent on disk)
    for i = 1:numel(FieldsToRead)
        if ~isfield(EigenMat, FieldsToRead{i})
            EigenMat.(FieldsToRead{i}) = [];
        end
    end
end

% Backward-compat shim: legacy eigen_ files store the per-hemisphere mode count
% as 'K'; the canonical field is now 'nModes'. Map it on read so old files load
% unchanged (no file rewrite needed).
if (~isfield(EigenMat,'nModes') || isempty(EigenMat.nModes)) && isfield(EigenMat,'K') && ~isempty(EigenMat.K)
    EigenMat.nModes = EigenMat.K;
end
