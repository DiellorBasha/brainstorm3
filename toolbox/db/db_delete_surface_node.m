function iModifiedSubjects = db_delete_surface_node(FileNames, isForced)
% DB_DELETE_SURFACE_NODE: Delete derived-anatomy nodes nested under a surface.
%
% USAGE:  iModifiedSubjects = db_delete_surface_node(FileNames, isForced)
%
% DESCRIPTION:
%     Deletes manifold_/operator_/eigen_ derived-anatomy nodes: removes each node file
%     from disk AND its entry from the DB cache (the parent surface's Manifold/Operator/
%     Eigen list). Each node is resolved by filename via bst_get (file_gettype picks the
%     type), so the caller only needs the relative/absolute file path(s).
%
%     This is the MATLAB-side deletion primitive shared by node_delete (tree "Delete")
%     and the interactive overwrite flow in tess_eigen/tess_operators/tess_manifold.
%     The caller is responsible for refreshing the tree (e.g.
%     panel_protocols('UpdateNode', 'Subject', iModifiedSubjects)).
%
% INPUT:
%    - FileNames : char or cell array of derived-node file paths (relative or absolute)
%    - isForced  : if 1, delete the files without the confirmation dialog (default 0)
% OUTPUT:
%    - iModifiedSubjects : unique subject indices whose Surface lists changed (or [])
%
% SEE ALSO: node_delete, db_add_eigen, db_add_operator, db_add_manifold, file_gettype

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

    if (nargin < 2) || isempty(isForced); isForced = 0; end
    if ischar(FileNames); FileNames = {FileNames}; end
    iModifiedSubjects = [];

    ProtocolInfo     = bst_get('ProtocolInfo');
    ProtocolSubjects = bst_get('ProtocolSubjects');

    % --- Resolve each node by filename: (subject, surface, node, type field, full path) ---
    FullFiles = {};
    R         = zeros(0, 3);   % [iSubject, iSurface, iNode] per resolved node
    fields    = {};            % 'Manifold' | 'Operator' | 'Eigen' per resolved node
    for i = 1:numel(FileNames)
        relFile = file_short(FileNames{i});
        switch lower(file_gettype(relFile))
            case 'manifold', field = 'Manifold'; getCase = 'ManifoldFile';
            case 'operator', field = 'Operator'; getCase = 'OperatorFile';
            case 'eigen',    field = 'Eigen';    getCase = 'EigenFile';
            otherwise,       continue;   % not a derived surface node
        end
        [~, iSubj, iSurf, iNode] = bst_get(getCase, relFile);
        if isempty(iSubj) || isempty(iSurf) || isempty(iNode); continue; end
        FullFiles{end+1} = bst_fullfile(ProtocolInfo.SUBJECTS, relFile); %#ok<AGROW>
        R(end+1, :)      = [iSubj, iSurf, iNode];                         %#ok<AGROW>
        fields{end+1}    = field;                                         %#ok<AGROW>
    end
    if isempty(FullFiles); return; end

    % --- Delete the files that still exist (one batched confirmation). Missing files
    %     (stale entries / DB drift) are tolerated: we still drop their DB entry below,
    %     otherwise a dangling registration can never be cleaned up. Only an interactive
    %     cancel (isForced=0 and file_delete reports failure) aborts the removal. ---
    existing = FullFiles(logical(cellfun(@(f) double(file_exist(f)), FullFiles)));
    if ~isempty(existing) && (file_delete(existing, isForced) ~= 1) && ~isForced
        return;   % user cancelled the confirmation dialog
    end

    % --- Remove the DB entries, grouped by (subject, surface, field), in descending
    %     index order so earlier removals do not shift the indices of later ones ---
    key = arrayfun(@(k) sprintf('%d|%d|%s', R(k,1), R(k,2), fields{k}), 1:size(R,1), 'UniformOutput', false);
    [uKey, ~, grp] = unique(key, 'stable');
    for g = 1:numel(uKey)
        idx    = find(grp == g);
        iSubj  = R(idx(1), 1);
        iSurf  = R(idx(1), 2);
        field  = fields{idx(1)};
        iNodes = sort(R(idx, 3), 'descend')';
        if (iSubj == 0)
            ProtocolSubjects.DefaultSubject.Surface(iSurf).(field)(iNodes) = [];
        else
            ProtocolSubjects.Subject(iSubj).Surface(iSurf).(field)(iNodes) = [];
        end
        iModifiedSubjects = [iModifiedSubjects, iSubj]; %#ok<AGROW>
    end
    bst_set('ProtocolSubjects', ProtocolSubjects);
    iModifiedSubjects = unique(iModifiedSubjects);
end
