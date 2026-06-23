function Op = bst_get_operator_node(SurfaceFile, Variant)
% BST_GET_OPERATOR_NODE: Find-or-create the operator_ node of a given Variant under a
% surface and return the loaded node.
%
% USAGE:  Op = bst_get_operator_node(SurfaceFile, Variant)
%
% Resolves the surface's Operator children via bst_get, loads the one whose Variant
% matches (in_bst_operator), and returns it. If no such operator exists yet, it is created
% with tess_operators(SurfaceFile, Variant) and re-resolved. This is the single shared
% resolver used by the Helmholtz callers (view_helmholtz, process_vortex_track), replacing
% the per-caller copies.
%
% INPUTS:
%    - SurfaceFile : protocol-relative path to the parent cortical surface
%    - Variant     : operator variant, e.g. 'Dirac', 'Laplace-Beltrami', 'Hodge-Face'
% OUTPUT:
%    - Op          : the loaded operator node (db_template('operatormat'))
%
% SEE ALSO: tess_operators, in_bst_operator, bst_helmholtz

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

    [sSubject, ~, iSurf] = bst_get('SurfaceFile', SurfaceFile);
    Op = [];
    if ~isempty(iSurf) && isfield(sSubject.Surface(iSurf), 'Operator')
        for k = 1:numel(sSubject.Surface(iSurf).Operator)
            S = in_bst_operator(sSubject.Surface(iSurf).Operator(k).FileName);
            if strcmpi(S.Variant, Variant), Op = S; break; end
        end
    end
    if isempty(Op)
        tess_operators(SurfaceFile, Variant);
        Op = bst_get_operator_node(SurfaceFile, Variant);
    end
end
