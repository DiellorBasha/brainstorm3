function numElems = node_create_subject(nodeSubject, nodeRoot, sSubject, iSubject, iSearch)
% NODE_CREATE_SUBJECT: Create subject node from subject structure.
%
% USAGE:  node_create_subject(nodeSubject, nodeRoot, sSubject, iSubject)
%
% INPUT: 
%     - nodeSubject : BstNode object with Type 'subject' => Root of the subject subtree
%     - nodeRoot    : BstNode object, root of the whole database tree
%     - sSubject    : Brainstorm subject structure
%     - iSubject    : indice of the subject node in Brainstorm subjects list
%     - iSearch     : ID of the active DB search, or empty/0 if none
% OUTPUT:
%     - numElems    : Number of node children elements (including self) that
%                     pass the active search filter. If 0, node should be hidden

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
% Authors: Francois Tadel, 2008-2016
%          Martin Cousineau, 2019-2020

% If iSubject=0 => default subject
import org.brainstorm.tree.*;

% Parse inputs
if nargin < 4 || isempty(iSearch) || iSearch == 0
    iSearch = 0;
    % No search applied: ensure the node is added to the database
    numElems = 1;
else
    numElems = 0;
end
showParentNodes = node_show_parents(iSearch);

% Update node fields
nodeSubject.setFileName(sSubject.FileName);
nodeSubject.setItemIndex(0);
nodeSubject.setStudyIndex(iSubject);
if (iSubject ~= 0)
    nodeSubject.setComment(sSubject.Name);
else
    nodeSubject.setComment('(Default anatomy)');
end

% Anatomy files to use : Individual or Protocol defaults
% ==== Default anatomy ====
if sSubject.UseDefaultAnat && (iSubject ~= 0)
    nodeLink = BstNode('defaultanat', '(Default anatomy)', '', 0, 0);
    nodeSubject.add(nodeLink);

% ==== Individual anatomy ====
else
    % Create list of anat files (put the default at the top)
    iAnatList = 1:length(sSubject.Anatomy);
    iAtlas = find(~cellfun(@(c)(isempty(strfind(char(c), '_volatlas')) && isempty(strfind(char(c), '_tissues'))), {sSubject.Anatomy.FileName}));
    iCt    = find(cellfun(@(c)(~isempty(strfind(char(c), '_volct'))), {sSubject.Anatomy.FileName}));
    iPet   = find(cellfun(@(c)(~isempty(strfind(char(c), '_volpet'))), {sSubject.Anatomy.FileName}));
    if (length(sSubject.Anatomy) > 1)
        iAnatList = [sSubject.iAnatomy, setdiff(iAnatList,[iAtlas,sSubject.iAnatomy]), setdiff(iAtlas,sSubject.iAnatomy)];
    end
    % Create and add anatomy nodes
    for iAnatomy = iAnatList
        if ismember(iAnatomy, iAtlas)
            nodeType = 'volatlas';
        elseif ismember(iAnatomy, iCt)
            nodeType = 'volct';
        elseif ismember(iAnatomy, iPet)
            nodeType = 'volpet';
        else
            nodeType = 'anatomy';
        end
        [nodeCreated, nodeAnatomy] = CreateNode(nodeType, ...
            char(sSubject.Anatomy(iAnatomy).Comment), ...
            char(sSubject.Anatomy(iAnatomy).FileName), ...
            iAnatomy, iSubject, iSearch);

        if nodeCreated
            % If current item is default one
            if ismember(iAnatomy, sSubject.iAnatomy)
                nodeAnatomy.setMarked(1);
            end
            if showParentNodes
                nodeSubject.add(nodeAnatomy);
            else
                nodeRoot.add(nodeAnatomy);
            end
            numElems = numElems + 1;
        end
    end

    % Sort surfaces by category
    SortedSurfaces = db_surface_sort(sSubject.Surface);
    iSorted = [SortedSurfaces.IndexScalp, SortedSurfaces.IndexOuterSkull, SortedSurfaces.IndexInnerSkull, ...
               SortedSurfaces.IndexCortex, SortedSurfaces.IndexOther, SortedSurfaces.IndexFibers, SortedSurfaces.IndexFEM];
    % Process all the surfaces
    for i = 1:length(iSorted)
        iSurface = iSorted(i);
        SurfaceType = sSubject.Surface(iSurface).SurfaceType;
        % Create a node adapted to represent this surface
        [nodeCreated, nodeSurface] = CreateNode(lower(SurfaceType), ...
            char(sSubject.Surface(iSurface).Comment), ...
            char(sSubject.Surface(iSurface).FileName), ...
            iSurface, iSubject, iSearch);
        if nodeCreated
            % If current item is default one
            if ismember(iSurface, sSubject.(['i' SurfaceType]))
                nodeSurface.setMarked(1);
            end
            if showParentNodes
                nodeSubject.add(nodeSurface);
            else
                nodeRoot.add(nodeSurface);
            end
            numElems = numElems + 1;
            % Nest manifold child nodes under this surface
            if isfield(sSubject.Surface(iSurface), 'Manifold')
                for iM = 1:numel(sSubject.Surface(iSurface).Manifold)
                    [chCreated, chNode] = CreateNode('manifold', ...
                        char(sSubject.Surface(iSurface).Manifold(iM).Comment), ...
                        char(sSubject.Surface(iSurface).Manifold(iM).FileName), ...
                        iM, iSubject, iSearch);
                    if chCreated
                        nodeSurface.add(chNode);
                    end
                end
            end
            % Nest operator child nodes under this surface. With more than one operator,
            % collapse them under an 'operatorlist' container (a pure display grouping, like
            % the recordings 'datalist': no file, no DB entry); a single operator stays direct.
            if isfield(sSubject.Surface(iSurface), 'Operator')
                nOp = numel(sSubject.Surface(iSurface).Operator);
                opParent = nodeSurface;  opGroup = [];
                if nOp > 1
                    opGroup = org.brainstorm.tree.BstNode('operatorlist', 'Operators', '', 0, iSubject);
                    opParent = opGroup;
                end
                nOpAdded = 0;
                for iO = 1:nOp
                    % Child label = the operator Variant ('Dirac', 'Laplace-Beltrami', ...): the
                    % 'Operators' container conveys the kind, so the bare variant avoids the
                    % redundant "... operator" suffix. Fall back to the Comment if Variant is empty.
                    opLabel = char(sSubject.Surface(iSurface).Operator(iO).Variant);
                    if isempty(opLabel)
                        opLabel = char(sSubject.Surface(iSurface).Operator(iO).Comment);
                    end
                    [chCreated, chNode] = CreateNode('operator', opLabel, ...
                        char(sSubject.Surface(iSurface).Operator(iO).FileName), ...
                        iO, iSubject, iSearch);
                    if chCreated
                        opParent.add(chNode);
                        nOpAdded = nOpAdded + 1;
                    end
                end
                if ~isempty(opGroup) && (nOpAdded > 0)
                    nodeSurface.add(opGroup);   % attach the container once its children are built
                end
            end
            % Nest eigen child nodes under this surface. With more than one eigenbasis,
            % collapse them under an 'eigenlist' container (same pure-display grouping as the
            % operators); a single eigenbasis stays direct. Wavelet children still nest under
            % their own eigen node, wherever that eigen node sits.
            if isfield(sSubject.Surface(iSurface), 'Eigen')
                nEig = numel(sSubject.Surface(iSurface).Eigen);
                eigParent = nodeSurface;  eigGroup = [];
                if nEig > 1
                    eigGroup = org.brainstorm.tree.BstNode('eigenlist', 'Eigen', '', 0, iSubject);
                    eigParent = eigGroup;
                end
                nEigAdded = 0;
                for iE = 1:nEig
                    % Child label = the eigenbasis Variant ('Dirac', 'Laplace-Beltrami', ...): the
                    % 'Eigen' container conveys the kind, so the bare variant avoids the redundant
                    % "... eigenmodes (K=...)" suffix. Fall back to the Comment if Variant is empty.
                    eigLabel = char(sSubject.Surface(iSurface).Eigen(iE).Variant);
                    if isempty(eigLabel)
                        eigLabel = char(sSubject.Surface(iSurface).Eigen(iE).Comment);
                    end
                    [chCreated, chNode] = CreateNode('eigen', eigLabel, ...
                        char(sSubject.Surface(iSurface).Eigen(iE).FileName), ...
                        iE, iSubject, iSearch);
                    if chCreated
                        eigParent.add(chNode);
                        nEigAdded = nEigAdded + 1;
                        % Nest wavelet child nodes under THIS eigen node
                        if isfield(sSubject.Surface(iSurface), 'Wavelet')
                            eigName = char(sSubject.Surface(iSurface).Eigen(iE).FileName);
                            ws = sSubject.Surface(iSurface).Wavelet;
                            for iW = 1:numel(ws)
                                if file_compare(ws(iW).ParentEigen, eigName)
                                    [wCreated, wNode] = CreateNode('wavelet', ...
                                        char(ws(iW).Comment), char(ws(iW).FileName), ...
                                        iW, iSubject, iSearch);
                                    if wCreated
                                        chNode.add(wNode);
                                    end
                                end
                            end
                        end
                    end
                end
                if ~isempty(eigGroup) && (nEigAdded > 0)
                    nodeSurface.add(eigGroup);   % attach the container once its children are built
                end
            end
        end
    end
end
end

% Create a Java object for a database node if it passes the active search
%
% Inputs:
%  - nodeType to iStudy: See BstJava's constructor
%  - iSearch: ID of the active search filter (or 0 if none)
%
% Outputs:
%  - isCreated: Whether the node was actually created (1 or 0)
%  - node: Newly created Java object for the node
function [isCreated, node] = CreateNode(nodeType, nodeComment, ...
        nodeFileName, iItem, iStudy, iSearch)
    import org.brainstorm.tree.BstNode;
    % Only create Java object is required
    [isCreated, filteredComment] = node_apply_search(iSearch, nodeType, nodeComment, nodeFileName);
    if isCreated
        node = BstNode(nodeType, filteredComment, nodeFileName, iItem, iStudy);
    else
        node = [];
    end
end
