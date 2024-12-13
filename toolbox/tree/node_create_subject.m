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

        % Handle volpet nodes with frame children
        if strcmp(nodeType, 'volpet')
            % Check for frames in the corresponding @subjectimage directory
            [volumeDir, baseFileName, ~] = bst_fileparts(file_fullpath(char(sSubject.Anatomy(iAnatomy).FileName)));
            % Locate the frame directory (starts with @ and matches volume name)
            frameDir = bst_fullfile(volumeDir, ['@' baseFileName]);
            if isfolder(frameDir)
                % Multi-frame PET: Create the parent PET node
                [nodeCreated, nodePet] = CreateNode('volpet', ...
                    [char(sSubject.Anatomy(iAnatomy).Comment) ' (PET)'], ...
                    char(sSubject.Anatomy(iAnatomy).FileName), ...
                    iAnatomy, iSubject, iSearch);
                if nodeCreated
                    % Add frame children under PET node
                    petFrames = get_pet_frames(char(sSubject.Anatomy(iAnatomy).FileName)); % Retrieve frame file paths
                    for iFrame = 1:length(petFrames)
                        frameComment = sprintf('Frame %d', iFrame);
                        frameFileName = petFrames{iFrame};
                        [frameCreated, nodeFrame] = CreateNode('volpet', ...
                            frameComment, ...
                            frameFileName, ...
                            iFrame, iSubject, iSearch);
                        if frameCreated
                            nodePet.add(nodeFrame); % Add frame node to PET node
                            numElems = numElems + 1;
                        end
                    end
                    % Add PET node to subject
                    if showParentNodes
                        nodeSubject.add(nodePet);
                    else
                        nodeRoot.add(nodePet);
                    end
                    numElems = numElems + 1;
                end
            else
                % Single-frame PET: Create node directly
                [nodeCreated, nodePet] = CreateNode('volpet', ...
                    char(sSubject.Anatomy(iAnatomy).Comment), ...
                    char(sSubject.Anatomy(iAnatomy).FileName), ...
                    iAnatomy, iSubject, iSearch);
                if nodeCreated
                    if showParentNodes
                        nodeSubject.add(nodePet);
                    else
                        nodeRoot.add(nodePet);
                    end
                    numElems = numElems + 1;
                end
            end
        else
            % Create other anatomy nodes
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

function petFrames = get_pet_frames(petFileName)
% GET_PET_FRAMES: Retrieve list of PET frames from a PET volume file.
% INPUT:
%    - petFileName : Path to the PET volume file.
% OUTPUT:
%    - petFrames   : Cell array of frame file paths.

petFileFull = file_fullpath(petFileName);
% Extract directory and filename information
[volumeDir, baseFileName, ~] = bst_fileparts(petFileFull);
% Locate the frame directory (starts with @ and matches volume name)
frameDir = bst_fullfile(volumeDir, ['@' baseFileName]);
% Check if frame directory exists
if ~isfolder(frameDir)
    warning('Frame directory not found: %s', frameDir);
    petFrames = {};
    return;
end

% Look for frame files matching pattern `frame_*_volpet`
frameFiles = dir(fullfile(frameDir, 'frame_*_volpet.mat'));

% Get full file paths for each frame
petFrames = fullfile({frameFiles.folder}, {frameFiles.name});
end

