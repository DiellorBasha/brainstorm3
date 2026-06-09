function [U, V, N] = tess_frame(SurfaceFile, varargin)
% TESS_FRAME: Derived per-element intrinsic frame {U,V,N} from the stored bundle.
%
% USAGE:  [U,V,N] = tess_frame(SurfaceFile)                    % 'vertex' domain
%         [U,V,N] = tess_frame(SurfaceFile, 'Domain','face')
%
% DESCRIPTION:
%     Reads the per-hemisphere TessMat.Geometry/.Gauge (1x2), applies the gauge
%     rotation to the complex grid (c = e1 + i*e2), and scatters to full-mesh
%     order via GlobalVertices/GlobalFaces. Returns U=real, V=imag, N=cross(U,V).
%     Computes nothing new; persists nothing.
%
%     Vertex domain works for every gauge. Face domain works only for
%     euclidean/levi-civita: the trivial gauge's Gauge.face.rotation is deferred
%     in nxr, so face+trivial errors clearly.
%
% SEE ALSO: tess_bundle, tess_tangents, tess_normals

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

    Domain = 'vertex';
    for i = 1:2:numel(varargin)
        if strcmpi(varargin{i}, 'domain'); Domain = lower(varargin{i+1}); end
    end
    if ~ismember(Domain, {'vertex','face'})
        error('tess_frame:badDomain', 'Domain must be ''vertex'' or ''face''.');
    end

    TessMat = in_tess_bst(SurfaceFile, 0);
    if ~isfield(TessMat,'Geometry') || isempty(TessMat.Geometry) ...
       || ~isfield(TessMat,'Gauge') || isempty(TessMat.Gauge)
        error('tess_frame:noBundle', 'Surface has no stored bundle; run tess_bundle first.');
    end
    Geo = TessMat.Geometry; Ga = TessMat.Gauge;

    if strcmp(Domain,'vertex')
        nElem = size(TessMat.Vertices,1);
    else
        nElem = size(TessMat.Faces,1);
    end
    U = zeros(nElem,3); V = zeros(nElem,3);

    for hh = 1:numel(Geo)
        gridH = Geo(hh).(Domain).grid;          % nElemH x 3 complex
        if strcmpi(Ga(hh).type, 'trivial')
            if strcmp(Domain,'face')
                % nxr ships Gauge.face.rotation as a deferred (empty) placeholder.
                % Future-proof: support it the moment the build populates it.
                if ~isfield(Ga(hh).face,'rotation') || isempty(Ga(hh).face.rotation)
                    error('tess_frame:faceTrivialDeferred', ...
                        'Face-domain trivial frame needs Gauge.face.rotation (empty/deferred in nxr).');
                end
                rot = Ga(hh).face.rotation;      % nElemH x 1 complex (when populated)
            else
                rot = Ga(hh).vertex.rotation;    % nElemH x 1 complex
            end
        else
            rot = ones(size(gridH,1),1);
        end
        cRot = gridH .* rot;                      % broadcast over the 3 columns
        if strcmp(Domain,'vertex')
            idx = Geo(hh).GlobalVertices;
        else
            idx = Geo(hh).GlobalFaces;
        end
        U(idx,:) = real(cRot);
        V(idx,:) = imag(cRot);
    end

    N = cross(U, V, 2);
end
