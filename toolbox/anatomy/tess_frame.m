function [U, V, N] = tess_frame(SurfaceFile, varargin)
% TESS_FRAME: Compute/store the nxr facet bundle per hemisphere; return the 3D frame.
%
% USAGE:  [U,V,N] = tess_frame(SurfaceFile)
%         [U,V,N] = tess_frame(SurfaceFile, 'Gauge','trivial', 'Domain','vertex', ...
%                              'ForceRecompute',1, 'NoSave',1)
%
% DESCRIPTION:
%     Loads the surface, splits hemispheres with tess_hemisplit (atlas L/R, never
%     conncomp), runs nxr_compute('facets', ...) on each hemisphere submesh (with
%     FreeSurfer-pole singularities for the trivial gauge), and stores the result
%     as five top-level 1x2 per-hemisphere struct arrays on the surface file:
%     TessMat.{Topology, Embedded, Intrinsic, Extrinsic, Gauge}. Each element is
%     the verbatim nxr facet struct (LOCAL indexing) plus GlobalVertices/
%     GlobalFaces/Hemisphere/Provenance scatter maps to the full-mesh order.
%
%     Returns the full-mesh intrinsic frame {U,V,N}: U=real(grid.*rot),
%     V=imag(grid.*rot), N=cross(U,V), where grid=Embedded.vertex.grid and
%     rot=Gauge.vertex.rotation (rot==1 for euclidean/levi-civita).
%
%     If the five fields are already present and ForceRecompute is false, the
%     frame is derived from the stored bundle without recomputing.
%
%     Vertex domain works for every gauge. Face domain works only for
%     euclidean/levi-civita: the trivial gauge's Gauge.face.rotation is deferred
%     in nxr, so face+trivial errors clearly.
%
% Requires the nxr-compute plugin; the trivial gauge needs a FreeSurfer reg sphere.
%
% SEE ALSO: tess_hemisplit, tess_tangents, tess_normals

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

    % --- options ---
    Gauge='trivial'; Domain='vertex'; NoSave=false; ForceRecompute=false;
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'gauge',          Gauge=lower(varargin{i+1});
            case 'domain',         Domain=lower(varargin{i+1});
            case 'nosave',         NoSave=varargin{i+1};
            case 'forcerecompute', ForceRecompute=varargin{i+1};
        end
    end
    if ~ismember(Domain, {'vertex','face'})
        error('tess_frame:badDomain', 'Domain must be ''vertex'' or ''face''.');
    end
    if ~ismember(Gauge, {'trivial','levi-civita','euclidean'})
        error('tess_frame:badGauge', 'Gauge must be ''trivial'', ''levi-civita'', or ''euclidean''.');
    end

    TessFile = file_fullpath(SurfaceFile);
    TessMat  = in_tess_bst(SurfaceFile, 0);

    Groups = {'Topology','Embedded','Intrinsic','Extrinsic','Gauge'};
    haveAll = all(cellfun(@(f) isfield(TessMat,f) && ~isempty(TessMat.(f)), Groups));

    if ForceRecompute || ~haveAll
        TessMat = local_compute_store(TessFile, TessMat, Gauge, NoSave);
    end

    [U, V, N] = local_derive_frame(TessMat, Domain);
end

% ----------------------------------------------------------------------------
function TessMat = local_compute_store(TessFile, TessMat, Gauge, NoSave)
    Groups = {'Topology','Embedded','Intrinsic','Extrinsic','Gauge'};

    % require nxr-compute (no MATLAB fallback)
    [isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
    if ~isOk
        error('tess_frame:nxrUnavailable', 'tess_frame requires nxr-compute: %s', errMsg);
    end

    % trivial gauge needs a FreeSurfer registration sphere
    if strcmpi(Gauge,'trivial')
        if ~isfield(TessMat,'Reg') || ~isstruct(TessMat.Reg) || ~isfield(TessMat.Reg,'Sphere') ...
           || ~isfield(TessMat.Reg.Sphere,'Vertices') || isempty(TessMat.Reg.Sphere.Vertices)
            error('tess_frame:noRegSphere', ...
                'Trivial gauge needs a FreeSurfer registration sphere (Reg.Sphere.Vertices).');
        end
    end

    % require a Structures atlas with L/R labels so tess_hemisplit uses the
    % atlas split — never the geometric/connectivity fallback
    hasLabels = false;
    if isfield(TessMat,'Atlas') && ~isempty(TessMat.Atlas)
        iStruct = find(strcmpi({TessMat.Atlas.Name}, 'Structures'), 1);
        if ~isempty(iStruct) && ~isempty(TessMat.Atlas(iStruct).Scouts)
            scouts = TessMat.Atlas(iStruct).Scouts;
            labels = {scouts.Label}; regions = {scouts.Region};
            reg1 = cellfun(@(c) c(1), regions(~cellfun(@isempty, regions)), 'UniformOutput', false);
            hasL = any(strcmpi(labels,'lh')) || any(strcmpi(reg1,'L'));
            hasR = any(strcmpi(labels,'rh')) || any(strcmpi(reg1,'R'));
            hasLabels = hasL && hasR;
        end
    end
    if ~hasLabels
        error('tess_frame:noHemisphereLabels', ...
            ['Surface has no Structures atlas with left/right hemisphere labels ' ...
             '(required for the atlas-based hemisphere split; the geometric fallback is not allowed).']);
    end

    % hemisphere split from atlas labels (never conncomp)
    [rH, lH, isConn] = tess_hemisplit(TessMat);
    if isConn
        error('tess_frame:connectedHemispheres', ...
            'Hemispheres are connected; nxr bundles each as an independent component.');
    end
    hemis = {lH(:), rH(:)}; tags = {'L','R'};
    Vtx = double(TessMat.Vertices); Fcs = double(TessMat.Faces); nVtot = size(Vtx,1);

    nxrVer = '';
    try, nxrVer = nxr_compute('version'); catch, end   %#ok<CTCH>
    prov = struct('Backend','nxr', 'Package','facets', 'NxrVersion',nxrVer, 'Gauge',Gauge, ...
                  'ComputeDate',datestr(now,'yyyy-mm-dd HH:MM:SS'));

    Arr = struct();
    for hh = 1:2
        vH = hemis{hh};
        if isempty(vH)
            error('tess_frame:emptyHemisphere', 'Hemisphere %s has no vertices.', tags{hh});
        end
        isV = false(nVtot,1); isV(vH) = true;
        fMask = all(isV(Fcs), 2);
        map = zeros(nVtot,1); map(vH) = 1:numel(vH);
        Vloc = Vtx(vH,:);
        Floc = map(Fcs(fMask,:));

        opts = struct();
        if strcmpi(Gauge,'trivial')
            sph = TessMat.Reg.Sphere.Vertices(vH,:);
            [~, iN] = max(sph(:,3));   % north pole (local index)
            [~, iS] = min(sph(:,3));   % south pole (local index)
            opts.singVerts  = [iN; iS];
            opts.singValues = [1; 1];
        end

        h = nxr_compute('create', Vloc, Floc);
        S = nxr_compute('facets', h, Gauge, opts);
        nxr_compute('destroy', h);

        for f = 1:numel(Groups)
            s = S.(Groups{f});
            s.GlobalVertices = vH;
            s.GlobalFaces    = find(fMask);
            s.Hemisphere     = tags{hh};
            s.Provenance     = prov;
            Arr.(Groups{f})(hh) = s;   %#ok<AGROW>
        end
    end

    for f = 1:numel(Groups), TessMat.(Groups{f}) = Arr.(Groups{f}); end

    if ~NoSave
        TessMat_full = load(TessFile);
        for f = 1:numel(Groups), TessMat_full.(Groups{f}) = TessMat.(Groups{f}); end
        TessMat_full = bst_history('add', TessMat_full, 'facets', ...
            sprintf(['Stored nxr facet bundle (gauge=%s) as per-hemisphere ' ...
                     'Topology/Embedded/Intrinsic/Extrinsic/Gauge.'], Gauge));
        bst_save(TessFile, TessMat_full, 'v7');
    end
end

% ----------------------------------------------------------------------------
function [U, V, N] = local_derive_frame(TessMat, Domain)
    Emb = TessMat.Embedded; Ga = TessMat.Gauge;
    if strcmp(Domain,'vertex')
        nElem = size(TessMat.Vertices,1);
    else
        nElem = size(TessMat.Faces,1);
    end
    U = zeros(nElem,3); V = zeros(nElem,3);

    for hh = 1:numel(Emb)
        gridH = Emb(hh).(Domain).grid;            % nElemH x 3 complex
        if strcmpi(Ga(hh).type, 'trivial')
            if strcmp(Domain,'face')
                if ~isfield(Ga(hh).face,'rotation') || isempty(Ga(hh).face.rotation)
                    error('tess_frame:faceTrivialDeferred', ...
                        'Face-domain trivial frame needs Gauge.face.rotation (empty/deferred in nxr).');
                end
                rot = Ga(hh).face.rotation;
            else
                rot = Ga(hh).vertex.rotation;
            end
        else
            rot = ones(size(gridH,1),1);
        end
        cRot = gridH .* rot;                       % broadcast over the 3 columns
        if strcmp(Domain,'vertex')
            idx = Emb(hh).GlobalVertices;
        else
            idx = Emb(hh).GlobalFaces;
        end
        U(idx,:) = real(cRot);
        V(idx,:) = imag(cRot);
    end

    % N is unit-norm because nxr returns an orthonormal tangent grid (|e1|=|e2|=1, e1⊥e2).
    N = cross(U, V, 2);
end
