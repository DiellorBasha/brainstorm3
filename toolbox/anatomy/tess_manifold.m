function ManifoldMat = tess_manifold(SurfaceFile, varargin)
% TESS_MANIFOLD: Assemble the nxr geometry backbone as a manifold_ DB node.
%
% USAGE:  ManifoldMat = tess_manifold(SurfaceFile)
%         ManifoldMat = tess_manifold(SurfaceFile, 'Gauge','trivial', ...
%                                    'NoSave',false, 'ForceRecompute',false)
%
% DESCRIPTION:
%     Loads the surface, splits hemispheres with tess_hemisplit (atlas L/R,
%     never conncomp), runs nxr_compute('facets', ...) on each hemisphere
%     submesh (with FreeSurfer-pole singularities for the trivial gauge),
%     and assembles a ManifoldMat structure (db_template('manifoldmat'))
%     containing 1x2 per-hemisphere Topology/Embedded/Intrinsic/Extrinsic/Gauge
%     arrays with GlobalVertices/GlobalFaces/Hemisphere scatter maps, plus a
%     Provenance record.
%
%     Unless NoSave is true, the result is saved as a manifold_*.mat file
%     alongside the parent surface and registered in the Brainstorm DB via
%     db_add_manifold (creating a child node under the surface).
%
%     If ForceRecompute is false and the surface already has a manifold child
%     registered, the existing node is LOADED and returned (no recompute), so
%     tess_manifold is a single find-or-load-or-create entry point that always
%     hands back the complete ManifoldMat. It is the canonical (and only) handler
%     for manifold-file data; other code should obtain the manifold via this call
%     rather than loading the node ad hoc.
%
% OPTIONS:
%     'Gauge'          : 'trivial' (default), 'levi-civita', or 'euclidean'
%     'NoSave'         : true/false (default false) — compute but do not write
%                        to disk or register in the DB
%     'ForceRecompute' : true/false (default false) — recompute even if a
%                        manifold node is already registered
%
% OUTPUT:
%     ManifoldMat : struct matching db_template('manifoldmat'), with fields:
%                   Comment, ParentSurface, Topology(1x2), Embedded(1x2),
%                   Intrinsic(1x2), Extrinsic(1x2), Gauge(1x2), Provenance
%
% Requires the nxr-compute plugin.  The trivial gauge also needs a FreeSurfer
% registration sphere (TessMat.Reg.Sphere.Vertices).
%
% SEE ALSO: tess_hemisplit, db_add_manifold, db_template, view_manifold

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

    % --- parse options ---
    Gauge          = 'trivial';
    NoSave         = false;
    ForceRecompute = false;
    Interactive    = false;  % GUI: prompt Overwrite/Cancel when a same-gauge manifold exists
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'gauge',          Gauge          = lower(varargin{i+1});
            case 'nosave',         NoSave         = logical(varargin{i+1});
            case 'forcerecompute', ForceRecompute = logical(varargin{i+1});
            case 'interactive',    Interactive    = logical(varargin{i+1});
            otherwise
                error('tess_manifold:badOption', 'Unknown option: %s', varargin{i});
        end
    end
    if ~ismember(Gauge, {'trivial','levi-civita','euclidean'})
        error('tess_manifold:badGauge', ...
            'Gauge must be ''trivial'', ''levi-civita'', or ''euclidean''.');
    end

    % --- cache / skip check ---
    % If not forcing, and a manifold of the requested Gauge is already registered, load
    % and return it (tess_manifold is the single find-or-load-or-create entry point, so it
    % always hands back the complete ManifoldMat). bst_get resolves the match from the cache.
    %
    % Schema gate: a node whose Embedded group predates the current schema (i.e. lacks the
    % per-face fields the consumers now rely on, e.g. face.area added in schemaVersion 2) is
    % treated as a cache MISS and recomputed -- so callers always receive a complete node.
    REQUIRED_EMBEDDED_SCHEMA = 2;   % nxr facets schemaVersion that first carries Embedded.face.area
    if ~ForceRecompute
        [sMan, ~, iSurfMan, iMan] = bst_get('ManifoldFileForSurface', SurfaceFile, Gauge);
        if ~isempty(iMan)
            existFile = sMan.Surface(iSurfMan).Manifold(iMan).FileName;
            Mcand = in_bst_manifold(existFile);
            isStale = isempty(Mcand.Embedded) || ~isfield(Mcand.Embedded, 'schemaVersion') ...
                      || isempty(Mcand.Embedded(1).schemaVersion) ...
                      || (Mcand.Embedded(1).schemaVersion < REQUIRED_EMBEDDED_SCHEMA) ...
                      || ~isfield(Mcand.Embedded(1).face, 'area') ...
                      || ~isfield(Mcand, 'DEC') || isempty(Mcand.DEC) ...     % node predates the DEC group
                      || ~isfield(Mcand.DEC, 'sharp') || isfield(Mcand.DEC, 'flat');   % 'flat' = old unstable build -> recompute
            if ~isStale
                if Interactive
                    % Prompt Overwrite / Cancel. Cancel reuses; Overwrite deletes and recomputes.
                    msg = sprintf(['A manifold (gauge=%s) already exists for this surface.\n\n' ...
                        'Overwrite it (delete and recompute)?\n[No keeps and reuses the existing one.]'], Gauge);
                    if ~java_dialog('confirm', msg, 'Compute manifold')
                        ManifoldMat = Mcand;
                        return;
                    end
                    db_delete_surface_node(existFile, 1);   % overwrite: drop the old node
                else
                    ManifoldMat = Mcand;
                    return;
                end
            end
            % stale schema -> fall through to recompute (db_add_manifold replaces the same-gauge node)
        end
    end

    % --- load surface ---
    TessMat = in_tess_bst(SurfaceFile, 0);

    % --- guard: nxr-compute plugin ---
    [isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
    if ~isOk
        error('tess_manifold:nxrUnavailable', ...
            'tess_manifold requires nxr-compute: %s', errMsg);
    end

    % --- guard: trivial gauge needs Reg.Sphere ---
    if strcmpi(Gauge, 'trivial')
        if ~isfield(TessMat,'Reg') || ~isstruct(TessMat.Reg) || ...
           ~isfield(TessMat.Reg,'Sphere') || ...
           ~isfield(TessMat.Reg.Sphere,'Vertices') || ...
           isempty(TessMat.Reg.Sphere.Vertices)
            error('tess_manifold:noRegSphere', ...
                'Trivial gauge needs a FreeSurfer registration sphere (Reg.Sphere.Vertices).');
        end
    end

    % --- guard: Structures atlas with L/R labels ---
    hasLabels = false;
    if isfield(TessMat,'Atlas') && ~isempty(TessMat.Atlas)
        iStruct = find(strcmpi({TessMat.Atlas.Name}, 'Structures'), 1);
        if ~isempty(iStruct) && ~isempty(TessMat.Atlas(iStruct).Scouts)
            scouts  = TessMat.Atlas(iStruct).Scouts;
            labels  = {scouts.Label};
            regions = {scouts.Region};
            reg1 = cellfun(@(c) c(1), regions(~cellfun(@isempty, regions)), 'UniformOutput', false);
            hasL = any(strcmpi(labels,'lh')) || any(strcmpi(reg1,'L'));
            hasR = any(strcmpi(labels,'rh')) || any(strcmpi(reg1,'R'));
            hasLabels = hasL && hasR;
        end
    end
    if ~hasLabels
        error('tess_manifold:noHemisphereLabels', ...
            ['Surface has no Structures atlas with left/right hemisphere labels ' ...
             '(required for the atlas-based hemisphere split; the geometric fallback is not allowed).']);
    end

    % --- hemisphere split (atlas-based, never conncomp) ---
    [rH, lH, isConn] = tess_hemisplit(TessMat);
    if isConn
        error('tess_manifold:connectedHemispheres', ...
            'Hemispheres are connected; nxr requires each hemisphere as an independent component.');
    end
    hemis = {lH(:), rH(:)};
    tags  = {'L','R'};
    Vtx   = double(TessMat.Vertices);
    Fcs   = double(TessMat.Faces);
    nVtot = size(Vtx, 1);

    % nxr version for provenance
    nxrVer = '';
    try, nxrVer = nxr_compute('version'); catch, end  %#ok<CTCH>

    prov = struct('Backend','nxr', 'NxrVersion',nxrVer, 'Gauge',Gauge, ...
                  'ComputeDate',datestr(now,'yyyy-mm-dd HH:MM:SS'));

    Groups = {'Topology','Embedded','Intrinsic','Extrinsic','Gauge'};
    Arr    = struct();
    decArr = repmat(struct('d0',[],'d1',[],'h0',[],'h1',[],'h2',[],'sharp',[], ...
                           'GlobalVertices',[],'GlobalFaces',[],'Hemisphere',[],'Provenance',[]), 1, 2);

    for hh = 1:2
        vH = hemis{hh};
        if isempty(vH)
            error('tess_manifold:emptyHemisphere', ...
                'Hemisphere %s has no vertices.', tags{hh});
        end

        % build local submesh (local indices)
        isV   = false(nVtot, 1);  isV(vH) = true;
        fMask = all(isV(Fcs), 2);
        mapV  = zeros(nVtot, 1);  mapV(vH) = 1:numel(vH);
        Vloc  = Vtx(vH, :);
        Floc  = mapV(Fcs(fMask, :));

        % gauge-specific singularity options (trivial: FreeSurfer poles)
        opts = struct();
        if strcmpi(Gauge, 'trivial')
            sph = TessMat.Reg.Sphere.Vertices(vH, :);
            [~, iN] = max(sph(:,3));   % north pole (local index)
            [~, iS] = min(sph(:,3));   % south pole (local index)
            opts.singVerts  = [iN; iS];
            opts.singValues = [1; 1];
        end

        % run nxr facets + the discrete exterior calculus (DEC) operators of this hemisphere.
        % DEC is the manifold's calculus structure (not a parametrized analysis operator), so it
        % lives here next to the geometry it derives from, NOT as a standalone operator_ node.
        h = nxr_compute('create', Vloc, Floc);
        S = nxr_compute('facets', h, Gauge, opts);
        DC = nxr_compute('operators', h, 'dec');                 % d0 [E x V], d1 [F x E]
        dh0 = nxr_compute('operators', h, 'hodge', 'h0');        % Hodge stars (diagonal)
        dh1 = nxr_compute('operators', h, 'hodge', 'h1');
        dh2 = nxr_compute('operators', h, 'hodge', 'h2');
        nxr_compute('destroy', h);
        % Whitney sharp ♯ [3F x E] (musical isomorphism): 1-form -> tangent face vector. nxr only
        % APPLIES Whitney, so assemble the matrix here, using the WINDING normal (Whitney's
        % face-normal convention). The flat ♭ is NOT stored: the metric flat ⋆1^{-1}♯'W_F is
        % numerically unstable (⋆1 = cotan has negative/near-zero entries on obtuse triangles),
        % and toolbox/differential builds div/curl as STABLE composites that avoid ⋆1^{-1}
        % (div = -⋆0^{-1} d0' ♯' W_F ; curl v = -div(N x v), the surface Hodge rotation).
        e1f = Vloc(Floc(:,2),:) - Vloc(Floc(:,1),:);
        e2f = Vloc(Floc(:,3),:) - Vloc(Floc(:,1),:);
        Ncr = cross(e1f, e2f, 2);  twoA = sqrt(sum(Ncr.^2, 2));
        Nf  = Ncr ./ max(twoA, eps);  fArea = twoA / 2;
        sharp = local_whitney_sharp(Vloc, Floc, Nf, fArea, DC.d0);

        % scatter-map attachment + store (facets)
        for f = 1:numel(Groups)
            s = S.(Groups{f});
            s.GlobalVertices = vH;
            s.GlobalFaces    = find(fMask);
            s.Hemisphere     = tags{hh};
            s.Provenance     = prov;
            Arr.(Groups{f})(hh) = s;  %#ok<AGROW>
        end
        % store the DEC group (same scatter-map attachment as the facets)
        decArr(hh) = struct('d0', DC.d0, 'd1', DC.d1, 'h0', dh0, 'h1', dh1, 'h2', dh2, ...
                            'sharp', sharp, ...
                            'GlobalVertices', vH, 'GlobalFaces', find(fMask), ...
                            'Hemisphere', tags{hh}, 'Provenance', prov);
    end

    % --- assemble ManifoldMat ---
    ManifoldMat = db_template('manifoldmat');
    for f = 1:numel(Groups)
        ManifoldMat.(Groups{f}) = Arr.(Groups{f});
    end
    ManifoldMat.DEC        = decArr;
    ManifoldMat.Provenance = prov;

    % --- save / register in DB ---
    if ~NoSave
        % Resolve subject index (needed by db_add_manifold)
        [sSubjectSave, iSubjectSave] = bst_get('SurfaceFile', SurfaceFile);
        if isempty(sSubjectSave)
            error('tess_manifold:subjectNotFound', ...
                'Could not resolve subject for surface: %s', SurfaceFile);
        end
        Comment = sprintf('Manifold (gauge=%s)', Gauge);
        db_add_manifold(iSubjectSave, SurfaceFile, ManifoldMat, Comment);
    end
end


%% ===== Whitney sharp (musical isomorphism #) operator [3F x E] =====
function sharp = local_whitney_sharp(Vloc, Floc, Nf, Af, d0)
% Discrete sharp ♯: edge 1-form -> per-face ambient vector, exactly geometry-central's Whitney
% reconstruction (nxr only APPLIES it). For a face (i,j,k) in winding order:
%   V_f = N_f x ( c_ij(e_ki - e_jk) + c_jk(e_ij - e_ki) + c_ki(e_jk - e_ij) ) / (6 A_f),
% where c_XY is the 1-form value on that edge signed by the edge's CANONICAL orientation (from
% d0). Rows interleaved so face k holds (x,y,z) in rows 3k-2:3k. Composing ♯ * d0 recovers the
% FEM gradient (the correctness oracle); ♯ matches nxr's applied Whitney to machine precision.
    nFh = size(Floc, 1);  nVh = size(Vloc, 1);  nE = size(d0, 1);
    % per-edge tail/head from d0 ((d0 f)_e = f(head) - f(tail) => d0(e,head)=+1, d0(e,tail)=-1)
    [er, vc, vv] = find(d0);
    headV = zeros(nE,1);  tailV = zeros(nE,1);
    headV(er(vv > 0)) = vc(vv > 0);  tailV(er(vv < 0)) = vc(vv < 0);
    Elookup = sparse(min(tailV,headV), max(tailV,headV), (1:nE)', nVh, nVh);
    vi = Floc(:,1);  vj = Floc(:,2);  vk = Floc(:,3);
    eij = Vloc(vj,:) - Vloc(vi,:);
    ejk = Vloc(vk,:) - Vloc(vj,:);
    eki = Vloc(vi,:) - Vloc(vk,:);
    [eID_ij, s_ij] = i_face_edge(vi, vj, Elookup, tailV, headV);
    [eID_jk, s_jk] = i_face_edge(vj, vk, Elookup, tailV, headV);
    [eID_ki, s_ki] = i_face_edge(vk, vi, Elookup, tailV, headV);
    coef = 1 ./ (6 * Af);
    Vij = cross(Nf, eki - ejk, 2) .* (coef .* s_ij);   % dV/dc_ij per face [nFh x 3]
    Vjk = cross(Nf, eij - eki, 2) .* (coef .* s_jk);
    Vki = cross(Nf, ejk - eij, 2) .* (coef .* s_ki);
    rx = (3*(1:nFh)-2)';  ry = (3*(1:nFh)-1)';  rz = (3*(1:nFh))';
    rows = [rx; ry; rz;  rx; ry; rz;  rx; ry; rz];
    cols = [repmat(eID_ij,3,1); repmat(eID_jk,3,1); repmat(eID_ki,3,1)];
    vals = [Vij(:,1); Vij(:,2); Vij(:,3);  Vjk(:,1); Vjk(:,2); Vjk(:,3);  Vki(:,1); Vki(:,2); Vki(:,3)];
    sharp = sparse(rows, cols, vals, 3*nFh, nE);
end

% global edge index + orientation sign for a face's directed edge a->b (sign +1 iff canonical)
function [eID, s] = i_face_edge(a, b, Elookup, tailV, headV)
    eID = full(Elookup(sub2ind(size(Elookup), min(a,b), max(a,b))));
    s   = 2*(tailV(eID) == a & headV(eID) == b) - 1;
end
