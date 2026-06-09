function B = tess_store_perhemi(SurfaceFile, FieldNames, NeedSphere, Provenance, UseCache, NoSave, ComputeFn)
% TESS_STORE_PERHEMI: Shared per-hemisphere nxr-bundle orchestration.
%
% Splits the cortex by hemisphere, runs ComputeFn on each hemisphere submesh,
% attaches scatter maps, assembles 1x2 struct arrays, and stores the requested
% fields on the surface file. Backs the tess_topology/geometry/gauge/bundle writers.
%
% INPUTS:
%   SurfaceFile  Brainstorm surface file name.
%   FieldNames   cellstr subset of {'Topology','Geometry','Gauge'} to store.
%   NeedSphere   logical; require a FreeSurfer reg sphere + compute poles (gauge).
%   Provenance   struct attached as .Provenance on every stored element.
%   UseCache     logical; if true and all FieldNames already present, return cached.
%   NoSave       logical; skip writing to disk.
%   ComputeFn    @(h, polesLocal) -> struct whose fields are (a superset of)
%                FieldNames, each the nxr struct for that hemisphere submesh.
%                polesLocal is [iN; iS] local indices (empty when ~NeedSphere).
%
% OUTPUT: B, struct with the FieldNames as 1x2 per-hemisphere struct arrays.
%
% Authors: Diellor Basha, 2026

    TessFile = file_fullpath(SurfaceFile);
    TessMat  = in_tess_bst(SurfaceFile, 0);

    % --- cache return: all requested fields already present ---
    if UseCache
        haveAll = true;
        for f = 1:numel(FieldNames)
            if ~isfield(TessMat, FieldNames{f}) || isempty(TessMat.(FieldNames{f}))
                haveAll = false; break;
            end
        end
        if haveAll
            B = struct();
            for f = 1:numel(FieldNames), B.(FieldNames{f}) = TessMat.(FieldNames{f}); end
            return;
        end
    end

    % --- require nxr-compute ---
    [isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
    if ~isOk
        error('tess_store_perhemi:nxrUnavailable', 'requires nxr-compute: %s', errMsg);
    end

    % --- require FreeSurfer registration sphere (gauge only) ---
    if NeedSphere
        if ~isfield(TessMat,'Reg') || ~isstruct(TessMat.Reg) || ~isfield(TessMat.Reg,'Sphere') ...
           || ~isfield(TessMat.Reg.Sphere,'Vertices') || isempty(TessMat.Reg.Sphere.Vertices)
            error('tess_store_perhemi:noRegSphere', ...
                'Trivial gauge needs a FreeSurfer registration sphere (Reg.Sphere.Vertices).');
        end
    end

    % --- require Structures atlas with L/R labels ---
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
        error('tess_store_perhemi:noHemisphereLabels', ...
            'Surface has no Structures atlas with left/right hemisphere labels.');
    end

    % --- hemisphere split (import labels) ---
    [rH, lH, isConn] = tess_hemisplit(TessMat);
    if isConn
        error('tess_store_perhemi:connectedHemispheres', ...
            'Hemispheres are connected; nxr bundles each as an independent component.');
    end
    hemis = {lH(:), rH(:)}; tags = {'L','R'};
    Vtx = double(TessMat.Vertices); Fcs = double(TessMat.Faces); nVtot = size(Vtx,1);

    Arr = struct();
    for hh = 1:2
        vH = hemis{hh};
        if isempty(vH)
            error('tess_store_perhemi:emptyHemisphere', 'Hemisphere %s has no vertices.', tags{hh});
        end
        isV = false(nVtot,1); isV(vH) = true;
        fMask = all(isV(Fcs), 2);
        map = zeros(nVtot,1); map(vH) = 1:numel(vH);
        Vloc = Vtx(vH,:);
        Floc = map(Fcs(fMask,:));

        poles = [];
        if NeedSphere
            sph = TessMat.Reg.Sphere.Vertices(vH,:);
            [~, iN] = max(sph(:,3));
            [~, iSp] = min(sph(:,3));
            poles = [iN; iSp];
        end

        h = nxr_compute('create', Vloc, Floc);
        S = ComputeFn(h, poles);
        nxr_compute('destroy', h);

        for f = 1:numel(FieldNames)
            s = S.(FieldNames{f});
            s.GlobalVertices = vH;
            s.GlobalFaces    = find(fMask);
            s.Hemisphere     = tags{hh};
            s.Provenance     = Provenance;
            Arr.(FieldNames{f})(hh) = s;   %#ok<AGROW>
        end
    end

    B = struct();
    for f = 1:numel(FieldNames), B.(FieldNames{f}) = Arr.(FieldNames{f}); end

    % --- save ---
    if ~NoSave
        TessMat_full = load(TessFile);
        for f = 1:numel(FieldNames), TessMat_full.(FieldNames{f}) = B.(FieldNames{f}); end
        TessMat_full = bst_history('add', TessMat_full, 'bundle', ...
            sprintf('Stored %s (per-hemisphere nxr).', strjoin(FieldNames, '/')));
        bst_save(TessFile, TessMat_full, 'v7');
    end
end
