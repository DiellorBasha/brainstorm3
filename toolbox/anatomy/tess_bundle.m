function B = tess_bundle(SurfaceFile, varargin)
% TESS_BUNDLE: Store the nxr-compute coordinate bundle as per-hemisphere fields.
%
% USAGE:  B = tess_bundle(SurfaceFile)
%         B = tess_bundle(SurfaceFile, 'Gauge','trivial', 'NoSave',1, 'ForceRecompute',1)
%
% DESCRIPTION:
%     nxr-compute never integrates the two hemispheres, so each cortex is bundled
%     one connected genus-0 hemisphere at a time. The three results are stored as
%     1x2 struct arrays ((1)=left, (2)=right). Each element is the verbatim nxr
%     submesh bundle (LOCAL indexing) plus GlobalVertices/GlobalFaces/Hemisphere
%     scatter maps to TessMat.Vertices/.Faces global order.
%
% OPTIONS:
%     'Gauge'           'trivial' (default) | 'levi-civita' | 'euclidean'
%     'Operators'       false (default) | true   - heavy .operators (Task 2)
%     'Coupling'        'ambient' (default) | 'product'   - covariantLaplacian (heavy)
%     'Mass'            'lumped' (default) | 'galerkin'    - mass variant (heavy)
%     'NoSave'          false (default) | true
%     'ForceRecompute'  false (default) | true
%
% OUTPUT (also stored as TessMat.Topology/.Geometry/.Gauge, each 1x2):
%     B.Topology, B.Geometry, B.Gauge   - 1x2 struct arrays
%
% Requires the nxr-compute plugin; trivial gauge needs a FreeSurfer reg sphere.
%
% Authors: Diellor Basha, 2026

    % --- options ---
    Gauge='trivial'; Operators=false; Coupling='ambient'; Mass='lumped';
    NoSave=false; ForceRecompute=false;
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'gauge',          Gauge=varargin{i+1};
            case 'operators',      Operators=varargin{i+1};
            case 'coupling',       Coupling=varargin{i+1};
            case 'mass',           Mass=varargin{i+1};
            case 'nosave',         NoSave=varargin{i+1};
            case 'forcerecompute', ForceRecompute=varargin{i+1};
        end
    end

    TessFile = file_fullpath(SurfaceFile);
    TessMat  = in_tess_bst(SurfaceFile, 0);

    % --- cache return (light only; heavy is opt-in and always recomputes) ---
    if ~ForceRecompute && ~Operators && isfield(TessMat,'Topology') && ~isempty(TessMat.Topology)
        B = struct();                          % explicit assignment: struct() ctor
        B.Topology = TessMat.Topology;         % mishandles struct-array field values
        B.Geometry = TessMat.Geometry;
        B.Gauge    = TessMat.Gauge;
        return;
    end

    % --- trivial gauge needs a FreeSurfer registration sphere ---
    if strcmpi(Gauge,'trivial')
        if ~isfield(TessMat,'Reg') || ~isstruct(TessMat.Reg) || ~isfield(TessMat.Reg,'Sphere') ...
           || ~isfield(TessMat.Reg.Sphere,'Vertices') || isempty(TessMat.Reg.Sphere.Vertices)
            error('tess_bundle:noRegSphere', ...
                'Trivial gauge needs a FreeSurfer registration sphere (Reg.Sphere.Vertices).');
        end
    end

    % --- require nxr-compute (no MATLAB fallback) ---
    [isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
    if ~isOk
        error('tess_bundle:nxrUnavailable', 'tess_bundle requires nxr-compute: %s', errMsg);
    end
    nxrVer = '';
    try, nxrVer = nxr_compute('version'); catch, end   %#ok<CTCH>
    prov = struct('Backend','nxr', 'NxrVersion',nxrVer, 'Gauge',Gauge, ...
                  'Operators',logical(Operators), 'Coupling',Coupling, 'Mass',Mass, ...
                  'ComputeDate',datestr(now,'yyyy-mm-dd HH:MM:SS'));

    % --- hemisphere split from import labels (not geometry) ---
    [rH, lH, isConn] = tess_hemisplit(TessMat);
    if isConn
        error('tess_bundle:connectedHemispheres', ...
            'Hemispheres are connected; nxr bundles each as an independent component.');
    end
    hemis = {lH(:), rH(:)}; tags = {'L','R'};

    Vtx = double(TessMat.Vertices);
    Fcs = double(TessMat.Faces);
    nVtot = size(Vtx,1);

    for hh = 1:2
        vH = hemis{hh};
        if isempty(vH)
            error('tess_bundle:emptyHemisphere', 'Hemisphere %s has no vertices.', tags{hh});
        end
        isV = false(nVtot,1); isV(vH) = true;
        fMask = all(isV(Fcs), 2);
        map = zeros(nVtot,1); map(vH) = 1:numel(vH);
        Vloc = Vtx(vH,:);
        Floc = map(Fcs(fMask,:));

        % nxr submesh bundle
        h = nxr_compute('create', Vloc, Floc);
        opts = struct();
        if strcmpi(Gauge,'trivial')
            sph = TessMat.Reg.Sphere.Vertices(vH,:);
            [~, iN] = max(sph(:,3));   % north pole (local index)
            [~, iS] = min(sph(:,3));   % south pole (local index)
            opts.singVerts  = [iN; iS];
            opts.singValues = [1; 1];
        end
        Bh = nxr_compute('bundle', h, Gauge, opts);
        nxr_compute('destroy', h);

        % attach scatter maps; build the 1x2 struct arrays
        sT = Bh.Topology; sT.GlobalVertices = vH; sT.GlobalFaces = find(fMask); sT.Hemisphere = tags{hh}; sT.Provenance = prov;
        sG = Bh.Geometry; sG.GlobalVertices = vH; sG.GlobalFaces = find(fMask); sG.Hemisphere = tags{hh}; sG.Provenance = prov;
        sA = Bh.Gauge;    sA.GlobalVertices = vH; sA.GlobalFaces = find(fMask); sA.Hemisphere = tags{hh}; sA.Provenance = prov;
        TopoArr(hh) = sT;  GeoArr(hh) = sG;  GaArr(hh) = sA;  %#ok<AGROW>
    end

    B = struct();                  % explicit assignment (struct() ctor mishandles
    B.Topology = TopoArr;          % struct-array field values)
    B.Geometry = GeoArr;
    B.Gauge    = GaArr;

    % --- save ---
    if ~NoSave
        TessMat_full = load(TessFile);
        TessMat_full.Topology = B.Topology;
        TessMat_full.Geometry = B.Geometry;
        TessMat_full.Gauge    = B.Gauge;
        TessMat_full = bst_history('add', TessMat_full, 'bundle', ...
            sprintf('Stored nxr bundle (gauge=%s) as per-hemisphere Topology/Geometry/Gauge.', Gauge));
        bst_save(TessFile, TessMat_full, 'v7');
    end
end
