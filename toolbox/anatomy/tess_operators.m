function OperatorMat = tess_operators(SurfaceFile, OperatorName, varargin)
% TESS_OPERATORS: Assemble per-hemisphere discrete operators as an operator_ DB node.
%
% USAGE:  OperatorMat = tess_operators(SurfaceFile, OperatorName)
%         OperatorMat = tess_operators(SurfaceFile, OperatorName, ...
%                                      'Tau',0.5, 'NoSave',false, 'ForceRecompute',false)
%
% DESCRIPTION:
%     Loads the surface, splits hemispheres with tess_hemisplit (atlas L/R,
%     never conncomp), builds each hemisphere as an independent nxr submesh
%     (mask + reindex to local indices), and computes the requested operator
%     pencil (A, B) per hemisphere via nxr_compute('operators', ...).
%
%     Three operator variants are supported (case-insensitive OperatorName):
%       'Laplace-Beltrami'     A = laplacian/cotan      [nVh x nVh] real symmetric
%                              B = mass/galerkin        [nVh x nVh]
%       'Connection Laplacian' A = laplacian/connection [nVh x nVh] Hermitian
%                              B = mass/galerkin        [nVh x nVh]
%       'Dirac'                A = dirac(Tau)           [4nVh x 4nVh]
%                              B = kron(mass/galerkin, I4)
%
%     The result is assembled into an OperatorMat structure
%     (db_template('operatormat')) holding 1x2 per-hemisphere Operator/Mass
%     arrays (cell-wrapped sparse matrices), a 1x2 GlobalVertices scatter map,
%     the Variant, and a Provenance record.
%
%     Unless NoSave is true, the result is saved as an operator_*.mat file
%     alongside the parent surface and registered in the Brainstorm DB via
%     db_add_operator (creating a child node under the surface).
%
% OPTIONS:
%     'Tau'            : scalar in [0,1] (default 0.5) — relative-Dirac mixing
%                        parameter (only used by the 'Dirac' variant)
%     'NoSave'         : true/false (default false) — compute but do not write
%                        to disk or register in the DB
%     'ForceRecompute' : true/false (default false) — currently accepted for
%                        API symmetry with the other assemblers; tess_operators
%                        always recomputes (operator nodes are not de-duplicated)
%
% OUTPUT:
%     OperatorMat : struct matching db_template('operatormat'), with fields:
%                   Comment, ParentSurface, Variant, Operator(1x2), Mass(1x2),
%                   GlobalVertices(1x2), Provenance
%
% Requires the nxr-compute plugin.  The hemisphere split requires a Structures
% atlas with left/right labels.  Unlike tess_manifold, no FreeSurfer
% registration sphere is needed: the operators act directly on the discrete
% per-hemisphere submesh and the connection Laplacian uses the intrinsic
% Levi-Civita connection (no trivial gauge / FS-pole singularities).
%
% SEE ALSO: tess_manifold, tess_eigen, tess_hemisplit, db_add_operator

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
    Tau            = 0.5;
    NoSave         = false;
    ForceRecompute = false;  %#ok<NASGU> % accepted for API symmetry; see help
    Interactive    = false;  % GUI: prompt Overwrite/Cancel (+ dependent-eigen cascade)
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'tau',            Tau            = varargin{i+1};
            case 'nosave',         NoSave         = logical(varargin{i+1});
            case 'forcerecompute', ForceRecompute = logical(varargin{i+1}); %#ok<NASGU>
            case 'interactive',    Interactive    = logical(varargin{i+1});
            otherwise
                error('tess_operators:badOption', 'Unknown option: %s', varargin{i});
        end
    end
    if ~(isscalar(Tau) && Tau>=0 && Tau<=1)
        error('tess_operators:badTau', 'Tau must be a scalar in [0,1].');
    end

    % --- map OperatorName -> Variant (case-insensitive) ---
    switch lower(strrep(OperatorName, ' ', '-'))
        case {'laplace-beltrami','lbo','laplacian'}
            Variant = 'Laplace-Beltrami';
        case {'connection-laplacian','connection'}
            Variant = 'Connection Laplacian';
        case {'dirac'}
            Variant = 'Dirac';
        case {'dirac-face','diracface'}
            Variant = 'Dirac-Face';
        case {'hodge-face','hodgeface'}
            Variant = 'Hodge-Face';
        otherwise
            error('tess_operators:badVariant', ...
                ['Unknown operator ''%s''. Valid options: ' ...
                 '''Laplace-Beltrami'', ''Connection Laplacian'', ''Dirac'', ''Dirac-Face'', ''Hodge-Face''.'], OperatorName);
    end

    % --- interactive overwrite: a matching operator already exists (GUI only) ---
    % Non-interactive callers (incl. tess_eigen's find-or-create) keep the historical
    % always-recompute behaviour. Interactive overwrite cascades to the dependent eigen
    % nodes, which reference this operator and would be orphaned by the replacement.
    if Interactive
        [sOp, ~, iSurfOp, iOp] = bst_get('OperatorFileForSurface', SurfaceFile, Variant, Tau);
        if ~isempty(iOp)
            existFile = sOp.Surface(iSurfOp).Operator(iOp).FileName;
            iDep      = bst_get('EigenFilesForOperator', SurfaceFile, existFile);
            depFiles  = {};
            depMsg    = '';
            if ~isempty(iDep)
                depFiles = arrayfun(@(k) sOp.Surface(iSurfOp).Eigen(k).FileName, iDep, 'UniformOutput', false);
                depMsg   = sprintf('\n\nThis will ALSO remove %d dependent eigenbasis node(s).', numel(iDep));
            end
            tauStr = '';
            if any(strcmpi(Variant, {'Dirac','Dirac-Face'})); tauStr = sprintf(', tau=%.3g', Tau); end
            msg = sprintf(['A %s operator already exists for this surface%s.%s\n\n' ...
                'Overwrite it (delete and recompute)?\n[No keeps and reuses the existing one.]'], ...
                Variant, tauStr, depMsg);
            if ~java_dialog('confirm', msg, 'Compute operator')
                OperatorMat = in_bst_operator(existFile);
                return;
            end
            % Overwrite: cascade-delete dependent eigen nodes first, then the operator.
            if ~isempty(depFiles); db_delete_surface_node(depFiles, 1); end
            db_delete_surface_node(existFile, 1);
        end
    end

    % --- load surface ---
    TessMat = in_tess_bst(SurfaceFile, 0);

    % --- guard: nxr-compute plugin ---
    [isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
    if ~isOk
        error('tess_operators:nxrUnavailable', ...
            'tess_operators requires nxr-compute: %s', errMsg);
    end

    % NOTE: No registration-sphere guard. The LBO/connection/Dirac operators act
    % directly on the discrete per-hemisphere submesh; the connection Laplacian
    % builds the intrinsic Levi-Civita connection internally, so no trivial-gauge
    % or FreeSurfer-pole singularity configuration is required (verified in SP2
    % Step 0: bit-identical with/without a gauge/facets call). This is the
    % deliberate difference from tess_manifold, which does call nxr 'facets'.

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
        error('tess_operators:noHemisphereLabels', ...
            ['Surface has no Structures atlas with left/right hemisphere labels ' ...
             '(required for the atlas-based hemisphere split; the geometric fallback is not allowed).']);
    end

    % --- hemisphere split (atlas-based, never conncomp) ---
    [rH, lH, isConn] = tess_hemisplit(TessMat);
    if isConn
        error('tess_operators:connectedHemispheres', ...
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

    prov = struct('Backend','nxr', 'NxrVersion',nxrVer, 'Variant',Variant, ...
                  'ComputeDate',datestr(now,'yyyy-mm-dd HH:MM:SS'));
    if ismember(Variant, {'Dirac','Dirac-Face'})
        prov.Tau = Tau;
    end

    Operator       = cell(1, 2);
    Mass           = cell(1, 2);
    GlobalVertices = cell(1, 2);
    GlobalFaces    = cell(1, 2);   % face-domain variants (e.g. 'Dirac-Face')
    diracScales    = cell(1, 2);   % [sL sE] per hemisphere (Dirac co-normalization)
    FirstOrderInt  = cell(1, 2);   % Dirac: intrinsic first-order D_int [4F x 4V] per hemisphere
    FirstOrderExt  = cell(1, 2);   % Dirac: extrinsic first-order D     [4F x 4V] per hemisphere
    FaceMass       = cell(1, 2);   % Dirac: face-area mass W_F [4F x 4F] per hemisphere
    FaceAux        = cell(1, 2);   % Hodge-Face: struct(ScalarMass,GradFace,FaceNormal) for the Hodge lift

    for hh = 1:2
        vH = hemis{hh};
        if isempty(vH)
            error('tess_operators:emptyHemisphere', ...
                'Hemisphere %s has no vertices.', tags{hh});
        end

        % build local submesh (local indices)
        isV   = false(nVtot, 1);  isV(vH) = true;
        fMask = all(isV(Fcs), 2);
        fGlob = find(fMask);                      % global face indices for this hemisphere
        mapV  = zeros(nVtot, 1);  mapV(vH) = 1:numel(vH);
        Vloc  = Vtx(vH, :);
        Floc  = mapV(Fcs(fMask, :));

        % create nxr submesh (validated; never a hand-built mesh)
        h = nxr_compute('create', Vloc, Floc);   % create validates the mesh (clean error, no segfault)
        try
            switch Variant
                case 'Laplace-Beltrami'
                    A = nxr_compute('operators', h, 'laplacian', 'cotan');
                    B = nxr_compute('operators', h, 'mass', 'galerkin');
                case 'Connection Laplacian'
                    % Hermitian [nVh x nVh]; the Levi-Civita connection is
                    % built internally from the discrete geometry (no separate
                    % gauge/facets call required, verified on the canonical cortex).
                    A = nxr_compute('operators', h, 'laplacian', 'connection');
                    B = nxr_compute('operators', h, 'mass', 'galerkin');
                case 'Dirac'
                    % Full frame-transport Dirac operator ("f and N together"):
                    % (1-Tau)*D_int^2 + Tau*E. BOTH blocks couple the quaternion
                    % components (the rotating cortical frame), so neither presupposes a
                    % tangent/normal split:
                    %   D_int^2 -- INTRINSIC Dirac squared (built from the immersion f,
                    %     edge vectors): the spin-connection / tangent-frame transport
                    %     (its scalar part is the cotan Laplacian, plus a genuine
                    %     quaternionic coupling that cotanL(x)I4 alone discards).
                    %   E       -- EXTRINSIC (relative) Dirac squared (built from the
                    %     Gauss map N): the shape operator / normal tilt.
                    % Each carries 1/length^k with different k, so normalize each to unit
                    % largest generalized eigenvalue vs the mass B; Tau is then a true
                    % dimensionless dial (Tau=0.5 == equal intrinsic/extrinsic frame
                    % weight), portable across mesh size / units.
                    L4 = local_dirac_intrinsic_sq(Vloc, Floc);            % intrinsic Dirac^2 (immersion f)
                    E  = nxr_compute('operators', h, 'dirac', 1);         % extrinsic Dirac^2 (Gauss map N)
                    Mg = nxr_compute('operators', h, 'mass', 'galerkin'); % [nVh x nVh]
                    B  = kron(Mg, speye(4));                              % [4nVh x 4nVh]
                    sL = local_lambda_max(L4, B);
                    sE = local_lambda_max(E,  B);
                    A  = (1 - Tau) * (L4 / sL) + Tau * (E / sE);          % [4nVh x 4nVh]
                    diracScales{hh} = [sL, sE];   % record for provenance / eigenvalue scale

                    % First-order Dirac ROOTS (nxr), stored for the signed/propagation
                    % branch: D_int (immersion/edge, root of L4) and D_ext (Gauss map,
                    % root of E). Each [4F x 4V]; their W_F-Galerkin square reproduces
                    % L4 / E (verified). W_F = kron(face-area, I4) pairs with them.
                    Dint1 = nxr_compute('operators', h, 'diracIntrinsicD');   % intrinsic first-order
                    Dext1 = nxr_compute('operators', h, 'diracD');            % extrinsic first-order
                    nFh   = size(Floc, 1);
                    e1f   = Vloc(Floc(:,2),:) - Vloc(Floc(:,1),:);
                    e2f   = Vloc(Floc(:,3),:) - Vloc(Floc(:,1),:);
                    fArea = 0.5 * sqrt(sum(cross(e1f, e2f, 2).^2, 2));
                    FirstOrderInt{hh} = Dint1;
                    FirstOrderExt{hh} = Dext1;
                    FaceMass{hh}      = kron(spdiags(fArea, 0, nFh, nFh), speye(4));   % [4F x 4F]

                case 'Dirac-Face'
                    % Face-domain Dirac operator (the dual of the vertex 'Dirac'):
                    %   (1-Tau)*E~_int + Tau*E~_ext  [4F x 4F], mode mass = W_F.
                    %   E~_int = D~_int' W_V D~_int  -- intrinsic (centroid-immersion dual),
                    %            assembled in MATLAB (mirrors local_dirac_intrinsic_sq on the
                    %            vertex side); its scalar block is the face cotan Laplacian.
                    %   E~_ext = extrinsicBlockFace = D~_ext' W_V D~_ext  (nxr, Gauss map).
                    % BOTH use the SAME vertex dual-area mass W_V (so E~_int is built exactly
                    % as nxr builds E~_ext); co-normalized vs W_F so Tau is dimensionless.
                    Dt_int = nxr_compute('operators', h, 'diracFaceIntrinsicD');   % [4V x 4F]
                    Mlump  = nxr_compute('operators', h, 'mass', 'lumped');        % vertex dual area [nVh x nVh]
                    WV     = kron(Mlump, speye(4));                                % [4V x 4V]
                    Eint   = Dt_int' * WV * Dt_int;  Eint = (Eint + Eint')/2;      % [4F x 4F]
                    Eext   = nxr_compute('operators', h, 'diracFace', 1);          % extrinsicBlockFace [4F x 4F]
                    nFh    = size(Floc, 1);
                    e1f    = Vloc(Floc(:,2),:) - Vloc(Floc(:,1),:);
                    e2f    = Vloc(Floc(:,3),:) - Vloc(Floc(:,1),:);
                    fArea  = 0.5 * sqrt(sum(cross(e1f, e2f, 2).^2, 2));
                    B      = kron(spdiags(fArea, 0, nFh, nFh), speye(4));          % W_F [4F x 4F]
                    sI = local_lambda_max(Eint, B);
                    sX = local_lambda_max(Eext, B);
                    A  = (1 - Tau) * (Eint / sI) + Tau * (Eext / sX);             % [4F x 4F]
                    diracScales{hh} = [sI, sX];

                case 'Hodge-Face'
                    % Face Hodge vector eigenbasis support: the SCALAR face Laplacian
                    % lapFace = gradFace' W_F gradFace [F x F] (full-rank, tall gradFace
                    % root -> smooth Weyl spectrum), its scalar eigenmodes lifted to vectors
                    % by tess_eigen. The node carries the scalar pencil (Operator=lapFace,
                    % FaceAux.ScalarMass=M_F) for the eigensolve, the lift operators
                    % (gradFace, face normals) in FaceAux, and the 4-component face mass
                    % Mass=W_F [4F x 4F] consumed by bst_dirac as B.
                    Gf = nxr_compute('operators', h, 'gradFace');                 % [3F x F]
                    Kf = nxr_compute('operators', h, 'lapFace');  Kf = (Kf+Kf')/2; % [F x F]
                    nFh = size(Floc, 1);
                    e1f = Vloc(Floc(:,2),:) - Vloc(Floc(:,1),:);
                    e2f = Vloc(Floc(:,3),:) - Vloc(Floc(:,1),:);
                    Nf  = cross(e1f, e2f, 2);  twoA = sqrt(sum(Nf.^2,2));
                    Nf  = Nf ./ max(twoA, eps);  fArea = twoA/2;
                    vn  = TessMat.VertNormals(vH(Floc(:,1)),:);                    % outward orientation
                    flip = sum(Nf.*vn,2) < 0;  Nf(flip,:) = -Nf(flip,:);
                    A = Kf;                                                        % scalar lapFace [F x F]
                    B = kron(spdiags(fArea,0,nFh,nFh), speye(4));                  % W_F [4F x 4F] (for bst_dirac)
                    FaceAux{hh} = struct('ScalarMass', spdiags(fArea,0,nFh,nFh), ...% M_F [F x F]
                                         'GradFace', Gf, 'FaceNormal', Nf);
            end
        catch ME
            nxr_compute('destroy', h);
            rethrow(ME);
        end
        nxr_compute('destroy', h);

        Operator{hh}       = A;
        Mass{hh}           = B;
        GlobalVertices{hh} = vH;
        GlobalFaces{hh}    = fGlob;
    end

    % Record the Dirac block co-normalization (makes Tau dimensionless / portable and
    % lets a consumer recover the unnormalized scale of the stored eigenvalues).
    if strcmpi(Variant, 'Dirac')
        prov.Blocks        = '(1-Tau)*intrinsic_Dirac^2 (f) + Tau*extrinsic_Dirac^2 (N)';
        prov.Normalization = 'lambda_max-vs-B (per-block, co-normalized)';
        prov.DiracScale    = diracScales;   % 1x2 cell, [sL sE] per hemisphere
    elseif strcmpi(Variant, 'Dirac-Face')
        prov.Blocks        = '(1-Tau)*intrinsic_faceDirac^2 (D~_int'' W_V D~_int) + Tau*extrinsic_faceDirac^2 (extrinsicBlockFace)';
        prov.Normalization = 'lambda_max-vs-W_F (per-block, co-normalized)';
        prov.DiracScale    = diracScales;   % 1x2 cell, [sI sX] per hemisphere
    end

    % --- assemble OperatorMat ---
    OperatorMat                = db_template('operatormat');
    OperatorMat.Variant        = Variant;
    OperatorMat.ParentSurface  = SurfaceFile;
    OperatorMat.Operator       = Operator;        % 1x2 cell of sparse matrices
    OperatorMat.Mass           = Mass;            % 1x2 cell of sparse matrices
    OperatorMat.GlobalVertices = GlobalVertices;  % 1x2 cell of global vertex indices
    OperatorMat.GlobalFaces    = GlobalFaces;     % 1x2 cell of global face indices (face-domain variants)
    if strcmpi(Variant, 'Hodge-Face')
        OperatorMat.FaceAux = FaceAux;            % 1x2 struct: ScalarMass / GradFace / FaceNormal (Hodge lift)
    end
    if strcmpi(Variant, 'Dirac')
        OperatorMat.FirstOrder = struct('Intrinsic', {FirstOrderInt}, 'Extrinsic', {FirstOrderExt});
        OperatorMat.FaceMass   = FaceMass;         % 1x2 cell of W_F [4F x 4F]
    end
    OperatorMat.Provenance     = prov;

    % --- save / register in DB ---
    if ~NoSave
        [sSubjectSave, iSubjectSave] = bst_get('SurfaceFile', SurfaceFile);
        if isempty(sSubjectSave)
            error('tess_operators:subjectNotFound', ...
                'Could not resolve subject for surface: %s', SurfaceFile);
        end
        Comment = sprintf('%s operator', Variant);
        db_add_operator(iSubjectSave, SurfaceFile, OperatorMat, Comment);
    end
end

% ----------------------------------------------------------------------------
function lmax = local_lambda_max(A, B)
% Largest generalized eigenvalue of a symmetric/PSD pencil (A, B), B SPD. Used to
% co-normalize the Dirac blocks. Factorization-free: 'largestabs' forms B^{-1}A and
% factorizes only the well-conditioned mass B (never the possibly-singular A), so it
% cannot trip the ill-conditioning warning. A coarse tolerance suffices -- only the
% scale (one figure) matters for the normalization.
    A = (A + A') / 2;
    B = (B + B') / 2;
    opts = struct('tol', 1e-4, 'maxit', 300, 'disp', 0);
    lmax = abs(eigs(A, B, 1, 'largestabs', opts));
    if ~isfinite(lmax) || (lmax <= 0)
        error('tess_operators:badBlockScale', ...
            'Could not estimate a positive largest eigenvalue for Dirac block normalization.');
    end
end

% ----------------------------------------------------------------------------
function L = local_dirac_intrinsic_sq(V, F)
% Intrinsic (immersion/edge-based) quaternionic Dirac operator squared,
% L = D_int' * MF * D_int  [4nV x 4nV], where MF is the face-area 2-form mass.
% Ported verbatim from gptoolbox dirac_operator.m (Crane et al., "Spin
% Transformations of Discrete Surfaces"): each 4x4 block is the quaternion
% left-multiplication by the opposite EDGE VECTOR (the immersion f) over twice the
% area. Unlike cotanL(x)I4, its square couples the quaternion components -- it
% carries the intrinsic spin-connection (tangent-frame transport); its scalar part
% equals the cotan Laplacian. Pairs with the extrinsic (Gauss-map) Dirac to give
% the full rotating-frame operator.
    nF = size(F, 1);  nV = size(V, 1);
    e1 = V(F(:,2),:) - V(F(:,1),:);
    e2 = V(F(:,3),:) - V(F(:,1),:);
    dblA = sqrt(sum(cross(e1, e2, 2).^2, 2));               % doublearea [nF x 1]
    EV = [zeros(numel(F),1), V(F(:,[2 3 1]),:) - V(F(:,[3 1 2]),:)];   % edge vectors (Im quaternion)
    Q = [1 0 0 0;0 -1  0 0;0 0 -1  0;0  0 0 -1; ...
         0 1 0 0;1  0  0 0;0 0  0 -1;0  0 1  0; ...
         0 0 1 0;0  0  0 1;1 0  0  0;0 -1 0  0; ...
         0 0 0 1;0  0 -1 0;0 1  0  0;1  0 0  0]';
    II = repmat(repmat((0:nF-1)'*4 + (1:4), 1, 4), 3, 1);
    JJ = (repmat(F(:),1,16)-1)*4 + reshape(repmat(1:4,4,1),1,[]);
    D  = sparse(II, JJ, -EV*Q ./ [dblA;dblA;dblA], 4*nF, 4*nV);
    MF = kron(spdiags(dblA/2, 0, nF, nF), speye(4));        % face-area 2-form mass
    L  = D' * MF * D;
    L  = (L + L') / 2;                                      % symmetrize (1e-16 noise)
end
