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
% SEE ALSO: tess_manifold, tess_dirac_eigenmodes, tess_hemisplit, db_add_operator

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
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'tau',            Tau            = varargin{i+1};
            case 'nosave',         NoSave         = logical(varargin{i+1});
            case 'forcerecompute', ForceRecompute = logical(varargin{i+1}); %#ok<NASGU>
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
        otherwise
            error('tess_operators:badVariant', ...
                ['Unknown operator ''%s''. Valid options: ' ...
                 '''Laplace-Beltrami'', ''Connection Laplacian'', ''Dirac''.'], OperatorName);
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
    if strcmpi(Variant, 'Dirac')
        prov.Tau = Tau;
    end

    Operator       = cell(1, 2);
    Mass           = cell(1, 2);
    GlobalVertices = cell(1, 2);
    diracScales    = cell(1, 2);   % [sL sE] per hemisphere (Dirac co-normalization)

    for hh = 1:2
        vH = hemis{hh};
        if isempty(vH)
            error('tess_operators:emptyHemisphere', ...
                'Hemisphere %s has no vertices.', tags{hh});
        end

        % build local submesh (local indices)
        isV   = false(nVtot, 1);  isV(vH) = true;
        fMask = all(isV(Fcs), 2);
        mapV  = zeros(nVtot, 1);  mapV(vH) = 1:numel(vH);
        Vloc  = Vtx(vH, :);
        Floc  = mapV(Fcs(fMask, :));

        % create nxr submesh (validated; never a hand-built mesh)
        h = nxr_safe_create(Vloc, Floc);
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
                    % Co-normalized relative-Dirac family. The intrinsic block
                    % (cotanL ⊗ I4) is dimensionless while the extrinsic curvature
                    % block E carries 1/length^2, so a raw (1-Tau)*cotanL + Tau*E mixes
                    % blocks whose generalized spectra scale as 1/s^2 vs 1/s^4 -- Tau is
                    % then unit-dependent (in meters, Tau=0.5 is ~99.999% extrinsic).
                    % Normalize each block to unit largest generalized eigenvalue vs the
                    % mass B so Tau is a true dimensionless dial: Tau=0.5 == equal spectral
                    % weight, portable across mesh size / units. (Single-scalar block
                    % normalization is eigenvector-safe within each block; it only sets
                    % their relative weighting in the combination.)
                    L4 = nxr_compute('operators', h, 'dirac', 0);         % (cotanL ⊗ I4), intrinsic
                    E  = nxr_compute('operators', h, 'dirac', 1);         % extrinsic curvature block
                    Mg = nxr_compute('operators', h, 'mass', 'galerkin'); % [nVh x nVh]
                    B  = kron(Mg, speye(4));                              % [4nVh x 4nVh]
                    sL = local_lambda_max(L4, B);
                    sE = local_lambda_max(E,  B);
                    A  = (1 - Tau) * (L4 / sL) + Tau * (E / sE);          % [4nVh x 4nVh]
                    diracScales{hh} = [sL, sE];   % record for provenance / eigenvalue scale
            end
        catch ME
            nxr_compute('destroy', h);
            rethrow(ME);
        end
        nxr_compute('destroy', h);

        Operator{hh}       = A;
        Mass{hh}           = B;
        GlobalVertices{hh} = vH;
    end

    % Record the Dirac block co-normalization (makes Tau dimensionless / portable and
    % lets a consumer recover the unnormalized scale of the stored eigenvalues).
    if strcmpi(Variant, 'Dirac')
        prov.Normalization = 'lambda_max-vs-B (per-block, co-normalized)';
        prov.DiracScale    = diracScales;   % 1x2 cell, [sL sE] per hemisphere
    end

    % --- assemble OperatorMat ---
    OperatorMat                = db_template('operatormat');
    OperatorMat.Variant        = Variant;
    OperatorMat.ParentSurface  = SurfaceFile;
    OperatorMat.Operator       = Operator;        % 1x2 cell of sparse matrices
    OperatorMat.Mass           = Mass;            % 1x2 cell of sparse matrices
    OperatorMat.GlobalVertices = GlobalVertices;  % 1x2 cell of global vertex indices
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
