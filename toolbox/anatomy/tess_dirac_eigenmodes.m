function DiracEigen = tess_dirac_eigenmodes(SurfaceFile, varargin)
% TESS_DIRAC_EIGENMODES: Per-hemisphere Dirac eigenbasis stored on the surface.
%
% USAGE:  DiracEigen = tess_dirac_eigenmodes(SurfaceFile)
%         DiracEigen = tess_dirac_eigenmodes(SurfaceFile, 'Tau',0.5, 'K',400, ...
%                                            'NoSave',1, 'ForceRecompute',1)
%
% DESCRIPTION:
%     Splits the cortex by hemisphere (tess_hemisplit, atlas L/R; never conncomp)
%     and, per hemisphere, solves the generalized Dirac eigenproblem
%         L(Tau) * phi = lambda * B * phi,   B = kron(Mg, I4)
%     where L(Tau)=operators(h,'dirac',Tau) is the [4Vh x 4Vh] relative-Dirac
%     family and Mg=operators(h,'mass','galerkin') is the vertex Galerkin mass.
%     The K smallest-magnitude eigenpairs are B-orthonormalized and stored as
%     TessMat.DiracEigen, a 1x2 per-hemisphere struct array ((1)=L,(2)=R), each:
%       .Vectors        [4Vh x K]  B-orthonormal (vertex-interleaved 4v+c, [w,x,y,z])
%       .Values         [K x 1]    ascending, >= 0
%       .Mass           [Vh x Vh]  galerkin vertex mass (so B = kron(Mass,I4))
%       .nModes, .Order, .Tau, .GlobalVertices, .Hemisphere, .Provenance
%
%     Eigenvalues come in 4-fold quaternionic multiplets; only the B-orthonormal
%     subspace matters (no multiplet canonicalization).
%
% OPTIONS:  'Tau' (0.5) | 'K' (400) | 'NoSave' (false) | 'ForceRecompute' (false)
%
% Requires the nxr-compute plugin (operators 'dirac'/'mass').
%
% SEE ALSO: tess_eigenmodes, tess_frame, tess_hemisplit

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
    Tau=0.5; K=400; NoSave=false; ForceRecompute=false;
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'tau',            Tau=varargin{i+1};
            case 'k',              K=varargin{i+1};
            case 'nosave',         NoSave=varargin{i+1};
            case 'forcerecompute', ForceRecompute=varargin{i+1};
        end
    end
    if ~(isscalar(Tau) && Tau>=0 && Tau<=1)
        error('tess_dirac_eigenmodes:badTau', 'Tau must be a scalar in [0,1].');
    end

    TessFile = file_fullpath(SurfaceFile);
    % VertConn not requested (0): the Structures-atlas L/R guard below ensures
    % tess_hemisplit returns via its atlas path and never reaches the VertConn fallback.
    TessMat  = in_tess_bst(SurfaceFile, 0);

    % --- cache return: both hemispheres must match Tau AND nModes ---
    if ~ForceRecompute && isfield(TessMat,'DiracEigen') && ~isempty(TessMat.DiracEigen) ...
       && isequal(size(TessMat.DiracEigen),[1 2]) ...
       && isequal([TessMat.DiracEigen.Tau], [Tau Tau]) ...
       && all([TessMat.DiracEigen.nModes] == K)
        DiracEigen = TessMat.DiracEigen;
        return;
    end

    % --- require nxr-compute ---
    [isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
    if ~isOk
        error('tess_dirac_eigenmodes:nxrUnavailable', 'requires nxr-compute: %s', errMsg);
    end

    % --- require Structures atlas with L/R labels (atlas split, never conncomp) ---
    hasLabels = false;
    if isfield(TessMat,'Atlas') && ~isempty(TessMat.Atlas)
        iStruct = find(strcmpi({TessMat.Atlas.Name}, 'Structures'), 1);
        if ~isempty(iStruct) && ~isempty(TessMat.Atlas(iStruct).Scouts)
            scouts = TessMat.Atlas(iStruct).Scouts;
            labels = {scouts.Label}; regions = {scouts.Region};
            reg1 = cellfun(@(c) c(1), regions(~cellfun(@isempty, regions)), 'UniformOutput', false);
            hasLabels = (any(strcmpi(labels,'lh')) || any(strcmpi(reg1,'L'))) && ...
                        (any(strcmpi(labels,'rh')) || any(strcmpi(reg1,'R')));
        end
    end
    if ~hasLabels
        error('tess_dirac_eigenmodes:noHemisphereLabels', ...
            'Surface has no Structures atlas with left/right hemisphere labels.');
    end

    [rH, lH, isConn] = tess_hemisplit(TessMat);
    if isConn
        error('tess_dirac_eigenmodes:connectedHemispheres', ...
            'Cannot split the surface into two disconnected hemispheres (they share vertices/faces); check the Structures atlas and mesh.');
    end
    hemis = {lH(:), rH(:)}; tags = {'L','R'};
    Vtx = double(TessMat.Vertices); Fcs = double(TessMat.Faces); nVtot = size(Vtx,1);

    nxrVer = '';
    try, nxrVer = nxr_compute('version'); catch, end   %#ok<CTCH>
    prov = struct('Backend','nxr', 'Package','dirac-eigenmodes', 'NxrVersion',nxrVer, ...
                  'Tau',Tau, 'K',K, 'ComputeDate',datestr(now,'yyyy-mm-dd HH:MM:SS'));

    clear DE
    for hh = 1:2
        vH = hemis{hh};
        if isempty(vH)
            error('tess_dirac_eigenmodes:emptyHemisphere', 'Hemisphere %s has no vertices.', tags{hh});
        end
        isV = false(nVtot,1); isV(vH) = true;
        fMask = all(isV(Fcs), 2);
        map = zeros(nVtot,1); map(vH) = 1:numel(vH);
        Vloc = Vtx(vH,:);
        Floc = map(Fcs(fMask,:));
        nVh  = numel(vH);

        h  = nxr_compute('create', Vloc, Floc);
        L  = nxr_compute('operators', h, 'dirac', Tau);     % [4nVh x 4nVh]
        Mg = nxr_compute('operators', h, 'mass', 'galerkin'); % [nVh x nVh]
        nxr_compute('destroy', h);

        B = kron(Mg, speye(4));                              % [4nVh x 4nVh]
        if K > 4*nVh - 2
            error('tess_dirac_eigenmodes:tooManyModes', ...
                'K=%d exceeds 4*nV-2=%d on hemisphere %s.', K, 4*nVh-2, tags{hh});
        end
        % over-fetch generously: the 4-fold multiplets make eigs return a
        % rank-deficient spanning set (~12-15% deficient), so K+8 is not enough
        % for large K; ~30% over-fetch keeps rank >= K after the rank drop.
        % When K is near 4*nVh-2 the cap shrinks the over-fetch and
        % local_ritz_basis may then throw :rankDeficient (caught, not silent).
        nRequest = min(K + max(8, ceil(0.3*K)), 4*nVh - 2);
        opts = struct('tol', 1e-6, 'maxit', 1000, 'disp', 0);
        [Uraw, ~] = eigs(L, B, nRequest, 'smallestabs', opts);
        Uraw = real(Uraw);
        [U, lam] = local_ritz_basis(L, B, Uraw, K);   % matched (Vectors, Values)

        s = struct();
        s.Vectors        = U;
        s.Values         = lam;
        s.Mass           = Mg;
        s.nModes         = K;
        s.Order          = (1:K)';
        s.Tau            = Tau;
        s.GlobalVertices = vH;
        s.Hemisphere     = tags{hh};
        s.Provenance     = prov;
        DE(hh) = s;   %#ok<AGROW>
    end

    DiracEigen = DE;

    if ~NoSave
        TessMat_full = load(TessFile);
        TessMat_full.DiracEigen = DiracEigen;
        TessMat_full = bst_history('add', TessMat_full, 'dirac-eigenmodes', ...
            sprintf('Stored Dirac eigenbasis (tau=%.3g, K=%d) as per-hemisphere DiracEigen.', Tau, K));
        bst_save(TessFile, TessMat_full, 'v7');
    end
end

% ----------------------------------------------------------------------------
function [Phi, mu] = local_ritz_basis(L, B, U, K)
% Rank-revealing B-orthonormalization of the eigs spanning set U, then a
% Rayleigh-Ritz diagonalization of L on that span, returning K genuine
% B-orthonormal eigenpairs (ascending), matched (Phi, mu). No randomness.
%
% Why: eigs(L,B,'smallestabs') returns linearly-dependent vectors within the
% 4-fold quaternionic multiplets, so the raw Gram U'BU is rank-deficient. We
% keep only the well-conditioned directions, build a B-orthonormal basis W of
% span(U), then diagonalize the small Lr = W'LW so every output column is a
% true eigenvector with its matched Ritz value.
    G = U' * (B * U); G = (G + G')/2;
    [Vg, Dg] = eig(full(G));
    dg = real(diag(Dg));
    keep = dg > 1e-10 * max(dg);                 % drop multiplet-induced rank deficiency
    W = U * (Vg(:, keep) * diag(1 ./ sqrt(dg(keep))));   % W'BW = I (B-orthonormal span)
    Lr = W' * (L * W); Lr = (Lr + Lr')/2;        % Rayleigh-Ritz on span(W)
    [Vr, Dr] = eig(full(Lr));
    mu = real(diag(Dr));
    [mu, idx] = sort(mu, 'ascend');
    Phi = W * Vr(:, idx);                         % genuine eigenvectors, B-orthonormal
    if size(Phi, 2) < K
        error('tess_dirac_eigenmodes:rankDeficient', ...
            'Only %d independent modes recovered (< K=%d); increase the eigs over-fetch.', ...
            size(Phi, 2), K);
    end
    Phi = Phi(:, 1:K);
    mu  = mu(1:K);
    negThresh = -1e-6 * max(abs(mu));
    if any(mu < negThresh)
        warning('tess_dirac_eigenmodes:nonPSD', ...
            ['Dirac eigenproblem returned %d eigenvalue(s) below %.2e (most negative %.3e); ' ...
             'clamping to 0. The nxr mass/operator may be non-PSD.'], ...
            sum(mu < negThresh), negThresh, min(mu));
    end
    mu(mu < 0) = 0;
end
