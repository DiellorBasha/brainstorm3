function EigenMat = tess_eigen(SurfaceFile, OperatorName, varargin)
% TESS_EIGEN: Assemble a per-hemisphere eigenbasis as an eigen_ DB node.
%
% USAGE:  EigenMat = tess_eigen(SurfaceFile, OperatorName)
%         EigenMat = tess_eigen(SurfaceFile, OperatorName, ...
%                               'K',400, 'Tau',0.5, 'NoSave',false, 'ForceRecompute',false)
%
% DESCRIPTION:
%     Finds (or creates) the operator node of the requested variant under the
%     parent surface, loads its per-hemisphere operator pencil (A, B), and
%     solves the generalized eigenproblem
%         A * phi = lambda * B * phi
%     per hemisphere for the K smallest-magnitude eigenpairs. The pairs are
%     B-orthonormalized (Rayleigh-Ritz for Dirac, to handle the 4-fold
%     quaternionic multiplets; standard column normalization otherwise) and the
%     result is assembled into an EigenMat structure (db_template('eigenmat')).
%
%     Three operator variants are supported (case-insensitive OperatorName):
%       'Laplace-Beltrami'     A = laplacian/cotan      [nVh x nVh]
%       'Connection Laplacian' A = laplacian/connection [nVh x nVh]
%       'Dirac'                A = dirac(Tau)           [4nVh x 4nVh]
%
%     The operator node is found-or-created via tess_operators: if the surface
%     already has an Operator child of the requested Variant it is reused;
%     otherwise tess_operators(SurfaceFile, OperatorName, 'Tau', Tau) is called
%     to create one. tess_operators is the ONLY path that reaches
%     nxr_compute('create') (transitively, via nxr_safe_create); tess_eigen
%     never builds a mesh itself.
%
%     Unless NoSave is true, the result is saved as an eigen_*.mat file alongside
%     the parent surface and registered in the Brainstorm DB via db_add_eigen
%     (creating a child node under the surface).
%
% OPTIONS:
%     'K'              : number of modes to keep per hemisphere (default 400)
%     'Tau'            : scalar in [0,1] (default 0.5) — relative-Dirac mixing
%                        parameter, forwarded to operator creation (Dirac only)
%     'NoSave'         : true/false (default false) — compute but do not write
%                        to disk or register in the DB
%     'ForceRecompute' : true/false (default false) — when false, an existing
%                        eigen_*.mat child of the surface whose Variant matches,
%                        whose stored K is at least the requested K, and (for Dirac)
%                        whose Tau matches, is loaded and truncated to K modes
%                        instead of re-solving. Set true to force a fresh eigensolve.
%
% OUTPUT:
%     EigenMat : struct matching db_template('eigenmat'), with fields:
%                Comment, ParentSurface, OperatorFile, Variant, Phi(1x2),
%                Lambda(1x2), K, GlobalVertices(1x2), Provenance
%
% Requires the nxr-compute plugin (transitively, via tess_operators).
%
% SEE ALSO: tess_operators, bst_dirac, db_add_eigen

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
    K              = 400;
    Tau            = 0.5;
    NoSave         = false;
    ForceRecompute = false;  % when false, reuse a cached eigen node; see help
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'k',              K              = varargin{i+1};
            case 'tau',            Tau            = varargin{i+1};
            case 'nosave',         NoSave         = logical(varargin{i+1});
            case 'forcerecompute', ForceRecompute = logical(varargin{i+1});
            otherwise
                error('tess_eigen:badOption', 'Unknown option: %s', varargin{i});
        end
    end
    if ~(isscalar(K) && K == round(K) && K >= 1)
        error('tess_eigen:badK', 'K must be a positive integer.');
    end
    if ~(isscalar(Tau) && Tau>=0 && Tau<=1)
        error('tess_eigen:badTau', 'Tau must be a scalar in [0,1].');
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
        otherwise
            error('tess_eigen:badVariant', ...
                ['Unknown operator ''%s''. Valid options: ' ...
                 '''Laplace-Beltrami'', ''Connection Laplacian'', ''Dirac'', ''Dirac-Face''.'], OperatorName);
    end
    isFace  = strcmpi(Variant, 'Dirac-Face');               % face-domain (modes on faces)
    isDirac = strcmpi(Variant, 'Dirac') || isFace;          % quaternion Dirac-type: over-fetch + Rayleigh-Ritz + Tau
    isLBO   = strcmpi(Variant, 'Laplace-Beltrami');

    % --- guard: nxr-compute plugin (operators reach nxr transitively) ---
    [isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
    if ~isOk
        error('tess_eigen:nxrUnavailable', ...
            'tess_eigen requires nxr-compute: %s', errMsg);
    end

    % --- resolve subject / surface for the parent ---
    [sSubject, iSubject, iSurface] = bst_get('SurfaceFile', SurfaceFile);
    if isempty(sSubject) || isempty(iSurface)
        error('tess_eigen:surfaceNotFound', ...
            'Could not resolve surface in current protocol: %s', SurfaceFile);
    end

    % --- find-or-load a cached eigen node before the (expensive) eigensolve ---
    %     Reuse an existing eigen_*.mat child of this surface whose Variant matches,
    %     whose stored K >= the requested K, and (for Dirac) whose Tau matches. The
    %     eigenbasis is nested (ascending eigenvalues, B-orthonormal), so a node with
    %     more modes is truncated to its first K columns. ForceRecompute skips this.
    if ~ForceRecompute
        cached = local_find_eigen(sSubject, iSurface, Variant, K, Tau, isDirac);
        if ~isempty(cached)
            EigenMat = local_truncate_eigen(cached, K);
            return;
        end
    end

    % --- progress feedback: indeterminate bar so the user can see the
    %     (potentially long) eigensolve is running. Only open our own bar if
    %     one is not already visible (so we don't disrupt a parent process'
    %     bar). onCleanup guarantees it stops on normal return AND on error. ---
    eigenBar = [];
    if ~bst_progress('isVisible')
        bst_progress('start', 'Eigenmodes', sprintf('Computing %s eigenmodes...', Variant));
        eigenBar = onCleanup(@() bst_progress('stop'));  %#ok<NASGU>
    end

    % --- find-or-create the operator node of the requested Variant ---
    OperatorFile = local_find_operator(sSubject, iSurface, Variant);
    if isempty(OperatorFile)
        % Create it (tess_operators -> nxr_safe_create -> nxr_compute('create')).
        bst_progress('text', sprintf('Computing %s operator...', Variant));
        tess_operators(SurfaceFile, OperatorName, 'Tau', Tau);
        % Re-read the surface's Operator list from a fresh subject struct.
        sSubject = bst_get('Subject', iSubject);
        OperatorFile = local_find_operator(sSubject, iSurface, Variant);
        if isempty(OperatorFile)
            error('tess_eigen:operatorCreateFailed', ...
                'tess_operators did not produce a ''%s'' operator node under %s.', ...
                Variant, SurfaceFile);
        end
    end

    % --- verify the operator node resolves in the DB ---
    [sOpSubject, ~, ~, iOperator] = bst_get('OperatorFile', OperatorFile);
    if isempty(sOpSubject) || isempty(iOperator)
        error('tess_eigen:operatorUnresolved', ...
            'Operator node does not resolve via bst_get(''OperatorFile''): %s', OperatorFile);
    end

    % --- load operator pencil (per-hemisphere cells) ---
    Op = load(file_fullpath(OperatorFile));
    if ~isfield(Op,'Operator') || ~isfield(Op,'Mass') || ~isfield(Op,'GlobalVertices')
        error('tess_eigen:badOperatorFile', ...
            'Operator file is missing Operator/Mass/GlobalVertices: %s', OperatorFile);
    end

    % --- nxr version for provenance ---
    nxrVer = '';
    try, nxrVer = nxr_compute('version'); catch, end  %#ok<CTCH>

    prov = struct('Backend','nxr', 'NxrVersion',nxrVer, 'Variant',Variant, 'K',K);
    if isDirac
        prov.Tau = Tau;
    end
    if isLBO
        prov.Ortho = 'B-orthonormal';
    else
        % Connection Laplacian and Dirac both use Rayleigh-Ritz to resolve their
        % degenerate multiplets into genuine B-orthonormal eigenpairs.
        prov.Ortho = 'Rayleigh-Ritz';
    end
    prov.ComputeDate = datestr(now,'yyyy-mm-dd HH:MM:SS');

    Phi            = cell(1, 2);
    Lambda         = cell(1, 2);
    GlobalVertices = cell(1, 2);
    GlobalFaces    = cell(1, 2);
    tags           = {'L','R'};
    if isFace && (~isfield(Op,'GlobalFaces') || isempty(Op.GlobalFaces) || isempty(Op.GlobalFaces{1}))
        error('tess_eigen:noGlobalFaces', ...
            'Dirac-Face operator node is missing GlobalFaces (recompute the operator).');
    end

    for hh = 1:2
        A  = Op.Operator{hh};
        B  = Op.Mass{hh};
        gv = Op.GlobalVertices{hh};
        n  = size(A, 1);

        if isDirac
            % quaternion domain count: faces for 'Dirac-Face', vertices for 'Dirac'
            if isFace, nDom = numel(Op.GlobalFaces{hh}); else, nDom = numel(gv); end
            if size(A,1) ~= 4*nDom
                error('tess_eigen:diracSizeMismatch', ...
                    'Dirac operator on hemisphere %s is %dx%d, expected 4*nDom=%d.', ...
                    tags{hh}, size(A,1), size(A,2), 4*nDom);
            end
            nVh = nDom;   % reuse the over-fetch / size logic below
            if K > 4*nVh - 2
                error('tess_eigen:tooManyModes', ...
                    'K=%d exceeds 4*nV-2=%d on hemisphere %s.', K, 4*nVh-2, tags{hh});
            end
            % Over-fetch generously: the 4-fold quaternionic multiplets make eigs
            % return a rank-deficient spanning set; ~30%% over-fetch keeps rank>=K.
            nRequest = min(K + max(8, ceil(0.3*K)), 4*nVh - 2);
        else
            if K > n - 2
                error('tess_eigen:tooManyModes', ...
                    'K=%d exceeds nV-2=%d on hemisphere %s.', K, n-2, tags{hh});
            end
            nRequest = min(K + 8, n - 2);
        end

        hemiName = 'left'; if hh == 2, hemiName = 'right'; end
        bst_progress('text', sprintf('Eigensolve: %s hemisphere (%s, K=%d)...', hemiName, Variant, K));
        opts = struct('tol', 1e-6, 'maxit', 1000, 'disp', 0);
        % Robust smallest-eigenpair solve: A has a near-zero kernel (constant /
        % constant-quaternion modes), so eigs(...,'smallestabs') would shift-invert at
        % sigma=0 and factorize the singular A (RCOND ~ 1e-18 warning). bst_eigs_smallest
        % symmetrizes A,B (forces the Lanczos path) and uses a small negative sigma shift.
        [V, D] = bst_eigs_smallest(A, B, nRequest, opts);
        lam = real(diag(D));
        [lam, idx] = sort(lam, 'ascend');
        V = V(:, idx);

        if isLBO
            % LBO is real symmetric with a simple spectrum on the canonical
            % cortex; standard B-normalization (scale each column so phi'Bphi=1)
            % suffices and keeps the modes real.
            V   = real(V);
            nrm = sqrt(real(diag(V' * (B * V))));
            nrm(nrm < eps) = eps;
            V = V * spdiags(1 ./ nrm, 0, numel(nrm), numel(nrm));
            Vk   = V(:, 1:K);
            lamk = lam(1:K);
        else
            % Connection Laplacian (complex Hermitian) and Dirac both have
            % degenerate multiplets, so eigs returns a B-non-orthonormal /
            % rank-deficient spanning set. Rank-revealing B-orthonormalization +
            % Rayleigh-Ritz recovers K genuine B-orthonormal eigenpairs. The
            % helper is complex-safe (conjugate transpose throughout), so the
            % complex connection modes are preserved.
            [Vk, lamk] = local_ritz_basis(A, B, V, K);
        end

        Phi{hh}            = Vk;
        Lambda{hh}         = lamk;
        GlobalVertices{hh} = gv;
        if isFace, GlobalFaces{hh} = Op.GlobalFaces{hh}; end
    end

    % --- assemble EigenMat ---
    EigenMat                = db_template('eigenmat');
    EigenMat.OperatorFile   = file_short(OperatorFile);
    EigenMat.Variant        = Variant;
    EigenMat.ParentSurface  = SurfaceFile;
    EigenMat.Phi            = Phi;             % 1x2 cell of eigenvector matrices
    EigenMat.Lambda         = Lambda;          % 1x2 cell of eigenvalue vectors [K x 1]
    EigenMat.K              = K;
    EigenMat.GlobalVertices = GlobalVertices;  % 1x2 cell of global vertex indices
    EigenMat.GlobalFaces    = GlobalFaces;     % 1x2 cell of global face indices (face-domain variants)
    EigenMat.Provenance     = prov;

    % --- save / register in DB ---
    if ~NoSave
        bst_progress('text', 'Saving eigenmodes to database...');
        Comment = sprintf('%s eigenmodes (K=%d)', Variant, K);
        db_add_eigen(iSubject, SurfaceFile, EigenMat, Comment);
    end
end

% ----------------------------------------------------------------------------
function OperatorFile = local_find_operator(sSubject, iSurface, Variant)
% Find the first Operator child of the parent surface whose variant matches.
% The DB child entry's .Variant field is not reliably populated (db_add_operator
% only sets FileName/Comment), so fall back to loading the operator file and
% reading its stored Variant field.
    OperatorFile = '';
    if isempty(sSubject) || ~isfield(sSubject,'Surface') || isempty(sSubject.Surface) ...
       || iSurface > numel(sSubject.Surface) ...
       || ~isfield(sSubject.Surface(iSurface),'Operator') ...
       || isempty(sSubject.Surface(iSurface).Operator)
        return;
    end
    ops = sSubject.Surface(iSurface).Operator;
    for k = 1:numel(ops)
        % Prefer the entry's Variant if present and non-empty.
        if isfield(ops(k),'Variant') && ~isempty(ops(k).Variant) ...
           && strcmpi(ops(k).Variant, Variant)
            OperatorFile = ops(k).FileName;
            return;
        end
    end
    % Fallback: load each operator file and read its Variant.
    for k = 1:numel(ops)
        try
            f = file_fullpath(ops(k).FileName);
            S = load(f, 'Variant');
            if isfield(S,'Variant') && strcmpi(S.Variant, Variant)
                OperatorFile = ops(k).FileName;
                return;
            end
        catch
            % skip unreadable/missing entries
        end
    end
end

% ----------------------------------------------------------------------------
function EigenMat = local_find_eigen(sSubject, iSurface, Variant, K, Tau, isDirac)
% Scan the parent surface's Eigen child nodes for a cached eigenbasis of the
% requested Variant carrying at least K modes (and, for Dirac, matching Tau).
% Returns the loaded EigenMat (with its stored OperatorFile reference) or [] if
% none qualifies. The node may hold more than K modes; the caller truncates.
    EigenMat = [];
    if isempty(sSubject) || ~isfield(sSubject,'Surface') || isempty(sSubject.Surface) ...
       || iSurface > numel(sSubject.Surface) ...
       || ~isfield(sSubject.Surface(iSurface),'Eigen') ...
       || isempty(sSubject.Surface(iSurface).Eigen)
        return;
    end
    nodes = sSubject.Surface(iSurface).Eigen;
    for k = 1:numel(nodes)
        % Skip nodes whose registered Variant is set and does not match.
        if isfield(nodes(k),'Variant') && ~isempty(nodes(k).Variant) ...
           && ~strcmpi(nodes(k).Variant, Variant)
            continue;
        end
        try
            E = load(file_fullpath(nodes(k).FileName));
        catch
            continue;   % unreadable / missing entry
        end
        if ~isfield(E,'Variant') || ~strcmpi(E.Variant, Variant);  continue; end
        if ~isfield(E,'K') || isempty(E.K) || E.K < K;             continue; end
        if isDirac && abs(local_eigen_tau(E, NaN) - Tau) > 1e-12;  continue; end
        % The basis must actually carry usable eigen data.
        if ~isfield(E,'Phi') || isempty(E.Phi) || ~isfield(E,'Lambda') || isempty(E.Lambda)
            continue;
        end
        EigenMat = E;
        return;
    end
end

% ----------------------------------------------------------------------------
function E = local_truncate_eigen(E, K)
% Truncate a (possibly larger) cached eigenbasis to its first K modes per
% hemisphere. The basis is nested (ascending eigenvalues, B-orthonormal), so the
% leading K columns are exactly the K-mode basis. No-op when E.K == K.
    if ~isfield(E,'K') || isempty(E.K) || E.K <= K
        return;
    end
    for hh = 1:numel(E.Phi)
        if ~isempty(E.Phi{hh}) && size(E.Phi{hh},2) > K
            E.Phi{hh} = E.Phi{hh}(:, 1:K);
        end
        if isfield(E,'Lambda') && hh <= numel(E.Lambda) ...
           && ~isempty(E.Lambda{hh}) && numel(E.Lambda{hh}) > K
            E.Lambda{hh} = E.Lambda{hh}(1:K);
        end
    end
    E.K = K;
end

% ----------------------------------------------------------------------------
function tau = local_eigen_tau(EigenMat, default)
% Read Tau from an eigen node's Provenance, falling back to default if absent.
    tau = default;
    if isfield(EigenMat,'Provenance') && isstruct(EigenMat.Provenance) ...
       && isfield(EigenMat.Provenance,'Tau') && ~isempty(EigenMat.Provenance.Tau)
        tau = EigenMat.Provenance.Tau;
    end
end

% ----------------------------------------------------------------------------
function [Phi, mu] = local_ritz_basis(L, B, U, K)
% Rank-revealing B-orthonormalization of the eigs spanning set U, then a
% Rayleigh-Ritz diagonalization of L on that span, returning K genuine
% B-orthonormal eigenpairs (ascending), matched (Phi, mu). No randomness.
%
% Why: eigs(L,B,'smallestabs') returns vectors that are not B-orthonormal across
% degenerate eigenvalues (the Dirac 4-fold quaternionic multiplets and the
% connection Laplacian's parallel-transport degeneracies), so the raw Gram U'BU
% is non-identity / rank-deficient. We keep only the well-conditioned directions,
% build a B-orthonormal basis W of span(U), then diagonalize the small Lr = W'LW
% so every output column is a true eigenvector with its matched Ritz value.
% Complex-safe (conjugate transpose throughout): preserves the complex connection
% modes (a real-valued Rayleigh-Ritz routine would collapse them).
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
        error('tess_eigen:rankDeficient', ...
            'Only %d independent modes recovered (< K=%d); increase the eigs over-fetch.', ...
            size(Phi, 2), K);
    end
    Phi = Phi(:, 1:K);
    mu  = mu(1:K);
    negThresh = -1e-6 * max(abs(mu));
    if any(mu < negThresh)
        warning('tess_eigen:nonPSD', ...
            ['Eigenproblem returned %d eigenvalue(s) below %.2e (most negative %.3e); ' ...
             'clamping to 0. The nxr mass/operator may be non-PSD.'], ...
            sum(mu < negThresh), negThresh, min(mu));
    end
    mu(mu < 0) = 0;
end
