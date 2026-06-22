function EigenMat = tess_eigen(SurfaceFile, OperatorName, varargin)
% TESS_EIGEN: Assemble a per-hemisphere eigenbasis as an eigen_ DB node.
%
% USAGE:  EigenMat = tess_eigen(SurfaceFile, OperatorName)
%         EigenMat = tess_eigen(SurfaceFile, OperatorName, ...
%                               'nModes',400, 'Tau',0.5, 'NoSave',false, 'ForceRecompute',false)
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
%     'nModes'         : number of modes to keep per hemisphere (default 400).
%                        ('K' is accepted as a legacy alias.)
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
            case {'k','nmodes'},   K              = varargin{i+1};
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
        case {'hodge-face','hodgeface'}
            Variant = 'Hodge-Face';
        otherwise
            error('tess_eigen:badVariant', ...
                ['Unknown operator ''%s''. Valid options: ' ...
                 '''Laplace-Beltrami'', ''Connection Laplacian'', ''Dirac'', ''Dirac-Face'', ''Hodge-Face''.'], OperatorName);
    end
    isDiracFace = strcmpi(Variant, 'Dirac-Face');
    isHodgeFace = strcmpi(Variant, 'Hodge-Face');           % scalar lapFace eigensolve + Hodge vector lift
    isFace  = isDiracFace || isHodgeFace;                   % face-domain (modes on faces)
    isDirac = strcmpi(Variant, 'Dirac') || isDiracFace;     % quaternion Dirac-type: over-fetch + Rayleigh-Ritz + Tau
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
    %     whose stored nModes >= the requested K, and (for Dirac) whose Tau matches. The
    %     eigenbasis is nested (ascending eigenvalues, B-orthonormal), so a node with
    %     more modes is truncated to its first K columns. ForceRecompute skips this.
    %     bst_get resolves the match from the cache (no file load); we then load the file.
    if ~ForceRecompute
        [sCached, ~, iSurfCached, iEigCached] = bst_get('EigenFileForSurface', SurfaceFile, Variant, K, Tau);
        if ~isempty(iEigCached)
            cached   = in_bst_eigen(sCached.Surface(iSurfCached).Eigen(iEigCached).FileName);
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

    % --- find-or-create the operator node of the requested Variant (Tau-aware) ---
    %     bst_get matches Variant and, for the Dirac-type variants, Tau -- so a Dirac
    %     operator computed at a different Tau is not wrongly reused.
    OperatorFile = local_operator_file(SurfaceFile, Variant, Tau);
    if isempty(OperatorFile)
        % Create it (tess_operators -> nxr_safe_create -> nxr_compute('create')).
        bst_progress('text', sprintf('Computing %s operator...', Variant));
        tess_operators(SurfaceFile, OperatorName, 'Tau', Tau);
        % Re-resolve from the refreshed cache (db_add_operator updated GlobalData).
        OperatorFile = local_operator_file(SurfaceFile, Variant, Tau);
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
    Op = in_bst_operator(OperatorFile);
    if ~isfield(Op,'Operator') || ~isfield(Op,'Mass') || ~isfield(Op,'GlobalVertices')
        error('tess_eigen:badOperatorFile', ...
            'Operator file is missing Operator/Mass/GlobalVertices: %s', OperatorFile);
    end

    % --- nxr version for provenance ---
    nxrVer = '';
    try, nxrVer = nxr_compute('version'); catch, end  %#ok<CTCH>

    prov = struct('Backend','nxr', 'NxrVersion',nxrVer, 'Variant',Variant, 'nModes',K);
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

        if isHodgeFace
            % Scalar face Laplacian eigensolve (A=lapFace [F x F]) + Hodge vector lift
            % {grad psi_k, n x grad psi_k}/sqrt(lambda_k), W_F-orthonormalized, embedded as
            % pure-imaginary quaternions [4F x 2K]. (Not a quaternion Dirac eigensolve.)
            [Vk, lamk] = local_hodge_face_modes(A, Op.FaceAux{hh}, K);
            Phi{hh} = Vk;  Lambda{hh} = lamk;
            GlobalVertices{hh} = gv;  GlobalFaces{hh} = Op.GlobalFaces{hh};
            continue;
        end

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
    EigenMat.nModes         = K;
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
function [Phi, lam] = local_hodge_face_modes(lapFace, aux, Kvec)
% Face Hodge vector eigenbasis: the scalar face Laplacian eigenmodes psi_k (smooth,
% Weyl spectrum) lifted to vector fields {grad psi_k, n x grad psi_k}/sqrt(lambda_k),
% Gram-Cholesky W_F-orthonormalized, embedded as pure-imaginary quaternions [4F x Kvec]
% (w=0) so the existing bst_dirac face Transform/Reconstruct consumes them unchanged.
    G  = aux.GradFace;            % [3F x F] gradFace
    Nf = aux.FaceNormal;          % [F x 3] outward face normals
    Mf = aux.ScalarMass;          % [F x F] diag(faceArea)
    nF = size(Mf,1);
    Ks = min(ceil(Kvec/2), nF-2);
    % smallest Ks+1 scalar lapFace modes (drop the near-zero constant), M_f-orthonormal
    [Psi, D] = bst_eigs_smallest(lapFace, Mf, Ks+1, struct('tol',1e-8,'maxit',1000,'disp',0));
    lam_s = real(diag(D));  [lam_s, idx] = sort(lam_s,'ascend');  Psi = real(Psi(:,idx));
    Psi = Psi(:, 2:Ks+1);  lam_s = lam_s(2:Ks+1);                 % [F x Ks]
    nrm = sqrt(max(real(diag(Psi'*(Mf*Psi))), eps));  Psi = Psi * spdiags(1./nrm,0,Ks,Ks);
    % Hodge lift: grad and skew-grad of each scalar mode, scaled by 1/sqrt(lambda)
    fArea = full(diag(Mf));  WF3 = spdiags(repelem(fArea,3),0,3*nF,3*nF);
    sc = 1 ./ sqrt(max(lam_s, eps));
    U = zeros(3*nF, 2*Ks);
    for k = 1:Ks
        gk = G*Psi(:,k);                              % [3F] grad psi_k (x,y,z per face)
        gm = reshape(gk, 3, [])';                     % [F x 3]
        sk = reshape(cross(Nf, gm, 2)', [], 1);       % [3F] n x grad
        U(:, 2*k-1) = gk * sc(k);
        U(:, 2*k)   = sk * sc(k);
    end
    % W_F-orthonormalize the stacked [3F x 2Ks] set (rank-revealing)
    Gram = U'*WF3*U;  Gram = (Gram+Gram')/2;
    [Rc,p] = chol(Gram);
    if p == 0
        U = U / Rc;
    else
        [Vg,Dg] = eig(full(Gram));  dg = real(diag(Dg));  keep = dg > 1e-10*max(dg);
        U = U*Vg(:,keep)*spdiags(1./sqrt(dg(keep)),0,sum(keep),sum(keep));
    end
    nC = min(size(U,2), Kvec);  U = U(:,1:nC);
    % embed as pure-imaginary quaternions [4F x nC]: rows 2,3,4 of each face block = (x,y,z)
    Phi = zeros(4*nF, nC);
    Phi(2:4:end,:) = U(1:3:end,:);
    Phi(3:4:end,:) = U(2:3:end,:);
    Phi(4:4:end,:) = U(3:3:end,:);
    lam = repelem(lam_s(:),2);  lam = lam(1:nC);       % paired scalar eigenvalue per vector mode
end

% ----------------------------------------------------------------------------
function OperatorFile = local_operator_file(SurfaceFile, Variant, Tau)
% Resolve the operator node of the requested Variant (and, for the Dirac-type variants,
% Tau) from the DB cache via bst_get. Returns the node's relative FileName, or '' if none.
    OperatorFile = '';
    [sOp, ~, iSurf, iOp] = bst_get('OperatorFileForSurface', SurfaceFile, Variant, Tau);
    if ~isempty(iOp)
        OperatorFile = sOp.Surface(iSurf).Operator(iOp).FileName;
    end
end

% ----------------------------------------------------------------------------
function E = local_truncate_eigen(E, K)
% Truncate a (possibly larger) cached eigenbasis to its first K modes per
% hemisphere. The basis is nested (ascending eigenvalues, B-orthonormal), so the
% leading K columns are exactly the K-mode basis. No-op when E.nModes == K.
    if ~isfield(E,'nModes') || isempty(E.nModes) || E.nModes <= K
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
    E.nModes = K;
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
