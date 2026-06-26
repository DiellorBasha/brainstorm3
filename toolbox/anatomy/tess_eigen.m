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
%     nxr_compute('create'); tess_eigen never builds a mesh itself.
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
%                        whose stored nModes is at least the requested nModes, and (for
%                        Dirac) whose Tau matches, is loaded and truncated to nModes
%                        instead of re-solving. Set true to force a fresh eigensolve.
%     'Interactive'    : true/false (default false) — GUI only. When true and an
%                        exact-spec match already exists, prompt Overwrite/Cancel:
%                        Cancel reuses the existing node; Overwrite deletes it and
%                        recomputes (no duplicate). Programmatic callers leave it false.
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
    Interactive    = false;  % GUI: prompt Overwrite/Cancel when an exact-spec match exists
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case {'k','nmodes'},   K              = varargin{i+1};
            case 'tau',            Tau            = varargin{i+1};
            case 'nosave',         NoSave         = logical(varargin{i+1});
            case 'forcerecompute', ForceRecompute = logical(varargin{i+1});
            case 'interactive',    Interactive    = logical(varargin{i+1});
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
            existEntry = sCached.Surface(iSurfCached).Eigen(iEigCached);
            if Interactive
                % Prompt Overwrite / Cancel. Cancel reuses the existing node; Overwrite
                % deletes it and falls through to recompute (no duplicate accumulates).
                msg = sprintf(['A %s eigenbasis already exists for this surface ' ...
                    '(%d modes%s).\n\nOverwrite it (delete and recompute)?\n' ...
                    '[No keeps and reuses the existing one.]'], ...
                    Variant, existEntry.nModes, local_tau_str(existEntry.Tau));
                if ~java_dialog('confirm', msg, 'Compute eigenmodes')
                    EigenMat = local_truncate_eigen(in_bst_eigen(existEntry.FileName), K);
                    return;
                end
                db_delete_surface_node(existEntry.FileName, 1);   % overwrite: drop the old node
            else
                EigenMat = local_truncate_eigen(in_bst_eigen(existEntry.FileName), K);
                return;
            end
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
        % Create it (tess_operators -> nxr_compute('create')).
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
    local_registry_consistency(Op, isFace, isDirac);

    % --- nxr version for provenance ---
    nxrVer = '';
    try, nxrVer = nxr_compute('version'); catch, end  %#ok<CTCH>

    prov = struct('Backend','nxr', 'NxrVersion',nxrVer, 'Variant',Variant, 'nModes',K);
    if isDirac
        prov.Tau = Tau;
    end
    % Orthonormalization / solve strategy is per-operator (see the local_solve_* helpers).
    switch Variant
        case 'Laplace-Beltrami',     prov.Ortho = 'B-orthonormal (real)';
        case 'Connection Laplacian', prov.Ortho = 'M-orthonormal (smallest-positive)';
        case {'Dirac','Dirac-Face'}, prov.Ortho = 'Rayleigh-Ritz';
        case 'Hodge-Face',           prov.Ortho = 'Hodge lift (W_F-orthonormal)';
        otherwise,                   prov.Ortho = '';
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

        hemiName = 'left'; if hh == 2, hemiName = 'right'; end
        bst_progress('text', sprintf('Eigensolve: %s hemisphere (%s, K=%d)...', hemiName, Variant, K));

        % Per-operator eigensolve. Each operator is a different mathematical object and
        % needs its own tuned solver -- a single shared path produces wrong/spurious modes
        % for at least one of them (see the local_solve_* helpers):
        %   Laplace-Beltrami     real symmetric PSD      -> scalar fields
        %   Connection Laplacian complex Hermitian, ~PSD -> 2D tangent-vector fields
        %   Dirac / Dirac-Face   quaternionic, 4-fold    -> 3D embedded vector fields
        %   Hodge-Face           scalar lapFace + lift   -> 3D face vector fields
        switch Variant
            case 'Laplace-Beltrami'
                [Vk, lamk] = local_solve_lbo(A, B, K);
            case 'Connection Laplacian'
                [Vk, lamk] = local_solve_connection(A, B, K);
            case {'Dirac', 'Dirac-Face'}
                if isFace, nDom = numel(Op.GlobalFaces{hh}); else, nDom = numel(gv); end
                [Vk, lamk] = local_solve_dirac(A, B, K, nDom, tags{hh});
            case 'Hodge-Face'
                [Vk, lamk] = local_hodge_face_modes(A, Op.FaceAux{hh}, K);
            otherwise
                error('tess_eigen:badVariant', 'No eigensolve defined for variant ''%s''.', Variant);
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
function [Phi, lam] = local_solve_lbo(A, B, K)
% Laplace-Beltrami eigensolve: A is real symmetric and PSD with a simple spectrum on the
% canonical cortex. Take the K smallest eigenpairs (smoothest scalar modes, for scalar
% cortical data analysis), keep them REAL, and B-normalize each column so phi'*B*phi = 1.
    n = size(A, 1);
    if K > n - 2
        error('tess_eigen:tooManyModes', 'K=%d exceeds nV-2=%d.', K, n - 2);
    end
    nReq = min(K + 8, n - 2);
    [V, D] = bst_eigs_smallest(A, B, nReq, struct('tol', 1e-6, 'maxit', 1000, 'disp', 0));
    lam = real(diag(D));
    [lam, idx] = sort(lam, 'ascend');
    V = real(V(:, idx));
    nrm = sqrt(real(diag(V' * (B * V))));
    nrm(nrm < eps) = eps;
    V = V * spdiags(1 ./ nrm, 0, numel(nrm), numel(nrm));
    Phi = V(:, 1:K);
    lam = lam(1:K);
end

% ----------------------------------------------------------------------------
function [Phi, lam] = local_solve_connection(A, B, K)
% Connection-Laplacian eigensolve: A is COMPLEX HERMITIAN and only ~PSD -- the discrete
% connection Laplacian carries a few small spurious NEGATIVE eigenvalues (a connection
% Laplacian on a closed curved surface has no zero/negative mode, so anything <= 0 is a
% discretization artifact). The genuine modes are the smooth low-frequency POSITIVE ones
% (the 2D tangent-vector / n-RoSy basis). So: over-fetch the smallest-magnitude eigenpairs,
% DROP the non-positive spurious modes, keep the smallest K positives, M-orthonormalize.
% No Rayleigh-Ritz: the connection spectrum is generically simple, and re-diagonalizing on
% an over-fetched span re-introduces the spurious negatives. Complex modes preserved.
    n = size(A, 1);
    if K > n - 2
        error('tess_eigen:tooManyModes', 'K=%d exceeds nV-2=%d.', K, n - 2);
    end
    % Over-fetch so that, after dropping the few spurious non-positive modes, at least K
    % genuine positive eigenpairs remain.
    nReq = min(K + max(16, ceil(0.5 * K)), n - 2);
    [V, D] = bst_eigs_smallest(A, B, nReq, struct('tol', 1e-8, 'maxit', 2000, 'disp', 0));
    lam = real(diag(D));
    keep = (lam > 0);                       % drop the spurious non-positive artifacts
    V = V(:, keep);  lam = lam(keep);
    [lam, idx] = sort(lam, 'ascend');
    V = V(:, idx);
    if numel(lam) < K
        error('tess_eigen:connTooFewPositive', ...
            ['Only %d positive connection eigenpairs recovered (< K=%d); increase the ' ...
             'over-fetch in local_solve_connection.'], numel(lam), K);
    end
    V   = V(:, 1:K);
    lam = lam(1:K);
    % M-orthonormal magnitude (eigs returns B-normalized; guard it, complex-safe).
    nrm = sqrt(real(sum(conj(V) .* (B * V), 1)));
    nrm(nrm < eps) = eps;
    Phi = V ./ nrm;
    lam = lam(:);
end

% ----------------------------------------------------------------------------
function [Phi, lam] = local_solve_dirac(A, B, K, nDom, tag)
% Dirac / Dirac-Face eigensolve: A is quaternionic [4*nDom] with 4-fold quaternionic
% multiplets, so eigs returns a B-non-orthonormal / rank-deficient spanning set. Over-fetch
% generously, then Rayleigh-Ritz (rank-revealing B-orthonormalization + re-diagonalization)
% recovers K genuine B-orthonormal eigenpairs (for 3D embedded vector-field analysis).
% Complex-safe throughout.
    if size(A, 1) ~= 4 * nDom
        error('tess_eigen:diracSizeMismatch', ...
            'Dirac operator on hemisphere %s is %dx%d, expected 4*nDom=%d.', ...
            tag, size(A, 1), size(A, 2), 4 * nDom);
    end
    if K > 4 * nDom - 2
        error('tess_eigen:tooManyModes', ...
            'K=%d exceeds 4*nV-2=%d on hemisphere %s.', K, 4 * nDom - 2, tag);
    end
    % ~30% over-fetch keeps rank >= K despite the 4-fold multiplets.
    nReq = min(K + max(8, ceil(0.3 * K)), 4 * nDom - 2);
    [V, ~] = bst_eigs_smallest(A, B, nReq, struct('tol', 1e-6, 'maxit', 1000, 'disp', 0));
    [Phi, lam] = local_ritz_basis(A, B, V, K);
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
function s = local_tau_str(tau)
% Format a Tau value for a dialog message: ', tau=0.5' or '' when empty (non-Dirac).
    if isempty(tau); s = ''; else; s = sprintf(', tau=%.3g', tau); end
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

% ----------------------------------------------------------------------------
function [V, D] = bst_eigs_smallest(A, B, k, opts)
% The k smallest generalized eigenpairs of a (near-)singular symmetric/Hermitian
% pencil (A, B) with B SPD, avoiding the sigma=0 shift-invert of the singular A that
% triggers MATLAB's RCOND "matrix close to singular" warning. Returns V, D ascending
% and real-valued (the pencil is symmetric/Hermitian), matching
% eigs(A,B,k,'smallestabs'); V may be complex for a Hermitian A (e.g. the connection
% Laplacian). Forces the symmetric/Hermitian Lanczos path (real spectrum, faster, clean
% degenerate multiplets) and selects the smallest modes via a small negative sigma shift
% so (A - sigma*B) = A + |sigma|*B is SPD/well-conditioned. The shared eigensolve backend
% of every local_solve_* helper above; verified end-to-end by test_bst_eigs_smallest.
    if nargin < 4 || isempty(opts)
        opts = struct('tol', 1e-6, 'maxit', 1000, 'disp', 0);
    end
    % Symmetrize / Hermitian-ize: no-op to ~1e-16 for already-symmetric operators,
    % but makes issymmetric()/ishermitian() true so eigs uses Lanczos, not Arnoldi.
    A = (A + A') / 2;
    B = (B + B') / 2;
    % Factorization-free spectrum-scale estimate. Only an order of magnitude is needed
    % (to place the shift), so use a coarse tolerance, not the caller's tight opts.
    % 'largestabs' factorizes only the SPD, well-conditioned mass B (never the singular
    % A), so it cannot warn. Its own try/catch makes a non-converged estimate fall
    % through to the legacy path instead of hard-erroring before the shift fallbacks.
    estOpts = struct('tol', 1e-2, 'maxit', 300, 'disp', 0);
    try
        lmax = abs(eigs(A, B, 1, 'largestabs', estOpts));
    catch
        lmax = 0;
    end
    if ~isfinite(lmax) || (lmax <= 0)
        [V, D] = eigs(A, B, k, 'smallestabs', opts);     % degenerate estimate: legacy path
    else
        % Small negative shift below the (PSD) spectrum bottom: (A - sigma*B) =
        % A + |sigma|*B is SPD/well-conditioned, yet 'nearest sigma' still returns the
        % k smallest modes (including the near-zero kernel). The -1e-7 -> -1e-4
        % escalation lifts the factorization further if the first shift lands too close
        % to the kernel for eigs to converge; 'smallestabs' is the last-resort path.
        try
            [V, D] = eigs(A, B, k, -1e-7 * lmax, opts);
        catch
            try
                [V, D] = eigs(A, B, k, -1e-4 * lmax, opts);
            catch
                [V, D] = eigs(A, B, k, 'smallestabs', opts);
            end
        end
    end
    % Enforce the documented contract: real eigenvalues (symmetric/Hermitian pencil) in
    % ascending order, matching eigs(...,'smallestabs'). Eigenvectors are NOT realified
    % -- they are genuinely complex for a Hermitian A (e.g. the connection Laplacian).
    d = real(diag(D));
    [d, idx] = sort(d, 'ascend');
    V = V(:, idx);
    D = diag(d);
end

% ----------------------------------------------------------------------------
function local_registry_consistency(Op, isFace, isDirac)
% Guard: if the operator node carries nxr v0.2.0 registry metadata, verify its
% domain / field-type agrees with the Variant-derived flags. Warn (never error)
% on drift so a future nxr id repurpose is caught loudly without breaking the
% solve. No-op when Registry is absent (old nodes / pre-registry binary).
    if ~isfield(Op,'Registry') || isempty(Op.Registry) || ~isfield(Op.Registry,'Primary') ...
            || isempty(Op.Registry.Primary)
        return;
    end
    P = Op.Registry.Primary;
    regFace  = isfield(P,'domain')     && strcmpi(P.domain, 'face');
    regQuat  = isfield(P,'field_type') && strcmpi(P.field_type, 'quaternion');
    if (regFace ~= logical(isFace)) || (regQuat ~= logical(isDirac))
        id_str  = '?'; if isfield(P,'id'),         id_str  = P.id;         end
        dom_str = '?'; if isfield(P,'domain'),     dom_str = P.domain;     end
        ft_str  = '?'; if isfield(P,'field_type'), ft_str  = P.field_type; end
        warning('tess_eigen:registryMismatch', ...
            ['Operator registry (%s: domain=%s field_type=%s) disagrees with the ' ...
             'Variant-derived flags (isFace=%d isDirac=%d). nxr ids may have drifted.'], ...
            id_str, dom_str, ft_str, isFace, isDirac);
    end
end
