function varargout = bst_dirac(HeadModel, varargin)
% BST_DIRAC: Spectral transform of the UNCONSTRAINED leadfield in the Dirac eigenbasis.
%
% USAGE:
%   % TRANSFORM (default): base unconstrained head model -> composed Dirac eigenmode head model
%     CompHM = bst_dirac(HeadModel)
%     CompHM = bst_dirac(HeadModel, 'nModes', K, 'Tau', tau)
%
%   % RECONSTRUCT: composed Dirac head model + coefficient rows -> per-vertex 3-vector field
%     J = bst_dirac(CompHM, 'Reconstruct', GainRows)
%
% DESCRIPTION:
%     Forward and inverse spectral change-of-basis for Dirac eigenmode source
%     mapping. (Spectral "transform" / "reconstruct", NOT the MEG/EEG forward and
%     inverse source-mapping solutions.)
%
%     TRANSFORM keeps the full unconstrained 3-vector gain and expands it in the
%     curvature-aware Dirac eigenbasis. Each per-vertex gain 3-vector (ambient
%     coords) is embedded as a pure-imaginary quaternion psi = [0, gx, gy, gz] and
%     projected onto the per-hemisphere B-orthonormal Dirac eigenvectors:
%         L~_h = Psi_h' * B_h * Phi_h      [nCh x K],   B_h = kron(Mass_h, I4)
%     The hemispheres are stacked into a composed head model
%         Gain = [L~_L, L~_R]              [nCh x 2K]
%     consumed unchanged by bst_inverse_eigenmodes.
%
%     RECONSTRUCT is the inverse spectral transform: since Phi is B-orthonormal,
%     a coefficient row c [1 x 2K] maps back to a cortical field by f_h = Phi_h * c_h,
%     keeping the quaternion vector part (drop the w slot) as the per-vertex 3-vector.
%     Used by view_leadfield_vectors to draw each channel's band-limited field.
%
%     The Dirac eigenbasis is fetched from the surface's eigen_ DB node (not passed
%     in). bst_dirac resolves HeadModel.SurfaceFile, finds a 'Dirac' eigen node
%     matching (nModes, Tau), and creates one via tess_eigen if absent. The mass B_h
%     is read from the operator_ node referenced by the eigen node (Op.Mass{h} already
%     stores the full kron(Mass_h, I4)); reconstruction needs only Phi (no mass load).
%
% INPUTS / OUTPUTS:
%   TRANSFORM:
%     HeadModel : base UNCONSTRAINED surface head model (in_bst_headmodel, ApplyOrient=0):
%                 .Gain [nCh x 3*nVert], .HeadModelType='surface', .SurfaceFile, .Comment
%     CompHM    : composed head model (Gain [nCh x 2K], Eigenvalues [2K x 1], nModes=2K,
%                 ModeHemisphere [2K x 1], HemiGlobalVertices {L,R}, DiracEigenFile,
%                 DiracTau, isEigenmode/isDiracEigenmode flags).
%   RECONSTRUCT:
%     HeadModel : composed Dirac head model (CompHM, as produced above)
%     GainRows  : [m x 2K] eigenmode-leadfield rows (columns stacked L then R)
%     J         : [m x 3*nVert] reconstructed cortical leadfield (x,y,z per vertex)
%
% OPTIONS (TRANSFORM):
%     'nModes' : modes per hemisphere (default 400). The Dirac eigen node is
%                found-or-created with this many modes; an existing node with at
%                least nModes (and matching Tau) is reused and truncated.
%     'Tau'    : relative-Dirac curvature weighting in [0,1] (default 0.5),
%                forwarded to tess_eigen when a node must be created.
%
% SEE ALSO: tess_eigen, tess_operators, bst_inverse_eigenmodes, view_leadfield_vectors
%
% Authors: Diellor Basha, 2026

    % --- options / mode selection ---
    nModes        = 400;
    Tau           = 0.5;
    doReconstruct = false;
    GainRows      = [];
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'nmodes',      nModes = varargin{i+1};
            case 'tau',         Tau    = varargin{i+1};
            case 'reconstruct', doReconstruct = true; GainRows = varargin{i+1};
            otherwise
                error('bst_dirac:badOption', 'Unknown option: %s', varargin{i});
        end
    end

    if doReconstruct
        varargout{1} = local_reconstruct(HeadModel, GainRows);
        return;
    end

    % ===== TRANSFORM =====
    if isempty(nModes) || nModes <= 0
        nModes = 400;
    end
    if ~(isscalar(Tau) && Tau >= 0 && Tau <= 1)
        error('bst_dirac:badTau', 'Tau must be a scalar in [0,1].');
    end

    % surface head model only
    if isfield(HeadModel,'HeadModelType') && ~isempty(HeadModel.HeadModelType) ...
            && ~strcmpi(HeadModel.HeadModelType,'surface')
        error('bst_dirac:NotSurface', ...
            'Dirac eigenmode leadfield requires a surface head model (got ''%s'').', HeadModel.HeadModelType);
    end

    G = double(HeadModel.Gain);          % [nCh x 3*nVert] unconstrained
    nCh = size(G,1);
    if mod(size(G,2), 3) ~= 0
        error('bst_dirac:NotUnconstrained', ...
            'Dirac leadfield requires an unconstrained head model: Gain must be [nCh x 3*nVert].');
    end

    if ~isfield(HeadModel,'SurfaceFile') || isempty(HeadModel.SurfaceFile)
        error('bst_dirac:noSurfaceFile', ...
            'HeadModel has no SurfaceFile; cannot resolve the Dirac eigen node.');
    end

    % FACE-domain vs vertex-domain: a face leadfield (isFaceBased=1) is expressed in
    % the face Dirac eigenbasis ('Dirac-Face', modes on faces); everything else (the
    % quaternion embed, Psi'*B*Phi projection, reconstruct) is domain-agnostic.
    isFace  = isfield(HeadModel,'isFaceBased') && isequal(HeadModel.isFaceBased,1);
    Variant = 'Dirac';
    if isFace
        Variant = 'Hodge-Face';                         % default face basis (smooth, full-rank, localizes)
        if isfield(HeadModel,'FaceBasis') && ~isempty(HeadModel.FaceBasis)
            switch lower(HeadModel.FaceBasis)
                case {'dirac','dirac-face'}, Variant = 'Dirac-Face';   % wide-root, non-localizing (kept for study)
                case {'hodge','hodge-face'}, Variant = 'Hodge-Face';
            end
        end
    end

    % --- find-or-create the Dirac eigenbasis as an eigen_ DB node ---
    [EigenMat, OperatorMat, EigenFile] = local_get_dirac_eigen(HeadModel.SurfaceFile, nModes, Tau, Variant);

    % Truncate to nModes (a reused node may carry more modes than requested).
    K = min(nModes, EigenMat.nModes);
    Tau = local_eigen_tau(EigenMat, Tau);

    Lblk = cell(1,2); vblk = cell(1,2); hblk = cell(1,2); gblk = cell(1,2);
    for hh = 1:2
        if isFace, vH = EigenMat.GlobalFaces{hh}(:); else, vH = EigenMat.GlobalVertices{hh}(:); end
        nVh = numel(vH);                          % domain elements (faces or vertices) this hemi
        Phi = double(EigenMat.Phi{hh});
        if size(Phi,1) ~= 4*nVh
            error('bst_dirac:shapeMismatch', ...
                'Hemisphere %d: Phi has %d rows, expected 4*nV=%d.', hh, size(Phi,1), 4*nVh);
        end
        if size(Phi,2) < K
            error('bst_dirac:shapeMismatch', ...
                'Hemisphere %d: Phi has %d columns, fewer than K=%d.', hh, size(Phi,2), K);
        end
        Phi  = Phi(:, 1:K);
        Vals = double(EigenMat.Lambda{hh}(:)); Vals = Vals(1:K);

        B = OperatorMat.Mass{hh};                 % [4Vh x 4Vh] = kron(Mass_h, I4)
        if size(B,1) ~= 4*nVh || size(B,2) ~= 4*nVh
            error('bst_dirac:shapeMismatch', ...
                'Hemisphere %d: operator Mass is %dx%d, expected 4*nV=%d.', hh, size(B,1), size(B,2), 4*nVh);
        end

        % embed unconstrained gain as a pure-imaginary quaternion field (w=0)
        Psi = zeros(4*nVh, nCh);
        Psi(2:4:end, :) = G(:, (vH-1)*3 + 1).';   % i (x)
        Psi(3:4:end, :) = G(:, (vH-1)*3 + 2).';   % j (y)
        Psi(4:4:end, :) = G(:, (vH-1)*3 + 3).';   % k (z)   (w rows 1:4:end stay 0)

        Lblk{hh} = Psi' * (B * Phi);              % [nCh x K]
        vblk{hh} = Vals;
        hblk{hh} = hh * ones(K,1);
        gblk{hh} = vH;                            % domain scatter (faces or vertices)
    end

    CompHM = HeadModel;
    CompHM.Gain               = [Lblk{1}, Lblk{2}];     % [nCh x 2K]
    CompHM.GridLoc            = [];
    CompHM.GridOrient         = [];
    CompHM.GridAtlas          = [];
    CompHM.isEigenmode        = 1;
    CompHM.isDiracEigenmode   = 1;
    CompHM.isFaceBased        = isFace;
    CompHM.FaceBasisVariant   = Variant;                % which face basis was used (for Reconstruct fallback)
    CompHM.nModes             = 2*K;
    CompHM.Eigenvalues        = [vblk{1}; vblk{2}];     % [2K x 1]
    CompHM.ModeHemisphere     = [hblk{1}; hblk{2}];     % [2K x 1] hemisphere index per mode
    CompHM.HemiGlobalVertices = {EigenMat.GlobalVertices{1}(:), EigenMat.GlobalVertices{2}(:)};
    if isFace
        CompHM.HemiGlobalFaces = {gblk{1}, gblk{2}};    % face scatter for Reconstruct
    end
    CompHM.HeadModelType      = 'surface';
    % Stamp the eigen-node provenance so Reconstruct re-fetches the exact basis.
    CompHM.DiracEigenFile     = '';
    if ~isempty(EigenFile); CompHM.DiracEigenFile = file_short(EigenFile); end
    CompHM.DiracTau           = Tau;
    CompHM.Comment            = sprintf('Dirac eigenmode leadfield (%d modes, tau=%.3g) | %s', ...
        2*K, Tau, local_default(HeadModel,'Comment',''));

    varargout{1} = CompHM;
end

% ----------------------------------------------------------------------------
function J = local_reconstruct(CompHM, GainRows)
% Inverse spectral transform: map Dirac eigenmode coefficient rows back to a
% per-vertex 3-vector cortical field. f_h = Phi_h * c_h, keeping the quaternion
% vector part (rows 2:4 of each vertex block; the w slot 1 is dropped).
    if ~isfield(CompHM,'isDiracEigenmode') || isempty(CompHM.isDiracEigenmode) || ~CompHM.isDiracEigenmode
        error('bst_dirac:NotDiracEigenmode', ...
            'Reconstruct requires a composed Dirac eigenmode head model (isDiracEigenmode=1).');
    end
    if isempty(GainRows)
        error('bst_dirac:noGainRows', 'Reconstruct requires a [m x 2K] GainRows matrix.');
    end

    isFace   = isfield(CompHM,'isFaceBased') && isequal(CompHM.isFaceBased,1);
    EigenMat = local_eigen_for_reconstruct(CompHM);

    Kh    = CompHM.nModes / 2;
    if isFace, hemiIdx = CompHM.HemiGlobalFaces; else, hemiIdx = CompHM.HemiGlobalVertices; end
    nVert = sum(cellfun(@numel, hemiIdx));               % total domain elements (faces or vertices)
    m     = size(GainRows, 1);
    J     = zeros(m, 3*nVert);
    for hh = 1:2
        vH   = hemiIdx{hh}(:);
        Phi  = double(EigenMat.Phi{hh});
        if size(Phi,2) < Kh
            error('bst_dirac:shapeMismatch', ...
                'Hemisphere %d: eigen node has %d modes, fewer than the %d needed.', hh, size(Phi,2), Kh);
        end
        Phi  = Phi(:, 1:Kh);                                % [4Vh x Kh]
        cols = (CompHM.ModeHemisphere(:) == hh);            % [2K x 1] logical
        R    = Phi * GainRows(:, cols).';                   % [4Vh x m]
        J(:, (vH-1)*3 + 1) = R(2:4:end, :).';               % x  (drop w rows 1:4:end)
        J(:, (vH-1)*3 + 2) = R(3:4:end, :).';               % y
        J(:, (vH-1)*3 + 3) = R(4:4:end, :).';               % z
    end
end

% ----------------------------------------------------------------------------
function EigenMat = local_eigen_for_reconstruct(CompHM)
% Resolve the Dirac eigen node for reconstruction. Prefer the stamped
% DiracEigenFile (exact basis used by Transform); fall back to a find-by
% (Kh, DiracTau) on the surface. Never creates (the leadfield already exists).
    EigenMat = [];
    if isfield(CompHM,'DiracEigenFile') && ~isempty(CompHM.DiracEigenFile)
        try
            EigenMat = in_bst_eigen(CompHM.DiracEigenFile);
        catch
            EigenMat = [];
        end
    end
    if isempty(EigenMat)
        if ~isfield(CompHM,'SurfaceFile') || isempty(CompHM.SurfaceFile)
            error('bst_dirac:noSurfaceFile', ...
                'Composed head model has no DiracEigenFile and no SurfaceFile; cannot fetch the eigenbasis.');
        end
        Kh  = CompHM.nModes / 2;
        Tau = local_default(CompHM, 'DiracTau', 0.5);
        Variant = 'Dirac';
        if isfield(CompHM,'isFaceBased') && isequal(CompHM.isFaceBased,1)
            Variant = local_default(CompHM, 'FaceBasisVariant', 'Hodge-Face');
        end
        EigenMat = local_find_dirac_eigen(CompHM.SurfaceFile, Kh, Tau, Variant);
    end
    if isempty(EigenMat) || ~isfield(EigenMat,'Phi') || isempty(EigenMat.Phi)
        error('bst_dirac:badBasis', ...
            'Could not resolve a Dirac eigenbasis to reconstruct the field.');
    end
end

% ----------------------------------------------------------------------------
function [EigenMat, OperatorMat, EigenFile] = local_get_dirac_eigen(SurfaceFile, nModes, Tau, Variant)
% Find a Dirac eigen node under the surface matching (Variant, nModes, Tau), else
% create one via tess_eigen. Returns the loaded EigenMat, its referenced OperatorMat
% (per-hemisphere mass matrices), and the eigen node's file path. Variant is 'Dirac'
% (vertex, default) or 'Dirac-Face' (face-domain).
    if nargin < 4 || isempty(Variant), Variant = 'Dirac'; end
    [EigenMat, EigenFile] = local_find_dirac_eigen(SurfaceFile, nModes, Tau, Variant);
    if isempty(EigenMat)
        EigenMat = tess_eigen(SurfaceFile, Variant, 'nModes', nModes, 'Tau', Tau);
        % Re-find to recover the saved node's filename for provenance stamping.
        [found, EigenFile] = local_find_dirac_eigen(SurfaceFile, nModes, Tau, Variant);
        if ~isempty(found); EigenMat = found; end
    end
    if isempty(EigenMat) || ~isfield(EigenMat,'Phi') || isempty(EigenMat.Phi)
        error('bst_dirac:badBasis', ...
            'Could not obtain a Dirac eigenbasis for surface: %s', SurfaceFile);
    end
    if ~isfield(EigenMat,'OperatorFile') || isempty(EigenMat.OperatorFile)
        error('bst_dirac:badBasis', ...
            'Dirac eigen node has no OperatorFile reference; cannot load the mass matrix.');
    end
    OperatorMat = in_bst_operator(EigenMat.OperatorFile);
    if ~isfield(OperatorMat,'Mass') || isempty(OperatorMat.Mass) || numel(OperatorMat.Mass) ~= 2
        error('bst_dirac:badBasis', ...
            'Operator node is missing a 1x2 Mass; cannot project: %s', EigenMat.OperatorFile);
    end
end

% ----------------------------------------------------------------------------
function [EigenMat, EigenFile] = local_find_dirac_eigen(SurfaceFile, nModes, Tau, Variant)
% Scan the surface's Eigen child nodes for a Variant ('Dirac' or 'Dirac-Face')
% eigenbasis with matching Tau and at least nModes modes. Returns the loaded
% EigenMat and its file path, or ([], '') if none found.
    if nargin < 4 || isempty(Variant), Variant = 'Dirac'; end
    EigenMat = [];
    EigenFile = '';
    [sSubject, ~, iSurface] = bst_get('SurfaceFile', SurfaceFile);
    if isempty(sSubject) || isempty(iSurface) ...
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
            E = in_bst_eigen(nodes(k).FileName);
        catch
            continue;   % unreadable / missing entry
        end
        if ~isfield(E,'Variant') || ~strcmpi(E.Variant, Variant); continue; end
        if ~isfield(E,'nModes') || isempty(E.nModes) || E.nModes < nModes; continue; end
        if abs(local_eigen_tau(E, NaN) - Tau) > 1e-12;            continue; end
        EigenMat  = E;
        EigenFile = nodes(k).FileName;
        return;
    end
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
function v = local_default(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
