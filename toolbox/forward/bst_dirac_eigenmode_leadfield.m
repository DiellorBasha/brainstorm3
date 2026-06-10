function CompHM = bst_dirac_eigenmode_leadfield(HeadModel, DiracEigen, varargin)
% BST_DIRAC_EIGENMODE_LEADFIELD: Compose the UNCONSTRAINED leadfield into the Dirac eigenbasis.
%
% USAGE:  CompHM = bst_dirac_eigenmode_leadfield(HeadModel, DiracEigen, 'nModes', K)
%
% DESCRIPTION:
%     Forward step for Dirac eigenmode source mapping. Unlike the scalar LBO
%     composer (which constrains the leadfield to the surface normal), this keeps
%     the full unconstrained 3-vector gain and expands it in the curvature-aware
%     Dirac eigenbasis. Each per-vertex gain 3-vector (ambient/world coords) is
%     embedded as a pure-imaginary quaternion psi = [0, gx, gy, gz] and projected
%     onto the per-hemisphere Dirac eigenvectors:
%         L~_h = Psi_h' * B_h * Phi_h      [nCh x K],   B_h = kron(Mass_h, I4)
%     The two hemispheres are stacked into a composed head model
%         Gain = [L~_L, L~_R]              [nCh x 2K]
%     consumed unchanged by bst_inverse_eigenmodes.
%
% INPUTS:
%     HeadModel  : base UNCONSTRAINED surface head model (in_bst_headmodel, ApplyOrient=0):
%                  .Gain [nCh x 3*nVert], .HeadModelType, .SurfaceFile, .Comment
%     DiracEigen : 1x2 per-hemisphere struct (TessMat.DiracEigen, from tess_dirac_eigenmodes):
%                  .Vectors [4Vh x K], .Values [K x 1], .Mass [Vh x Vh],
%                  .nModes, .Tau, .GlobalVertices
% OPTIONS:
%     'nModes' : modes per hemisphere to keep (default all; clamped to available)
%
% OUTPUT:
%     CompHM : composed head-model struct (Gain [nCh x 2K], Eigenvalues [2K x 1],
%              nModes=2K, ModeHemisphere [2K x 1], HemiGlobalVertices {L,R},
%              isEigenmode/isDiracEigenmode flags) ready to save / feed the inverse.
%
% Authors: Diellor Basha, 2026

    nModes = [];
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'nmodes', nModes = varargin{i+1};
        end
    end

    % surface head model only
    if isfield(HeadModel,'HeadModelType') && ~isempty(HeadModel.HeadModelType) ...
            && ~strcmpi(HeadModel.HeadModelType,'surface')
        error('bst_dirac_eigenmode_leadfield:NotSurface', ...
            'Dirac eigenmode leadfield requires a surface head model (got ''%s'').', HeadModel.HeadModelType);
    end

    G = double(HeadModel.Gain);          % [nCh x 3*nVert] unconstrained
    nCh = size(G,1);
    if mod(size(G,2), 3) ~= 0
        error('bst_dirac_eigenmode_leadfield:NotUnconstrained', ...
            'Dirac leadfield requires an unconstrained head model: Gain must be [nCh x 3*nVert].');
    end

    if numel(DiracEigen) ~= 2
        error('bst_dirac_eigenmode_leadfield:badBasis', ...
            'DiracEigen must be a 1x2 per-hemisphere struct array (from tess_dirac_eigenmodes).');
    end
    if DiracEigen(1).Tau ~= DiracEigen(2).Tau
        error('bst_dirac_eigenmode_leadfield:tauMismatch', ...
            'Hemispheres have different Tau (%.3g vs %.3g); recompute DiracEigen with a single Tau.', ...
            DiracEigen(1).Tau, DiracEigen(2).Tau);
    end

    Kfull = min([DiracEigen.nModes]);
    if isempty(nModes) || nModes <= 0
        K = Kfull;
    else
        K = min(nModes, Kfull);
    end

    Lblk = cell(1,2); vblk = cell(1,2); hblk = cell(1,2);
    for hh = 1:2
        vH  = DiracEigen(hh).GlobalVertices(:);
        nVh = numel(vH);
        Phi = double(DiracEigen(hh).Vectors);
        if size(Phi,1) ~= 4*nVh
            error('bst_dirac_eigenmode_leadfield:shapeMismatch', ...
                'Hemisphere %d: Vectors has %d rows, expected 4*nV=%d.', hh, size(Phi,1), 4*nVh);
        end
        if size(Phi,2) < K
            error('bst_dirac_eigenmode_leadfield:shapeMismatch', ...
                'Hemisphere %d: Vectors has %d columns, fewer than K=%d.', hh, size(Phi,2), K);
        end
        Phi  = Phi(:, 1:K);
        Vals = double(DiracEigen(hh).Values(:)); Vals = Vals(1:K);

        % embed unconstrained gain as a pure-imaginary quaternion field (w=0)
        Psi = zeros(4*nVh, nCh);
        Psi(2:4:end, :) = G(:, (vH-1)*3 + 1).';   % i (x)
        Psi(3:4:end, :) = G(:, (vH-1)*3 + 2).';   % j (y)
        Psi(4:4:end, :) = G(:, (vH-1)*3 + 3).';   % k (z)   (w rows 1:4:end stay 0)

        B = kron(DiracEigen(hh).Mass, speye(4));   % [4Vh x 4Vh]
        Lblk{hh} = Psi' * (B * Phi);               % [nCh x K]
        vblk{hh} = Vals;
        hblk{hh} = hh * ones(K,1);
    end

    CompHM = HeadModel;
    CompHM.Gain               = [Lblk{1}, Lblk{2}];     % [nCh x 2K]
    CompHM.GridLoc            = [];
    CompHM.GridOrient         = [];
    CompHM.GridAtlas          = [];
    CompHM.isEigenmode        = 1;
    CompHM.isDiracEigenmode   = 1;
    CompHM.nModes             = 2*K;
    CompHM.Eigenvalues        = [vblk{1}; vblk{2}];     % [2K x 1]
    CompHM.ModeHemisphere     = [hblk{1}; hblk{2}];     % [2K x 1] hemisphere index per mode
    CompHM.HemiGlobalVertices = {DiracEigen(1).GlobalVertices(:), DiracEigen(2).GlobalVertices(:)};
    CompHM.HeadModelType      = 'surface';
    CompHM.Comment            = sprintf('Dirac eigenmode leadfield (%d modes, tau=%.3g) | %s', ...
        2*K, DiracEigen(1).Tau, local_default(HeadModel,'Comment',''));
end

function v = local_default(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
