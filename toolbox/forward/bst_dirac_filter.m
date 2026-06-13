function Out = bst_dirac_filter(SurfaceFile, iVertex, Direction, Filter, Chirality, varargin)
% BST_DIRAC_FILTER: Dirac eigenmode vector wavelet (scale x direction x chirality).
%
% USAGE:
%   Out = bst_dirac_filter(SurfaceFile, iVertex, Direction, Filter, Chirality)
%   Out = bst_dirac_filter(..., 'nModes', K, 'Tau', tau)
%
% DESCRIPTION:
%   Builds a vector wavelet on the cortex by reconstructing a QUATERNION vertex
%   delta through a kernel defined on the Dirac eigenvalue axis, with an optional
%   chirality (helicity) projector. The Dirac structure exposes three independent
%   filter axes (see bst_dirac, tess_operators):
%     1. SCALE     - a kernel g(lambda) on the eigenvalues (heat, mexhat, ...),
%                    exactly as in the scalar Laplace-Beltrami filterbank.
%     2. DIRECTION - the seed is a pure-imaginary quaternion (0, Direction) at the
%                    delta vertex; the spin-connection in the Dirac operator
%                    parallel-transports it along the curved cortical frame.
%     3. CHIRALITY - right-multiplication R_n by a unit-imaginary axis n commutes
%                    with the Dirac operator EXACTLY (verified) and squares to -I,
%                    so P_pm = (I -/+ i*R_n)/2 splits every eigen-quartet into its
%                    two helicities. P_pm yields a circularly-polarized vector
%                    wavelet rotating about n, with handedness set by the sign.
%
%   The pipeline (per hemisphere h holding the seed):
%       psi = delta_v (x) (0, Direction)            % pure-imaginary quaternion delta
%       c   = Phi_h' * B_h * psi                    % project onto Dirac eigenbasis
%       c   = diag(g(lambda)) * c                   % scale kernel
%       c   = P_pm(n) * c                           % optional chirality (complex)
%       F   = Phi_h * c                             % reconstruct quaternion field
%       V   = imag-quaternion part of F (drop w)    % per-vertex 3-vector (complex)
%
% INPUTS:
%   SurfaceFile : cortex surface; the 'Dirac' eigen node is found-or-loaded (tess_eigen)
%   iVertex     : GLOBAL vertex index of the delta seed
%   Direction   : [1x3] launch direction (ambient coords); normalized internally
%   Filter      : struct('Type', name, 'Params', paramStruct) for the eigfilter library
%                 (e.g. struct('Type','heat','Params',struct('t',8e4))). [] => flat (g=1).
%   Chirality   : struct('Axis',[1x3], 'Sign', +1|-1) selecting a helicity about a
%                 unit-imaginary axis. [] or Sign==0 => no projector (real field).
% OPTIONS:
%   'nModes' (default 400), 'Tau' (default 0.5) : forwarded to tess_eigen if a Dirac
%                                                 node must be created.
%
% OUTPUT struct Out:
%   .V        [nVert x 3] complex per-vertex 3-vector wavelet (zeros off the seed hemi)
%   .W        [nVert x 1] complex scalar (w) channel (spin-connection leakage)
%   .Helicity [nVert x 1] real local helicity density (Re V x Im V) . Axis
%   .Coeffs   [K x 1]     final (complex) mode coefficients
%   .Lambda   [K x 1]     Dirac eigenvalues of the seed hemisphere
%   .iHemi    1|2         hemisphere holding the seed
%   .Gain     scalar      filter gain vector summary (g at lambda)
%
% SEE ALSO: bst_dirac, tess_eigen, tess_operators, bst_eigfilter_kernel
%
% Authors: Diellor Basha, 2026

    % --- options ---
    nModes = 400;  Tau = 0.5;
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'nmodes', nModes = varargin{i+1};
            case 'tau',    Tau    = varargin{i+1};
            otherwise, error('bst_dirac_filter:badOption', 'Unknown option: %s', varargin{i});
        end
    end
    if (nargin < 5), Chirality = []; end
    if (nargin < 4), Filter    = []; end

    % --- Dirac eigenbasis + mass (found-or-loaded) ---
    Eig = tess_eigen(SurfaceFile, 'Dirac', 'K', nModes, 'Tau', Tau);
    Op  = load(file_fullpath(Eig.OperatorFile));

    % --- locate the seed vertex's hemisphere ---
    iHemi = 0;  vloc = 0;
    for h = 1:numel(Eig.GlobalVertices)
        loc = find(Eig.GlobalVertices{h}(:) == iVertex, 1);
        if ~isempty(loc), iHemi = h; vloc = loc; break; end
    end
    if iHemi == 0
        error('bst_dirac_filter:vertexNotFound', 'Vertex %d is not in the Dirac eigen support.', iVertex);
    end

    Phi = double(Eig.Phi{iHemi});           % [4nVh x K]
    B   = Op.Mass{iHemi};                    % [4nVh x 4nVh] = kron(Mass_h, I4)
    lam = Eig.Lambda{iHemi}(:);              % [K x 1]
    gv  = Eig.GlobalVertices{iHemi}(:);
    K   = size(Phi,2);  nVh = size(Phi,1)/4;

    % --- pure-imaginary quaternion delta at the seed vertex ---
    d = Direction(:).';  nd = norm(d);  if nd > 0, d = d/nd; end
    psi = zeros(4*nVh, 1);
    psi(4*(vloc-1) + (2:4)) = d(:);          % rows [x;y;z] of block vloc (w stays 0)
    c = Phi' * (B * psi);                     % [K x 1] mode coefficients (real)

    % --- scale kernel g(lambda) ---
    if isempty(Filter)
        g = ones(K,1);
    else
        params = struct();
        if isfield(Filter,'Params') && ~isempty(Filter.Params), params = Filter.Params; end
        gfun = bst_eigfilter_kernel(Filter.Type, params);
        g    = bst_eigfilter_evaluate(gfun, lam);
    end
    c = g .* c;

    % --- optional chirality (helicity) projector P_pm(n) ---
    hAxis = [1 0 0];
    if ~isempty(Chirality) && isfield(Chirality,'Sign') && ~isempty(Chirality.Sign) && (Chirality.Sign ~= 0)
        if isfield(Chirality,'Axis') && ~isempty(Chirality.Axis), hAxis = Chirality.Axis(:).'; end
        a = hAxis / norm(hAxis);
        rhoR = [0 -a(1) -a(2) -a(3); a(1) 0 a(3) -a(2); a(2) -a(3) 0 a(1); a(3) a(2) -a(1) 0];
        Rfield = kron(speye(nVh), sparse(rhoR));          % [4nVh x 4nVh] right-mult by a
        Rmode  = Phi' * (B * (Rfield * Phi));             % [K x K] in the eigenbasis
        c = 0.5 * (c - 1i * sign(Chirality.Sign) * (Rmode * c));
    end

    % --- reconstruct quaternion field, split into vector (x,y,z) and scalar (w) ---
    F  = Phi * c;                             % [4nVh x 1] complex
    Vh = [F(2:4:end), F(3:4:end), F(4:4:end)];% [nVh x 3] imaginary part (the 3-vector)
    Wh = F(1:4:end);                          % [nVh x 1] scalar (w) channel

    % --- scatter to full-surface arrays ---
    nVertTotal = double(max(cellfun(@(x) max(x(:)), Eig.GlobalVertices)));
    V = zeros(nVertTotal, 3);  W = zeros(nVertTotal, 1);
    V(gv, :) = Vh;  W(gv) = Wh;
    Hel = zeros(nVertTotal, 1);
    Hel(gv) = sum(cross(real(Vh), imag(Vh), 2) .* repmat(hAxis/norm(hAxis), nVh, 1), 2);

    Out = struct('V', V, 'W', W, 'Helicity', Hel, 'Coeffs', c, 'Lambda', lam, ...
                 'iHemi', iHemi, 'iVertexLocal', vloc, 'Gain', g, 'ChiralityAxis', hAxis/norm(hAxis));
end
