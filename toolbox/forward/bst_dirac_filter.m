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
%   .V             [nVert x 3] complex per-vertex 3-vector wavelet (zeros off the seed hemi)
%   .Helicity      [nVert x 1] real local helicity density (Re V x Im V) . Axis
%   .ChiralityAxis [1 x 3]     the (unit) helicity axis used
%
% The core (project -> scale -> chirality -> reconstruct) is delegated to
% bst_dirac_eigenmodes_filter: this synthesis is that filter applied to a
% directional vertex delta.
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
    Op  = in_bst_operator(Eig.OperatorFile);

    % --- check the seed vertex is in the eigen support ---
    if ~any(cellfun(@(gv) any(gv(:) == iVertex), Eig.GlobalVertices))
        error('bst_dirac_filter:vertexNotFound', 'Vertex %d is not in the Dirac eigen support.', iVertex);
    end

    % --- build the seed: a pure-imaginary quaternion vertex delta = the launch
    %     direction at iVertex (full-cortex 3-vector field, zero elsewhere) ---
    d = Direction(:).';  nd = norm(d);  if nd > 0, d = d/nd; end
    nVertTotal = double(max(cellfun(@(x) max(x(:)), Eig.GlobalVertices)));
    Jdelta = zeros(3*nVertTotal, 1);
    Jdelta(3*(iVertex-1) + (1:3)) = d(:);

    % --- scale kernel g(lambda) as a handle (full eigfilter library) ---
    if isempty(Filter)
        gfun = @(l) ones(size(l));
    else
        params = struct();
        if isfield(Filter,'Params') && ~isempty(Filter.Params), params = Filter.Params; end
        gfun = bst_eigfilter_kernel(Filter.Type, params);
    end

    % --- delegate project -> scale -> chirality -> reconstruct to the shared filter
    %     (single source of truth: bst_dirac_eigenmodes_filter; the wavelet is just
    %     that filter applied to a directional vertex delta) ---
    Vcol = bst_dirac_eigenmodes_filter(Eig, Op.Mass, Jdelta, 'custom', 'TransferFn', gfun, 'Chirality', Chirality);
    V = reshape(Vcol, 3, []).';               % [nVertTotal x 3] complex (imag-quaternion part)

    % --- helicity density about the chirality axis ---
    hAxis = [1 0 0];
    if ~isempty(Chirality) && isfield(Chirality,'Axis') && ~isempty(Chirality.Axis), hAxis = Chirality.Axis(:).'; end
    hAxis = hAxis / norm(hAxis);
    Hel = sum(cross(real(V), imag(V), 2) .* repmat(hAxis, nVertTotal, 1), 2);

    Out = struct('V', V, 'Helicity', Hel, 'ChiralityAxis', hAxis);
end
