function [JFilt, h, Coeffs] = bst_dirac_eigenmodes_filter(EigenMat, MassCell, J, FilterType, varargin)
% BST_DIRAC_EIGENMODES_FILTER: Spectral filtering of a cortical VECTOR field in the
% Dirac eigenmode domain (the vector analogue of bst_eigenmodes_filter).
%
% USAGE:  JFilt = bst_dirac_eigenmodes_filter(EigenMat, MassCell, J, 'heat', 'DiffusionTime', t)
%         JFilt = bst_dirac_eigenmodes_filter(EigenMat, MassCell, J, 'bandpass', 'ModeRange', [m1 m2])
%         JFilt = bst_dirac_eigenmodes_filter(..., 'Chirality', struct('Axis',[1 0 0],'Sign',+1))
%
% INPUTS:
%     EigenMat   : Dirac eigen node (from tess_eigen, Variant 'Dirac'). Uses the
%                  per-hemisphere cells .Phi{h} [4nVh x K] (B-orthonormal Dirac
%                  eigenvectors), .Lambda{h} [K x 1], .GlobalVertices{h}.
%     MassCell   : {B1, B2} per-hemisphere mass B_h = kron(Mass_h, I4) [4nVh x 4nVh]
%                  (Operator node .Mass; the inner product Phi is orthonormal in).
%     J          : [3*nVert x nTime] unconstrained source vector field (x,y,z per
%                  vertex), or [nVert x 3] for a single time sample.
%     FilterType : SCALE kernel on the Dirac eigenvalues — 'lowpass','highpass',
%                  'bandpass','heat','inverse_heat','tikhonov','custom'. The gain is
%                  computed by bst_eigenmodes_filter_gain (the shared eigfilter
%                  library), so every scalar-LBO kernel is available unchanged.
%
% OPTIONS:
%     Scale-kernel options (CutoffMode, ModeRange, DiffusionTime, MaxGain, RegBeta,
%        TransferFn) are forwarded verbatim to bst_eigenmodes_filter_gain.
%     'Chirality' : struct('Axis',[1x3], 'Sign',+1|-1) — apply the helicity projector
%        P_pm(n) = (I -/+ i*R_n)/2 about unit-imaginary axis Axis. R_n is right-
%        multiplication by the quaternion axis; it commutes with the Dirac operator
%        and squares to -I, so P_pm selects one circular-polarization (helicity).
%        Output is then COMPLEX (the analytic helicity field). Default [] = none.
%
% OUTPUT:
%     JFilt : [3*nVert x nTime] filtered vector field (same shape as J; complex if a
%             chirality projector was applied), zeros on vertices outside the eigen
%             support.
%     h     : {h1, h2} the per-hemisphere scale-gain vectors h(lambda) actually used.
%
% MATHEMATICAL NOTES (per hemisphere h):
%     psi = embed(J)                         pure-imaginary quaternion  [4nVh x nT]
%     c   = Phi_h' * B_h * psi               project onto the Dirac eigenbasis  [K x nT]
%     c   = diag(h(lambda_h)) * c            SCALE filter (h via the eigfilter library)
%     c   = P_pm(n) * c                      optional CHIRALITY/helicity projector
%     F   = Phi_h * c                        reconstruct  [4nVh x nT]
%     JFilt(hemi) = imaginary-quaternion part of F   (the w slot is dropped)
%
%     The scale gain h(lambda) is constant within each 4-fold quaternionic multiplet
%     (all four share lambda), so it commutes with the chirality projector. With
%     FilterType 'flat'/identity and no chirality, JFilt == J (B-orthonormal round trip).
%
% SEE ALSO: bst_eigenmodes_filter, bst_eigenmodes_filter_gain, bst_dirac_filter,
%           bst_dirac, tess_eigen, bst_eigfilter_kernel
%
% Authors: Diellor Basha, 2026

    % --- parse Chirality / ReturnCoeffs / Coeffs options; the rest go to the gain ---
    Chirality = [];  ReturnCoeffs = false;  CoeffsIn = [];  %#ok<NASGU>
    keep = true(1, numel(varargin));
    for i = 1:2:numel(varargin)
        if ischar(varargin{i})
            switch lower(varargin{i})
                case 'chirality',    Chirality    = varargin{i+1};          keep(i:i+1) = false;
                case 'returncoeffs', ReturnCoeffs = logical(varargin{i+1}); keep(i:i+1) = false; %#ok<NASGU>
                case 'coeffs',       CoeffsIn     = varargin{i+1};          keep(i:i+1) = false;
            end
        end
    end
    gainArgs = varargin(keep);

    % --- normalize J to [3*nVert x nTime] (J may be empty when Coeffs are supplied) ---
    nVert = double(max(cellfun(@(x) max(x(:)), EigenMat.GlobalVertices)));
    if (size(J,2) == 3) && (size(J,1) == nVert)
        J = reshape(J.', [], 1);                 % [nVert x 3] single frame -> [3nVert x 1]
    end
    if ~isempty(CoeffsIn)
        nTime = size(CoeffsIn{1}, 2);            % number of columns is set by the cached coeffs
    else
        nTime = size(J, 2);
    end
    JFilt = zeros(3*nVert, nTime);

    nHemi = numel(EigenMat.Phi);
    h = cell(1, nHemi);
    Coeffs = cell(1, nHemi);
    for hh = 1:nHemi
        Phi = double(EigenMat.Phi{hh});          % [4nVh x K]
        B   = MassCell{hh};                       % [4nVh x 4nVh]
        lam = EigenMat.Lambda{hh}(:);             % [K x 1]
        gv  = EigenMat.GlobalVertices{hh}(:);
        nVh = size(Phi, 1) / 4;

        % --- project (or use cached coeffs), SCALE-filter, optional CHIRALITY, reconstruct ---
        if ~isempty(CoeffsIn)
            c = CoeffsIn{hh};                              % cached projection (skip embed + project)
        else
            % embed the hemisphere's 3-vectors as pure-imaginary quaternions, then project
            Jx = J(3*(gv-1)+1, :);  Jy = J(3*(gv-1)+2, :);  Jz = J(3*(gv-1)+3, :);
            Psi = zeros(4*nVh, nTime);
            Psi(2:4:end, :) = Jx;  Psi(3:4:end, :) = Jy;  Psi(4:4:end, :) = Jz;   % w rows = 0
            c = Phi' * (B * Psi);                          % [K x nTime]
        end
        Coeffs{hh} = c;
        h{hh} = bst_eigenmodes_filter_gain(lam, FilterType, gainArgs{:});   % [K x 1]
        c = bsxfun(@times, h{hh}, c);

        if ~isempty(Chirality) && isfield(Chirality,'Sign') && ~isempty(Chirality.Sign) && (Chirality.Sign ~= 0)
            a = [1 0 0];
            if isfield(Chirality,'Axis') && ~isempty(Chirality.Axis), a = Chirality.Axis(:).'; end
            a = a / norm(a);
            rhoR = [0 -a(1) -a(2) -a(3); a(1) 0 a(3) -a(2); a(2) -a(3) 0 a(1); a(3) a(2) -a(1) 0];
            Rfield = kron(speye(nVh), sparse(rhoR));        % right-mult by axis a
            Rmode  = Phi' * (B * (Rfield * Phi));           % [K x K], block-diag per multiplet
            c = 0.5 * (c - 1i * sign(Chirality.Sign) * (Rmode * c));
        end

        F = Phi * c;                                        % [4nVh x nTime]
        JFilt(3*(gv-1)+1, :) = F(2:4:end, :);               % x  (drop w = rows 1:4:end)
        JFilt(3*(gv-1)+2, :) = F(3:4:end, :);               % y
        JFilt(3*(gv-1)+3, :) = F(4:4:end, :);               % z
    end
end
