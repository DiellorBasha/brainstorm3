function varargout = bst_eigenwavelet(varargin)
% BST_EIGENWAVELET: Spectral graph wavelet orchestrator on an operator eigenbasis.
%
% USAGE:
%   frame = bst_eigenwavelet('Design', family, Nf, lrange)   % family: 'mexhat'|'heat'|'itersine'
%   H     = bst_eigenwavelet('Evaluate', frame, Lambda)      % [K x M] gains
%   b     = bst_eigenwavelet('Bounds', frame, Lambda)        % struct A,B,Tightness
%   [W,Messages,isError]    = bst_eigenwavelet('Analysis', F, EigenMat, OperatorMat, frame)
%   [Frec,Messages,isError] = bst_eigenwavelet('Synthesis', W, EigenMat, OperatorMat, frame)
%   [A,Messages,isError]    = bst_eigenwavelet('Atom', EigenMat, OperatorMat, frame, seedVert, seedDir)
%   Wq    = bst_eigenwavelet('Steer', W, q0, EigenMat)       % right-quaternion steer (Dirac family)
%   V     = bst_eigenwavelet('ToVec',  W, EigenMat)          % full-quaternion -> physical 3-vector
%   W     = bst_eigenwavelet('ToQuat', V, EigenMat)          % physical 3-vector -> pure quaternion (w=0)
%
% DESCRIPTION:
%     GSPBox-style spectral graph wavelet frame (Perraudin et al., GSPBOX, arXiv:1408.5781,
%     2014; Hammond, Vandergheynst & Gribonval, "Wavelets on graphs via spectral graph
%     theory", Applied and Computational Harmonic Analysis 30(2):129-150, 2011). Adapted to
%     Brainstorm's exact precomputed eigenbasis: filters are applied as Phi*diag(g(lambda))*Phi'
%     (NOT the Chebyshev polynomial approximation GSPBox uses to avoid eigendecomposition).
%
%     A wavelet transform is a multi-member FRAME {g_1(lambda)..g_M(lambda)}; a single filter
%     (bst_eigenfilter) is the 1-member case of the same machinery. The per-variant row mapping
%     and gain evaluation are shared with bst_eigenfilter ('RowMap'/'Design'/'Evaluate').
%
%     Families:
%       'itersine' : tight half-cosine frame (sum_m g_m^2 = const -> exact reconstruction)
%       'mexhat'   : Nf log-scaled band-pass Mexican-hat wavelets + 1 low-pass scaling function
%       'heat'     : Nf low-pass heat scales (diffusion scale-space, coarse -> fine)
%
%     FULL-QUATERNION DIRAC: for the quaternionic Dirac family ('Dirac','Dirac-Face',
%     'Hodge-Face') the transform carries the FULL quaternion field [4*nV x nT x M]
%     end-to-end -- the real (w) component is NOT discarded -- so Analysis->Synthesis is
%     exact and atoms steer EXACTLY by RIGHT quaternion multiplication (the 4-fold mode
%     multiplicity = the orientation freedom). Physical input is the 3-vector dipole map
%     [3*nV x nT] (embedded as a pure quaternion, w=0); use ToVec to project the
%     full-quaternion result back to a 3-vector for display. Scalar (Laplace-Beltrami) and
%     complex (Connection Laplacian) variants keep their native field layout unchanged.
%
% SEE ALSO: bst_eigen, bst_eigenfilter, bst_eigfilter_kernel, manifold_ft, manifold_ift

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

eval(macro_method);
end


%% ===== DESIGN: assemble a frame (cell of handles) =====
function frame = Design(family, Nf, lrange) %#ok<DEFNU>
    if (nargin < 2) || isempty(Nf);     Nf = 6; end
    if (nargin < 3) || isempty(lrange); error('bst_eigenwavelet:Design','lrange = [lmin lmax] is required.'); end
    lmin = lrange(1); lmax = lrange(2);
    if ~(lmax > 0) || (lmax <= lmin);   error('bst_eigenwavelet:Design','invalid lrange = [%g %g].', lmin, lmax); end
    switch lower(family)
        case 'itersine'
            overlap = 2;
            scale = lmax / (Nf - overlap + 1) * overlap;
            kf = @(x) sin(0.5*pi*(cos(pi*x)).^2) .* (x>=-0.5 & x<=0.5);
            g = cell(1, Nf);
            for ii = 1:Nf
                g{ii} = @(l) kf(double(l(:))/scale - (ii-overlap/2)/overlap) ./ sqrt(overlap) .* sqrt(2);
            end
        case 'mexhat'
            % Nf log-spaced band-pass mexhat wavelets + 1 low-pass scaling function
            tmin = 1/lmax;
            tmax = 2/max(lmin, lmax/(4*Nf));
            tw   = logspace(log10(tmin), log10(tmax), Nf);
            g = cell(1, Nf+1);
            g{1} = bst_eigfilter_kernel('heat', struct('t', 2*tmax));   % scaling fn (low-pass)
            for ii = 1:Nf
                g{ii+1} = bst_eigfilter_kernel('mexhat', struct('t', tw(ii)));
            end
        case 'heat'
            % Nf low-pass heat scales (diffusion scale-space), coarse -> fine
            tmin = 1/lmax;
            tmax = 4/max(lmin, lmax/(4*Nf));
            tw   = logspace(log10(tmax), log10(tmin), Nf);   % coarse (large t) first
            g = cell(1, Nf);
            for ii = 1:Nf
                g{ii} = bst_eigfilter_kernel('heat', struct('t', tw(ii)));
            end
        otherwise
            error('bst_eigenwavelet:Design','unknown family ''%s'' (use mexhat|heat|itersine).', family);
    end
    % Per-member characteristic spatial frequency = gain-weighted centroid of sqrt(lambda)
    lg = linspace(max(lmin, eps), lmax, 512)';
    M  = numel(g); cen = zeros(1, M);
    for m = 1:M
        gm = abs(g{m}(lg));
        if sum(gm) > 0; cen(m) = sqrt( sum(lg .* gm) / sum(gm) ); end
    end
    frame = struct('Family', lower(family), 'Nf', Nf, 'g', {g}, 'Lrange', [lmin lmax], 'Centers', cen);
end


%% ===== EVALUATE: frame -> gains on Lambda =====
function H = Evaluate(frame, Lambda) %#ok<DEFNU>
    Lambda = double(Lambda(:));
    M = numel(frame.g);
    H = zeros(numel(Lambda), M);
    for m = 1:M
        v = frame.g{m}(Lambda);
        H(:, m) = v(:);
    end
end


%% ===== BOUNDS: frame bounds A,B and tightness B/A =====
function b = Bounds(frame, Lambda) %#ok<DEFNU>
    H = Evaluate(frame, Lambda);
    S = sum(H.^2, 2);
    A = min(S); B = max(S);
    b = struct('A', A, 'B', B, 'Tightness', B / max(A, eps));
end


%% ===== ANALYSIS: map -> scalogram =====
% Scalar/complex variants: W has the SAME row layout as the input F.
% Dirac family: input F is the physical 3-vector map [3*nV x nT]; W is the FULL
% quaternion scalogram [4*nV x nT x M] (real w retained -> exact, steerable).
function [W, Messages, isError] = Analysis(F, EigenMat, OperatorMat, frame) %#ok<DEFNU>
    Messages = ''; isError = 0;
    M  = numel(frame.g);
    nT = size(F, 2);
    isQuat = i_isQuat(EigenMat);
    if isQuat
        W = zeros(4 * i_maxidx(EigenMat), nT, M);   % full-quaternion global field
    else
        W = zeros(size(F, 1), nT, M);
        if ~isreal(F); W = complex(W); end
    end
    for h = 1:numel(EigenMat.Phi)
        Phi = EigenMat.Phi{h};
        if isempty(Phi); continue; end
        Lam = EigenMat.Lambda{h}(:);
        B   = OperatorMat.Mass{h};
        map = i_hemimap(EigenMat, h);
        if max(map.gIn) > size(F, 1)
            Messages = sprintf('bst_eigenwavelet:Analysis hemi %d needs input row %d but F has %d rows.', h, max(map.gIn), size(F,1));
            isError = 1; return;
        end
        if ~isreal(Phi); W = complex(W); end        % connection -> complex storage
        H = Evaluate(frame, Lam);                    % [K x M]
        U = zeros(map.nrows, nT);  if ~isreal(Phi); U = complex(U); end
        U(map.lIn, :) = F(map.gIn, :);
        C = manifold_ft(Phi, B, U);                  % [K x nT]
        for m = 1:M
            Uf = manifold_ift(Phi, H(:, m) .* C);
            W(map.gOut, :, m) = Uf(map.lOut, :);
        end
    end
end


%% ===== SYNTHESIS: scalogram -> map (sum of g_m(L) c_m) =====
% Inverse of Analysis. For the Dirac family W and Frec are FULL quaternion [4*nV x .];
% apply ToVec to obtain the physical 3-vector for display.
function [Frec, Messages, isError] = Synthesis(W, EigenMat, OperatorMat, frame) %#ok<DEFNU>
    Messages = ''; isError = 0;
    M  = numel(frame.g);
    nT = size(W, 2);
    Frec = zeros(size(W, 1), nT);
    if ~isreal(W); Frec = complex(Frec); end
    for h = 1:numel(EigenMat.Phi)
        Phi = EigenMat.Phi{h};
        if isempty(Phi); continue; end
        Lam = EigenMat.Lambda{h}(:);
        B   = OperatorMat.Mass{h};
        map = i_hemimap(EigenMat, h);
        if max(map.gOut) > size(W, 1)
            Messages = sprintf('bst_eigenwavelet:Synthesis hemi %d needs row %d but W has %d rows.', h, max(map.gOut), size(W,1));
            isError = 1; return;
        end
        H = Evaluate(frame, Lam);
        acc = zeros(map.nrows, nT);  if ~isreal(Phi); acc = complex(acc); end
        for m = 1:M
            U = zeros(map.nrows, nT);  if ~isreal(Phi); U = complex(U); end
            U(map.lOut, :) = W(map.gOut, :, m);
            acc = acc + manifold_ift(Phi, H(:, m) .* manifold_ft(Phi, B, U));
        end
        Frec(map.gOut, :) = acc(map.lOut, :);
    end
end


%% ===== ATOM: wavelet atoms at a seed (filtered delta) =====
% seedDir: scalar (Laplace-Beltrami), complex phase (Connection Laplacian), or 3-vector
% dipole (Dirac family). Returns A = the per-member atoms [fieldrows x 1 x M], full
% quaternion for the Dirac family. Steer/ToVec post-process for display.
function [A, Messages, isError] = Atom(EigenMat, OperatorMat, frame, seedVert, seedDir) %#ok<DEFNU>
    A = [];   % Messages/isError are set by the 'otherwise' branch or by Analysis below
    nV = i_maxidx(EigenMat);
    switch EigenMat.Variant
        case 'Laplace-Beltrami'
            if nargin < 5 || isempty(seedDir); seedDir = 1; end
            F = zeros(nV, 1);  F(seedVert) = seedDir;
        case 'Connection Laplacian'
            if nargin < 5 || isempty(seedDir); seedDir = 1; end
            F = zeros(nV, 1);  F(seedVert) = seedDir;  F = complex(F);
        case {'Dirac', 'Dirac-Face', 'Hodge-Face'}
            if nargin < 5 || isempty(seedDir); seedDir = [1; 0; 0]; end
            seedDir = seedDir(:);
            F = zeros(3 * nV, 1);  F((seedVert-1)*3 + (1:3)) = seedDir(1:3);
        otherwise
            Messages = sprintf('bst_eigenwavelet:Atom unsupported variant ''%s''.', EigenMat.Variant);
            isError = 1; return;
    end
    [A, Messages, isError] = Analysis(F, EigenMat, OperatorMat, frame);
end


%% ===== STEER: right-quaternion multiply a full-quaternion field/scalogram =====
function Wq = Steer(W, q0, EigenMat) %#ok<DEFNU>
    if ~i_isQuat(EigenMat)
        error('bst_eigenwavelet:Steer', 'steering applies to the Dirac family only.');
    end
    q0 = q0(:);  q0 = q0 / max(norm(q0), eps);
    sz = size(W);
    Rr = i_qrightmat(q0);                 % 4x4 right-multiplication matrix
    Wq = reshape(Rr * reshape(W, 4, []), sz);
end


%% ===== ToVec / ToQuat: full quaternion <-> physical 3-vector =====
function V = ToVec(W, EigenMat) %#ok<DEFNU>
    if ~i_isQuat(EigenMat); V = W; return; end
    sz = size(W);  nV = sz(1) / 4;
    Wr = reshape(W, 4, nV, []);
    V  = reshape(Wr(2:4, :, :), [3*nV, sz(2:end)]);
end

function W = ToQuat(V, EigenMat) %#ok<DEFNU>
    if ~i_isQuat(EigenMat); W = V; return; end
    sz = size(V);  nV = sz(1) / 3;
    Vr = reshape(V, 3, nV, []);
    Wr = zeros([4, nV, size(Vr, 3)]);
    Wr(2:4, :, :) = Vr;
    W  = reshape(Wr, [4*nV, sz(2:end)]);
end


%% ===== internal: per-hemi row mapping (variant-aware) =====
function tf = i_isQuat(EigenMat)
    tf = any(strcmp(EigenMat.Variant, {'Dirac', 'Dirac-Face', 'Hodge-Face'}));
end

function n = i_maxidx(EigenMat)
    if any(strcmp(EigenMat.Variant, {'Dirac-Face', 'Hodge-Face'}))
        idxc = EigenMat.GlobalFaces;
    else
        idxc = EigenMat.GlobalVertices;
    end
    n = 0;
    for h = 1:numel(idxc)
        if ~isempty(idxc{h}); n = max(n, max(idxc{h}(:))); end
    end
end

function map = i_hemimap(EigenMat, h)
% Per-hemisphere index sets bridging the global field and the local (Phi-row) ordering.
%   gIn/lIn  : embed the INPUT (3-vector for Dirac, native otherwise) into the local field
%   gOut/lOut: scatter the OUTPUT (FULL quaternion for Dirac, native otherwise) to global
    switch EigenMat.Variant
        case {'Laplace-Beltrami', 'Connection Laplacian'}
            gv = EigenMat.GlobalVertices{h}(:);
            map.gIn = gv;  map.lIn = (1:numel(gv))';
            map.gOut = gv; map.lOut = (1:numel(gv))';
            map.nrows = numel(gv);
        case {'Dirac', 'Dirac-Face', 'Hodge-Face'}
            if strcmp(EigenMat.Variant, 'Dirac')
                idx = EigenMat.GlobalVertices{h}(:);
            else
                idx = EigenMat.GlobalFaces{h}(:);
            end
            n = numel(idx);
            map.gIn  = reshape([(idx-1)*3+1, (idx-1)*3+2, (idx-1)*3+3].', [], 1);          % global 3-vec rows
            map.lIn  = reshape([(0:n-1)*4+2; (0:n-1)*4+3; (0:n-1)*4+4], [], 1);            % local quat imag slots
            map.gOut = reshape([(idx-1)*4+1, (idx-1)*4+2, (idx-1)*4+3, (idx-1)*4+4].', [], 1); % global quat rows
            map.lOut = (1:4*n)';                                                           % local quat (full)
            map.nrows = 4 * n;
        otherwise
            error('bst_eigenwavelet:i_hemimap', 'unsupported variant ''%s''.', EigenMat.Variant);
    end
end

%% ===== internal: quaternion right-multiplication matrix =====
function R = i_qrightmat(p)
% 4x4 matrix R such that R*q == q (x) p  (Hamilton product), q,p,result = [w x y z].
    pw = p(1); px = p(2); py = p(3); pz = p(4);
    R = [ pw, -px, -py, -pz; ...
          px,  pw,  pz, -py; ...
          py, -pz,  pw,  px; ...
          pz,  py, -px,  pw ];
end
