function varargout = bst_eigenfilter(varargin)
% BST_EIGENFILTER: Spatial filter orchestrator in an operator eigenbasis (1-member frame).
%
% USAGE:
%   g     = bst_eigenfilter('Design', KernelName, KernelParams)
%   h     = bst_eigenfilter('Evaluate', g, Lambda)
%   [Ffilt,Messages,isError] = bst_eigenfilter('Analysis', F, EigenMat, OperatorMat, KernelName, KernelParams)
%   [Frec, Messages,isError] = bst_eigenfilter('Synthesis', C, EigenMat, OperatorMat, KernelName, KernelParams)
%   [srcRows,dstRows,nrows,msg] = bst_eigenfilter('RowMap', F, EigenMat, h)
%
% DESCRIPTION:
%     The spatial-domain analogue of a bandpass filter, applied in an operator eigenbasis.
%     'Analysis' projects a surface-mapped source field onto the modes, scales each mode by a
%     gain h(lambda) from the eigfilter library, and reconstructs:
%         C = Phi'*(B*U);  C <- h(lambda).*C;  U_filt = Phi*C
%     General over the operator family carried in the eigen_ node (EigenMat.Variant): Laplace-
%     Beltrami (scalar), Connection Laplacian (complex tangent), Dirac (3D vector as a pure-
%     imaginary quaternion). LH and RH are filtered SEPARATELY (the basis is B-orthonormal per
%     hemisphere). The per-variant source<->basis row mapping is exposed as the 'RowMap' verb
%     so bst_eigenwavelet reuses it.
%
%     GSPBox-style frame framework (Perraudin et al., GSPBOX, arXiv:1408.5781, 2014; Hammond,
%     Vandergheynst & Gribonval, "Wavelets on graphs via spectral graph theory", ACHA
%     30(2):129-150, 2011), adapted to Brainstorm's exact precomputed eigenbasis (filters
%     applied as Phi*diag(g(lambda))*Phi', not Chebyshev). A single filter is the 1-member case
%     of the bst_eigenwavelet multi-member frame.
%
% SEE ALSO: bst_eigen, bst_eigenwavelet, bst_eigenspectrum, bst_eigfilter_kernel, manifold_ft, manifold_ift

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


%% ===== DESIGN: name/params -> handle =====
function g = Design(KernelName, KernelParams) %#ok<DEFNU>
    if (nargin < 2) || isempty(KernelParams); KernelParams = struct(); end
    if isa(KernelName, 'function_handle') || isnumeric(KernelName)
        g = KernelName;   % already a handle / precomputed gain
    else
        g = bst_eigfilter_kernel(KernelName, KernelParams);
    end
end


%% ===== ATOM: realise ONE localised atom (static / ts / js), domain-aware =====
function [W, gv] = Atom(ax, KernelName, KernelParams, seedVert, seedDir) %#ok<DEFNU>
    % Drop a unit delta at seedVert and propagate it through the kernel over the eigenbasis.
    % Domain-aware (kernel registry 'domain' tag): static g(lambda) (constant in time) | ts g(lambda,t)
    % | js g(lambda,omega) (realised by inverse JOINT time-vertex transform). ax = bst_eigen('Axes', ...).
    % seedDir (optional): scalar amplitude | tangent complex amplitude | Dirac 3-vector dipole direction.
    if (nargin < 4) || isempty(seedVert), error('bst_eigenfilter(''Atom''): seedVert required.'); end
    if (nargin < 5); seedDir = []; end
    blk = 0;
    for h = 1:numel(ax.GlobalVertices)
        if ~isempty(ax.GlobalVertices{h}) && any(ax.GlobalVertices{h} == seedVert), blk = h; break; end
    end
    if blk == 0, error('bst_eigenfilter(''Atom''): seed vertex %d not in the eigenbasis support.', seedVert); end
    Phi = ax.Phi{blk};  Lam = ax.Lambda{blk};  M = ax.Mass{blk};  gv = ax.GlobalVertices{blk};
    [C, kind] = i_fiber(ax); %#ok<ASGLU>
    nSrc = 0; for hh=1:numel(ax.GlobalVertices), nSrc = max(nSrc, max(ax.GlobalVertices{hh}(:))); end
    switch kind
        case 'scalar'
            if isempty(seedDir), seedDir = 1; end
            F = zeros(nSrc, 1);       F(seedVert) = seedDir(1);
        case 'tangent'
            if isempty(seedDir), seedDir = 1; end
            F = complex(zeros(nSrc,1)); F(seedVert) = seedDir(1);           % complex tangent in the frame
        case 'quaternion'
            if isempty(seedDir), seedDir = [1;0;0]; end
            seedDir = seedDir(:);
            F = zeros(3*nSrc, 1);     F((seedVert-1)*3 + (1:3)) = seedDir(1:3);   % 3-vector dipole
    end
    [srcRows, dstRows, nrows, mapMsg] = RowMap(F, ax, blk);
    if ~isempty(mapMsg), error('bst_eigenfilter:Atom', '%s', mapMsg); end
    U = zeros(nrows, 1);  if ~isreal(F), U = complex(U); end
    U(dstRows) = F(srcRows);
    c0 = manifold_ft(Phi, M, U);                                          % seed in the eigenbasis [K x 1]
    kp  = KernelParams;  if ~isfield(kp,'lmax') || isempty(kp.lmax), kp.lmax = max(Lam); end
    g    = bst_eigfilter_kernel(KernelName, kp);
    meta = bst_eigfilter_kernel('info', KernelName);
    dom  = 'static';  if isfield(meta,'domain') && ~isempty(meta.domain), dom = meta.domain; end
    switch lower(dom)
        case 'static'
            W = repmat(manifold_ift(Phi, g(Lam) .* c0), 1, ax.nT);
        case 'ts'
            W = manifold_ift(Phi, g(Lam, ax.tlag) .* c0);
        case 'js'
            % js kernel = temporal frequency response g(lambda,omega). The seed (spatial delta, temporal
            % impulse) has a flat temporal spectrum, so the atom is the inverse temporal FFT of g(lambda,omega)
            % scaled by the spatial seed coeffs -> the exact temporal-Fourier dual of the ts path.
            G  = bst_eigfilter_jtv_evaluate(g, 'js', Lam, ax.omega, ax.NFFT);   % [K x NFFT] = g(lambda,omega)
            Gt = real(ifft(G, ax.NFFT, 2));                                     % -> temporal impulse response
            W  = manifold_ift(Phi, Gt(:, 1:ax.nT) .* c0);
        otherwise
            error('bst_eigenfilter(''Atom''): unknown kernel domain ''%s''.', dom);
    end
end


%% ===== EVALUATE: handle -> gains on Lambda =====
function h = Evaluate(g, Lambda) %#ok<DEFNU>
    if isnumeric(g)
        h = g(:);
    else
        h = bst_eigfilter_evaluate(g, Lambda);
    end
end


%% ===== ROWMAP: per-variant source<->basis row mapping =====
function [srcRows, dstRows, nrows, msg] = RowMap(F, EigenMat, h) %#ok<DEFNU>
    msg = ''; srcRows = []; dstRows = []; nrows = 0;
    switch EigenMat.Variant
        case {'Laplace-Beltrami', 'Connection Laplacian'}
            gv = EigenMat.GlobalVertices{h}(:);
            if max(gv) > size(F, 1)
                msg = sprintf('bst_eigenfilter: hemisphere %d indexes row %d but the map has %d rows.', h, max(gv), size(F,1));
                return;
            end
            if strcmp(EigenMat.Variant, 'Connection Laplacian') && isreal(F)
                msg = ['bst_eigenfilter: Connection Laplacian needs a COMPLEX tangent field in the ' ...
                       'operator frame (the real-3D embedding needs the not-yet-persisted per-vertex frame).'];
                return;
            end
            srcRows = gv;  dstRows = (1:numel(gv))';  nrows = numel(gv);
        case {'Dirac', 'Dirac-Connectome', 'Dirac-Face', 'Hodge-Face'}
            if strcmp(EigenMat.Variant, 'Dirac') || strcmp(EigenMat.Variant, 'Dirac-Connectome')
                idx = EigenMat.GlobalVertices{h}(:);
            else
                idx = EigenMat.GlobalFaces{h}(:);
            end
            if isempty(idx) || (3*max(idx) > size(F, 1))
                msg = sprintf(['bst_eigenfilter: %s hemisphere %d needs a 3-vector map [3*N x nTime] ' ...
                    '(row %d) but the map has %d rows.'], EigenMat.Variant, h, 3*max([idx;0]), size(F,1));
                return;
            end
            n = numel(idx);
            srcRows = reshape([(idx-1)*3+1, (idx-1)*3+2, (idx-1)*3+3].', [], 1);
            dstRows = reshape([(0:n-1)*4+2; (0:n-1)*4+3; (0:n-1)*4+4], [], 1);
            nrows   = 4*n;
        otherwise
            msg = sprintf('bst_eigenfilter: unsupported eigen variant ''%s''.', EigenMat.Variant);
    end
end


%% ===== ANALYSIS: filter a map (project -> gain -> reconstruct) =====
function [Ffilt, Messages, isError] = Analysis(F, EigenMat, OperatorMat, KernelName, KernelParams) %#ok<DEFNU>
    if (nargin < 5) || isempty(KernelParams); KernelParams = struct(); end
    Messages = ''; isError = 0;
    Ffilt = zeros(size(F)); nT = size(F, 2);
    g = Design(KernelName, KernelParams);
    if iscell(g)
        Messages = 'bst_eigenfilter: kernel returned a filterbank (vector param); pass a single-scale kernel.';
        isError = 1; return;
    end
    for h = 1:numel(EigenMat.Phi)
        Phi = EigenMat.Phi{h};
        if isempty(Phi); continue; end
        Lam = EigenMat.Lambda{h}(:);
        B   = OperatorMat.Mass{h};
        [srcRows, dstRows, nrows, msg] = RowMap(F, EigenMat, h);
        if ~isempty(msg); Messages = msg; isError = 1; return; end
        hgain = Evaluate(g, Lam);
        if numel(hgain) ~= numel(Lam)
            Messages = sprintf('bst_eigenfilter: gain has %d entries but the basis has %d modes.', numel(hgain), numel(Lam));
            isError = 1; return;
        end
        U = zeros(nrows, nT);
        U(dstRows, :)  = F(srcRows, :);
        C  = manifold_ft(Phi, B, U);
        Uf = manifold_ift(Phi, hgain .* C);
        Ffilt(srcRows, :) = Uf(dstRows, :);
    end
end


%% ===== SYNTHESIS: re-apply the (symmetric) filter =====
function [Frec, Messages, isError] = Synthesis(C, EigenMat, OperatorMat, KernelName, KernelParams) %#ok<DEFNU>
    if (nargin < 5) || isempty(KernelParams); KernelParams = struct(); end
    [Frec, Messages, isError] = Analysis(C, EigenMat, OperatorMat, KernelName, KernelParams);
end


%% ===== FIBER: operator dimensionality from field_type =====
% Fiber dimensionality of the operator: C components/vertex + kind tag. Source of truth is the nxr
% registry field_type (real/complex/quaternion); falls back to the Phi row/vertex ratio.
function [C, kind] = Fiber(ax) %#ok<DEFNU>
    [C, kind] = i_fiber(ax);
end
function [C, kind] = i_fiber(ax)
    C = [];  ft = '';
    if isfield(ax,'Operator') && isstruct(ax.Operator) && isfield(ax.Operator,'Registry') ...
            && ~isempty(ax.Operator.Registry) && isfield(ax.Operator.Registry,'Primary') ...
            && ~isempty(ax.Operator.Registry.Primary) && isfield(ax.Operator.Registry.Primary,'field_type')
        ft = ax.Operator.Registry.Primary.field_type;
    end
    switch lower(ft)
        case 'real',       C = 1;
        case 'complex',    C = 2;
        case 'quaternion', C = 4;
    end
    if isempty(C)                                   % pre-registry binary -> derive from the Phi layout
        nV = numel(ax.GlobalVertices{1});
        C  = round(size(ax.Phi{1},1) / max(nV,1));
    end
    switch C
        case 1, kind = 'scalar';
        case 2, kind = 'tangent';
        case 4, kind = 'quaternion';
        otherwise, kind = 'scalar';
    end
end
