function out = tess_cholesky(varargin)
% TESS_CHOLESKY: Lazy, cached Cholesky factor of an assembled operator matrix.
%
% The factor of a pinned SPD operator (A = OperatorNode.Operator{hh}) is the factorization
% of ONE specific assembled matrix, so it is persisted ON the operator_ node (next to A),
% not on the geometry-only manifold_ node. Three modes:
%
%   dF   = tess_cholesky(OperatorNode, hh, pin)            % PURE getter, I/O-free
%          -> dF = struct('L',L,'p',p,'free',free,'n',n), [L,~,p]=chol(A(free,free),'lower','vector')
%          Resolution: node.Cholesky{hh} -> in-session memory cache -> compute (+cache in memory)
%   x    = tess_cholesky('solve', dF, rhs)                 % A x = rhs on the free block; pinned rows 0
%   Node = tess_cholesky('attach', OperatorNode, File, pin)% compute both hemis, save onto the node (I/O)
%
% The opaque `decomposition` object is NOT serialized; L (sparse lower) + p (permutation) are
% plain data that reload cleanly and reconstruct the two triangular solves.
%
% SEE ALSO: bst_poisson, bst_get_operator_node, tess_operators

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

    persistent MEM
    if isempty(MEM), MEM = containers.Map('KeyType','char','ValueType','any'); end

    if ischar(varargin{1})
        switch lower(varargin{1})
            case 'solve'
                out = i_solve(varargin{2}, varargin{3});
            case 'attach'
                out = i_attach(varargin{2}, varargin{3}, i_pin(varargin{4:end}));
            otherwise
                error('tess_cholesky:badMode', 'Unknown mode: %s', varargin{1});
        end
        return;
    end

    % pure getter: tess_cholesky(OperatorNode, hh, pin)
    Node = varargin{1};  hh = varargin{2};  pin = i_pin(varargin{3:end});
    % 1) factor already on the node?
    if isfield(Node,'Cholesky') && iscell(Node.Cholesky) && numel(Node.Cholesky) >= hh && ~isempty(Node.Cholesky{hh})
        out = Node.Cholesky{hh};  return;
    end
    % 2) in-session memory cache?
    key = i_key(Node, hh, pin);
    if isKey(MEM, key), out = MEM(key);  return; end
    % 3) compute + cache in memory (no disk write on the pure path)
    out = i_factor(Node.Operator{hh}, pin);
    MEM(key) = out;
end

function pin = i_pin(varargin)
    if isempty(varargin) || isempty(varargin{1}), pin = 1; else, pin = varargin{1}(:)'; end
end

function key = i_key(Node, hh, pin)
    key = sprintf('%s|%s|%d|%s', Node.Variant, Node.ParentSurface, hh, mat2str(pin));
end

function dF = i_factor(A, pin)
    A = (A + A')/2;
    n = size(A,1);
    free = setdiff((1:n)', pin(:));
    [L, flag, p] = chol(A(free,free), 'lower', 'vector');
    if flag ~= 0
        error('tess_cholesky:notSPD', 'Pinned operator block is not SPD (chol flag=%d).', flag);
    end
    dF = struct('L', L, 'p', p(:), 'free', free, 'n', n);
end

function x = i_solve(dF, rhs)
    % A(free,free)(p,p) = L L'  =>  solve on the free block, pinned rows stay 0
    k = size(rhs, 2);
    x = zeros(dF.n, k);
    bf = rhs(dF.free, :);
    bp = bf(dF.p, :);
    y  = dF.L \ bp;
    z  = dF.L' \ y;
    xf = zeros(numel(dF.free), k);
    xf(dF.p, :) = z;
    x(dF.free, :) = xf;
end

function Node = i_attach(Node, OperatorFile, pin)
    nH = numel(Node.Operator);
    if ~iscell(Node.Cholesky), Node.Cholesky = cell(1, nH); end
    for hh = 1:nH
        if isempty(Node.Operator{hh}), continue; end
        Node.Cholesky{hh} = i_factor(Node.Operator{hh}, pin);
    end
    bst_save(file_fullpath(OperatorFile), Node, 'v6');
end
