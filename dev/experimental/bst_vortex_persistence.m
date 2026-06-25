function Cores = bst_vortex_persistence(field, faces, varargin)
% BST_VORTEX_PERSISTENCE  Vortex cores as persistent extrema of a scalar field.
%
% Cores are local extrema of FIELD on a triangular mesh, ranked by topological
% persistence (peak-minus-saddle contrast) from one union-find pass over the
% superlevel-set filtration (the elder rule). Maxima of +FIELD are positive-
% chirality cores; maxima of -FIELD (minima of FIELD) are negative-chirality.
% Every connected component contributes one never-merging (Inf-persistence)
% global core, so the result is correct on disconnected meshes (e.g. 2 hemispheres).
%
% USAGE:
%   Cores = bst_vortex_persistence(field, faces)
%   Cores = bst_vortex_persistence(field, faces, 'MinPersistence', p, 'Neighbors', nb)
%
% INPUT:
%   field  : [nV x 1] scalar field on mesh vertices (e.g. the stream function psi)
%   faces  : [nF x 3] 1-based triangles (used only to build the 1-ring if Neighbors
%            is not supplied; pass [] when supplying Neighbors)
%   'MinPersistence' : keep cores with persistence >= this (default 0; Inf always kept)
%   'Neighbors'      : {nV x 1} precomputed 1-ring adjacency (reuse across frames)
%
% OUTPUT: Cores (columnar struct, sorted by persistence descending):
%   .vertex .value .persistence .chirality(+1/-1) .isGlobal .birth .death
%
% Author: Diellor Basha, 2026

    field = field(:);
    nV = numel(field);

    MinPersistence = 0;  neighbors = {};
    for k = 1:2:numel(varargin)
        switch lower(varargin{k})
            case 'minpersistence', MinPersistence = varargin{k+1};
            case 'neighbors',      neighbors      = varargin{k+1};
            otherwise, error('bst_vortex_persistence: unknown option %s', varargin{k});
        end
    end
    if isempty(neighbors)
        neighbors = i_one_ring(faces, nV);
    elseif numel(neighbors) ~= nV
        error('bst_vortex_persistence: Neighbors must have one cell per vertex (%d vs %d).', ...
              numel(neighbors), nV);
    end

    pos = i_superlevel(field,  neighbors, nV);     % +field
    neg = i_superlevel(-field, neighbors, nV);     % -field

    vertex      = [pos.peak;        neg.peak];
    persistence = [pos.persistence; neg.persistence];
    birth       = [pos.birth;       neg.birth];
    death       = [pos.death;       neg.death];
    chirality   = [ ones(numel(pos.peak),1); -ones(numel(neg.peak),1)];
    isGlobal    = isinf(persistence);
    value       = field(vertex);

    keep = persistence >= MinPersistence;          % Inf always passes
    idx  = find(keep);
    [persistence, ord] = sort(persistence(idx), 'descend');
    idx  = idx(ord);

    Cores = struct('vertex',vertex(idx), 'value',value(idx), 'persistence',persistence, ...
                   'chirality',chirality(idx), 'isGlobal',isGlobal(idx), ...
                   'birth',birth(idx), 'death',death(idx));
end

% ---- 0-D persistence of the superlevel-set filtration of FIELD ----
function F = i_superlevel(field, neighbors, nV)
    [~, order] = sort(field, 'descend');
    parent     = zeros(nV,1);     % union-find parent; 0 = not yet added
    peakVertex = zeros(nV,1);     % peak vertex carried by each component root
    isAdded    = false(nV,1);
    featPeak   = zeros(nV,1);
    featDeath  = zeros(nV,1);
    nFeat      = 0;

    for k = 1:nV
        v = order(k);  isAdded(v) = true;
        nb = neighbors{v};  nb = nb(isAdded(nb));
        if isempty(nb)
            parent(v) = v;  peakVertex(v) = v;            % new component born (a maximum)
        else
            roots = unique(arrayfun(@findRoot, nb(:)));
            [~, hi]  = max(field(peakVertex(roots)));      % elder = highest peak
            survivor = roots(hi);
            survivorPeak = peakVertex(survivor);
            parent(v) = survivor;
            for r = roots(:)'
                if r ~= survivor
                    nFeat = nFeat + 1;
                    featPeak(nFeat)  = peakVertex(r);       % younger dies here
                    featDeath(nFeat) = field(v);            % at saddle value
                    parent(r) = survivor;
                end
            end
            peakVertex(survivor) = survivorPeak;
        end
    end
    % every surviving root never merges -> a global (Inf-persistence) feature
    for v = 1:nV
        if parent(v) == v
            nFeat = nFeat + 1;
            featPeak(nFeat)  = peakVertex(v);
            featDeath(nFeat) = -Inf;
        end
    end

    featPeak  = featPeak(1:nFeat);
    featDeath = featDeath(1:nFeat);
    F.peak        = featPeak;
    F.birth       = field(featPeak);
    F.death       = featDeath;
    F.persistence = F.birth - featDeath;                    % Inf for globals

    function r = findRoot(x)
        r = x;
        while parent(r) ~= r, r = parent(r); end
        while parent(x) ~= r, nx = parent(x); parent(x) = r; x = nx; end
    end
end

% ---- 1-ring vertex adjacency from triangle faces ----
function neighbors = i_one_ring(faces, nV)
    if isempty(faces)
        error('bst_vortex_persistence: need faces or Neighbors.');
    end
    e = [faces(:,[1 2]); faces(:,[2 3]); faces(:,[3 1])];
    A = sparse([e(:,1);e(:,2)], [e(:,2);e(:,1)], true, nV, nV);
    neighbors = cell(nV,1);
    for v = 1:nV, neighbors{v} = find(A(:,v)); end
end
