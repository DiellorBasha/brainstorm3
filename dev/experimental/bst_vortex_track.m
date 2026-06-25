function Tracks = bst_vortex_track(coresPerFrame, varargin)
% BST_VORTEX_TRACK  Link per-frame cores into trajectories (greedy, chirality-consistent).
%
% USAGE: Tracks = bst_vortex_track(coresPerFrame, 'MinPersistence', p, 'MaxJump', d)
% INPUT:
%   coresPerFrame : {1 x nT} cell; coresPerFrame{t} = struct array of cores with
%                   fields .iVertex .chirality(+/-1) .persistence .pos (1x3, meters)
%   'MinPersistence' (default 0)  : drop cores below this persistence (Inf always kept)
%   'MaxJump' (default 0.010 m)   : max core displacement between consecutive frames
% OUTPUT: Tracks struct array (one per trajectory):
%   .frames [1xL] .iVertex [1xL] .pos [Lx3] .persistence [1xL]
%   .chirality (+/-1) .birthFrame .deathFrame
% Author: Diellor Basha, 2026

    MinPersistence = 0;  MaxJump = 0.010;
    for k = 1:2:numel(varargin)
        switch lower(varargin{k})
            case 'minpersistence', MinPersistence = varargin{k+1};
            case 'maxjump',        MaxJump        = varargin{k+1};
            otherwise, error('bst_vortex_track: unknown option %s', varargin{k});
        end
    end

    nT = numel(coresPerFrame);
    Tracks = i_empty_tracks();
    openIdx = [];                         % indices into Tracks still open
    for t = 1:nT
        cur = coresPerFrame{t};
        if ~isempty(cur)
            cur = cur([cur.persistence] >= MinPersistence);   % Inf passes
        end
        nc = numel(cur);
        matched = false(1, nc);

        if ~isempty(openIdx) && nc > 0
            H = numel(openIdx);
            pairs = zeros(0,3);            % [headSlot, curIdx, dist]
            for hs = 1:H
                tr = Tracks(openIdx(hs));
                hp = tr.pos(end,:);  hc = tr.chirality;
                for c = 1:nc
                    if cur(c).chirality ~= hc, continue; end
                    d = norm(cur(c).pos - hp);
                    if d <= MaxJump, pairs(end+1,:) = [hs, c, d]; end %#ok<AGROW>
                end
            end
            usedH = false(1,H);
            if ~isempty(pairs)
                pairs = sortrows(pairs, 3);
                for p = 1:size(pairs,1)
                    hs = pairs(p,1);  c = pairs(p,2);
                    if usedH(hs) || matched(c), continue; end
                    usedH(hs) = true;  matched(c) = true;
                    ti = openIdx(hs);
                    Tracks(ti) = i_extend(Tracks(ti), t, cur(c));
                end
            end
            newOpen = openIdx(usedH);
            for hs = find(~usedH), Tracks(openIdx(hs)).deathFrame = t-1; end
            openIdx = newOpen;
        elseif ~isempty(openIdx)          % nc == 0: close everything open
            for hs = 1:numel(openIdx), Tracks(openIdx(hs)).deathFrame = t-1; end
            openIdx = [];
        end

        for c = 1:nc                       % births
            if ~matched(c)
                Tracks(end+1) = i_new(t, cur(c)); %#ok<AGROW>
                openIdx(end+1) = numel(Tracks); %#ok<AGROW>
            end
        end
    end
    for hs = 1:numel(openIdx), Tracks(openIdx(hs)).deathFrame = nT; end
end

function s = i_empty_tracks()
    s = struct('frames',{},'iVertex',{},'pos',{},'persistence',{}, ...
               'chirality',{},'birthFrame',{},'deathFrame',{});
end
function tr = i_new(t, c)
    tr = struct('frames',t, 'iVertex',c.iVertex, 'pos',c.pos, 'persistence',c.persistence, ...
                'chirality',c.chirality, 'birthFrame',t, 'deathFrame',t);
end
function tr = i_extend(tr, t, c)
    tr.frames(end+1)      = t;
    tr.iVertex(end+1)     = c.iVertex;
    tr.pos(end+1,:)       = c.pos;
    tr.persistence(end+1) = c.persistence;
    tr.deathFrame         = t;
end
