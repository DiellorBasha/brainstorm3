function c = bst_git_commit(bstDir)
% BST_GIT_COMMIT  Commit of the Brainstorm checkout at bstDir, for provenance: $NSP_BST_COMMIT if set,
% else read from bstDir/.git without a git binary (the container may not have one). Handles a
% detached HEAD (a pinned checkout: HEAD holds the sha), a branch (loose ref or packed-refs) and a
% worktree (.git is a file). 'unknown' when none resolves.
%
% USAGE:  c = bst_git_commit(bstDir)
%
% Authors: Diellor Basha, 2026 (nsp brainstorm-pet pathway)
    c = getenv('NSP_BST_COMMIT');
    if ~isempty(c), return; end
    c = 'unknown';
    g = fullfile(bstDir, '.git');
    if exist(g, 'file') == 2                       % worktree: ".git" is a file "gitdir: <path>"
        t = regexp(fileread(g), 'gitdir:\s*(\S+)', 'tokens', 'once');
        if isempty(t), return; end
        g = t{1};
    end
    h = fullfile(g, 'HEAD');
    if exist(h, 'file') ~= 2, return; end
    head = strtrim(fileread(h));
    ref = regexp(head, '^ref:\s*(\S+)', 'tokens', 'once');
    if isempty(ref), c = head; return; end
    common = g;                                    % a worktree's refs live in the common dir
    if exist(fullfile(g, 'commondir'), 'file') == 2
        common = fullfile(g, strtrim(fileread(fullfile(g, 'commondir'))));
    end
    r = fullfile(common, ref{1});
    if exist(r, 'file') == 2
        c = strtrim(fileread(r));
    elseif exist(fullfile(common, 'packed-refs'), 'file') == 2
        t = regexp(fileread(fullfile(common, 'packed-refs')), ['(\w{40}) ' regexptranslate('escape', ref{1})], 'tokens', 'once');
        if ~isempty(t), c = t{1}; end
    end
end
