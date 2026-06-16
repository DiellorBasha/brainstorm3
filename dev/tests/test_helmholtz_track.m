function test_helmholtz_track()
% Unit tests for the pure helpers behind the Helmholtz trajectory overlay.
% Author: Diellor Basha, 2026
    nFail = 0;
    mkc = @(v,ch,p,h,xyz) struct('iVertex',v,'charge',ch,'chirality',sign(ch), ...
                                 'persistence',p,'hemi',h,'pos',xyz);

    % --- bst_persistence_gate: per-hemisphere threshold ---
    % hemi1 max finite=5; hemi2 max finite=0.5. A global gate (max=5) at frac .5 -> thr 2.5
    % would drop hemi2's real core (0.5); per-hemi keeps it.
    mk = [ mkc(1, 1, inf, 1,[0 0 0]), mkc(2, 1, 5,   1,[0 0 0]), ...
           mkc(3,-1, inf, 2,[0 0 0]), mkc(4,-1, 0.5, 2,[0 0 0]), mkc(5,-1,0.05,2,[0 0 0]) ];
    g = bst_persistence_gate(mk, 0.5);
    kv = [g.iVertex];
    nFail = nFail + chk('globals kept (1,3)', all(ismember([1 3], kv)));
    nFail = nFail + chk('hemi1 core 2 kept (>=2.5)', ismember(2, kv));
    nFail = nFail + chk('hemi2 core 4 kept per-hemi (>=0.25)', ismember(4, kv));
    nFail = nFail + chk('hemi2 core 5 dropped (<0.25)', ~ismember(5, kv));

    % --- bst_vortex_link_step: chirality + hemi gated greedy match ---
    prev = [ mkc(10, 1, 1, 1,[0 0 0]), mkc(11,-1,1,1,[0.10 0 0]) ];
    cur  = [ mkc(20, 1, 1, 1,[0.001 0 0]), ...   % matches prev(1)
             mkc(21,-1, 1, 2,[0.10 0 0]), ...     % hemi mismatch -> no match to prev(2)
             mkc(22, 1, 1, 1,[0.20 0 0]) ];       % too far -> birth
    m = bst_vortex_link_step(prev, cur, 0.012);
    nFail = nFail + chk('prev1 -> cur1', m(1)==1);
    nFail = nFail + chk('prev2 unmatched (hemi)', m(2)==0);
    nFail = nFail + chk('chirality mismatch blocks', bst_vortex_link_step(mkc(1,1,1,1,[0 0 0]), mkc(2,-1,1,1,[0 0 0]), 0.012)==0);
    nFail = nFail + chk('within radius same chir+hemi matches', bst_vortex_link_step(mkc(1,1,1,1,[0 0 0]), mkc(2,1,1,1,[0.005 0 0]), 0.012)==1);

    fprintf('\n==== test_helmholtz_track: %d failed ====\n', nFail);
    if nFail > 0, error('test_helmholtz_track FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
