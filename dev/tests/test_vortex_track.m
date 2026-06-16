function test_vortex_track()
% Pure unit tests for bst_vortex_track on hand-built core sequences.
% Author: Diellor Basha, 2026
    nFail = 0;
    mk = @(v,ch,p,xyz) struct('iVertex',v,'chirality',ch,'persistence',p,'pos',xyz);

    % A) one core drifting in small steps over 5 frames -> one length-5 track
    A = cell(1,5);
    for t=1:5, A{t} = mk(t, 1, 1, [0.001*t 0 0]); end
    T = bst_vortex_track(A, 'MaxJump', 0.010);
    nFail = nFail + chk('A: single track',        numel(T)==1);
    nFail = nFail + chk('A: length 5',            numel(T(1).frames)==5);
    nFail = nFail + chk('A: birth1 death5',       T(1).birthFrame==1 && T(1).deathFrame==5);

    % B) two separated stationary cores -> two length-5 tracks
    B = cell(1,5);
    for t=1:5, B{t} = [mk(1,1,1,[0 0 0]), mk(2,1,1,[0.05 0 0])]; end
    Tb = bst_vortex_track(B, 'MaxJump', 0.010);
    nFail = nFail + chk('B: two tracks',          numel(Tb)==2);
    nFail = nFail + chk('B: both length 5',       all(arrayfun(@(x)numel(x.frames)==5, Tb)));

    % C) birth + death: core A (frames 1-2), core B far away (frames 3-5)
    C = {mk(1,1,1,[0 0 0]), mk(1,1,1,[0 0 0]), mk(2,1,1,[0.05 0 0]), mk(2,1,1,[0.05 0 0]), mk(2,1,1,[0.05 0 0])};
    Tc = bst_vortex_track(C, 'MaxJump', 0.010);
    nFail = nFail + chk('C: two tracks',          numel(Tc)==2);
    nFail = nFail + chk('C: a death at 2',        any([Tc.deathFrame]==2));
    nFail = nFail + chk('C: a birth at 3',        any([Tc.birthFrame]==3));

    % D) chirality mismatch must not link
    D = {mk(1,1,1,[0 0 0]), mk(2,-1,1,[0.001 0 0])};
    Td = bst_vortex_track(D, 'MaxJump', 0.010);
    nFail = nFail + chk('D: chirality blocks link', numel(Td)==2);

    % E) jump beyond MaxJump must not link
    E = {mk(1,1,1,[0 0 0]), mk(2,1,1,[0.05 0 0])};
    Te = bst_vortex_track(E, 'MaxJump', 0.010);
    nFail = nFail + chk('E: long jump splits',    numel(Te)==2);

    fprintf('\n==== test_vortex_track: %d failed ====\n', nFail);
    if nFail > 0, error('test_vortex_track FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
