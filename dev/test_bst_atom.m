function test_bst_atom()
% TEST_BST_ATOM: unit + round-trip tests for the (center,extent) localization accessor.
%
% USAGE:  test_bst_atom   % Brainstorm running
%
% Authors: Diellor Basha, 2026

    PF = {'FAIL','PASS'};  pass = true;

    % ---------- T1: Axes metadata + Get on a synthetic group (4 axes, 3 states) ----------
    A = bst_atom('Axes');
    axesOK = (numel(A)==4) && strcmp(A(1).name,'time') && strcmp(A(2).name,'freq') ...
          && strcmp(A(3).name,'source') && strcmp(A(4).name,'scale') ...
          && (A(1).perOcc && ~A(2).perOcc && A(3).perOcc && ~A(4).perOcc);

    G = bst_dynamics('NewGroup', 'g');
    G.times = [0.10 0.20; 0.30 0.40];      % extended (2xN): occ1 window [.10 .30]->c=.20 w=.10
    G.band  = [8 12];                       % freq window: c=10 w=2 (alpha)
    G.bandName = 'alpha';
    G.scale = [];                           % scale unlocalized
    G.vertices = [100 200];  G.pos = [1 2 3; 4 5 6];  G.radius = [0.005 0];  % occ1 window r=5mm, occ2 point

    lt = bst_atom('Get', G, 'time', 1);
    okT = strcmp(lt.state,'window') && abs(lt.center-0.20)<1e-9 && abs(lt.extent-0.10)<1e-9;
    lf = bst_atom('Get', G, 'freq');        % group-level -> occ ignored
    okF = strcmp(lf.state,'window') && (lf.center==10) && (lf.extent==2) && strcmp(lf.label,'alpha');
    ls1 = bst_atom('Get', G, 'source', 1);
    okS1 = strcmp(ls1.state,'window') && (ls1.center==100) && abs(ls1.extent-0.005)<1e-12 && isequal(ls1.pos,[1 2 3]);
    ls2 = bst_atom('Get', G, 'source', 2);
    okS2 = strcmp(ls2.state,'point') && (ls2.center==200) && (ls2.extent==0);
    lk = bst_atom('Get', G, 'scale');
    okK = strcmp(lk.state,'unlocalized') && ~isfinite(lk.center);
    % weighting default + radius field present
    okW = strcmp(lt.weighting,'hard') && isfield(db_template('atomgroup'),'radius');

    ok1 = axesOK && okT && okF && okS1 && okS2 && okK && okW;
    fprintf('T1 Get: axes=%d time=%d freq=%d src1=%d src2=%d scale=%d weight/field=%d => %s\n', ...
        axesOK, okT, okF, okS1, okS2, okK, okW, PF{ok1+1});
    pass = pass && ok1;

    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end
