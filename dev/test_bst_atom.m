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

    % ---------- T2: Set then Get round-trips on every axis; writes the right fields ----------
    G2 = bst_dynamics('NewGroup', 'g2');
    G2.times = [0.5];   % simple, 1 occurrence
    % time: set a window on occ 1 -> promotes to extended
    G2 = bst_atom('Set', G2, 'time', 1, i_loc('time', 0.40, 0.10));
    rtT = (size(G2.times,1)==2) && strcmp(G2.type,'extended');
    gt  = bst_atom('Get', G2, 'time', 1);
    rtT = rtT && abs(gt.center-0.40)<1e-9 && abs(gt.extent-0.10)<1e-9;
    % freq: set alpha band (group-level)
    G2 = bst_atom('Set', G2, 'freq', [], i_loc_lbl('freq', 10, 2, 'alpha'));
    gf  = bst_atom('Get', G2, 'freq');
    rtF = isequal(G2.band,[8 12]) && abs(gf.center-10)<1e-9 && abs(gf.extent-2)<1e-9 && strcmp(gf.label,'alpha');
    % source: set seed + radius on occ 1 (with pos)
    lcS = i_loc('source', 250, 0.006);  lcS.pos = [7 8 9];
    G2 = bst_atom('Set', G2, 'source', 1, lcS);
    gs  = bst_atom('Get', G2, 'source', 1);
    rtS = (G2.vertices(1)==250) && abs(G2.radius(1)-0.006)<1e-12 && isequal(G2.pos(1,:),[7 8 9]) ...
       && strcmp(gs.state,'window') && abs(gs.extent-0.006)<1e-12;
    % scale: set eigen-band (group-level)
    G2 = bst_atom('Set', G2, 'scale', [], i_loc_lbl('scale', 50, 10, 'gyrus'));
    gk  = bst_atom('Get', G2, 'scale');
    rtK = isequal(G2.scale,[40 60]) && abs(gk.center-50)<1e-9 && abs(gk.extent-10)<1e-9 && strcmp(gk.label,'gyrus');

    ok2 = rtT && rtF && rtS && rtK;
    fprintf('T2 Set round-trip: time=%d freq=%d source=%d scale=%d => %s\n', rtT, rtF, rtS, rtK, PF{ok2+1});
    pass = pass && ok2;

    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end

function loc = i_loc(axis, c, w)
    loc = bst_atom('NewLoc', axis);  loc.center = c;  loc.extent = w;
end
function loc = i_loc_lbl(axis, c, w, lbl)
    loc = i_loc(axis, c, w);  loc.label = lbl;
end
