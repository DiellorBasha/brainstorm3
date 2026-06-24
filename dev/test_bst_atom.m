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

    % ---------- T3: read a real refphase-detected band group through bst_atom ----------
    [linkFile, relData] = i_find_kernel_atom();
    if isempty(linkFile)
        fprintf('T3: SKIPPED (no unconstrained kernel link)\n');
        fprintf('\n==== SUITE: %s ====\n', PF{pass+1});  return;
    end
    % build a band-window group exactly as OnDetect does: refphase on alpha
    DataMat = in_bst_data(relData, 'F', 'Time');
    ChannelMat = in_bst_channel(bst_get('ChannelFileForStudy', relData));
    iMEG = channel_find(ChannelMat.Channel, 'MEG');
    if isstruct(DataMat.F)
        IO = db_template('ImportOptions');  IO.ImportMode='Time'; IO.UseCtfComp=1; IO.UseSsp=1;
        IO.EventsMode='ignore'; IO.DisplayMessages=0; IO.RemoveBaseline='no';
        [F, TimeVector] = in_fread(DataMat.F, ChannelMat, 1, [], iMEG, IO);
    else
        F = DataMat.F(iMEG,:);  TimeVector = DataMat.Time;
    end
    OPTIONS = process_evt_refphase('Compute');  OPTIONS.freqRange = [8 13];
    [evt, ~] = process_evt_refphase('Compute', F, TimeVector, OPTIONS);
    W = bst_dynamics('NewGroup', 'alpha (8-13 Hz)');
    W.times = evt;  W.band = [8 13];  W.bandName = 'alpha';
    % freq axis: center 10.5, extent 2.5, label alpha
    lf = bst_atom('Get', W, 'freq');
    okF = abs(lf.center-10.5)<1e-9 && abs(lf.extent-2.5)<1e-9 && strcmp(lf.label,'alpha') && strcmp(lf.state,'window');
    % time axis: occ 1 is an extended window with center = mean(onset,offset), extent>0
    lt = bst_atom('Get', W, 'time', 1);
    okT = strcmp(lt.state,'window') && (lt.extent>0) && abs(lt.center-mean(evt(:,1)))<1e-9;
    % source axis: a band-window group has no source -> unlocalized
    ls = bst_atom('Get', W, 'source', 1);
    okS = strcmp(ls.state,'unlocalized');
    ok3 = (size(evt,2)>0) && okF && okT && okS;
    fprintf('T3 real detect: nWin=%d freq=%d time=%d srcUnloc=%d => %s\n', size(evt,2), okF, okT, okS, PF{ok3+1});
    pass = pass && ok3;

    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end

function loc = i_loc(axis, c, w)
    loc = bst_atom('NewLoc', axis);  loc.center = c;  loc.extent = w;
end
function loc = i_loc_lbl(axis, c, w, lbl)
    loc = i_loc(axis, c, w);  loc.label = lbl;
end

function [linkFile, relData] = i_find_kernel_atom()
    linkFile = '';
    relData = 'Subject01/S01_AEF_20131218_01_notch/data_block001_02.mat';
    [sStudy, ~] = bst_get('DataFile', relData);
    if isempty(sStudy), return; end
    comments = {sStudy.Result.Comment};  fnames = {sStudy.Result.FileName};
    isMN = ~cellfun(@isempty, regexp(comments, 'MN: MEG\(Unconstr\)', 'once')) & ...
           ~cellfun(@isempty, regexp(fnames,   'KERNEL', 'once'));
    for j = find(isMN)
        try
            r = in_bst_results(fnames{j}, 0, 'nComponents','ImagingKernel');
            if (r.nComponents==3) && ~isempty(r.ImagingKernel)
                linkFile = ['link|' fnames{j} '|' relData];  return;
            end
        catch
        end
    end
end
