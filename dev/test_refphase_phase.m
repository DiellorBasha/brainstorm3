function test_refphase_phase()
% TEST_REFPHASE_PHASE: regression for process_evt_refphase phase markers.
%
% The phase markers (peak/trough/rising/falling) must, WITHIN each detected period,
% strictly alternate off a monotonic phase: the extremum family {peak,trough} and the
% crossing family {rising,falling} interleave (max,min,max,min), consecutive extrema
% alternate peak<->trough, and consecutive crossings alternate rising<->falling. The
% old per-event sign/slope resolution dropped ~4% of crossings (rising often missing);
% the phase-alternation rewrite makes skips impossible.
%
% USAGE:  test_refphase_phase   % Brainstorm running, TutorialAuditory loaded
%
% Authors: Diellor Basha, 2026

    PF = {'FAIL','PASS'};  pass = true;
    band = [8 13];

    % ---------- load an alpha-rich stretch of the auditory recording ----------
    rawLink = 'Subject01/@rawS01_AEF_20131218_01_notch/data_0raw_S01_AEF_20131218_01_notch.mat';
    [sStudy,~] = bst_get('DataFile', rawLink);
    if isempty(sStudy)
        fprintf('SKIPPED (no raw link %s)\n', rawLink);
        fprintf('\n==== SUITE: %s ====\n', PF{pass+1});  return;
    end
    DataMat = in_bst_data(rawLink, 'F');  sFile = DataMat.F;
    sfreq = sFile.prop.sfreq;  Tprop = sFile.prop.times;
    ChannelMat = in_bst_channel(bst_get('ChannelFileForStudy', rawLink));
    iMEG = channel_find(ChannelMat.Channel, 'MEG');
    IO = db_template('ImportOptions');
    IO.ImportMode='Time'; IO.UseCtfComp=1; IO.UseSsp=1; IO.EventsMode='ignore';
    IO.DisplayMessages=0; IO.RemoveBaseline='no';
    s0 = round((0-Tprop(1))*sfreq)+1;  s1 = round((min(Tprop(2),180)-Tprop(1))*sfreq)+1;
    [F, TimeVector] = in_fread(sFile, ChannelMat, 1, [s0 s1], iMEG, IO);

    OPTIONS = process_evt_refphase('Compute');  OPTIONS.freqRange = band;
    [evt, markers] = process_evt_refphase('Compute', F, TimeVector, OPTIONS);

    % ---------- T1: within-period alternation (no skips) ----------
    nViol = 0;  nMarkTot = 0;  nPer = size(evt,2);
    extFam = containers.Map({'peak','trough','rising','falling'}, {1,1,0,0});  % 1=extremum,0=crossing
    for w = 1:nPer
        [tw, lw] = i_period_markers(markers, evt(:,w));
        nMarkTot = nMarkTot + numel(lw);
        if numel(lw) < 2, continue; end
        fam = cellfun(@(x) extFam(x), lw);
        % (a) families must interleave: no two extrema or two crossings in a row
        nViol = nViol + sum(fam(2:end)==fam(1:end-1));
        % (b) consecutive extrema alternate peak<->trough; crossings alternate rising<->falling
        ext = lw(fam==1);  cro = lw(fam==0);
        nViol = nViol + sum(strcmp(ext(2:end), ext(1:end-1)));
        nViol = nViol + sum(strcmp(cro(2:end), cro(1:end-1)));
    end
    ok1 = (nPer > 0) && (nViol == 0);
    fprintf('T1 within-period alternation: periods=%d markers=%d violations=%d => %s\n', ...
        nPer, nMarkTot, nViol, PF{ok1+1});
    pass = pass && ok1;

    % ---------- T2: balanced counts (crossings ~= extrema; no systematic drop) ----------
    nP=numel(markers.peak); nT=numel(markers.trough); nR=numel(markers.rising); nF=numel(markers.falling);
    ext = nP+nT;  cro = nR+nF;
    ok2 = (ext>0) && (abs(ext-cro) <= nPer) && (abs(nP-nT) <= nPer) && (abs(nR-nF) <= nPer);
    fprintf('T2 balanced counts: peak/trough/rising/falling=%d/%d/%d/%d (ext=%d cross=%d, <=%d apart) => %s\n', ...
        nP,nT,nR,nF, ext, cro, nPer, PF{ok2+1});
    pass = pass && ok2;

    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end

%% ===== collect one period's phase markers, time-sorted, as {label} list =====
function [tw, lw] = i_period_markers(markers, period)
    on = period(1);  off = period(2);
    names = {'peak','trough','rising','falling'};
    tw = [];  lw = {};
    for k = 1:numel(names)
        tt = markers.(names{k});
        tt = tt(tt >= on & tt <= off);
        tw = [tw, tt];                                  %#ok<AGROW>
        lw = [lw, repmat(names(k), 1, numel(tt))];      %#ok<AGROW>
    end
    [tw, ord] = sort(tw);  lw = lw(ord);
end
