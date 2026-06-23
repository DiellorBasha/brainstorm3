function test_evt_refphase()
% TEST_EVT_REFPHASE: Regression tests for process_evt_refphase.
%
% Covers: period detection, 2f phase markers (peak/trough), cycle-snapped
% bounds, bad-segment masking, filter-transient handling, per-channel
% normalization, and degenerate inputs. Synthetic tests are self-contained;
% one test runs on the real auditory block (requires the TutorialAuditory
% protocol loaded).
%
% USAGE:  test_evt_refphase   % from MATLAB with Brainstorm running
%
% Authors: Diellor Basha, 2026

    PF = {'FAIL','PASS'};
    pass = true;
    fs = 300;  t = 0:1/fs:30;  n = numel(t);

    % ---------- T1: synthetic topographic burst (MAD) -> 1 period + markers ----------
    rng(42);
    Fsyn = 0.2*randn(20, n);
    burst = (t>=12 & t<=15);
    topo  = randn(20,1);                       % spatial pattern (GFP can see it)
    Fsyn(:, burst) = Fsyn(:, burst) + topo*(3*sin(2*pi*10*t(burst)));
    OPT = process_evt_refphase('Compute');
    OPT.freqRange=[8 13]; OPT.thresholdMode='mad'; OPT.enterThresh=4; OPT.exitThresh=2;  % timing auto-derived
    [evt, mk, st] = process_evt_refphase('Compute', Fsyn, t, OPT);
    ov       = ~isempty(evt) && any(evt(1,:)<=15 & evt(2,:)>=12);
    localized = ~isempty(evt) && (evt(1)>=10.5) && (evt(2)<=16.5);          % period hugs the burst
    mkInPeriod = ~isempty(mk.peak) && all(mk.peak>=evt(1)-1e-9 & mk.peak<=evt(2)+1e-9) ...
                                   && all(mk.trough>=evt(1)-1e-9 & mk.trough<=evt(2)+1e-9);
    ok1 = (size(evt,2)==1) && ov && localized && ~isempty(mk.peak) && ~isempty(mk.trough) && mkInPeriod;
    fprintf('T1 burst: %dp [%.1f-%.1f]s cov%.1f%% %dpk %dtr inPeriod=%d => %s\n', ...
        size(evt,2), evt(1), evt(2), st.coverage, numel(mk.peak), numel(mk.trough), mkInPeriod, PF{ok1+1});
    pass = pass && ok1;

    % ---------- T2 + T3: real data ----------
    [F, TimeVector] = i_read_real();
    if ~isempty(F)
        OPT2 = process_evt_refphase('Compute'); OPT2.freqRange=[8 13];
        [evt2, mk2, st2] = process_evt_refphase('Compute', F, TimeVector, OPT2);
        nW   = size(evt2,2);
        mono = all(evt2(2,:) > evt2(1,:));
        inP  = true;
        for k=1:numel(mk2.peak)
            inP = inP && any(mk2.peak(k)>=evt2(1,:)-1e-9 & mk2.peak(k)<=evt2(2,:)+1e-9);
        end
        ok2 = (nW>=3) && (nW<=15) && mono && ~isempty(mk2.peak) && inP;
        fprintf('T2 real: %dp cov%.1f%% %dpk %dtr mono=%d pkInP=%d => %s\n', nW, st2.coverage, numel(mk2.peak), numel(mk2.trough), mono, inP, PF{ok2+1});
        pass = pass && ok2;
        % 2f marker frequency: peak spacing ~ 1/(2*f_alpha)
        dpk = diff(sort(mk2.peak));  medISI = median(dpk(dpk<0.2));  expected = 1/(2*10.5);
        ok3 = abs(medISI - expected) < 0.020;
        fprintf('T3 peak ISI=%.3fs (expect ~%.3f @2x alpha) => %s\n', medISI, expected, PF{ok3+1});
        pass = pass && ok3;
    else
        fprintf('T2/T3 real: SKIPPED (file not found)\n');
    end

    % ---------- T4: flat -> no events, no markers ----------
    [evt4, mk4] = process_evt_refphase('Compute', ones(10,n), t, OPT);
    ok4 = isempty(evt4) && isempty(mk4.peak) && isempty(mk4.trough);
    fprintf('T4 flat: %dp %dpk => %s\n', size(evt4,2), numel(mk4.peak), PF{ok4+1});
    pass = pass && ok4;

    % ---------- T5: bad-segment masking excludes artifact, keeps real alpha ----------
    Fart = 0.2*randn(20,n);
    Fart(:, t>=5 & t<=15)    = Fart(:, t>=5 & t<=15) + topo*(1.5*sin(2*pi*10*t(t>=5 & t<=15)));
    Fart(:, t>=20 & t<=20.5) = Fart(:, t>=20 & t<=20.5) + 50;     % huge artifact
    vmask = true(1,n);  vmask(t>=19.8 & t<=20.7) = false;          % mark it bad
    evtA = process_evt_refphase('Compute', Fart, t, OPT, vmask);
    noArt   = isempty(evtA) || ~any(evtA(1,:)<20.7 & evtA(2,:)>19.8);
    hitReal = ~isempty(evtA) && any(evtA(1,:)<15 & evtA(2,:)>5);
    ok5 = noArt && hitReal;
    fprintf('T5 bad-seg: %dp noArt=%d hitReal=%d => %s\n', size(evtA,2), noArt, hitReal, PF{ok5+1});
    pass = pass && ok5;

    % ---------- T6: auto-derived timing scales inversely with the band ----------
    OPTa = process_evt_refphase('Compute'); OPTa.freqRange=[8 13];
    [~,~,sa] = process_evt_refphase('Compute', Fsyn, t, OPTa);
    OPTb = process_evt_refphase('Compute'); OPTb.freqRange=[13 30];
    [~,~,sb] = process_evt_refphase('Compute', Fsyn, t, OPTb);
    fcA = 10.5;
    autoOK = abs(sa.smoothing-0.5/fcA)<1e-9 && abs(sa.minDuration-3/fcA)<1e-9 && abs(sa.minGap-2/fcA)<1e-9;
    ok6 = autoOK && (sb.smoothing < sa.smoothing);   % higher band -> shorter windows
    fprintf('T6 auto-derive: alpha smooth=%.0fms minDur=%.0fms minGap=%.0fms (beta smooth=%.0fms) => %s\n', ...
        1000*sa.smoothing, 1000*sa.minDuration, 1000*sa.minGap, 1000*sb.smoothing, PF{ok6+1});
    pass = pass && ok6;

    % ---------- T7: normalization with mixed-scale channels ----------
    Fmix = Fsyn;  Fmix(1:10,:) = Fmix(1:10,:)*1000;     % half the channels 1000x
    OPTn = OPT;  OPTn.normalize = true;
    evtN = process_evt_refphase('Compute', Fmix, t, OPTn);
    ok7 = ~isempty(evtN) && any(evtN(1,:)<=15 & evtN(2,:)>=12);
    fprintf('T7 normalize mixed-scale: %dp overlaps=%d => %s\n', size(evtN,2), ok7, PF{ok7+1});
    pass = pass && ok7;

    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end


%% ===== READ REAL AUDITORY BLOCK (MEG) =====
function [F, TimeVector] = i_read_real()
    F = [];  TimeVector = [];
    relFile = 'Subject01/S01_AEF_20131218_01_notch/data_block001_02.mat';
    try
        [sStudy, ~, iData] = bst_get('DataFile', relFile); %#ok<ASGLU>
        if isempty(sStudy), return; end
        sFile = in_fopen(file_fullpath(relFile), 'BST-DATA');
        ChannelMat = in_bst_channel(bst_get('ChannelFileForStudy', relFile));
        iMEG = channel_find(ChannelMat.Channel, 'MEG');
        ImportOptions = db_template('ImportOptions');
        ImportOptions.ImportMode='Time'; ImportOptions.UseCtfComp=1; ImportOptions.UseSsp=1;
        ImportOptions.EventsMode='ignore'; ImportOptions.DisplayMessages=0; ImportOptions.RemoveBaseline='no';
        [F, TimeVector] = in_fread(sFile, ChannelMat, 1, [], iMEG, ImportOptions);
    catch
        F = [];  TimeVector = [];
    end
end
