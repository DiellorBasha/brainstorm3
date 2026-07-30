function out = omega_feature_track(subShort, moduleDs, mnFile, outFile)
% OMEGA_FEATURE_TRACK  Lagrangian critical-point tracking on the DIRAC EIGENBASIS flow, full 600 s.
%
% Streams the recording, reconstructs the ALPHA (8-13 Hz) current with the Dirac-eigenbasis
% currentKernel (J = currentKernel * bandpass/resample(b), 100 Hz), detects classified critical
% points per frame (detect.criticalPoints; keep top-15 PER TYPE to stop the strong vortices from
% starving source/sink detections), then associates them across frames with the ByteTrack tracker
% (detect.track, opts.frameFeatures). Reports the "vortex-gas" statistics: track counts by type,
% lifetime, straightness (1=traveling, 0=standing). Saves frameFeatures + tracks to outFile.
%
% FUTURE: band-matched (passband, rate) = a temporal level-of-detail. This driver HARDCODES alpha
% ([8 13] Hz, fs2=100) so it is alpha-only. Tracking a slow rhythm at a fixed high rate (or the raw
% 2400 Hz) tracks NOISE not motion: per-frame feature displacement falls far below the mesh spacing
% and the detection jitter, so the ByteTrack step follows the jitter. Generalize to a
% band -> (passband, rate) map at ~8-12 samples/cycle: delta 2-4->~30, theta 4-8->~60, alpha
% 8-13->~100 (current), beta 13-30->~250, gamma 30-45->~450 Hz. This is the temporal twin of the
% eigenbasis spatial LOD (mode count K) -- a joint spatio-temporal level-of-detail.
%
% Author: Diellor Basha, 2026

    addpath('/Users/diellorbasha/workspace/research/code/nxr-cortical-flow-matlab');
    ctx = flow.context(moduleDs, 400);
    S = ctx.S;  K = ctx.currentKernel;
    op = detect.operator(S, load.bases(moduleDs));
    types = {'vortex','source','sink','saddle'};  keepPerType = 15;

    R  = in_bst_results(mnFile, 0);  gc = R.GoodChannel;
    sF = in_bst_data(R.DataFile, 'F');  sFile = sF.F;
    CM = in_bst_channel(bst_get('ChannelFileForStudy', R.DataFile));
    fs = sFile.prop.sfreq;  t0 = sFile.prop.times(1);  t1 = sFile.prop.times(2);
    fs2 = 100;  [pp, qq] = rat(fs2 / fs);
    [b,a] = butter(4, [8 13]/(fs/2), 'bandpass');

    % ---- pass 1: stream -> reconstruct alpha J -> detect per frame ----
    frameFeatures = {};  blk = 20;  td = tic;
    for ts = t0:blk:t1
        te = min(ts+blk, t1);  if te - ts < 1, break; end
        F  = panel_record('ReadRawBlock', sFile, CM, 1, [ts te], 0, 1, 'all', 1);
        Fb = filtfilt(b, a, F(gc,:)')';               % alpha band
        Fr = resample(Fb', pp, qq)';                  % 100 Hz  [C x nF]
        J  = K * Fr;                                  % [3V x nF] reconstructed alpha current
        for t = 1:size(J,2)
            cp = detect.criticalPoints(J(:,t), S, types, op);
            frameFeatures{end+1} = i_keepPerType(cp, keepPerType, types); %#ok<AGROW>
        end
    end
    nFr = numel(frameFeatures);
    fprintf('[%s] detect done: %d frames in %.1f min\n', subShort, nFr, toc(td)/60);

    % ---- pass 2: associate (ByteTrack) ----
    ta = tic;
    T = detect.track([], S, struct('frameFeatures', {frameFeatures}, 'dt', 1/fs2, 'types', {types}));
    fprintf('[%s] associate done in %.1f min\n', subShort, toc(ta)/60);

    out = T.summary;
    out.subject = subShort;  out.nFrames = nFr;  out.sfreq = fs2;
    if nargin >= 4 && ~isempty(outFile)
        tracks = T.tracks; %#ok<NASGU>
        % save ONLY tracks + summary; frameFeatures (60k-cell, ~1 GB) corrupts the -v7.3 write.
        save(outFile, 'tracks', 'out', '-v7.3');
        fid = fopen([outFile '.done'], 'w'); fprintf(fid, '%d tracks\n', out.nTracks); fclose(fid);  % tiny atomic marker
    end
    s = T.summary;
    fprintf('[%s] TRACKS: n=%d (%d vortex / %d source / %d sink / %d saddle) | life med %.0fms | straight med %.2f\n', ...
        subShort, s.nTracks, i_cnt(T.tracks,'vortex'), i_cnt(T.tracks,'source'), i_cnt(T.tracks,'sink'), ...
        i_cnt(T.tracks,'saddle'), 1000*median(s.lifetime), median([T.tracks.straightness]));
end

% keep the top-K strongest detections PER TYPE (cp is sorted by descending strength)
function ff = i_keepPerType(cp, K, types)
    if isempty(cp.strength)
        ff = struct('pos', [], 'type', {{}}, 'strength', []);  return;
    end
    keep = false(numel(cp.strength), 1);
    for i = 1:numel(types)
        idx = find(strcmp(cp.type, types{i}));       % already strength-sorted
        keep(idx(1:min(K, numel(idx)))) = true;
    end
    ff = struct('pos', cp.pos(keep,:), 'type', {cp.type(keep)}, 'strength', cp.strength(keep));
end

function n = i_cnt(tracks, ty)
    n = sum(strcmp({tracks.type}, ty));
end
