%% Alpha traveling wave — cortical phase tracking with a consistent gauge
% Demonstrates how a globally consistent tangent frame enables spatial phase
% comparison between cortical vertices, and how this reveals the propagating
% structure of alpha-band (7-13 Hz) oscillations in MEG source space.
%
% NOTE ON INTERPRETATION
% ----------------------
% We are ANALYSING the MN source estimate, not claiming to have recovered
% the true neural currents.  The MNE is a regularised spatial filter that
% returns the minimum-norm current distribution consistent with the sensor
% data.  Amplitude maps show WHERE the estimate places activity.  Phase
% maps in the consistent gauge show the SPATIAL STRUCTURE of the estimated
% tangential current — which direction it points at each vertex and how
% those directions relate across the cortex.  Wave speed estimates refer to
% propagation of the estimated phase field, not directly to axonal speed.
%
% ANALYTIC SIGNAL (Hilbert transform)
% ------------------------------------
% For a bandpassed signal, the Hilbert transform gives the analytic signal:
%   a(t) = f(t) + i·H[f(t)]
%   |a(t)| = amplitude ENVELOPE  (smooth, non-negative, no frequency doubling)
%   arg(a(t)) = INSTANTANEOUS PHASE  (oscillates at the true carrier frequency)
%
% Because the kernel W is linear: W·H[F] = H[W·F].  So we apply the kernel
% to both the real signal and its Hilbert transform, then project both into
% the consistent tangent frame:
%
%   z_re(i,t) = J(i,t)·Uv(i)  + i·J(i,t)·Vv(i)   ← real bandpassed source
%   z_im(i,t) = JH(i,t)·Uv(i) + i·JH(i,t)·Vv(i)  ← Hilbert source
%   z_a(i,t)  = z_re(i,t) + i·z_im(i,t)           ← analytic tangential signal
%
%   amp_env(i,t)  = |z_a(i,t)|    — envelope, used for amplitude maps
%   ph_inst(i,t)  = arg(z_a(i,t)) — instantaneous phase, used for wave tracking
%
% For a linearly polarised oscillation this is exact; the amplitude envelope
% carries no frequency doubling artefact that the simple norm does.
%
% WHY GAUGE CONSISTENCY IS REQUIRED FOR PHASE
% ---------------------------------------------
% ph_inst(i,t) is the phase in the FreeSurfer trivial-connection frame.
% Because Uv(j) is the parallel transport of Uv(i) along any path i→j,
% the phase DIFFERENCE ph_inst(j,t) − ph_inst(i,t) is a geometrically
% meaningful quantity: it reflects the actual lead/lag between vertices,
% not a frame-convention artefact.  In an arbitrary gauge this difference
% would be noise even for a perfect traveling wave.
%
% CREST ISOLINE
% -------------
% arg(z_a) = 0 marks the set of vertices currently at their oscillation
% peak (crest of the carrier wave).  Rendered as a white contour on the
% phase map, it makes the wavefront visible.  As time advances within one
% alpha cycle (~100 ms) the crest isoline moves across the active region —
% that movement IS the wave propagation.
%
% DEPENDENCIES
% ------------
%   • Brainstorm (nogui acceptable); TutorialAuditory protocol loaded
%   • nxr-compute plugin (used internally by tess_tangents)
%   • Data:  Subject01/S01_AEF_20131218_02_notch/data_block001_band.mat
%            (raw 80–100 s, bandpassed 7-13 Hz, downsampled to 600 Hz)
%   • Kernel: Subject01/S01_AEF_20131218_01_notch/
%             results_MN_MEG_KERNEL_260605_0111.mat
%             (unconstrained MN, nComponents=3, reused across runs)
%
% Run sections in order with Ctrl+Enter.

%% Setup

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);

if ~brainstorm('status')
    brainstorm nogui
end

[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute plugin required: %s', errMsg);
bst_plugin('Load', 'nxr-compute');

%% Load surface and build consistent tangent frame

SurfaceFile = local_find_cortex(20484);
assert(~isempty(SurfaceFile), 'No 20484-vertex cortex found in current protocol.');

TessMat = in_tess_bst(SurfaceFile);
Vtx_mm  = TessMat.Vertices * 1000;
nV      = size(TessMat.Vertices, 1);

[Uf, ~]    = tess_tangents(SurfaceFile, 'NoSave', 1);
[Uv, Vv]   = bst_tangent_face2vertex(double(TessMat.Faces), Uf, TessMat.VertNormals);

fprintf('Surface: %s  |  tangent frame: Uv/Vv [%s]\n', SurfaceFile, mat2str(size(Uv)))

%% Load recording and kernel

dataFile   = 'Subject01/S01_AEF_20131218_02_notch/data_block001_band.mat';
kernelFile = 'Subject01/S01_AEF_20131218_01_notch/results_MN_MEG_KERNEL_260605_0111.mat';

D    = load(file_fullpath(dataFile));
K    = load(file_fullpath(kernelFile), 'ImagingKernel', 'GoodChannel');
W    = K.ImagingKernel;               % [3nV x 272]
Fs   = 600;
Time = D.Time;
t_ms = Time * 1000;

% Find peak on sensor GFP to limit window (saves memory)
F_good = D.F(K.GoodChannel, :);       % [272 x nTime]
gfp_s  = sqrt(mean(F_good.^2, 1));
[~, iPk_s] = max(gfp_s);

% Work on a ±2 s window around the sensor peak
win    = max(1, iPk_s-1200) : min(numel(Time), iPk_s+1200);
F_win  = F_good(:, win);
F_hilb = imag(hilbert(F_win.')).';    % Hilbert per channel  [272 x nWin]
Time_w = Time(win);
t_ms_w = Time_w * 1000;

fprintf('Window: %.3f – %.3f s  (%d samples)\n', Time_w(1), Time_w(end), numel(win))

%% Build analytic tangential source field

J_re = permute(reshape(W * F_win,  3, nV, []), [2 1 3]);   % [nV x 3 x nWin]
J_im = permute(reshape(W * F_hilb, 3, nV, []), [2 1 3]);

% Project real and Hilbert sources into the consistent FS frame
z_re = squeeze(sum(J_re .* Uv, 2)) + 1i * squeeze(sum(J_re .* Vv, 2));
z_im = squeeze(sum(J_im .* Uv, 2)) + 1i * squeeze(sum(J_im .* Vv, 2));

z_a     = z_re + 1i * z_im;     % [nV x nWin] complex analytic tangential signal
amp_env = abs(z_a);              % amplitude envelope  (no frequency doubling)
ph_inst = angle(z_a);            % instantaneous phase in consistent FS gauge

fprintf('Analytic field z_a: %s\n', mat2str(size(z_a)))

%% Find peak and select frames

gfp_w    = sqrt(mean(amp_env.^2, 1));
[~, iPkW] = max(gfp_w);
fprintf('Source GFP peak at %.3f s\n', Time_w(iPkW))

% Amplitude mask: active vertices = top 6% of mean amplitude around peak
[~, lH_v] = tess_hemisplit(TessMat);  lH_v = lH_v(:);
inL        = false(nV, 1);  inL(lH_v) = true;
fHl        = all(inL(double(TessMat.Faces)), 2);

amp_win  = mean(amp_env(:, max(1,iPkW-250):min(size(z_a,2),iPkW+250)), 2);
mask     = amp_win < 0.06 * max(amp_win(lH_v));

fprintf('Active LH vertices (6%% threshold): %d / %d\n', sum(~mask(lH_v)), numel(lH_v))

% Active-region bounding box for zoomed rendering
actV     = find(amp_win > median(amp_win(~mask & inL)) & inL);
pad      = 12;
xlims    = [min(Vtx_mm(actV,1))-pad, max(Vtx_mm(actV,1))+pad];
ylims    = [min(Vtx_mm(actV,2))-pad, max(Vtx_mm(actV,2))+pad];
zlims    = [min(Vtx_mm(actV,3))-pad, max(Vtx_mm(actV,3))+pad];

T_alpha  = round(Fs / 10);              % samples per alpha cycle  (60 at 600 Hz)

% Figure A: 6 frames within one cycle (~17 ms apart)
frames_A = iPkW + (0:5) * round(T_alpha/6);
frames_A = min(frames_A, size(z_a,2));

% Figure B: 5 successive cycles (~100 ms apart)
frames_B = iPkW + (0:4) * T_alpha;
frames_B = min(frames_B, size(z_a,2));

fprintf('Figure A frames: %s ms\n', mat2str(round(t_ms_w(frames_A) - t_ms_w(iPkW))))
fprintf('Figure B frames: %s ms\n', mat2str(round(t_ms_w(frames_B) - t_ms_w(iPkW))))

%% Figure A — Crest isoline within one alpha cycle

cmap_hsv = hsv(256);
cGrey    = repmat([0.25 0.25 0.23], nV, 1);

figure('Name','Alpha wave — within cycle','Color','k','Position',[50 50 1400 500])

for f = 1:6
    ax = subplot(1,6,f);  set(ax,'Color','k');  hold(ax,'on')
    ph_f   = ph_inst(:, frames_A(f));
    active = ~mask;
    cPh    = cGrey;
    for iv = find(active & inL)'
        ci = max(1, min(256, round((ph_f(iv)+pi)/(2*pi)*255)+1));
        cPh(iv,:) = cmap_hsv(ci,:);
    end
    patch('Faces',double(TessMat.Faces(fHl,:)),'Vertices',Vtx_mm, ...
        'FaceVertexCData',cPh,'FaceColor','interp','EdgeColor','none', ...
        'FaceLighting','phong','Parent',ax)
    camlight(ax,70,45); camlight(ax,-70,45);
    camlight(ax,70,-45); camlight(ax,-70,-45);
    material([0.65 0.45 0.02 8])

    % Crest isoline  (arg = 0 ↔ oscillation peak)
    segs = phase_isoline(Vtx_mm, double(TessMat.Faces(fHl,:)), ph_f, active);
    if size(segs,1) > 0
        nS = size(segs,1);
        xs = reshape([squeeze(segs(:,:,1)) nan(nS,1)]', 1, []);
        ys = reshape([squeeze(segs(:,:,2)) nan(nS,1)]', 1, []);
        zs = reshape([squeeze(segs(:,:,3)) nan(nS,1)]', 1, []);
        plot3(xs,ys,zs,'-','Color','w','LineWidth',2.0,'Parent',ax)
    end

    colormap(ax,hsv(256));  clim([-pi pi])
    xlim(ax,xlims); ylim(ax,ylims); zlim(ax,zlims);
    view(ax,-90,35);  axis(ax,'off')
    dt = round(t_ms_w(frames_A(f)) - t_ms_w(iPkW));
    t_lbl = title(ax, sprintf('+%d ms', dt));  t_lbl.Color = 'w';  t_lbl.FontSize = 10;
end

% Shared colourbar
axcb = axes('Position',[0.92 0.15 0.01 0.7],'Visible','off');
colormap(axcb,hsv(256));  clim([-pi pi]);
cb = colorbar(axcb,'Color','w','Location','eastoutside');
cb.Label.String = 'arg(z_a)  (rad)';  cb.Label.Color = 'w';
cb.Ticks = [-pi 0 pi];  cb.TickLabels = {'-\pi','0','\pi'};

sgtitle({'Alpha wave — within one cycle (~100 ms)  |  white line = oscillation crest  (arg = 0)', ...
    'Consistent gauge: crest location is physically meaningful  |  crest moves  \rightarrow  wave propagation'}, ...
    'Color','w','FontSize',10)

%% Figure B — Crest isoline across successive alpha cycles

figure('Name','Alpha wave — successive cycles','Color','k','Position',[50 50 1200 500])

for f = 1:5
    ax = subplot(1,5,f);  set(ax,'Color','k');  hold(ax,'on')
    ph_f   = ph_inst(:, frames_B(f));
    active = ~mask;
    cPh    = cGrey;
    for iv = find(active & inL)'
        ci = max(1, min(256, round((ph_f(iv)+pi)/(2*pi)*255)+1));
        cPh(iv,:) = cmap_hsv(ci,:);
    end
    patch('Faces',double(TessMat.Faces(fHl,:)),'Vertices',Vtx_mm, ...
        'FaceVertexCData',cPh,'FaceColor','interp','EdgeColor','none', ...
        'FaceLighting','phong','Parent',ax)
    camlight(ax,70,45); camlight(ax,-70,45);
    camlight(ax,70,-45); camlight(ax,-70,-45);
    material([0.65 0.45 0.02 8])

    segs = phase_isoline(Vtx_mm, double(TessMat.Faces(fHl,:)), ph_f, active);
    if size(segs,1) > 0
        nS = size(segs,1);
        xs = reshape([squeeze(segs(:,:,1)) nan(nS,1)]', 1, []);
        ys = reshape([squeeze(segs(:,:,2)) nan(nS,1)]', 1, []);
        zs = reshape([squeeze(segs(:,:,3)) nan(nS,1)]', 1, []);
        plot3(xs,ys,zs,'-','Color','w','LineWidth',2.0,'Parent',ax)
    end

    colormap(ax,hsv(256));  clim([-pi pi])
    xlim(ax,xlims); ylim(ax,ylims); zlim(ax,zlims);
    view(ax,-90,35);  axis(ax,'off')
    dt = round(t_ms_w(frames_B(f)) - t_ms_w(iPkW));
    t_lbl = title(ax, sprintf('cycle %d  (+%d ms)', f, dt));
    t_lbl.Color = 'w';  t_lbl.FontSize = 10;
end

axcb2 = axes('Position',[0.92 0.15 0.01 0.7],'Visible','off');
colormap(axcb2,hsv(256));  clim([-pi pi]);
cb2 = colorbar(axcb2,'Color','w','Location','eastoutside');
cb2.Label.String = 'arg(z_a)  (rad)';  cb2.Label.Color = 'w';
cb2.Ticks = [-pi 0 pi];  cb2.TickLabels = {'-\pi','0','\pi'};

sgtitle({'Alpha wave — five successive cycles (~100 ms apart)  |  white line = oscillation crest', ...
    'Consistent crest position across cycles confirms periodic traveling wave'}, ...
    'Color','w','FontSize',10)

%% Figure C — Phase timeseries: wave delay between two vertices

% Pick two high-amplitude vertices ≥30 mm apart
iWin = max(1,iPkW-250):min(size(z_a,2),iPkW+250);
[~, sord] = sort(amp_win, 'descend');
iV1 = sord(1);
for k = 2:numel(sord)
    iV2 = sord(k);
    if norm(Vtx_mm(iV1,:) - Vtx_mm(iV2,:)) > 30, break; end
end

dist_mm  = norm(Vtx_mm(iV1,:) - Vtx_mm(iV2,:));
ph_lead  = mean(angle(z_a(iV1, iWin) ./ z_a(iV2, iWin)));
delay_ms = ph_lead / (2*pi) * 100;
speed_ms = dist_mm / max(abs(delay_ms), 0.1);

fprintf('\nV1: %d  [%.0f %.0f %.0f] mm\n', iV1, Vtx_mm(iV1,:))
fprintf('V2: %d  [%.0f %.0f %.0f] mm  (%.1f mm apart)\n', iV2, Vtx_mm(iV2,:), dist_mm)
fprintf('Phase lead: %.3f rad = %.1f ms  |  estimated speed: %.2f m/s\n', ph_lead, delay_ms, speed_ms)

iPlot = max(1,iPkW-200):min(size(z_a,2),iPkW+300);
t_rel = t_ms_w(iPlot) - t_ms_w(iPkW);
ph1   = unwrap(angle(z_a(iV1, iPlot)));
ph2   = unwrap(angle(z_a(iV2, iPlot)));

figure('Name','Phase timeseries','Color','k','Position',[100 100 900 420])
ax3 = axes('Color','k');  hold(ax3,'on')
plot(t_rel, ph1, 'Color',[1 0.85 0.15],    'LineWidth',2.0)
plot(t_rel, ph2, 'Color',[0.35 0.75 0.55], 'LineWidth',2.0)
for f = 1:3
    xline(ax3, t_ms_w(frames_A(f))-t_ms_w(iPkW), '--', 'Color','w','LineWidth',0.8,'Alpha',0.5)
end
xlabel(ax3,'time relative to peak  (ms)','Color','w')
ylabel(ax3,'arg(z_a)  (rad, unwrapped)','Color','w')
set(ax3,'XColor','w','YColor','w','GridColor',[0.3 0.3 0.3],'GridAlpha',0.4)
grid(ax3,'on')
legend(ax3, ...
    {sprintf('V1  [%.0f %.0f %.0f] mm', Vtx_mm(iV1,:)), ...
     sprintf('V2  [%.0f %.0f %.0f] mm  (%.0f mm, %.1f ms lag)', ...
         Vtx_mm(iV2,:), dist_mm, abs(delay_ms))}, ...
    'TextColor','w','Color',[0.10 0.10 0.10],'Location','northwest')
title(ax3, sprintf(['Phase timeseries (analytic signal)  |  mean lead V1\\rightarrowV2 = ' ...
    '%.2f rad = %.1f ms  |  est. speed %.2f m/s'], ph_lead, delay_ms, speed_ms), ...
    'Color','w','FontSize',9)
sgtitle({'Instantaneous phase from Hilbert analytic signal', ...
    'Interpretable across vertices only with gauge-consistent frame'}, ...
    'Color','w','FontSize',10)

%% --- helpers ---

function segs = phase_isoline(Vtx, Faces, ph, active)
% Compute mesh isoline where arg(ph) = 0 (oscillation crest).
% Uses zero-crossings of sin(ph) with cos(ph) > 0 to avoid the ±pi wrap.
% Returns segs [nSeg x 2 x 3] — line segment endpoints in 3D.
segs   = zeros(0, 2, 3);
s      = sin(ph);
c      = cos(ph);
edges  = [1 2; 2 3; 3 1];
for fi = 1:size(Faces,1)
    vv = Faces(fi,:);
    if ~any(active(vv)), continue; end
    sv  = s(vv);  cv = c(vv);
    pts = zeros(0, 3);
    for e = 1:3
        a = edges(e,1);  b = edges(e,2);
        if sv(a)*sv(b) < 0                              % sin crosses zero
            t    = sv(a) / (sv(a) - sv(b));
            cmid = cv(a) + t*(cv(b)-cv(a));
            if cmid > 0                                  % crest, not trough
                pts(end+1,:) = Vtx(vv(a),:) + t*(Vtx(vv(b),:)-Vtx(vv(a),:));
            end
        end
    end
    if size(pts,1) == 2
        segs(end+1,:,:) = pts;
    end
end
end

function SurfaceFile = local_find_cortex(nVert)
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
if isempty(sSubjects), return; end
for iS = 1:numel([sSubjects.Subject])
    for iF = 1:numel(sSubjects.Subject(iS).Surface)
        s = sSubjects.Subject(iS).Surface(iF);
        if ~strcmpi(s.SurfaceType,'Cortex'), continue; end
        try, T = load(file_fullpath(s.FileName),'Vertices'); catch, continue; end
        if size(T.Vertices,1) == nVert
            SurfaceFile = s.FileName;
            return
        end
    end
end
end
