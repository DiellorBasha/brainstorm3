%% Alpha traveling wave — cortical phase tracking with a consistent gauge
% Demonstrates how a globally consistent tangent frame enables spatial phase
% comparison between cortical vertices, and how this reveals the propagating
% structure of alpha-band (7-13 Hz) oscillations in MEG source space.
%
% WHY GAUGE CONSISTENCY IS REQUIRED
% ----------------------------------
% The unconstrained MN inverse returns a 3D current dipole q(i,t) ∈ ℝ³ at
% each cortex vertex.  To track wave propagation we need a *direction* for
% the current, not just its amplitude.  Project onto a local tangent frame:
%
%   z(i,t)  =  q(i,t)·e1(i) + i·q(i,t)·e2(i)  ∈ ℂ
%
% The phase arg(z(i,t)) is the direction of tangential current flow at
% vertex i, measured against the local frame {e1(i), e2(i)}.
%
% Problem: if e1(i) and e1(j) at neighbouring vertices i,j point in
% unrelated directions (arbitrary local convention), then
%
%   arg(z_j) - arg(z_i)   ← conflates frame rotation with real phase shift
%
% Solution: use the FreeSurfer trivial-connection frame (tess_tangents).
% This frame is smooth — e1(j) is the parallel transport of e1(i) along
% any path — so the phase difference is purely due to current direction,
% not frame rotation.  The mean phase lead between two vertices then gives
% the wave delay, and distance ÷ delay = propagation speed.
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

%% Load cortical surface

SurfaceFile = local_find_cortex(20484);
assert(~isempty(SurfaceFile), 'No 20484-vertex cortex found in current protocol.');

TessMat = in_tess_bst(SurfaceFile);
Vtx_mm  = TessMat.Vertices * 1000;   % metres → mm
nV      = size(TessMat.Vertices, 1);

fprintf('Surface: %s  (%d vertices)\n', SurfaceFile, nV)

%% Build the FreeSurfer trivial-connection tangent frame
% tess_tangents computes a globally consistent per-face frame using the
% trivial connection (singularities at the FreeSurfer sphere poles).
% bst_tangent_face2vertex transfers it to per-vertex.
% The resulting Uv/Vv vectors are the reference frame for phase comparison.

[Uf, ~]    = tess_tangents(SurfaceFile, 'NoSave', 1);
[Uv, Vv]   = bst_tangent_face2vertex(double(TessMat.Faces), Uf, TessMat.VertNormals);

fprintf('Tangent frame built: Uv, Vv  [%s]\n', mat2str(size(Uv)))

%% Load alpha-band recording and apply unconstrained MN kernel
% The kernel maps MEG sensors → 3D current dipoles at each cortex vertex
% (nComponents=3, global Cartesian orientation).
% The recording is bandpassed to 7-13 Hz, downsampled to 600 Hz.

dataFile   = 'Subject01/S01_AEF_20131218_02_notch/data_block001_band.mat';
kernelFile = 'Subject01/S01_AEF_20131218_01_notch/results_MN_MEG_KERNEL_260605_0111.mat';

D = load(file_fullpath(dataFile));
K = load(file_fullpath(kernelFile), 'ImagingKernel', 'GoodChannel');

Fs     = 600;
Time   = D.Time;
t_ms   = Time * 1000;
F_good = D.F(K.GoodChannel, :);           % [272 x nTime]  good MEG channels

% Apply kernel: [61452 x nTime] → reshape to [nV x 3 x nTime]
J = K.ImagingKernel * F_good;             % [3nV x nTime]
J = permute(reshape(J, 3, nV, []), [2 1 3]);   % [nV x 3 x nTime]

fprintf('Source timeseries J: %s\n', mat2str(size(J)))

%% Project 3D dipoles into the consistent tangent frame
% For each vertex i and time t:
%   z(i,t) = J(i,t)·Uv(i) + i·J(i,t)·Vv(i)
%
% z(i,t) ∈ ℂ encodes the tangential current as a complex number whose
% argument arg(z) is the current direction in the FS frame.  Because
% Uv and Vv are consistent across vertices, phase differences between
% adjacent vertices are geometrically meaningful.

z_fs = squeeze(sum(J .* Uv, 2)) + 1i * squeeze(sum(J .* Vv, 2));  % [nV x nTime]
fprintf('Tangential complex field z_fs: %s\n', mat2str(size(z_fs)))

%% Identify high-power window and two tracking vertices

% Spatial GFP across all vertices
gfp    = sqrt(mean(abs(z_fs).^2, 1));
[~, iPk] = max(gfp);
iWin   = max(1, iPk-250) : min(numel(Time), iPk+250);

fprintf('Peak tangential power at t = %.3f s\n', Time(iPk))

% Mean amplitude over the window
amp_win = mean(abs(z_fs(:, iWin)), 2);

% Left hemisphere vertices (where alpha is typically dominant)
[~, lH] = tess_hemisplit(TessMat);
lH_v    = lH(:);

% Pick two vertices with high amplitude, separated by ≥30 mm
[~, sord] = sort(amp_win, 'descend');
iV1 = sord(1);
for k = 2:numel(sord)
    iV2 = sord(k);
    if norm(Vtx_mm(iV1,:) - Vtx_mm(iV2,:)) > 30
        break
    end
end

dist_mm  = norm(Vtx_mm(iV1,:) - Vtx_mm(iV2,:));
ph_lead  = mean(angle(z_fs(iV1, iWin) ./ z_fs(iV2, iWin)));
delay_ms = ph_lead / (2*pi) * 100;    % at 10 Hz: 1 cycle = 100 ms
speed_ms = dist_mm / max(abs(delay_ms), 0.1);   % m/s

fprintf('\nV1: vertex %d  at [%.0f %.0f %.0f] mm\n', iV1, Vtx_mm(iV1,:))
fprintf('V2: vertex %d  at [%.0f %.0f %.0f] mm\n', iV2, Vtx_mm(iV2,:))
fprintf('Distance:    %.1f mm\n', dist_mm)
fprintf('Phase lead:  %.3f rad = %.1f ms at 10 Hz\n', ph_lead, delay_ms)
fprintf('Speed est.:  %.2f m/s\n', speed_ms)

%% Figure 1 — Three phase maps spanning one half-cycle

% Amplitude mask: show phase only where source is strong enough
inL     = false(nV, 1);  inL(lH_v) = true;
fHl     = all(inL(double(TessMat.Faces)), 2);
mask    = amp_win < 0.08 * max(amp_win(lH_v));

% Three frames: 0, +33 ms, +67 ms  (≈ 0, T/3, 2T/3 at 10 Hz)
frame_step = round(Fs / (10 * 6));
frames3    = iPk + [0, 2, 4] * frame_step;
frames3    = min(frames3, numel(Time));

cmap_hsv = hsv(256);
cGrey    = repmat([0.25 0.25 0.23], nV, 1);

figure('Name','Alpha phase maps','Color','k','Position',[50 50 1100 420])

for f = 1:3
    ax = subplot(1,3,f);  set(ax,'Color','k');  hold(ax,'on')

    ph_f  = angle(z_fs(:, frames3(f)));
    cPh   = cGrey;
    for iv = lH_v'
        if mask(iv), continue; end
        ci = max(1, min(256, round((ph_f(iv)+pi)/(2*pi)*255)+1));
        cPh(iv,:) = cmap_hsv(ci,:);
    end

    patch('Faces', double(TessMat.Faces(fHl,:)), 'Vertices', Vtx_mm, ...
        'FaceVertexCData', cPh, 'FaceColor', 'interp', ...
        'EdgeColor', 'none', 'FaceLighting', 'phong', 'Parent', ax)
    camlight(ax,70,45);  camlight(ax,-70,45);
    camlight(ax,70,-45); camlight(ax,-70,-45);
    material([0.65 0.45 0.02 8])

    % Mark the two tracked vertices
    plot3(Vtx_mm(iV1,1),Vtx_mm(iV1,2),Vtx_mm(iV1,3), 'o', ...
        'MarkerSize',10,'MarkerFaceColor',[1 0.85 0.15],'MarkerEdgeColor','w','Parent',ax)
    plot3(Vtx_mm(iV2,1),Vtx_mm(iV2,2),Vtx_mm(iV2,3), 'o', ...
        'MarkerSize',10,'MarkerFaceColor',[0.35 0.75 0.55],'MarkerEdgeColor','w','Parent',ax)

    colormap(ax, hsv(256));  clim([-pi pi])
    cb = colorbar(ax,'Color','w','Ticks',[-pi 0 pi],'TickLabels',{'-\pi','0','\pi'});
    cb.Label.String = 'arg(z)  (rad)';  cb.Label.Color = 'w';

    view(ax,-90,10);  axis(ax,'equal','off')
    dt = round(t_ms(frames3(f)) - t_ms(frames3(1)));
    t_lbl = title(ax, sprintf('t = +%d ms', dt));  t_lbl.Color = 'w';
end

sgtitle({'Alpha (7-13 Hz) phase in consistent FS gauge  |  left hemisphere', ...
    'Phase pattern shifts in time — wave propagation  |  yellow/green = tracked vertices'}, ...
    'Color','w','FontSize',10)

%% Figure 2 — Phase timeseries at both vertices

iPlot  = max(1, iPk-200) : min(numel(Time), iPk+300);
t_rel  = t_ms(iPlot) - t_ms(iPk);
ph1    = unwrap(angle(z_fs(iV1, iPlot)));
ph2    = unwrap(angle(z_fs(iV2, iPlot)));

figure('Name','Alpha phase timeseries','Color','k','Position',[100 100 900 420])
ax2 = axes('Color','k');  hold(ax2,'on')

plot(t_rel, ph1, 'Color',[1 0.85 0.15], 'LineWidth',2.0)
plot(t_rel, ph2, 'Color',[0.35 0.75 0.55], 'LineWidth',2.0)

% Mark the three frame times
for f = 1:3
    xline(ax2, t_ms(frames3(f))-t_ms(iPk), '--', ...
        'Color','w','LineWidth',0.8,'Alpha',0.5)
end

xlabel(ax2,'time relative to peak  (ms)','Color','w')
ylabel(ax2,'arg(z)  (rad, unwrapped)','Color','w')
set(ax2,'XColor','w','YColor','w','GridColor',[0.3 0.3 0.3],'GridAlpha',0.4)
grid(ax2,'on')
legend(ax2, ...
    {sprintf('V1  [%.0f %.0f %.0f] mm', Vtx_mm(iV1,:)), ...
     sprintf('V2  [%.0f %.0f %.0f] mm  (%.0f mm away, %.1f ms lag)', ...
         Vtx_mm(iV2,:), dist_mm, abs(delay_ms))}, ...
    'TextColor','w','Color',[0.10 0.10 0.10],'Location','northwest')

title(ax2, sprintf(['Phase timeseries  |  mean lead V1\\rightarrowV2 = %.2f rad' ...
    ' = %.1f ms at 10 Hz  |  estimated speed %.2f m/s'], ph_lead, delay_ms, speed_ms), ...
    'Color','w','FontSize',9)

sgtitle({'Alpha traveling wave: phase advance between two cortex vertices', ...
    'Interpretable only with gauge-consistent frame  (arbitrary frame → phase diff is noise)'}, ...
    'Color','w','FontSize',10)

%% --- helper ---
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
