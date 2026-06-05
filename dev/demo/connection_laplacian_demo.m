%% Connection Laplacian on the Cortical Surface
% A step-by-step walkthrough of how we use the connection Laplacian and its
% Fiedler eigenvector to define a smooth, phase-bearing basis on the cortex.
%
% Prerequisites: Brainstorm running (nogui is fine), nxr-compute plugin
% installed. Run sections in order with Ctrl+Enter.
%
% Data: TutorialAnatomy protocol, Subject01, tess_cortex_pial_low (20484V).

%% Setup
% Start Brainstorm and make sure we're in the right working directory.

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);

if ~brainstorm('status')
    brainstorm nogui
end

%% Section 1 — Load the cortical surface
% Brainstorm stores surfaces as .mat files inside the protocol's anatomy
% folder. We locate the 20484-vertex pial-low surface and load it.
%
% in_tess_bst returns a struct with:
%   Vertices   [nV x 3]  — vertex coordinates in metres (MRI space)
%   Faces      [nF x 3]  — triangular face index list (1-based)
%   VertNormals [nV x 3] — pre-computed outward vertex normals

SurfaceFile = local_find_cortex(20484);

TessMat = in_tess_bst(SurfaceFile);

Vtx = TessMat.Vertices;   % [20484 x 3]  metres
Fcs = TessMat.Faces;      % [40960 x 3]  triangles

fprintf('Loaded: %d vertices, %d triangles\n', size(Vtx,1), size(Fcs,1))

%% Section 1 — Render the cortical surface
% We render the mesh as a triangulated patch.  Converting to mm keeps the
% axis labels readable; the brain is ~130 mm across.
%
% trisurf with FaceColor 'flat' and a single FaceVertexCData value gives a
% uniform grey shell.  We add a directional light and a specular highlight
% to give the gyral/sulcal structure some depth.

Vtx_mm = Vtx * 1000;   % metres -> millimetres for display

figure('Name','Cortical surface','Color','k')

% patch gives explicit color control — trisurf maps Z to the colormap by
% default. We supply per-vertex normals so Phong shading looks smooth.
patch('Faces', Fcs, 'Vertices', Vtx_mm, ...
    'FaceColor',  [0.78 0.72 0.65], ...
    'EdgeColor',  'none', ...
    'FaceNormals', TessMat.VertNormals, ...
    'FaceLighting','phong')

% Four camera-relative lights in a ring (camlight positions are always
% relative to the camera, so they follow the view on rotation).
% This covers all faces regardless of orientation with no dark sides.
camlight( 70,  45)   % upper-right
camlight(-70,  45)   % upper-left
camlight( 70, -45)   % lower-right
camlight(-70, -45)   % lower-left

% Moderate ambient prevents fully-dark faces; low specular stays matte.
material([0.72 0.38 0.03 8])  % [ambient diffuse specular shininess]

% Orient to a standard left-lateral view.
view(-90, 10)
axis equal off
t = title('tess\_cortex\_pial\_low  (20484 vertices)');
t.Color = 'w';

%% Section 2 — FreeSurfer registration sphere and trivial-connection poles
% Every cortex vertex has a counterpart on the FreeSurfer registration sphere
% (TessMat.Reg.Sphere.Vertices).  The sphere is the shared coordinate system
% used for inter-subject alignment — it is what makes anatomically homologous
% vertices correspond across subjects.
%
% To compute a *trivial connection* (a globally consistent tangent frame) on
% each hemisphere we need to place exactly two +1 singularities whose indices
% sum to the Euler characteristic chi = 2 (Poincaré–Hopf theorem for a
% topological sphere).  Brainstorm places them at the geographic north and
% south poles of each hemisphere's registration sphere — the vertices with the
% maximum and minimum Z coordinate in sphere space.
%
% Both hemispheres share the same sphere atlas, so L and R pole positions
% coincide at the same sphere coordinates.  Each hemisphere is solved
% independently and gets its own pair of singularities.

Sphere = TessMat.Reg.Sphere.Vertices;   % [nV x 3] unit sphere

% Split vertices into left and right hemispheres using the Structures atlas
% (same split used by tess_tangents — never re-split geometrically).
[rH, lH] = tess_hemisplit(TessMat);

hemis = {lH(:), rH(:)};
tags  = {'Left hemisphere', 'Right hemisphere'};
clrs  = {[0.35 0.55 0.75], [0.75 0.55 0.35]};   % blue / orange

r       = 0.1;     % sphere radius (Brainstorm convention)
pushOut = 1.18;    % offset pole markers beyond surface to avoid z-fighting

figure('Name','Registration sphere — trivial-connection poles','Color','k', ...
       'Position',[200 200 900 480])

for h = 1:2
    vH  = hemis{h};
    sph = Sphere(vH,:);

    % North and south poles: max / min sphere Z.
    [~, iN] = max(sph(:,3));  pN = vH(iN);
    [~, iS] = min(sph(:,3));  pS = vH(iS);

    % Push marker positions outward so they sit above the sphere surface.
    mN = Sphere(pN,:) / norm(Sphere(pN,:)) * r * pushOut;
    mS = Sphere(pS,:) / norm(Sphere(pS,:)) * r * pushOut;

    % Per-vertex colour: this hemisphere's colour, other hemisphere dimmed.
    cAll = repmat([0.22 0.22 0.22], size(Sphere,1), 1);
    cAll(vH,:) = repmat(clrs{h}, numel(vH), 1);

    % Only faces entirely inside this hemisphere.
    inH = false(size(Sphere,1),1);  inH(vH) = true;
    fH  = all(inH(Fcs), 2);

    ax = subplot(1,2,h);
    patch('Faces', Fcs(fH,:), 'Vertices', Sphere, ...
        'FaceVertexCData', cAll, 'FaceColor', 'interp', ...
        'EdgeColor', 'none', 'FaceLighting', 'phong', 'Parent', ax)
    camlight(ax,  70,  45);  camlight(ax, -70,  45);
    camlight(ax,  70, -45);  camlight(ax, -70, -45);
    material([0.72 0.38 0.03 8])

    % Wireframe overlay — shows the triangulation density of the mesh.
    hold(ax,'on')
    patch('Faces', Fcs(fH,:), 'Vertices', Sphere, ...
        'FaceColor', 'none', ...
        'EdgeColor', [1 1 1], ...
        'EdgeAlpha', 0.50, ...
        'Parent', ax)
    scatter3(mN(1), mN(2), mN(3), 220, 'p', ...
        'MarkerFaceColor', [1 0.85 0], 'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.8, 'Parent', ax)
    scatter3(mS(1), mS(2), mS(3), 220, 'p', ...
        'MarkerFaceColor', [1 0.85 0], 'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.8, 'Parent', ax)

    text(0, 0,  r*1.30, 'north pole', 'Color','w','FontSize',9, ...
        'HorizontalAlignment','center','Parent',ax)
    text(0, 0, -r*1.30, 'south pole', 'Color','w','FontSize',9, ...
        'HorizontalAlignment','center','Parent',ax)

    % Side view (elevation 0) — both poles visible at top and bottom.
    view(ax, 0, 0);  axis(ax,'equal','off')
    tt = title(ax, tags{h});  tt.Color = 'w';  tt.FontWeight = 'bold';
end

sgtitle('Registration sphere — trivial-connection singularities', ...
    'Color','w','FontSize',11)

%% Section 3 — Analytic sphere: vertices as sampling points on S²
% Because every ico vertex lies on S² (to machine precision), each vertex
% maps to exact spherical coordinates (theta, phi).  Any analytic function
% defined on S² can therefore be evaluated at these vertex positions without
% approximation — we are just sampling an analytic function at known points.
%
% This establishes the bridge between the discrete ico mesh and the analytic
% sphere: the connection Laplacian eigenmodes we compute on the cortex are
% a discrete approximation to analytic objects on S². Later sections compare
% the discrete Fiedler vector against its analytic counterpart.
%
% Coordinates:
%   theta in [0, pi]    — colatitude  (0 = north pole, pi = south pole)
%   phi   in [-pi, pi]  — longitude   (atan2(y, x))

r   = 0.1;
th  = acos( max(min(Sphere(:,3)/r, 1), -1) );  % colatitude [0, pi]
phi = atan2( Sphere(:,2), Sphere(:,1) );         % longitude  [-pi, pi]

% Work with the left hemisphere only for a clean hemisphere view.
vH  = lH(:);
inH = false(size(Sphere,1),1);  inH(vH) = true;
fH  = all(inH(Fcs), 2);

figure('Name','Analytic sphere','Color','k','Position',[100 100 820 480])

% --- Left panel: longitude phi (cyclic, HSV colormap) ---
% phi completes one full cycle around the sphere — this is the analytic
% "location coordinate" that the Fiedler vector's phase approximates on
% the folded cortex.
phi_norm = (phi(vH) + pi) / (2*pi);
cidx     = max(1, min(256, round(phi_norm * 255) + 1));
cAll     = repmat([0.18 0.18 0.18], size(Sphere,1), 1);
cmap_hsv   = hsv(256);
cAll(vH,:) = cmap_hsv(cidx, :);

ax1 = subplot(1,2,1);
patch('Faces', Fcs(fH,:), 'Vertices', Sphere, ...
    'FaceVertexCData', cAll, 'FaceColor', 'interp', ...
    'EdgeColor', 'none', 'FaceLighting', 'phong', 'Parent', ax1)
camlight(ax1,  70,  45);  camlight(ax1, -70,  45);
camlight(ax1,  70, -45);  camlight(ax1, -70, -45);
material([0.72 0.38 0.03 8])
colormap(ax1, hsv(256));  clim([-pi pi])
cb1 = colorbar(ax1, 'Color', 'w', ...
    'Ticks', [-pi -pi/2 0 pi/2 pi], ...
    'TickLabels', {'-\pi','-\pi/2','0','\pi/2','\pi'});
cb1.Label.String = '\phi  (longitude)';  cb1.Label.Color = 'w';
view(ax1, 0, 0);  axis(ax1, 'equal', 'off')
t1 = title(ax1, 'Analytic \phi on S^2');  t1.Color = 'w';

% --- Right panel: colatitude theta (parula colormap) ---
th_norm = th(vH) / pi;
cidx2   = max(1, min(256, round(th_norm * 255) + 1));
cAll2   = repmat([0.18 0.18 0.18], size(Sphere,1), 1);
cmap_par    = parula(256);
cAll2(vH,:) = cmap_par(cidx2, :);

ax2 = subplot(1,2,2);
patch('Faces', Fcs(fH,:), 'Vertices', Sphere, ...
    'FaceVertexCData', cAll2, 'FaceColor', 'interp', ...
    'EdgeColor', 'none', 'FaceLighting', 'phong', 'Parent', ax2)
camlight(ax2,  70,  45);  camlight(ax2, -70,  45);
camlight(ax2,  70, -45);  camlight(ax2, -70, -45);
material([0.72 0.38 0.03 8])
colormap(ax2, parula(256));  clim([0 pi])
cb2 = colorbar(ax2, 'Color', 'w', ...
    'Ticks', [0 pi/4 pi/2 3*pi/4 pi], ...
    'TickLabels', {'0','\pi/4','\pi/2','3\pi/4','\pi'});
cb2.Label.String = '\theta  (colatitude)';  cb2.Label.Color = 'w';
view(ax2, 0, 0);  axis(ax2, 'equal', 'off')
t2 = title(ax2, 'Analytic \theta on S^2');  t2.Color = 'w';

sgtitle('Analytic sphere: vertices as exact sampling points on S^2', ...
    'Color', 'w', 'FontSize', 11)

%% Section 4 — Eigenmode leadfield
% We now build the sensor-space representation of each connection-Laplacian
% eigenmode. The goal is a real matrix Lmu [nCh x 2K] where each pair of
% columns (L~_k^R, L~_k^I) is the MEG sensor pattern produced by eigenmode k
% with unit real and unit imaginary coefficient respectively.
%
% The construction has four steps:
%
% (A) Project the free-orientation leadfield into the nxr tangent frame.
%     Gain [nCh x 3V] carries three Cartesian columns per vertex.
%     Projecting onto {e1(i), e2(i)} gives two tangent columns per vertex:
%       L1(:,i) = Gain(:, 3i-2:3i) * e1(i,:)'
%       L2(:,i) = Gain(:, 3i-2:3i) * e2(i,:)'
%     These are the sensor patterns of unit dipoles along the nxr frame axes
%     — the exact gauge that the complex eigenvectors are stored in.
%
% (B) Mass-weight the eigenmode components.
%     A_re(i,k) = Area(i) * Re(Psi_k(i))
%     A_im(i,k) = Area(i) * Im(Psi_k(i))
%
% (C) Combine into eigenmode leadfield columns.
%     Lmu_R = L1 * A_re + L2 * A_im   [nCh x K]   (driven by c_k^R)
%     Lmu_I = -L1 * A_im + L2 * A_re  [nCh x K]   (driven by c_k^I)
%
% (D) Interleave into Lmu [nCh x 2K]:
%     columns 1,3,5,... = Lmu_R;   columns 2,4,6,... = Lmu_I

%% Section 4 — What is a leadfield column?
% Before building the eigenmode leadfield it helps to see what a single
% leadfield column looks like geometrically.
%
% The forward model for one sensor s is:
%
%   b_s = sum_i  g(i) . J(i)
%
% where J(i) in R^3 is the current dipole moment at vertex i and
% g(i) = Gain(s, 3i-2 : 3i) in R^3 is the LEAD VECTOR (sensitivity vector)
% at vertex i for sensor s.  The contribution of vertex i to sensor s is the
% dot product g(i) . J(i), so:
%
%   - Direction of g(i): the direction a dipole at vertex i must point to
%     drive sensor s most strongly.  A dipole aligned with g(i) produces a
%     response of |g(i)|; a dipole perpendicular to g(i) contributes zero.
%
%   - Magnitude |g(i)|: how sensitive sensor s is to vertex i at all.
%     Large near the sensor (high sensitivity), small far away.
%
% The lead vector is therefore NOT the dipole moment field — it is the
% sensor's sensitivity field.  At each vertex independently, g(i) points in
% the direction of the locally optimal dipole.  The globally optimal current
% distribution for driving sensor s (a rank-1 problem over all vertices) is
% a different object (related to the beamformer spatial filter).
%
% MEG physics note: g(i) is approximately tangential to the cortical surface
% near the sensor — radial currents cancel in the spherical conductor model
% and contribute near-zero to MEG.  This is why the plot forms a blade shape:
% the high-magnitude arrows are roughly perpendicular to the radial direction
% from each vertex to the sensor, which is the direction of zero sensitivity.
% This tangential sensitivity is exactly why projecting the leadfield onto the
% nxr tangent frame {e1, e2} captures almost all of the measurable signal.

%% Section 4 — Load headmodel and select MEG channels
% The Brainstorm headmodel Gain contains all 340 channels — including EEG,
% EOG, reference channels, etc. — as NaN rows. We select only the 274 MEG
% gradiometers/magnetometers.

HMFile = ['/Users/diellorbasha/workspace/library/datasets/brainstorm_db/' ...
          'tmp_aggregate/TutorialAuditory/data/Subject01/' ...
          'S01_AEF_20131218_01_notch/headmodel_surf_os_meg.mat'];
hm   = load(HMFile, 'Gain');
Gain = double(hm.Gain);       % [340 x 3*nV]

% Identify MEG channels from the channel file.
ChannelMat = load(file_fullpath( ...
    'Subject01/S01_AEF_20131218_01_notch/channel_ctf_acc1.mat'), 'Channel');
iMEG = find(strcmp({ChannelMat.Channel.Type}, 'MEG'));
Gain_meg = Gain(iMEG, :);     % [274 x 3*nV]  — no NaNs
[nCh, ~] = size(Gain_meg);
fprintf('MEG channels: %d  |  Gain NaN: %d\n', nCh, any(isnan(Gain_meg(:))))

%% Section 4 — Visualise one leadfield column (gain vector field)
% Pick the most sensitive MEG sensor and plot its gain vectors.

nV     = size(Vtx, 1);
Vtx_mm = Vtx * 1000;   % metres -> mm

rms_per_sensor = sqrt(mean(Gain_meg.^2, 2));
[~, iSensor] = max(rms_per_sensor);

g      = reshape(Gain_meg(iSensor,:), 3, nV)';   % [nV x 3]  T/Am
sLoc   = ChannelMat.Channel(iMEG(iSensor)).Loc(:,1)' * 1000;  % mm

iSub   = 1:6:nV;
P      = Vtx_mm(iSub,:);
G      = g(iSub,:);
mag    = sqrt(sum(G.^2,2));
sc     = 6.0 / max(mag);                         % max arrow = 6 mm

mag_n  = (mag - min(mag)) / max(max(mag)-min(mag), eps);
cmap_h = hot(256);
cidx_h = max(1,min(256,round(mag_n*255)+1));
clrs_h = cmap_h(cidx_h,:);

figure('Name','Leadfield column','Color','k','Position',[100 100 760 580])
ax_g = axes('Color','k');  hold(ax_g,'on')

for k = 1:size(P,1)
    quiver3(P(k,1),P(k,2),P(k,3), G(k,1)*sc,G(k,2)*sc,G(k,3)*sc, 0, ...
        'Color',clrs_h(k,:),'LineWidth',0.6,'MaxHeadSize',0.3,'Parent',ax_g)
end

plot3(sLoc(1),sLoc(2),sLoc(3), 'o','MarkerSize',14, ...
    'MarkerFaceColor',[0.3 0.8 1],'MarkerEdgeColor','w','LineWidth',1.5,'Parent',ax_g)
text(sLoc(1),sLoc(2),sLoc(3)+10,sprintf('sensor %d',iSensor), ...
    'Color','w','FontSize',10,'HorizontalAlignment','center','Parent',ax_g)

colormap(ax_g,hot(256));  clim([min(mag) max(mag)])
cb_g = colorbar(ax_g,'Color','w');
cb_g.Label.String = 'Gain magnitude  (T/Am)';  cb_g.Label.Color = 'w';

axis(ax_g,'equal');  grid(ax_g,'on')
set(ax_g,'XColor','w','YColor','w','ZColor','w', ...
         'GridColor',[0.25 0.25 0.25],'GridAlpha',0.4)
xlabel(ax_g,'x (mm)');  ylabel(ax_g,'y (mm)');  zlabel(ax_g,'z (mm)')
view(ax_g,-70,20)
tg = title(ax_g,sprintf('Leadfield column — sensor %d  (one 3D gain vector per vertex)',iSensor));
tg.Color = 'w';

%% Section 4 — Gain vector decomposition at the most sensitive vertex
% Take the single largest gain vector from the column and decompose it into
% the three orthogonal directions of the nxr tangent frame {ê₁, ê₂, n̂}.
%
%   ĝ = (ĝ·n̂) n̂  +  (ĝ·ê₁) ê₁  +  (ĝ·ê₂) ê₂
%       \_______/    \____________________/
%       normal part       tangential part
%
% The tangential part is what L1(:,i) and L2(:,i) capture; the normal part
% is discarded when we project onto the tangent frame.

% Vertex with the largest gain magnitude for sensor iSensor.
g_all    = reshape(Gain_meg(iSensor,:), 3, nV)';
[~,iMax] = max(sqrt(sum(g_all.^2,2)));

pos    = Vtx_mm(iMax,:);
g_hat  = g_all(iMax,:) / norm(g_all(iMax,:));
n_hat  = TessMat.VertNormals(iMax,:);
e1_hat = e1_nxr(iMax,:);
e2_hat = e2_nxr(iMax,:);

cn = dot(g_hat, n_hat);
c1 = dot(g_hat, e1_hat);
c2 = dot(g_hat, e2_hat);
g_tang = c1*e1_hat + c2*e2_hat;
g_n    = cn * n_hat;

% Use a vertex near median gain magnitude with balanced Cartesian components.
mag_all   = sqrt(sum(reshape(Gain_meg(iSensor,:),3,nV)'.^2, 2));
med_mag   = median(mag_all);
inBand    = find(abs(mag_all - med_mag)/med_mag < 0.20);
gBandN    = reshape(Gain_meg(iSensor,:),3,nV)';
gBandN    = gBandN(inBand,:) ./ sqrt(sum(gBandN(inBand,:).^2,2));
[~,iBest] = min(max(abs(gBandN),[],2));
iSrc      = inBand(iBest);

pos_src = Vtx_mm(iSrc,:);
g_src   = reshape(Gain_meg(iSensor,:),3,nV)'; g_src = g_src(iSrc,:);
g_hat_s = g_src / norm(g_src);

sc       = 12;    % arrow length mm
comp_clr = [0.70 0.70 0.70];
sLoc     = ChannelMat.Channel(iMEG(iSensor)).Loc(:,1)' * 1000;

% Tangent frame at this vertex (nxr gauge).
e1_src = e1_nxr(iSrc,:);
e2_src = e2_nxr(iSrc,:);

% Axis limits: large enough to show both source and sensor comfortably.
pad_src = 22;  pad_sns = 15;
xl = [min(pos_src(1)-pad_src, sLoc(1)-pad_sns), max(pos_src(1)+pad_src, sLoc(1)+pad_sns)];
yl = [min(pos_src(2)-pad_src, sLoc(2)-pad_sns), max(pos_src(2)+pad_src, sLoc(2)+pad_sns)];
zl = [min(pos_src(3)-pad_src, sLoc(3)-pad_sns), max(pos_src(3)+pad_src, sLoc(3)+pad_sns)];

% Tangent plane at source — square spanned by e1 and e2.
tp_sc = 18;
tp = [pos_src + tp_sc*( e1_src + e2_src); ...
      pos_src + tp_sc*( e1_src - e2_src); ...
      pos_src + tp_sc*(-e1_src - e2_src); ...
      pos_src + tp_sc*(-e1_src + e2_src)];

figure('Name','Gain vector — Cartesian components','Color','k', ...
       'Position',[100 100 680 580])
ax5 = axes('Color','k'); hold(ax5,'on')

% Cortical tangent plane (grey, semi-transparent)
fill3(tp([1 2 3 4 1],1),tp([1 2 3 4 1],2),tp([1 2 3 4 1],3), ...
    [0.5 0.5 0.5],'FaceAlpha',0.20,'EdgeColor',[0.6 0.6 0.6], ...
    'LineWidth',0.8,'Parent',ax5)

% Three Cartesian component arrows — same colour, no text labels
quiver3(pos_src(1),pos_src(2),pos_src(3), g_hat_s(1)*sc,0,0, 0, ...
    'Color',comp_clr,'LineWidth',2,'MaxHeadSize',0.5,'Parent',ax5)
quiver3(pos_src(1),pos_src(2),pos_src(3), 0,g_hat_s(2)*sc,0, 0, ...
    'Color',comp_clr,'LineWidth',2,'MaxHeadSize',0.5,'Parent',ax5)
quiver3(pos_src(1),pos_src(2),pos_src(3), 0,0,g_hat_s(3)*sc, 0, ...
    'Color',comp_clr,'LineWidth',2,'MaxHeadSize',0.5,'Parent',ax5)

% Gain vector (yellow)
quiver3(pos_src(1),pos_src(2),pos_src(3), ...
    g_hat_s(1)*sc,g_hat_s(2)*sc,g_hat_s(3)*sc, 0, ...
    'Color',[1 0.85 0.15],'LineWidth',2.5,'MaxHeadSize',0.25,'Parent',ax5)
text(pos_src(1)+g_hat_s(1)*sc*1.12, pos_src(2)+g_hat_s(2)*sc*1.12, ...
    pos_src(3)+g_hat_s(3)*sc*1.12, ...
    'ĝ','Color',[1 0.85 0.15],'FontSize',13,'FontWeight','bold','Parent',ax5)

% Source marker
plot3(pos_src(1),pos_src(2),pos_src(3), 'o','MarkerSize',9, ...
    'MarkerFaceColor',[1 1 1],'MarkerEdgeColor','w','Parent',ax5)
text(pos_src(1)+1,pos_src(2)+1,pos_src(3)+2, ...
    'Source','Color','w','FontSize',10,'Parent',ax5)

% Sensor marker
plot3(sLoc(1),sLoc(2),sLoc(3), 'o','MarkerSize',10, ...
    'MarkerFaceColor',[0.3 0.8 1],'MarkerEdgeColor','w','LineWidth',1.2,'Parent',ax5)
text(sLoc(1)+1,sLoc(2),sLoc(3), sprintf('Sensor %d',iSensor), ...
    'Color','w','FontSize',9,'Parent',ax5)

xlim(ax5,xl);  ylim(ax5,yl);  zlim(ax5,zl)
grid(ax5,'on')
set(ax5,'XColor','w','YColor','w','ZColor','w', ...
        'GridColor',[0.25 0.25 0.25],'GridAlpha',0.5)
xlabel(ax5,'x (mm)');  ylabel(ax5,'y (mm)');  zlabel(ax5,'z (mm)')
view(ax5,-40,25)
tg2 = title(ax5, sprintf('Gain vector at vertex %d — Cartesian components', iSrc));
tg2.Color = 'w';

%% Section 4 — Build the tangent-frame leadfield (L1, L2)
% Project Gain from global Cartesian into the nxr per-vertex tangent frame.
% Vectorised: reshape Gain to [nCh x 3 x nV], broadcast dot with e1 and e2.

% nxr tangent frame — same gauge as the stored eigenmodes.
mctx   = nxr.manifold.context(Vtx, Fcs);
vFrame = nxr.manifold.measure.vertexFrame(mctx);
e1_nxr = vFrame.e1;   % [nV x 3]
e2_nxr = vFrame.e2;   % [nV x 3]

G3 = reshape(Gain_meg, nCh, 3, nV);   % [nCh x 3 x nV]  — no data copy
L1 = squeeze(sum(G3 .* reshape(e1_nxr', 1, 3, nV), 2));   % [nCh x nV]
L2 = squeeze(sum(G3 .* reshape(e2_nxr', 1, 3, nV), 2));   % [nCh x nV]
fprintf('L1, L2: %s each\n', mat2str(size(L1)))

%% Section 4 — Build the eigenmode leadfield Lmu
% Psi [nV x K] complex — eigenmodes in nxr gauge (complex single → double).
% Areas [nV x 1]       — vertex areas from the mass matrix diagonal.

Psi_d   = double(TessMat.ConnEigenmodes.Vectors);   % [nV x K] complex
Areas   = full(diag(TessMat.ConnEigenmodes.MassMatrix));
eigVals = TessMat.ConnEigenmodes.Values;
K       = size(Psi_d, 2);

Are = real(Psi_d) .* Areas;   % [nV x K]
Aim = imag(Psi_d) .* Areas;   % [nV x K]

Lmu_R = L1 * Are + L2 * Aim;     % [nCh x K]
Lmu_I = -L1 * Aim + L2 * Are;    % [nCh x K]

% Interleave real and imaginary columns.
Lmu = zeros(nCh, 2*K);
Lmu(:, 1:2:end) = Lmu_R;
Lmu(:, 2:2:end) = Lmu_I;
fprintf('Eigenmode leadfield Lmu: %s\n', mat2str(size(Lmu)))

%% Section 4 — Observability plot
% sigma_k = ||L~_k^R||^2 + ||L~_k^I||^2 is the total power a unit-amplitude
% eigenmode k injects into the sensors. It should decay with eigenvalue
% (higher spatial-frequency modes produce weaker — less observable — signals).
% The R/I ratio tests whether the two tangent components have equal sensor
% visibility; values near 1 indicate the frame is well-conditioned.

sigma_R = sum(Lmu_R.^2, 1);    % [1 x K]
sigma_I = sum(Lmu_I.^2, 1);    % [1 x K]
sigma   = sigma_R + sigma_I;

figure('Name','Eigenmode observability','Color','k','Position',[100 100 720 400])

yyaxis left
scatter(eigVals, sigma, 60, 'filled', ...
    'MarkerFaceColor',[0.30 0.70 1.0],'MarkerEdgeColor','none'); hold on
plot(eigVals, sigma, '-', 'Color',[0.30 0.70 1.0 0.35],'LineWidth',1.2)
ylabel('\sigma_k^2  (total observability)','Color','w')
set(gca,'YColor','w')

yyaxis right
scatter(eigVals, sigma_R ./ max(sigma_I, eps), 35, 'filled', ...
    'MarkerFaceColor',[1.0 0.70 0.25],'MarkerEdgeColor','none')
yline(1,'--','Color',[1.0 0.70 0.25 0.55],'LineWidth',1)
ylabel('\sigma_k^R / \sigma_k^I  (R/I symmetry)','Color',[1.0 0.70 0.25])
set(gca,'YColor',[1.0 0.70 0.25])

set(gca,'Color','k','XColor','w','GridColor',[0.3 0.3 0.3],'GridAlpha',0.3)
grid on
xlabel('Connection Laplacian eigenvalue \mu_k','Color','w')
t4 = title('Eigenmode observability (connection Laplacian tangent basis)');
t4.Color = 'w';
legend({'Total \sigma_k^2','','R/I ratio','R = I'},'TextColor','w', ...
    'Color',[0.10 0.10 0.10],'Location','northeast')

%% Section 5 — Parameterization choices for the dipole at each vertex
% The forward model doesn't know a cortical surface exists. At vertex i the
% Green's function G(r_i) [nCh x 3] maps any 3D dipole q in R^3 to sensor
% data b = G(r_i)*q. The dipole can point in any direction. Three degrees of
% freedom per vertex.
%
% The inverse must recover the unknown dipole moments q_i(t) from sensors.
% It needs to choose a basis for the 3D space of possible dipole directions.
% That choice is a PARAMETERIZATION — it affects the inverse, not the
% forward model.
%
% Choice 1 — Cartesian (unconstrained MNE in Brainstorm/MNE-Python):
%   unknowns are (q_x, q_y, q_z)_i along global X Y Z axes. Equal Cartesian
%   penalty. The parameterization carries no cortical geometry — all three
%   axes are treated identically. 3 unknowns/vertex.
%
% Choice 2 — Local frame + loose prior (existing MNE-Python loose):
%   unknowns are (q_n, q_1, q_2)_i along {n̂(i), e_1(i), e_2(i)}. Normal
%   component is dominant (smaller penalty α_n); tangential components are
%   penalised more: α_t = loose · α_n, default loose=0.2 → 5× more penalty.
%   This is physiologically motivated — cortical pyramidal neurons are
%   oriented along the column normal. The tangent frame {e_1, e_2} is
%   defined arbitrarily (e.g. first-halfedge), so q_1(i) at vertex i bears
%   no consistent relationship to q_1(j) at vertex j. 3 unknowns/vertex.
%
% Choice 3 — Normal only (standard constrained inverse):
%   set q_1=q_2=0, keep only q_n. 1 unknown/vertex. The tangential degrees
%   of freedom are discarded entirely.
%
% Choice 4 — Connection-Laplacian frame + loose prior (our approach):
%   Same penalty structure as Choice 2. The key difference: the tangent frame
%   {e_1(i), e_2(i)} is defined by the Levi-Civita connection (nxr) rather
%   than an arbitrary local convention. This establishes a GLOBALLY
%   CONSISTENT relationship: e_1(j) is the parallel transport of e_1(i)
%   along any path from i to j. The tangential components q_1(i) and q_1(j)
%   are now in the same gauge — they can be meaningfully compared,
%   subtracted, or decomposed in the connection-Laplacian eigenmode basis.
%   3 unknowns/vertex.

n_src  = TessMat.VertNormals(iSrc,:);
e1_src = e1_nxr(iSrc,:);
e2_src = e2_nxr(iSrc,:);

tp_sc5 = 14;
tp5 = [pos_src + tp_sc5*( e1_src + e2_src); ...
       pos_src + tp_sc5*( e1_src - e2_src); ...
       pos_src + tp_sc5*(-e1_src - e2_src); ...
       pos_src + tp_sc5*(-e1_src + e2_src)];

% Arrow lengths encode regularisation weight:
%   normal component — dominant (short penalty) → full length
%   tangential components — penalised ~5x more (loose=0.2) → shorter
sc_n  = sc * 1.0;
sc_t  = sc * 0.55;

panels5 = { ...
    'Choice 1: Cartesian  (unconstrained MNE)', ...
    '3 unknowns/vertex  |  equal Cartesian penalty  |  no geometry'; ...
    'Choice 2: Local frame + loose prior  (existing MNE)', ...
    '3 unknowns/vertex  |  \alpha_n : \alpha_t = 1 : 0.2  |  arbitrary tangent frame'; ...
    'Choice 3: Normal only  (constrained MNE)', ...
    '1 unknown/vertex  |  q_1 = q_2 = 0'; ...
    'Choice 4: Connection-Laplacian frame + loose prior  (our approach)', ...
    '3 unknowns/vertex  |  same penalty  |  gauge-consistent tangent frame'};

clr_n       = [0.92 0.92 0.92];   % normal — bright white
clr_tang    = [0.85 0.65 0.20];   % tangential — existing loose (arbitrary, orange)
clr_tang_cl = [0.35 0.80 0.55];   % tangential — connection-Laplacian (green)
clr_dim5    = [0.28 0.28 0.28];   % discarded
clr_cart    = [0.85 0.65 0.20];   % Cartesian (equal weight)

figure('Name','Parameterization choices','Color','k','Position',[50 50 1200 940])

for p = 1:4
    ax = subplot(2,2,p);
    set(ax,'Color','k'); hold(ax,'on')

    % Tangent plane
    if p == 4
        fc=[0.25 0.45 0.60]; fa=0.30; ec=[0.35 0.55 0.75];
    elseif p == 3
        fc=[0.30 0.30 0.30]; fa=0.08; ec=[0.30 0.30 0.30];
    else
        fc=[0.30 0.30 0.30]; fa=0.18; ec=[0.30 0.30 0.30];
    end
    fill3(tp5([1 2 3 4 1],1),tp5([1 2 3 4 1],2),tp5([1 2 3 4 1],3), ...
        fc,'FaceAlpha',fa,'EdgeColor',ec,'LineWidth',0.8,'Parent',ax)

    if p == 1      % Cartesian: X Y Z, all equal
        vecs5 = {[1 0 0],[0 1 0],[0 0 1]};
        labs5 = {'X̂','Ŷ','Ẑ'};
        clrs5 = {clr_cart, clr_cart, clr_cart};
        scs5  = [sc sc sc];
    elseif p == 2  % Existing loose MNE: n̂ dominant, ê₁ ê₂ penalised, arbitrary frame
        vecs5 = {n_src, e1_src, e2_src};
        labs5 = {'n̂','ê_1','ê_2'};
        clrs5 = {clr_n, clr_tang, clr_tang};
        scs5  = [sc_n sc_t sc_t];
    elseif p == 3  % Constrained: only n̂
        vecs5 = {n_src, e1_src, e2_src};
        labs5 = {'n̂','ê_1','ê_2'};
        clrs5 = {clr_n, clr_dim5, clr_dim5};
        scs5  = [sc_n sc_t sc_t];
    else           % Our approach: n̂ dominant, ê₁ ê₂ gauge-consistent (green)
        vecs5 = {n_src, e1_src, e2_src};
        labs5 = {'n̂','ê_1','ê_2'};
        clrs5 = {clr_n, clr_tang_cl, clr_tang_cl};
        scs5  = [sc_n sc_t sc_t];
    end

    for k = 1:3
        v = vecs5{k}; c = clrs5{k}; s = scs5(k);
        quiver3(pos_src(1),pos_src(2),pos_src(3), v(1)*s,v(2)*s,v(3)*s, 0, ...
            'Color',c,'LineWidth',1.8,'MaxHeadSize',0.4,'Parent',ax)
        text(pos_src(1)+v(1)*s*1.20, pos_src(2)+v(2)*s*1.20, pos_src(3)+v(3)*s*1.20, ...
            labs5{k},'Color',c,'FontSize',10,'FontWeight','bold','Parent',ax)
    end

    % Distinguishing annotation for panels 2 and 4
    if p == 2
        text(pos_src(1),pos_src(2)-tp_sc5*1.1,pos_src(3)+tp_sc5*0.4, ...
            'arbitrary frame','Color',[0.60 0.60 0.60],'FontSize',8, ...
            'HorizontalAlignment','center','FontAngle','italic','Parent',ax)
    elseif p == 4
        text(pos_src(1),pos_src(2)-tp_sc5*1.1,pos_src(3)+tp_sc5*0.4, ...
            'parallel transport','Color',clr_tang_cl,'FontSize',8, ...
            'HorizontalAlignment','center','FontAngle','italic','Parent',ax)
    end

    plot3(pos_src(1),pos_src(2),pos_src(3),'o','MarkerSize',8, ...
        'MarkerFaceColor','w','MarkerEdgeColor','w','Parent',ax)
    view(ax,-55,28); axis(ax,'equal','off')
    tt = title(ax, sprintf('%s\n%s', panels5{p,1}, panels5{p,2}));
    tt.Color = 'w'; tt.FontSize = 7.5;
end

sgtitle({'Source parameterization choices at each cortex vertex', ...
         'Arrow length encodes regularisation weight  (short = more penalised)'}, ...
    'Color','w','FontSize',10)

%% --- helper (keep at bottom) ---
function SurfaceFile = local_find_cortex(nVert)
% Return the first cortex surface with exactly nVert vertices.
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
if isempty(sSubjects), return; end
for iS = 1:numel([sSubjects.Subject])
    for iF = 1:numel(sSubjects.Subject(iS).Surface)
        s = sSubjects.Subject(iS).Surface(iF);
        if ~strcmpi(s.SurfaceType,'Cortex'), continue; end
        try
            T = load(file_fullpath(s.FileName), 'Vertices');
        catch
            continue
        end
        if size(T.Vertices,1) == nVert
            SurfaceFile = s.FileName;
            return
        end
    end
end
end
