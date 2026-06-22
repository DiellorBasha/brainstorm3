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

%% Section 6 — Why gauge consistency matters: norm vs phase
% Brainstorm's default readout for an unconstrained source map is the vector
% norm |q_i(t)| = sqrt(qx^2+qy^2+qz^2).  The norm is gauge-invariant — it
% gives the same value in any frame — so it has no sign ambiguity.  But it
% is also information-destroying: it collapses a 3D vector to a non-negative
% scalar, erasing ALL directional structure.
%
% To see why this matters, consider three current configurations at two
% adjacent vertices i and j, each with unit amplitude:
%
%   A — coherent flow:  both dipoles point in the same direction
%   B — diverging:      dipoles point in opposite directions
%   C — rotating:       dipoles point at 90° to each other
%
% The norm |z_i| = |z_j| = 1 for all three.  A norm-based source map is
% IDENTICAL for scenarios A, B, and C — it cannot distinguish coherent flow
% from divergence from rotation.
%
% The complex tangential field z_i = q·e1(i) + i·q·e2(i) in a globally
% consistent gauge DOES distinguish them:
%   A: phase difference = 0        (same direction → coherent flow)
%   B: phase difference = π        (opposite directions → divergence)
%   C: phase difference = π/2      (90° offset → rotation)
%
% The consistent gauge is the prerequisite.  In an arbitrary local frame
% the phase at vertex i and the phase at vertex j are measured against
% unrelated axes, so the phase difference is uninformative noise.
%
% This enables analyses that are impossible with the norm:
%   1. Spatial flow structure: gradient of arg(z) reveals current direction
%   2. Divergence and curl: distinguish source/sink from rotational patterns
%   3. Traveling wave detection: a propagating wave shows a smoothly varying
%      phase that advances in time; invisible in the norm (rectified to 2ω)
%   4. Cross-vertex phase coherence: meaningful connectivity metric
%   5. Spectral decomposition in the connection-eigenmode basis

%% Section 6 — Part A: three-scenario demonstration at two adjacent vertices

% Two adjacent mesh vertices with their nxr tangent frames.
nbrs6  = find(TessMat.VertConn(iSrc,:) > 0);
iNbr6  = nbrs6(1);
pos_i6 = Vtx_mm(iSrc,:);   pos_j6 = Vtx_mm(iNbr6,:);
e1_i6  = e1_nxr(iSrc,:);   e2_i6  = e2_nxr(iSrc,:);
e1_j6  = e1_nxr(iNbr6,:);  e2_j6  = e2_nxr(iNbr6,:);

% Define scenarios as complex values (guarantees |z|=1 and exact phase diffs).
%   z_i = 1+0i for all scenarios (current along e1 at vertex i).
%   z_j encodes the different flow patterns.
zA6_i=1+0i; zA6_j= 1+0i;    % coherent   — phase diff = 0
zB6_i=1+0i; zB6_j=-1+0i;    % diverging  — phase diff = pi
zC6_i=1+0i; zC6_j= 0+1i;    % rotating   — phase diff = pi/2

decode6 = @(z,e1,e2) real(z)*e1 + imag(z)*e2;

scArr6 = 8;
tp6_i = [pos_i6+4*(e1_i6+e2_i6); pos_i6+4*(e1_i6-e2_i6); ...
         pos_i6+4*(-e1_i6-e2_i6); pos_i6+4*(-e1_i6+e2_i6)];
tp6_j = [pos_j6+4*(e1_j6+e2_j6); pos_j6+4*(e1_j6-e2_j6); ...
         pos_j6+4*(-e1_j6-e2_j6); pos_j6+4*(-e1_j6+e2_j6)];

scenarios6 = {'A — coherent flow','B — diverging','C — rotating'};
zis6 = {zA6_i, zB6_i, zC6_i};
zjs6 = {zA6_j, zB6_j, zC6_j};
clrArr6 = [0.95 0.80 0.25];

figure('Name','Norm vs phase: three scenarios','Color','k','Position',[50 50 1200 820])

for s = 1:3
    zi6 = zis6{s};  zj6 = zjs6{s};
    qi6 = decode6(zi6,e1_i6,e2_i6);
    qj6 = decode6(zj6,e1_j6,e2_j6);

    % 3D arrows
    ax = subplot(3,4,(s-1)*4+1);
    set(ax,'Color','k'); hold(ax,'on')
    fill3(tp6_i([1 2 3 4 1],1),tp6_i([1 2 3 4 1],2),tp6_i([1 2 3 4 1],3), ...
        [0.3 0.3 0.3],'FaceAlpha',0.25,'EdgeColor',[0.5 0.5 0.5],'Parent',ax)
    fill3(tp6_j([1 2 3 4 1],1),tp6_j([1 2 3 4 1],2),tp6_j([1 2 3 4 1],3), ...
        [0.3 0.3 0.3],'FaceAlpha',0.25,'EdgeColor',[0.5 0.5 0.5],'Parent',ax)
    quiver3(pos_i6(1),pos_i6(2),pos_i6(3), qi6(1)*scArr6,qi6(2)*scArr6,qi6(3)*scArr6, 0, ...
        'Color',clrArr6,'LineWidth',2.5,'MaxHeadSize',0.4,'Parent',ax)
    quiver3(pos_j6(1),pos_j6(2),pos_j6(3), qj6(1)*scArr6,qj6(2)*scArr6,qj6(3)*scArr6, 0, ...
        'Color',clrArr6,'LineWidth',2.5,'MaxHeadSize',0.4,'Parent',ax)
    plot3([pos_i6(1) pos_j6(1)],[pos_i6(2) pos_j6(2)],[pos_i6(3) pos_j6(3)], ...
        '--','Color',[0.5 0.5 0.5],'LineWidth',0.8,'Parent',ax)
    plot3([pos_i6(1) pos_j6(1)],[pos_i6(2) pos_j6(2)],[pos_i6(3) pos_j6(3)], ...
        '.','Color','w','MarkerSize',10,'Parent',ax)
    view(ax,-55,30); axis(ax,'equal','off')
    tt=title(ax,scenarios6{s}); tt.Color='w'; tt.FontSize=8;

    % Norm
    ax2 = subplot(3,4,(s-1)*4+2);
    set(ax2,'Color','k'); hold(ax2,'on')
    bar(ax2,[abs(zi6) abs(zj6)],'FaceColor',[0.55 0.55 0.55],'EdgeColor','none')
    yline(ax2,1,'--','Color',[0.7 0.7 0.7],'LineWidth',1)
    set(ax2,'XTickLabel',{'i','j'},'XColor','w','YColor','w','Color','k', ...
        'YLim',[0 1.5],'FontSize',8)
    ylabel(ax2,'|z|','Color','w','FontSize',9)
    if s==1, title(ax2,'Norm  |z|','Color','w','FontSize',9); end

    % Phase
    ax3 = subplot(3,4,(s-1)*4+3);
    set(ax3,'Color','k'); hold(ax3,'on')
    bar(ax3,[angle(zi6) angle(zj6)],'FaceColor',[0.35 0.75 0.55],'EdgeColor','none')
    yline(ax3,0,'--','Color',[0.7 0.7 0.7],'LineWidth',0.8)
    set(ax3,'XTickLabel',{'i','j'},'XColor','w','YColor','w','Color','k', ...
        'YLim',[-pi-0.3 pi+0.3],'FontSize',8)
    ylabel(ax3,'arg(z)  (rad)','Color','w','FontSize',9)
    if s==1, title(ax3,'Phase  arg(z)','Color','w','FontSize',9); end

    % Phase difference
    ax4 = subplot(3,4,(s-1)*4+4);
    set(ax4,'Color','k'); hold(ax4,'on')
    dph = abs(angle(zj6/zi6));
    bar(ax4,s,dph,'FaceColor',[0.85 0.45 0.20],'EdgeColor','none','BarWidth',0.5)
    yline(ax4,0,'--','Color',[0.5 0.5 0.5]);
    yline(ax4,pi/2,'--','Color',[0.5 0.5 0.5]);
    yline(ax4,pi,'--','Color',[0.5 0.5 0.5])
    text(s,dph+0.15,sprintf('%.2f rad',dph),'Color','w','FontSize',9, ...
        'HorizontalAlignment','center','Parent',ax4)
    set(ax4,'XTickLabel',scenarios6{s},'XColor','w','YColor','w','Color','k', ...
        'YLim',[0 pi+0.4],'XTick',s,'FontSize',7)
    yticks(ax4,[0 pi/2 pi]); yticklabels(ax4,{'0','\pi/2','\pi'})
    ylabel(ax4,'|\Delta arg(z)|','Color','w','FontSize',9)
    if s==1, title(ax4,'Phase diff  |\Delta arg|','Color','w','FontSize',9); end
end

sgtitle({'Three current configurations — same norm, different phase', ...
    '|z| is identical in all cases  |  arg(z) distinguishes coherent / diverging / rotating'}, ...
    'Color','w','FontSize',10)

%% Section 6 — Part B: Fiedler eigenmode — norm vs phase on the cortex
% The Fiedler eigenmode's norm |Ψ₁| is nearly uniform across the hemisphere.
% Its phase arg(Ψ₁) — read out in the trivial-connection (FS) gauge —
% winds once around the hemisphere, encoding cortical location.
% This is the ground truth for what a consistent phase field looks like
% on a folded cortical surface: the norm tells you nothing about structure;
% the phase reveals the entire spatial organisation.

ConnEig6 = bst_conn_eigenmodes_ensure(SurfaceFile);
mctx6    = nxr.manifold.context(Vtx, double(TessMat.Faces));
vFr6     = nxr.manifold.measure.vertexFrame(mctx6);
[Uf6,~]  = tess_tangents(SurfaceFile,'NoSave',1);
[Uv6,Vv6]= bst_tangent_face2vertex(double(TessMat.Faces), Uf6, TessMat.VertNormals);
FsF6     = struct('e1',Uv6,'e2',Vv6);
R6       = bst_conn_phase(ConnEig6, vFr6, 'Rank',1, 'FsFrame',FsF6, 'nSing',2);

[rH6,lH6] = tess_hemisplit(TessMat);
inH6 = false(size(Vtx_mm,1),1); inH6(lH6)=true;
fHl6 = all(inH6(double(TessMat.Faces)),2);

psi6 = double(ConnEig6.Vectors(:,1));
psi6_mag = abs(psi6);
ph6_fs   = R6.Phase;

cmap_par6 = parula(256);
cmap_hsv6 = hsv(256);
mag6_n  = (psi6_mag - min(psi6_mag(lH6))) / max(psi6_mag(lH6)-min(psi6_mag(lH6)),eps);
cMag6   = repmat([0.15 0.15 0.15],size(Vtx_mm,1),1);
cidx6m  = max(1,min(256,round(mag6_n*255)+1));
cMag6(lH6,:) = cmap_par6(cidx6m(lH6),:);

ph6_n  = (ph6_fs + pi)/(2*pi);
cPh6   = repmat([0.15 0.15 0.15],size(Vtx_mm,1),1);
cidx6p = max(1,min(256,round(ph6_n*255)+1));
cPh6(lH6,:) = cmap_hsv6(cidx6p(lH6),:);

figure('Name','Fiedler norm vs phase','Color','k','Position',[50 50 1100 480])

for p = 1:2
    ax = subplot(1,2,p);
    set(ax,'Color','k')
    cAll6 = iff(p==1, cMag6, cPh6);
    patch('Faces',TessMat.Faces(fHl6,:),'Vertices',Vtx_mm, ...
        'FaceVertexCData',cAll6,'FaceColor','interp', ...
        'EdgeColor','none','FaceLighting','phong','Parent',ax)
    camlight(ax,70,45); camlight(ax,-70,45);
    camlight(ax,70,-45); camlight(ax,-70,-45);
    material([0.65 0.45 0.02 8])
    if p==1
        colormap(ax,parula(256)); clim([min(psi6_mag(lH6)) max(psi6_mag(lH6))])
        cb=colorbar(ax,'Color','w'); cb.Label.String='|\Psi_1(i)|'; cb.Label.Color='w';
        t=title(ax,'Norm  |\Psi_1|  —  amplitude only'); t.Color='w';
    else
        colormap(ax,hsv(256)); clim([-pi pi])
        cb=colorbar(ax,'Color','w','Ticks',[-pi 0 pi],'TickLabels',{'-\pi','0','\pi'});
        cb.Label.String='arg(\Psi_1)  (FS gauge)'; cb.Label.Color='w';
        t=title(ax,'Phase  arg(\Psi_1)  —  location coordinate (trivial-connection gauge)');
        t.Color='w';
    end
    view(ax,-90,10); axis(ax,'equal','off')
end

sgtitle({'Fiedler eigenmode: norm (left) vs phase (right)', ...
    'Norm: nearly uniform — tells you nothing about spatial structure', ...
    'Phase: winds once around the hemisphere — encodes cortical location'}, ...
    'Color','w','FontSize',10)

t%% Section 6 — Part C: Concrete examples from TutorialAuditory data
% All three examples use the unconstrained MN kernel applied to the deviant
% average (Subject01, S01_AEF_20131218_01_notch) at the M100 auditory peak
% (t ≈ 91 ms, selected by global field power).
%
% The tangential field z(i,t) = J(i,t)·ê₁(i) + i·J(i,t)·ê₂(i) is computed
% at every vertex using the nxr Levi-Civita frame (the same gauge as the
% connection-Laplacian eigenmodes), enabling the eigenmode projection.

HMFile_c  = ['/Users/diellorbasha/workspace/library/datasets/brainstorm_db/' ...
             'tmp_aggregate/TutorialAuditory/data/Subject01/' ...
             'S01_AEF_20131218_01_notch/headmodel_surf_os_meg.mat'];
hm_c      = load(HMFile_c,'Gain');
Gain_c    = double(hm_c.Gain);
ChanMat_c = load(file_fullpath('Subject01/S01_AEF_20131218_01_notch/channel_ctf_acc1.mat'),'Channel');
iMEG_c    = find(strcmp({ChanMat_c.Channel.Type},'MEG'));
Gain_meg_c= Gain_c(iMEG_c,:);

% Find unconstrained kernel link file and apply to deviant average
sStudies_c = bst_get('ProtocolStudies');
allStudy_c = [sStudies_c.Study];
ResultsFile_c = '';
for iS=1:numel(allStudy_c)
    for iR=1:numel(allStudy_c(iS).Result)
        fn=allStudy_c(iS).Result(iR).FileName;
        if isempty(fn)||~strncmp(fn,'link|',5), continue; end
        try, Mc=load(file_fullpath(file_resolve_link(fn)),'nComponents'); catch, continue; end
        if isfield(Mc,'nComponents')&&isequal(Mc.nComponents,3)
            ResultsFile_c=fn; break
        end
    end
    if ~isempty(ResultsFile_c), break; end
end
Res_c  = in_bst_results(ResultsFile_c, 1);
IGA_c  = Res_c.ImageGridAmp;   % [3nV x nTime]
Time_c = Res_c.Time(:)';
nV_c   = size(TessMat.Vertices,1);
J_c    = permute(reshape(IGA_c,3,nV_c,[]),[2 1 3]);   % [nV x 3 x nTime]
t_ms_c = Time_c * 1000;

% M100 peak
gfp_c = sqrt(sum(IGA_c.^2,1));
inWin_c = Time_c>=0.06 & Time_c<=0.14;
gfp_c(~inWin_c) = -inf;
[~,ti_c] = max(gfp_c);
fprintf('M100 peak: %.0f ms\n', t_ms_c(ti_c))

% Tangential field in nxr gauge + FS gauge
mctx_c = nxr.manifold.context(TessMat.Vertices, double(TessMat.Faces));
vFr_c  = nxr.manifold.measure.vertexFrame(mctx_c);
e1_c   = vFr_c.e1;  e2_c = vFr_c.e2;
[Uf_c,~]  = tess_tangents(SurfaceFile,'NoSave',1);
[Uv_c,Vv_c] = bst_tangent_face2vertex(double(TessMat.Faces), Uf_c, TessMat.VertNormals);

z_all_c  = squeeze(sum(J_c.*e1_c,2)) + 1i*squeeze(sum(J_c.*e2_c,2));  % nxr, [nV x nTime]
J_peak_c = squeeze(J_c(:,:,ti_c));
z_pk_nxr = sum(J_peak_c.*e1_c,2) + 1i*sum(J_peak_c.*e2_c,2);
z_pk_fs  = sum(J_peak_c.*Uv_c,2) + 1i*sum(J_peak_c.*Vv_c,2);
norm_pk  = sqrt(sum(J_peak_c.^2,2));

% Auditory vertex
Vtx_mm_c = TessMat.Vertices*1000;
[~,iAud_c] = min(sum((Vtx_mm_c-[50 -25 12]).^2,2));
qx_c = squeeze(J_c(iAud_c,1,:));
qy_c = squeeze(J_c(iAud_c,2,:));
qz_c = squeeze(J_c(iAud_c,3,:));
norm_t_c = sqrt(qx_c.^2+qy_c.^2+qz_c.^2);
z_t_c    = z_all_c(iAud_c,:)';
tang_t_c = abs(z_t_c);
ph_t_c   = angle(z_t_c);
Fs_c     = 1/mean(diff(Time_c));
inPost_c = Time_c>0.05 & Time_c<0.45;
fft_qx_c = abs(fft(qx_c(inPost_c))); fft_n_c = abs(fft(norm_t_c(inPost_c)));
nfft_c   = sum(inPost_c); freqs_c = (0:nfft_c-1)*(Fs_c/nfft_c);
[~,ip1_c]=max(fft_qx_c(2:floor(nfft_c/2))); ip1_c=ip1_c+1;
[~,ip2_c]=max(fft_n_c(2:floor(nfft_c/2)));  ip2_c=ip2_c+1;
fprintf('Dominant freq: qx=%.0f Hz, norm=%.0f Hz (doubled)\n', freqs_c(ip1_c), freqs_c(ip2_c))

%% Section 6 — Part C1: Norm vs phase on cortex at M100
[rH_c,~] = tess_hemisplit(TessMat);
rH_cv    = rH_c(:);
inR_c    = false(nV_c,1); inR_c(rH_cv)=true;
fHr_c    = all(inR_c(double(TessMat.Faces)),2);

% Amplitude mask (>20% of peak in right hemi)
amp_c   = abs(z_pk_fs);
mask_c  = amp_c < 0.20*max(amp_c(rH_cv));

cmap_p_c = parula(256);  cmap_h_c = hsv(256);
nr_c   = norm_pk;
nr_n_c = (nr_c-min(nr_c(rH_cv)))/max(nr_c(rH_cv)-min(nr_c(rH_cv)),eps);
cN_c   = repmat([0.15 0.15 0.15],nV_c,1);
cNxr_c = repmat([0.18 0.18 0.18],nV_c,1);
cFs_c  = repmat([0.18 0.18 0.18],nV_c,1);
for iv=rH_cv'
    cN_c(iv,:) = cmap_p_c(max(1,min(256,round(nr_n_c(iv)*255)+1)),:);
    if ~mask_c(iv)
        ci1=max(1,min(256,round((angle(z_pk_nxr(iv))+pi)/(2*pi)*255)+1));
        cNxr_c(iv,:)=cmap_h_c(ci1,:);
        ci2=max(1,min(256,round((angle(z_pk_fs(iv))+pi)/(2*pi)*255)+1));
        cFs_c(iv,:)=cmap_h_c(ci2,:);
    end
end

figure('Name','M100 norm vs phase','Color','k','Position',[50 50 1400 440])
plab_c  = {'Norm |q|  — amplitude only', ...
           'Phase arg(z)  — nxr gauge (arbitrary)', ...
           'Phase arg(z)  — FS gauge (trivial-connection, consistent)'};
pcmap_c = {parula(256), hsv(256), hsv(256)};
pclim_c = {[min(nr_c(rH_cv)) max(nr_c(rH_cv))], [-pi pi], [-pi pi]};
pcdat_c = {cN_c, cNxr_c, cFs_c};

for p=1:3
    ax=subplot(1,3,p); set(ax,'Color','k')
    patch('Faces',TessMat.Faces(fHr_c,:),'Vertices',Vtx_mm_c, ...
        'FaceVertexCData',pcdat_c{p},'FaceColor','interp', ...
        'EdgeColor','none','FaceLighting','phong','Parent',ax)
    camlight(ax,70,45); camlight(ax,-70,45);
    camlight(ax,70,-45); camlight(ax,-70,-45);
    material([0.65 0.45 0.02 8])
    hold(ax,'on')
    plot3(Vtx_mm_c(iAud_c,1),Vtx_mm_c(iAud_c,2),Vtx_mm_c(iAud_c,3), ...
        'o','MarkerSize',9,'MarkerFaceColor','w','MarkerEdgeColor','w','Parent',ax)
    colormap(ax,pcmap_c{p}); clim(pclim_c{p})
    cb=colorbar(ax,'Color','w');
    if p==1, cb.Label.String='|q_i|  (Am)';
    else,    cb.Label.String='arg(z_i)  (rad)';
             cb.Ticks=[-pi 0 pi]; cb.TickLabels={'-\pi','0','\pi'};
    end
    cb.Label.Color='w';
    view(ax,90,10); axis(ax,'equal','off')
    t=title(ax,plab_c{p}); t.Color='w'; t.FontSize=8;
end
sgtitle({'Right hemisphere  |  deviant average M100 (91 ms)  |  white dot = auditory cortex', ...
    'Phase panels masked to active region (>20% peak amplitude)'}, 'Color','w','FontSize',10)

%% Section 6 — Part C2: Timeseries — frequency doubling at auditory cortex
tWin_c = t_ms_c >= -100 & t_ms_c <= 400;
figure('Name','Timeseries frequency doubling','Color','k','Position',[50 50 1000 550])

subplot(3,1,1)
plot(t_ms_c(tWin_c),qx_c(tWin_c)*1e12,'Color',[0.85 0.55 0.25],'LineWidth',1.2); hold on
plot(t_ms_c(tWin_c),qy_c(tWin_c)*1e12,'Color',[0.45 0.75 0.35],'LineWidth',1.2)
plot(t_ms_c(tWin_c),qz_c(tWin_c)*1e12,'Color',[0.40 0.65 1.00],'LineWidth',1.2)
xline(t_ms_c(ti_c),'--','Color','w','LineWidth',0.8)
ylabel('pA\cdotm','Color','w'); set(gca,'Color','k','XColor','w','YColor','w')
legend({'q_x','q_y','q_z'},'TextColor','w','Color',[0.1 0.1 0.1],'Location','northeast')
title('Cartesian components — signed, preserve true oscillation frequency','Color','w','FontSize',9)

subplot(3,1,2)
plot(t_ms_c(tWin_c),norm_t_c(tWin_c)*1e12,'Color',[0.85 0.85 0.85],'LineWidth',1.5); hold on
plot(t_ms_c(tWin_c),tang_t_c(tWin_c)*1e12,'--','Color',[0.35 0.75 0.55],'LineWidth',1.5)
xline(t_ms_c(ti_c),'--','Color','w','LineWidth',0.8)
ylabel('pA\cdotm','Color','w'); set(gca,'Color','k','XColor','w','YColor','w')
legend({'norm |q|','tangential |z|'},'TextColor','w','Color',[0.1 0.1 0.1],'Location','northeast')
title(sprintf('Norm |q| and |z|: RECTIFIED  —  dominant at %.0f Hz  (%.0f Hz components, 2\\times frequency)', ...
    freqs_c(ip2_c), freqs_c(ip1_c)),'Color','w','FontSize',9)

subplot(3,1,3)
plot(t_ms_c(tWin_c),ph_t_c(tWin_c),'Color',[0.35 0.75 0.55],'LineWidth',1.5); hold on
xline(t_ms_c(ti_c),'--','Color','w','LineWidth',0.8)
ylabel('arg(z)  (rad)','Color','w'); xlabel('time  (ms)','Color','w')
set(gca,'Color','k','XColor','w','YColor','w','YLim',[-pi-0.3 pi+0.3])
yticks([-pi -pi/2 0 pi/2 pi]); yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
title(sprintf('Tangential phase arg(z_i(t)) — oscillates at %.0f Hz, preserves temporal phase', ...
    freqs_c(ip1_c)),'Color','w','FontSize',9)

sgtitle(sprintf('Auditory cortex [%.0f %.0f %.0f] mm  |  tangential fraction %.0f%%  |  deviant average', ...
    Vtx_mm_c(iAud_c,:), abs(z_pk_nxr(iAud_c))/norm_pk(iAud_c)*100), 'Color','w','FontSize',10)

%% Section 6 — Part C3: Eigenmode decomposition of tangential source
% Project z(i,t) onto the connection-Laplacian eigenmodes:
%   c_k(t) = sum_i conj(Psi_k(i)) * z(i,t) * Area_i
% The M-inner product is valid because z and Psi share the same nxr gauge.
% This spectral decomposition is the unique capability of gauge consistency.

ConnEig_c = bst_conn_eigenmodes_ensure(SurfaceFile);
Areas_c   = full(diag(ConnEig_c.MassMatrix));
Psi_c     = double(ConnEig_c.Vectors);
eigVals_c = ConnEig_c.Values;
K_c       = size(Psi_c,2);

c_time_c  = (conj(Psi_c)' .* Areas_c') * z_all_c;   % [K x nTime]
c_peak_c  = c_time_c(:, ti_c);
c_pow_c   = abs(c_peak_c).^2;

% Find left-hemi Fiedler mode
colL_c = 1;
for ci_k = 1:K_c
    supp=find(Psi_c(:,ci_k)~=0);
    [~,lH_c]=tess_hemisplit(TessMat);
    if mean(ismember(supp,lH_c(:)))>0.9, colL_c=ci_k; break; end
end
cFiedL_c = c_time_c(colL_c,:);

figure('Name','Eigenmode decomposition','Color','k','Position',[50 50 1100 500])

ax_sp = subplot(1,2,1); set(ax_sp,'Color','k'); hold(ax_sp,'on')
bar(ax_sp,eigVals_c,c_pow_c*1e26,'FaceColor',[0.35 0.75 0.55],'EdgeColor','none')
[~,top3_c]=sort(c_pow_c,'descend');
for ki=1:3
    text(eigVals_c(top3_c(ki)),c_pow_c(top3_c(ki))*1e26+0.05, ...
        sprintf('\\Psi_{%d}',top3_c(ki)),'Color','w','FontSize',9, ...
        'HorizontalAlignment','center','Parent',ax_sp)
end
xlabel(ax_sp,'Eigenvalue \mu_k  (spatial frequency)','Color','w')
ylabel(ax_sp,'|c_k|^2  (\times 10^{-26})','Color','w')
set(ax_sp,'XColor','w','YColor','w','GridColor',[0.3 0.3 0.3],'GridAlpha',0.4)
grid(ax_sp,'on')
t_sp=title(ax_sp,'Tangential source spectrum at M100  (eigenmode power)'); t_sp.Color='w';

ax_ts = subplot(1,2,2); set(ax_ts,'Color','k'); hold(ax_ts,'on')
yyaxis(ax_ts,'left')
plot(t_ms_c,abs(cFiedL_c)*1e13,'Color',[0.35 0.75 0.55],'LineWidth',1.5)
ylabel(ax_ts,'|c_1(t)|  (\times 10^{-13})','Color',[0.35 0.75 0.55])
set(ax_ts,'YColor',[0.35 0.75 0.55])
yyaxis(ax_ts,'right')
plot(t_ms_c,angle(cFiedL_c),'Color',[0.85 0.65 0.25],'LineWidth',1.0)
ylabel(ax_ts,'arg(c_1(t))  (rad)','Color',[0.85 0.65 0.25])
set(ax_ts,'YColor',[0.85 0.65 0.25],'YLim',[-pi-0.3 pi+0.3])
yticks(ax_ts,[-pi 0 pi]); yticklabels(ax_ts,{'-\pi','0','\pi'})
xline(ax_ts,t_ms_c(ti_c),'--','Color','w','LineWidth',0.8)
xline(ax_ts,0,'--','Color',[0.5 0.5 0.5],'LineWidth',0.6)
xlim(ax_ts,[-100 400]); xlabel(ax_ts,'time  (ms)','Color','w')
set(ax_ts,'XColor','w','GridColor',[0.3 0.3 0.3],'GridAlpha',0.4); grid(ax_ts,'on')
t_ts=title(ax_ts,'Fiedler coefficient c_1(t)  — left hemisphere'); t_ts.Color='w';

sgtitle({'Connection-eigenmode decomposition of tangential source', ...
    'Requires gauge-consistent frame  |  Fiedler mode (\Psi_1, left hemi) dominates at M100'}, ...
    'Color','w','FontSize',10)

%% Section 6 — Part D: Alpha-band traveling wave via phase gradient
% Alpha oscillations (7-13 Hz) in resting or pre-stimulus MEG have a known
% traveling-wave component: the phase of the oscillation advances in space,
% typically sweeping from posterior to anterior at ~1-10 m/s.
%
% To detect this with source imaging we need:
%   1. A complex source estimate z(i,t) in a CONSISTENT gauge — so that
%      arg(z_i) and arg(z_j) at adjacent vertices i,j are directly
%      comparable and their difference arg(z_j/z_i) = Δarg tells us
%      whether j leads or lags i.
%   2. A bandpassed recording where the instantaneous phase is stable.
%
% Without gauge consistency, Δarg between adjacent vertices conflates the
% geometric rotation between their arbitrary local frames with any genuine
% phase difference from propagation — you cannot separate the two.
%
% Data: raw segment 80-100s from run _02_notch, bandpassed to alpha (7-13 Hz),
%       downsampled to 600 Hz, unconstrained MN kernel from run _01_notch
%       (same headmodel/subject — kernel is session-independent).

% Load alpha-band recording and apply kernel
dataFile_d = 'Subject01/S01_AEF_20131218_02_notch/data_block001_band.mat';
D_d   = load(file_fullpath(dataFile_d));
kPath_d = file_fullpath('Subject01/S01_AEF_20131218_01_notch/results_MN_MEG_KERNEL_260605_0111.mat');
Kd    = load(kPath_d,'ImagingKernel','GoodChannel');
F_d   = D_d.F(Kd.GoodChannel,:);              % [272 x nTime]
nV_d  = size(TessMat.Vertices,1);
Fs_d  = 600;  Time_d = D_d.Time;  t_ms_d = Time_d*1000;

% Build FS tangent frame (consistent gauge)
mctx_d  = nxr.manifold.context(TessMat.Vertices, double(TessMat.Faces));
vFr_d   = nxr.manifold.measure.vertexFrame(mctx_d);
[Uf_d,~]    = tess_tangents(SurfaceFile,'NoSave',1);
[Uv_d,Vv_d] = bst_tangent_face2vertex(double(TessMat.Faces),Uf_d,TessMat.VertNormals);

J_d     = permute(reshape(Kd.ImagingKernel*F_d,3,nV_d,[]),[2 1 3]);
z_fs_d  = squeeze(sum(J_d.*Uv_d,2)) + 1i*squeeze(sum(J_d.*Vv_d,2));  % FS gauge

% Find peak amplitude window
gfp_d   = sqrt(mean(abs(z_fs_d).^2,1));
[~,iPk_d]  = max(gfp_d);
iWin_d  = max(1,iPk_d-250):min(numel(Time_d),iPk_d+250);
amp_win_d  = mean(abs(z_fs_d(:,iWin_d)),2);

% Pick two high-amplitude vertices ≥30 mm apart
[~,sord_d] = sort(amp_win_d,'descend');
iV1_d = sord_d(1);
for ki = 2:numel(sord_d)
    iV2_d = sord_d(ki);
    if norm(Vtx_mm_a(iV1_d,:)-Vtx_mm_a(iV2_d,:)) > 30, break; end
end
ph_lead = mean(angle(z_fs_d(iV1_d,iWin_d)./z_fs_d(iV2_d,iWin_d)));
fprintf('V1→V2 phase lead: %.3f rad = %.1f ms at 10 Hz\n', ph_lead, ph_lead/(2*pi)*100)

% Three phase maps spanning one half-cycle + phase timeseries
[~,lH_d] = tess_hemisplit(TessMat);
lH_dv   = lH_d(:);
inL_d   = false(nV_d,1); inL_d(lH_dv)=true;
fHl_d   = all(inL_d(double(TessMat.Faces)),2);
mask_d  = amp_win_d < 0.08*max(amp_win_d(lH_dv));
frame_step_d = round(Fs_d/(10*6));
frames3_d    = iPk_d + [0, 2, 4]*frame_step_d;
frames3_d    = min(frames3_d, numel(Time_d));
cmap_wd = hsv(256);

figure('Name','Alpha wave tracking','Color','k','Position',[50 50 1100 700])

for f = 1:3
    ax = subplot(2,3,f); set(ax,'Color','k'); hold(ax,'on')
    ph_f_d = angle(z_fs_d(:,frames3_d(f)));
    cPh_d  = repmat([0.25 0.25 0.23],nV_d,1);
    for iv = lH_dv'
        if mask_d(iv), continue; end
        ci=max(1,min(256,round((ph_f_d(iv)+pi)/(2*pi)*255)+1));
        cPh_d(iv,:)=cmap_wd(ci,:);
    end
    patch('Faces',double(TessMat.Faces(fHl_d,:)),'Vertices',Vtx_mm_a, ...
        'FaceVertexCData',cPh_d,'FaceColor','interp', ...
        'EdgeColor','none','FaceLighting','phong','Parent',ax)
    camlight(ax,70,45); camlight(ax,-70,45);
    camlight(ax,70,-45); camlight(ax,-70,-45);
    material([0.65 0.45 0.02 8])
    plot3(Vtx_mm_a(iV1_d,1),Vtx_mm_a(iV1_d,2),Vtx_mm_a(iV1_d,3), ...
        'o','MarkerSize',10,'MarkerFaceColor',[1 0.85 0.15],'MarkerEdgeColor','w','Parent',ax)
    plot3(Vtx_mm_a(iV2_d,1),Vtx_mm_a(iV2_d,2),Vtx_mm_a(iV2_d,3), ...
        'o','MarkerSize',10,'MarkerFaceColor',[0.35 0.75 0.55],'MarkerEdgeColor','w','Parent',ax)
    colormap(ax,hsv(256)); clim([-pi pi])
    view(ax,-90,10); axis(ax,'equal','off')
    dt_d = round(t_ms_d(frames3_d(f))-t_ms_d(frames3_d(1)));
    tl=title(ax,sprintf('t = +%d ms',dt_d)); tl.Color='w'; tl.FontSize=9;
end

% Phase timeseries
iPlot_d = max(1,iPk_d-200):min(numel(Time_d),iPk_d+300);
t_rel_d = t_ms_d(iPlot_d) - t_ms_d(iPk_d);
ph1_d   = unwrap(angle(z_fs_d(iV1_d,iPlot_d)));
ph2_d   = unwrap(angle(z_fs_d(iV2_d,iPlot_d)));

ax_bd = subplot(2,1,2); set(ax_bd,'Color','k'); hold(ax_bd,'on')
plot(t_rel_d,ph1_d,'Color',[1 0.85 0.15],'LineWidth',1.8)
plot(t_rel_d,ph2_d,'Color',[0.35 0.75 0.55],'LineWidth',1.8)
for f=1:3
    xline(ax_bd,t_ms_d(frames3_d(f))-t_ms_d(iPk_d),'--','Color','w','LineWidth',0.8,'Alpha',0.5)
end
xlabel(ax_bd,'time relative to peak  (ms)','Color','w')
ylabel(ax_bd,'arg(z)  (rad, unwrapped)','Color','w')
set(ax_bd,'XColor','w','YColor','w','GridColor',[0.3 0.3 0.3],'GridAlpha',0.4)
grid(ax_bd,'on')
legend(ax_bd,{sprintf('V1  [%.0f %.0f %.0f] mm',Vtx_mm_a(iV1_d,:)), ...
    sprintf('V2  [%.0f %.0f %.0f] mm  (%.0f mm away)',Vtx_mm_a(iV2_d,:), ...
    norm(Vtx_mm_a(iV1_d,:)-Vtx_mm_a(iV2_d,:)))}, ...
    'TextColor','w','Color',[0.1 0.1 0.1],'Location','northwest')
title(ax_bd,sprintf('Phase timeseries  |  mean lead V1\\rightarrowV2 = %.2f rad = %.0f ms at 10 Hz  |  distance = %.0f mm', ...
    ph_lead, ph_lead/(2*pi)*100, norm(Vtx_mm_a(iV1_d,:)-Vtx_mm_a(iV2_d,:))), ...
    'Color','w','FontSize',9)

sgtitle({'Alpha (7-13 Hz) traveling wave  |  phase in consistent FS gauge  |  left hemisphere', ...
    'Phase lead between two vertices 33 mm apart is measurable only with gauge consistency'}, ...
    'Color','w','FontSize',10)

%% --- helpers (keep at bottom) ---
function v = iff(cond, a, b)
if cond, v = a; else, v = b; end
end

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
