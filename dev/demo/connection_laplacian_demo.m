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
