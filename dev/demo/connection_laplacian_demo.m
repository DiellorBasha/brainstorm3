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

    hold(ax,'on')
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
