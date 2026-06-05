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

figure('Name','Cortical surface','Color','w')

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
title('tess\_cortex\_pial\_low  (20484 vertices)')

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
