function F = sa_frames(SurfaceFile)
% SA_FRAMES: Build the three per-vertex tangent frames used by the sign-
% ambiguity continuity test, plus the Fiedler decode and the connection
% Laplacian operator.
%
% USAGE:  F = sa_frames(SurfaceFile)
%
% OUTPUT F (struct):
%   .Vtx, .Fcs, .Nv      : geometry ([nV x 3], [nF x 3] double, [nV x 3])
%   .VertConn            : sparse adjacency [nV x nV]
%   .ConnEig             : connection eigenmodes (incl. .ConnLaplacian = K)
%   .R                   : bst_conn_phase decode (Rank 1 Fiedler), with FsFrame
%   .fiedler             : struct e1,e2 ([nV x 3])  -- Fiedler frame (e1 = field dir)
%   .tang                : struct e1,e2 ([nV x 3])  -- tess_tangents (FS) frame
%   .xyz                 : struct e1,e2 ([nV x 3])  -- global-xyz projected frame
%   .vFrame              : nxr per-vertex frame (e1,e2,normals)
%
% Frames are zero-filled where undefined (e.g. Fiedler near singularities).

TessMat = in_tess_bst(SurfaceFile);
Vtx = TessMat.Vertices;
Fcs = double(TessMat.Faces);
Nv  = TessMat.VertNormals;
nV  = size(Vtx, 1);
if isfield(TessMat, 'VertConn') && ~isempty(TessMat.VertConn)
    VertConn = TessMat.VertConn;
else
    VertConn = tess_vertconn(Vtx, Fcs);
end

% Connection eigenmodes + nxr gauge frame.
ConnEig = bst_conn_eigenmodes_ensure(SurfaceFile);
mctx    = nxr.manifold.context(Vtx, Fcs);
vFrame  = nxr.manifold.measure.vertexFrame(mctx);

% tess_tangents (FreeSurfer trivial-connection) frame at vertices.
[Uf, ~]  = tess_tangents(SurfaceFile, 'NoSave', 1);
[Uv, Vv] = bst_tangent_face2vertex(Fcs, Uf, Nv);
FsFrame  = struct('e1', Uv, 'e2', Vv);

% Fiedler decode (gauge-independent 3D field) + FS-gauge phase.
R = bst_conn_phase(ConnEig, vFrame, 'Rank', 1, 'FsFrame', FsFrame, 'nSing', 2);

% Fiedler frame: e1 = normalized field direction, e2 = n x e1 (zero where |field|~0).
fld = R.Field;
mag = sqrt(sum(fld.^2, 2));
def = mag > 1e-9;
fe1 = zeros(nV, 3);
fe1(def, :) = fld(def, :) ./ mag(def);
fe2 = cross(Nv, fe1, 2);
fiedler = struct('e1', fe1, 'e2', fe2);

% Global-xyz frame: world [1 0 0] projected to each tangent plane (fallback
% [0 1 0] where n is ~parallel to x), e2 = n x e1.
ref = repmat([1 0 0], nV, 1);
bad = abs(sum(Nv .* ref, 2)) > 0.95;        % n nearly parallel to x
ref(bad, :) = repmat([0 1 0], sum(bad), 1);
xe1 = ref - (sum(ref .* Nv, 2)) .* Nv;       % project to tangent plane
xn  = sqrt(sum(xe1.^2, 2));
xok = xn > 1e-9;
xe1(xok, :) = xe1(xok, :) ./ xn(xok);
xe1(~xok, :) = 0;
xe2 = cross(Nv, xe1, 2);
xe2(~xok, :) = 0;
xyz = struct('e1', xe1, 'e2', xe2);

F = struct('Vtx', Vtx, 'Fcs', Fcs, 'Nv', Nv, 'VertConn', VertConn, ...
           'ConnEig', ConnEig, 'R', R, 'fiedler', fiedler, ...
           'tang', FsFrame, 'xyz', xyz, 'vFrame', vFrame);
end
