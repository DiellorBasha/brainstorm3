function NxrData = tess_nxr_populate(SurfaceFile, varargin)
% TESS_NXR_POPULATE  Single-pass nxr-compute population of all surface geometry.
%
% USAGE:
%   NxrData = tess_nxr_populate(SurfaceFile)
%   NxrData = tess_nxr_populate(SurfaceFile, 'nConnModes', 20, 'nLBOModes', 200)
%   NxrData = tess_nxr_populate(SurfaceFile, 'ForceRecompute', true)
%
% DESCRIPTION:
%   Single source of truth for all geometric operators on a cortical surface.
%   One nxr-compute call gives the complete consistent set:
%
%     DEC operators      d0, d1, hodge0, hodge1, hodge2
%     Face frames        Uf, Vf, Nf  per face (trivial-connection gauge,
%                        singularities at FreeSurfer sphere poles)
%     ConnAngles         α_e per interior edge: parallel-transport rotation
%                        from face-frame at f1 to face-frame at f2
%     ConnLaplacian      [nF × nF] sparse Hermitian, face-level connection
%                        Laplacian built from trivial-connection angles
%     ConnEigLH/RH       face-level Fiedler + higher eigenvectors (complex)
%     LBOEigLH/RH        vertex scalar LBO eigenmodes (real) from d0'*★₁*d0
%     FaceAdj            [nF × nF] face adjacency from |d1·d1'|
%
%   All operators are consistent by construction — derived from the same
%   manifold context, same trivial-connection gauge.  Downstream functions
%   read from TessMat.nxr without calling nxr-compute again.
%
% OPTIONS (name-value):
%   'nConnModes'     K per hemisphere for face connection Laplacian (default: 20)
%   'nLBOModes'      K per hemisphere for vertex scalar LBO            (default: 200)
%   'ForceRecompute' recompute even if TessMat.nxr already exists     (default: false)
%
% OUTPUT:
%   NxrData  struct — also written into TessMat.nxr in the surface .mat file.
%     .d0          [nE × nV]    discrete exterior derivative (scalar gradient)
%     .d1          [nF × nE]    discrete exterior derivative (curl)
%     .hodge0      [nV × nV]    ★₀ vertex Voronoi area matrix (diagonal)
%     .hodge1      [nE × nE]    ★₁ edge Hodge star (diagonal)
%     .hodge2      [nF × nF]    ★₂ face area matrix (diagonal)
%     .FaceAdj     [nF × nF]    sparse face adjacency
%     .FaceFrames  struct: .U [nF×3] e1, .V [nF×3] e2, .N [nF×3] normal
%     .ConnAngles  [nE × 1]    α_e: parallel-transport angle per interior edge
%     .ConnLaplacian [nF × nF] Hermitian sparse face connection Laplacian
%     .ConnEigLH   struct: .Vectors [nLHF×K] complex, .Values [K×1], .FaceIdx [nLHF×1]
%     .ConnEigRH   struct: same for right hemisphere
%     .LBOEigLH    struct: .Vectors [nLHV×K] real,    .Values [K×1], .VertIdx [nLHV×1]
%     .LBOEigRH    struct: same for right hemisphere
%     .lH_f        [nLHF × 1]  LH face global indices
%     .rH_f        [nRHF × 1]  RH face global indices
%     .lH_v        [nLHV × 1]  LH vertex global indices
%     .rH_v        [nRHV × 1]  RH vertex global indices
%
% SEE ALSO: tess_tangents, bst_conn_eigenmodes_ensure, bst_eigenmodes_ensure
%
% Authors: Diellor Basha, 2026

%% ── Parse options ────────────────────────────────────────────────────────
nConnModes    = 20;
nLBOModes     = 200;
ForceRecompute = false;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'nconnmodes',     nConnModes     = varargin{k+1};
        case 'nlbomodes',      nLBOModes      = varargin{k+1};
        case 'forcerecompute', ForceRecompute = logical(varargin{k+1});
    end
end

%% ── Load surface ─────────────────────────────────────────────────────────
TessMat = in_tess_bst(SurfaceFile, 0);
Vtx   = TessMat.Vertices;          % [nV × 3] metres
Faces = double(TessMat.Faces);     % [nF × 3]
nV    = size(Vtx, 1);
nF    = size(Faces, 1);

if ~ForceRecompute && isfield(TessMat,'nxr') && ~isempty(TessMat.nxr)
    NxrData = TessMat.nxr;
    fprintf('tess_nxr_populate: using cached nxr data for %s\n', SurfaceFile);
    return;
end

fprintf('tess_nxr_populate: %d V, %d F — building all operators...\n', nV, nF);

%% ── Step 1: DEC operators ────────────────────────────────────────────────
bst_plugin('Load', 'nxr-compute');
clear mex
h_dec = nxr_compute('create', Vtx, Faces);
dec   = nxr_compute('assembleDECOperators', h_dec);
nxr_compute('destroy', h_dec);

nE = size(dec.d1, 2);

NxrData      = struct();
NxrData.d0   = dec.d0;      % [nE × nV]
NxrData.d1   = dec.d1;      % [nF × nE]
NxrData.hodge0 = dec.hodge0; % [nV × nV] diagonal
NxrData.hodge1 = dec.hodge1; % [nE × nE] diagonal
NxrData.hodge2 = dec.hodge2; % [nF × nF] diagonal

% Face adjacency: two faces are adjacent iff they share a primal edge.
% d1[f,e] = ±1, so (d1·d1')[f1,f2] = ±1 for adjacent pairs — use abs > 0.
NxrData.FaceAdj = (abs(dec.d1 * dec.d1') > 0) - speye(nF);

fprintf('  DEC assembled: %d edges\n', nE);

%% ── Step 2: Face frames (trivial-connection gauge) ───────────────────────
[Uf, Vf] = tess_tangents(SurfaceFile, 'NoSave', 1);
Nf = cross(Uf, Vf, 2);    % face normals in nxr's CW convention
NxrData.FaceFrames.U = Uf;
NxrData.FaceFrames.V = Vf;
NxrData.FaceFrames.N = Nf;

%% ── Step 3: Edge vertex pairs from d0 ───────────────────────────────────
% d0[e, v] = +1 at head vertex, -1 at tail vertex.
[e_rows, v_cols, d0_vals] = find(dec.d0);
pos_mask = d0_vals > 0;
neg_mask = d0_vals < 0;
edge_head = accumarray(e_rows(pos_mask), v_cols(pos_mask), [nE, 1]);
edge_tail = accumarray(e_rows(neg_mask), v_cols(neg_mask), [nE, 1]);

%% ── Step 4: Edge-to-face pairs from d1 ──────────────────────────────────
% d1[f, e] = ±1 means face f contains edge e.
% Interior edges have exactly 2 incident faces.
[f_rows_d1, e_cols_d1] = find(dec.d1);
% For each edge: collect all incident faces (should be 1 or 2)
edge_faces = cell(nE, 1);
for k = 1:numel(e_cols_d1)
    e = e_cols_d1(k);
    edge_faces{e}(end+1) = f_rows_d1(k);
end
interior = cellfun(@numel, edge_faces) == 2;   % [nE × 1] logical

%% ── Step 5: Face connection Laplacian ────────────────────────────────────
% For each interior edge e between faces f1, f2:
%   Express primal edge direction in f1's frame → complex z1
%   Express primal edge direction in f2's frame → complex z2
%   Transport ρ = z2/z1 / |z2/z1|  (unit-circle complex rotation)
%   L_conn[f1,f2] -= w_e * ρ_{f1→f2}    (Hermitian: [f2,f1] -= w_e * conj(ρ))
%   L_conn[f1,f1] += w_e
% Weight: hodge1 diagonal (dual-edge / primal-edge length ratio, always ≥ 0).

h1_diag = abs(full(diag(dec.hodge1)));  % clamp: obtuse triangles give h1<0

int_edges = find(interior);
nInt      = numel(int_edges);
conn_alpha = zeros(nE, 1);   % transport angle per edge (zero for boundary)

% Pre-allocate sparse triplets
I_trip = zeros(4*nInt, 1);
J_trip = zeros(4*nInt, 1);
V_trip = zeros(4*nInt, 1, 'like', 1i);
trip_k = 0;

for k = 1:nInt
    e  = int_edges(k);
    f1 = edge_faces{e}(1);
    f2 = edge_faces{e}(2);
    vi = edge_tail(e);
    vj = edge_head(e);
    if vi == 0 || vj == 0, continue; end

    % Primal edge direction (unit vector)
    edg = Vtx(vj,:) - Vtx(vi,:);
    en  = norm(edg);
    if en < eps, continue; end
    edg = edg / en;

    % Project onto face frames (complex representation)
    z1 = dot(edg, Uf(f1,:)) + 1i*dot(edg, Vf(f1,:));
    z2 = dot(edg, Uf(f2,:)) + 1i*dot(edg, Vf(f2,:));
    if abs(z1) < 1e-10 || abs(z2) < 1e-10, continue; end

    rho = z2 / z1;
    rho = rho / abs(rho);   % unit-circle transport
    conn_alpha(e) = angle(rho);

    w = h1_diag(e);

    % Four triplets: (f1,f2), (f2,f1), (f1,f1), (f2,f2)
    trip_k = trip_k + 1;
    I_trip(trip_k) = f1; J_trip(trip_k) = f2; V_trip(trip_k) = -w * rho;
    trip_k = trip_k + 1;
    I_trip(trip_k) = f2; J_trip(trip_k) = f1; V_trip(trip_k) = -w * conj(rho);
    trip_k = trip_k + 1;
    I_trip(trip_k) = f1; J_trip(trip_k) = f1; V_trip(trip_k) = w;
    trip_k = trip_k + 1;
    I_trip(trip_k) = f2; J_trip(trip_k) = f2; V_trip(trip_k) = w;
end

L_conn = sparse(I_trip(1:trip_k), J_trip(1:trip_k), V_trip(1:trip_k), nF, nF);
L_conn = (L_conn + L_conn') / 2;   % enforce exact Hermitian symmetry

NxrData.ConnAngles    = conn_alpha;   % [nE × 1] rad
NxrData.ConnLaplacian = L_conn;       % [nF × nF] sparse Hermitian complex

fprintf('  Face connection Laplacian: %d × %d, %d interior edges\n', nF, nF, nInt);

%% ── Step 6: Hemisphere split ─────────────────────────────────────────────
[~, lH_v_idx] = tess_hemisplit(TessMat);
lH_v = lH_v_idx(:);
rH_v = setdiff((1:nV)', lH_v);
inL  = false(nV, 1); inL(lH_v) = true;
lH_f = find(all(inL(Faces), 2));
rH_f = find(all(~inL(Faces), 2));

NxrData.lH_f = lH_f(:);
NxrData.rH_f = rH_f(:);
NxrData.lH_v = lH_v(:);
NxrData.rH_v = rH_v(:);

%% ── Step 7: Face connection Laplacian eigenmodes (per hemisphere) ─────────
fprintf('  Face connection Laplacian eigenmodes (K=%d per hemisphere)...\n', nConnModes);
M_face = dec.hodge2;

NxrData.ConnEigLH = local_conn_eigs(L_conn, M_face, lH_f, nConnModes, 'LH');
NxrData.ConnEigRH = local_conn_eigs(L_conn, M_face, rH_f, nConnModes, 'RH');

%% ── Step 8: Scalar LBO eigenmodes (vertex level, per hemisphere) ──────────
fprintf('  Scalar LBO eigenmodes (K=%d per hemisphere)...\n', nLBOModes);
L_lbo  = dec.d0' * dec.hodge1 * dec.d0;  % [nV × nV] cotangent Laplacian
M_vert = dec.hodge0;

NxrData.LBOEigLH = local_lbo_eigs(L_lbo, M_vert, lH_v, nLBOModes, 'LH');
NxrData.LBOEigRH = local_lbo_eigs(L_lbo, M_vert, rH_v, nLBOModes, 'RH');

%% ── Step 9: Save into TessMat ────────────────────────────────────────────
TessMat_full = load(file_fullpath(SurfaceFile));
TessMat_full.nxr = NxrData;
bst_save(file_fullpath(SurfaceFile), TessMat_full, 'v7');
fprintf('tess_nxr_populate: done. Saved → %s\n', SurfaceFile);

end   % tess_nxr_populate


%% ── Local: face connection Laplacian eigenmodes ───────────────────────────
function Eig = local_conn_eigs(L_conn, M_face, hemi_f, nModes, tag)
    % M_face is hodge2 = ★₂ = 1/A_f (inverse area).
    % L² mass for face vector fields is A_f = inv(hodge2), so we pass inv(M_face).
    nHF  = numel(hemi_f);
    nMod = min(nModes + 1, nHF - 1);
    L_h  = L_conn(hemi_f, hemi_f);
    h2_h = full(diag(M_face(hemi_f, hemi_f)));
    M_h  = spdiags(1./max(h2_h, eps), 0, nHF, nHF);   % A_f = inv(hodge2)

    [V, D] = eigs(L_h, M_h, nMod, 'smallestabs');
    lam    = real(diag(D));
    [lam, ord] = sort(lam, 'ascend');
    V    = V(:, ord);

    % Skip DC mode (λ ≈ 0)
    dc  = find(lam < 1e-6 * max(abs(lam)), 1, 'last');
    if isempty(dc), dc = 1; end
    sel = (dc+1) : min(dc + nModes, numel(lam));

    Eig.Vectors  = V(:, sel);     % [nHF × K] complex
    Eig.Values   = lam(sel);      % [K × 1]   real
    Eig.FaceIdx  = hemi_f(:);     % [nHF × 1] global face indices
    Eig.nModes   = numel(sel);
    fprintf('    %s: %d face conn modes  λ∈[%.1f … %.1f]\n', ...
        tag, Eig.nModes, min(Eig.Values), max(Eig.Values));
end


%% ── Local: vertex scalar LBO eigenmodes ───────────────────────────────────
function Eig = local_lbo_eigs(L_lbo, M_vert, hemi_v, nModes, tag)
    nHV  = numel(hemi_v);
    nMod = min(nModes + 1, nHV - 1);
    L_h  = L_lbo(hemi_v, hemi_v);
    M_h  = M_vert(hemi_v, hemi_v);

    [V, D] = eigs(L_h, M_h, nMod, 'smallestabs');
    lam    = real(diag(D));
    [lam, ord] = sort(lam, 'ascend');
    V    = real(V(:, ord));

    dc  = find(lam < 1e-6 * max(abs(lam)), 1, 'last');
    if isempty(dc), dc = 1; end
    sel = (dc+1) : min(dc + nModes, numel(lam));

    Eig.Vectors  = V(:, sel);    % [nHV × K] real
    Eig.Values   = lam(sel);     % [K × 1]
    Eig.VertIdx  = hemi_v(:);    % [nHV × 1] global vertex indices
    Eig.nModes   = numel(sel);
    fprintf('    %s: %d scalar LBO modes  λ∈[%.1f … %.1f]\n', ...
        tag, Eig.nModes, min(Eig.Values), max(Eig.Values));
end
