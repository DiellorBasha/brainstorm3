# M3 Plan A — nxr per-vertex tangent-frame export

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `vertexFrames` command to the `nxr-compute` MATLAB binding that exports geometry-central's per-vertex tangent basis (`e1`/`e2`/`normal`, `[nV×3]`) — the exact gauge in which the vertex connection-Laplacian's complex eigenvectors are expressed — so M3 can decode a stored eigenmode into a 3D tangent field.

**Architecture:** Mirror the existing per-face `frames` command end-to-end: a C++ `vertexFrames(Manifold&)` returning a `VertexFrames` struct from `geometry.vertexTangentBasis`, a MEX dispatch case `'vertexFrames'` + a `vertexFramesToStruct` marshaller, and a `+nxr/+manifold/+measure/vertexFrame.m` wrapper. Then rebuild the macOS plugin and restage it locally so the live Brainstorm session can use it.

**Tech Stack:** C++ (geometry-central, Eigen), the MEX binding, MATLAB `+nxr` wrappers, `scripts/build.sh`.

**Reference spec:** `dev/connection_phase_readout_integration.md` (§3.1).

**Repo:** this plan operates in `/Users/diellorbasha/workspace/research/code/nxr-compute` (a *separate* git repo from brainstorm3), except Task 2 which copies artifacts into the local Brainstorm plugin install.

**Gauge-correctness note (why `vertexTangentBasis` is exactly right):** the vertex connection Laplacian builds its transport from `halfedgeVectorsInVertex` (first-outgoing-halfedge = angle 0). `vertexTangentBasis` is geometry-central's 3D realization of that same space (`basisX` aligned to the angle-0 axis, `basisY = N × basisX`). So a connection-mode coordinate `z = (a,b)` at vertex `i` decodes to `a·e1_i + b·e2_i`. Geometric validity is tested here in Plan A; the *gauge*/semantic correctness (decoded Fiedler vs `smoothVertex`) is verified end-to-end in Plan B.

---

## File Structure

| Path (in `nxr-compute`) | Change |
|---|---|
| `include/nxr/compute.h` | Add `struct VertexFrames` + `VertexFrames vertexFrames(Manifold&)` decl + `using geometry::VertexFrames;`. |
| `src/vertex_frames.cpp` | New: compute per-vertex frame from `vertexTangentBasis`. |
| `CMakeLists.txt` | Add `src/vertex_frames.cpp` to the library sources. |
| `bindings/mex/src/marshal.h` | Add `vertexFramesToStruct`. |
| `bindings/mex/src/nxr_compute_mex.cpp` | Add `cmdVertexFrames` + register `"vertexFrames"`. |
| `bindings/mex/matlab/+nxr/+manifold/+measure/vertexFrame.m` | New wrapper. |
| `bindings/mex/test/test_vertex_frames.m` | New test (raw-mex, icosahedron fixture). |

---

## Task 1: Implement, build, and test the `vertexFrames` export

**Files:** all the `nxr-compute` paths above.

- [ ] **Step 1: Write the failing test**

Create `bindings/mex/test/test_vertex_frames.m`:

```matlab
function test_vertex_frames
% test_vertex_frames — the per-vertex tangent-frame export ('vertexFrames').
% Raw-mex test on the icosahedron fixture (no Brainstorm, no +nxr): asserts the
% frame is [nV x 3], unit-length, orthonormal, and right-handed (e1 x e2 = n).
fprintf('[test_vertex_frames] starting\n');

thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fullfile(thisDir, '..', '..', '..');
hits = dir(fullfile(repoRoot, 'build', '**', ['nxr_compute.' mexext]));
assert(~isempty(hits), 'nxr_compute.%s not found under %s/build', mexext, repoRoot);
addpath(hits(1).folder);
clear nxr_compute   % drop any previously-loaded copy so the fresh build resolves

[V, F] = local_icosahedron();
nV = size(V, 1);
h  = nxr_compute('create', V, F);

vf = nxr_compute('vertexFrames', h);

assert(isequal(size(vf.e1),      [nV 3]), 'e1 must be nV x 3');
assert(isequal(size(vf.e2),      [nV 3]), 'e2 must be nV x 3');
assert(isequal(size(vf.normals), [nV 3]), 'normals must be nV x 3');

n1 = sqrt(sum(vf.e1.^2, 2));
n2 = sqrt(sum(vf.e2.^2, 2));
nn = sqrt(sum(vf.normals.^2, 2));
assert(max(abs(n1 - 1)) < 1e-9, 'e1 must be unit length');
assert(max(abs(n2 - 1)) < 1e-9, 'e2 must be unit length');
assert(max(abs(nn - 1)) < 1e-9, 'normals must be unit length');

assert(max(abs(sum(vf.e1 .* vf.e2,      2))) < 1e-9, 'e1 . e2 must be 0');
assert(max(abs(sum(vf.e1 .* vf.normals, 2))) < 1e-9, 'e1 . n must be 0');
assert(max(abs(sum(vf.e2 .* vf.normals, 2))) < 1e-9, 'e2 . n must be 0');

cr = cross(vf.e1, vf.e2, 2);
assert(max(max(abs(cr - vf.normals))) < 1e-9, 'e1 x e2 must equal n (right-handed)');

% Determinism
vf2 = nxr_compute('vertexFrames', h);
assert(isequal(vf.e1, vf2.e1) && isequal(vf.e2, vf2.e2) && isequal(vf.normals, vf2.normals), ...
    'vertexFrames must be deterministic');

nxr_compute('destroy', h);
fprintf('PASSED: per-vertex frames [%d x 3], orthonormal, right-handed, deterministic.\n', nV);
fprintf('ALL TESTS PASSED: test_vertex_frames\n');
end


function [V, F] = local_icosahedron()
% Unit icosahedron (12 vertices, 20 faces).
t = (1 + sqrt(5)) / 2;
V = [-1  t  0;  1  t  0; -1 -t  0;  1 -t  0; ...
      0 -1  t;  0  1  t;  0 -1 -t;  0  1 -t; ...
      t  0 -1;  t  0  1; -t  0 -1; -t  0  1];
V = V ./ sqrt(sum(V.^2, 2));
F = [1 12 6; 1 6 2; 1 2 8; 1 8 11; 1 11 12; ...
     2 6 10; 6 12 5; 12 11 3; 11 8 7; 8 2 9; ...
     4 10 5; 4 5 3; 4 3 7; 4 7 9; 4 9 10; ...
     5 10 6; 3 5 12; 7 3 11; 9 7 8; 10 9 2];
end
```

- [ ] **Step 2: Run the test to verify it FAILS**

Run via the MATLAB MCP `run_matlab_file` on the absolute path `/Users/diellorbasha/workspace/research/code/nxr-compute/bindings/mex/test/test_vertex_frames.m`.
Expected: FAIL — the current binary throws an unknown-command error for `'vertexFrames'` (the build predates it). Confirm the failure is the unknown command, not a missing-binary/path error.

- [ ] **Step 3: Declare `VertexFrames` in the header**

In `include/nxr/compute.h`, immediately after the `FaceFrames` struct + `FaceFrames frames(Manifold& m);` declaration (the `} // namespace nxr::manifold::geometry` closing that block is right after), add — *before* that closing brace — the struct and declaration:

```cpp
/* Per-vertex orthonormal tangent frame = geometry-central's vertexTangentBasis,
 * the 3D realization of the per-vertex tangent space (halfedgeVectorsInVertex) in
 * which the vertex connection-Laplacian's complex eigenvector coordinates live:
 *   e1 = basisX (the angle-0 / first-outgoing-halfedge axis)
 *   e2 = basisY = n x e1
 *   n  = vertex normal
 * A connection-mode coordinate z=(a,b) at vertex i decodes to a*e1_i + b*e2_i. */
struct VertexFrames {
    Eigen::MatrixXd e1;       // [nV, 3] — first tangent (basisX)
    Eigen::MatrixXd e2;       // [nV, 3] — second tangent (n × e1)
    Eigen::MatrixXd normals;  // [nV, 3] — vertex normals
};

VertexFrames vertexFrames(Manifold& m);
```

Then, near `using geometry::FaceFrames;` (the line that re-exports FaceFrames, ~line 1365), add directly after it:

```cpp
    using geometry::VertexFrames;
```

- [ ] **Step 4: Implement the compute in `src/vertex_frames.cpp`**

Create `src/vertex_frames.cpp`:

```cpp
#include "nxr/compute.h"

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <iostream>

namespace nxr::manifold::geometry {

using namespace geometrycentral;
using namespace geometrycentral::surface;

// See VertexFrames in compute.h for the gauge convention. This is exactly
// geometry-central's vertexTangentBasis (basisX, basisY = n × basisX) plus the
// vertex normal, so it matches the tangent space the connection Laplacian uses.
VertexFrames vertexFrames(Manifold& m) {
    auto& mesh = m.mesh();
    auto& geom = m.geometry();

    geom.requireVertexTangentBasis();
    geom.requireVertexNormals();

    int nV = m.nV();
    VertexFrames out;
    out.e1.resize(nV, 3);
    out.e2.resize(nV, 3);
    out.normals.resize(nV, 3);

    for (Vertex v : mesh.vertices()) {
        int vi = static_cast<int>(v.getIndex());
        Vector3 b0 = geom.vertexTangentBasis[v][0];
        Vector3 b1 = geom.vertexTangentBasis[v][1];
        Vector3 N  = geom.vertexNormals[v];
        out.e1     (vi, 0) = b0.x; out.e1     (vi, 1) = b0.y; out.e1     (vi, 2) = b0.z;
        out.e2     (vi, 0) = b1.x; out.e2     (vi, 1) = b1.y; out.e2     (vi, 2) = b1.z;
        out.normals(vi, 0) = N.x;  out.normals(vi, 1) = N.y;  out.normals(vi, 2) = N.z;
    }

    std::cout << "[vertex_frames] " << nV << " vertex frames computed" << std::endl;
    return out;
}

} // namespace nxr::manifold::geometry
```

- [ ] **Step 5: Register the source in `CMakeLists.txt`**

In `CMakeLists.txt`, in the library source list, add a line immediately after `src/face_frames.cpp`:

```cmake
  src/vertex_frames.cpp
```

- [ ] **Step 6: Add the marshaller in `bindings/mex/src/marshal.h`**

Immediately after the `faceFramesToStruct` function, add:

```cpp
inline mxArray* vertexFramesToStruct(const nxr::manifold::geometry::VertexFrames& f) {
    const char* fields[] = {"e1", "e2", "normals"};
    mxArray* s = mxCreateStructMatrix(1, 1, 3, fields);
    mxSetField(s, 0, "e1",      eigenMatrixToMx(f.e1));
    mxSetField(s, 0, "e2",      eigenMatrixToMx(f.e2));
    mxSetField(s, 0, "normals", eigenMatrixToMx(f.normals));
    return s;
}
```

- [ ] **Step 7: Add the MEX command + register it in `bindings/mex/src/nxr_compute_mex.cpp`**

Immediately after the `cmdFrames` function definition, add:

```cpp
void cmdVertexFrames(int /*nlhs*/, mxArray** plhs, int nrhs, const mxArray** prhs) {
    if (nrhs != 2) {
        throw std::invalid_argument("nxr_compute('vertexFrames', handle) takes exactly 1 argument");
    }
    ContextHolder& h = getHolder(prhs[1]);
    plhs[0] = vertexFramesToStruct(nxr::manifold::geometry::vertexFrames(*h.ctx));
}
```

Then, in the dispatch chain, immediately after the line `else if (cmd == "frames") cmdFrames(nlhs, plhs, nrhs, prhs);`, add:

```cpp
        else if (cmd == "vertexFrames")                cmdVertexFrames(nlhs, plhs, nrhs, prhs);
```

- [ ] **Step 8: Add the `+nxr` wrapper**

Create `bindings/mex/matlab/+nxr/+manifold/+measure/vertexFrame.m`:

```matlab
function out = vertexFrame(mctx)
%VERTEXFRAME  Per-vertex orthonormal tangent frame (the connection-Laplacian gauge).
%   out = nxr.manifold.measure.vertexFrame(mctx)
%
%   Returns a struct {e1, e2, normals}, each [nV x 3] — the 3D realization of the
%   per-vertex tangent space in which the vertex connection-Laplacian's complex
%   eigenvector coordinates are expressed: a coordinate z=(a,b) at vertex i is the
%   3D tangent vector a*e1(i,:) + b*e2(i,:). e1 ⟂ e2, both unit, e2 = n × e1.
%
%   See also: nxr.manifold.measure.frame (the per-FACE frame).
    out = nxr.manifold.impl.withHandle(mctx, @(h) nxr_compute('vertexFrames', h));
end
```

- [ ] **Step 9: Build the plugin**

Run (Bash tool, from the nxr-compute repo root):

```bash
cd /Users/diellorbasha/workspace/research/code/nxr-compute && bash scripts/build.sh Release
```

Expected: a clean build producing `build/Release/nxr_compute.mexmaca64` (re-configures CMake to pick up the new source). If the build errors, fix the C++ (do not touch the test) and rebuild.

- [ ] **Step 10: Run the test to verify it PASSES**

Run via the MATLAB MCP `run_matlab_file` on `/Users/diellorbasha/workspace/research/code/nxr-compute/bindings/mex/test/test_vertex_frames.m`.
Expected: prints `ALL TESTS PASSED: test_vertex_frames`. The `clear nxr_compute` in the test forces MATLAB to resolve the freshly-built binary. Paste the output.

- [ ] **Step 11: Commit (in the nxr-compute repo)**

```bash
cd /Users/diellorbasha/workspace/research/code/nxr-compute
git add include/nxr/compute.h src/vertex_frames.cpp CMakeLists.txt \
        bindings/mex/src/marshal.h bindings/mex/src/nxr_compute_mex.cpp \
        bindings/mex/matlab/+nxr/+manifold/+measure/vertexFrame.m \
        bindings/mex/test/test_vertex_frames.m
git commit -m "feat(mex): export per-vertex tangent frame (vertexFrames)"
```

---

## Task 2: Restage the rebuilt plugin into the local Brainstorm install

So the live Brainstorm session (and Plans B/C) can call `nxr.manifold.measure.vertexFrame`, copy the new binary + wrapper into the installed plugin and reload.

**Files:** the installed plugin dir (under `~/.brainstorm/plugins/nxr-compute/`).

- [ ] **Step 1: Write the verification test**

Create `dev/tests/test_vertex_frame_plugin.m` (in the **brainstorm3** repo):

```matlab
function test_vertex_frame_plugin
% Verifies nxr.manifold.measure.vertexFrame is callable through the INSTALLED
% Brainstorm plugin (after restaging the rebuilt binary), on the real cortex.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute plugin required: %s', errMsg);
bst_plugin('Load', 'nxr-compute');

SurfaceFile = find_cortex_20484V();
if isempty(SurfaceFile)
    fprintf('SKIP: no 20484-vertex cortex in the current protocol.\n');
    return;
end
TessMat = in_tess_bst(SurfaceFile);
V = TessMat.Vertices;
F = double(TessMat.Faces);
nV = size(V, 1);

mctx = nxr.manifold.context(V, F);
vf = nxr.manifold.measure.vertexFrame(mctx);

assert(isequal(size(vf.e1), [nV 3]) && isequal(size(vf.e2), [nV 3]) && isequal(size(vf.normals), [nV 3]), ...
    'vertexFrame fields must be nV x 3.');
assert(max(abs(sum(vf.e1 .* vf.e2, 2))) < 1e-6, 'e1 . e2 must be ~0.');
assert(max(abs(sqrt(sum(vf.e1.^2,2)) - 1)) < 1e-6, 'e1 must be unit.');
cr = cross(vf.e1, vf.e2, 2);
assert(max(max(abs(cr - vf.normals))) < 1e-6, 'e1 x e2 must equal n.');

fprintf('PASSED: nxr.manifold.measure.vertexFrame works through the installed plugin (nV=%d).\n', nV);
fprintf('ALL TESTS PASSED: test_vertex_frame_plugin\n');
end


function SurfaceFile = find_cortex_20484V()
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
if isempty(sSubjects)
    return;
end
allSubj = [sSubjects.Subject];
fallback = '';
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex')
            continue;
        end
        try
            T = load(file_fullpath(surf(iF).FileName), 'Vertices', 'Reg');
        catch
            continue;
        end
        if size(T.Vertices, 1) ~= 20484
            continue;
        end
        hasReg = isfield(T,'Reg') && isstruct(T.Reg) && isfield(T.Reg,'Sphere') ...
                 && isfield(T.Reg.Sphere,'Vertices') && ~isempty(T.Reg.Sphere.Vertices);
        if hasReg
            SurfaceFile = surf(iF).FileName;
            return;
        elseif isempty(fallback)
            fallback = surf(iF).FileName;
        end
    end
end
if isempty(SurfaceFile)
    SurfaceFile = fallback;
end
end
```

- [ ] **Step 2: Run it to verify it FAILS**

Run `dev/tests/test_vertex_frame_plugin.m` via the MATLAB MCP.
Expected: FAIL — `Unrecognized function 'vertexFrame'` (the installed plugin still has the old binary + no `vertexFrame.m` wrapper).

- [ ] **Step 3: Restage the rebuilt binary + wrapper into the installed plugin**

Find the installed plugin dir and copy the new artifacts (Bash tool):

```bash
NXR=/Users/diellorbasha/workspace/research/code/nxr-compute
DEST=$(dirname "$(ls ~/.brainstorm/plugins/nxr-compute/*/nxr_compute.mexmaca64)")
echo "Installed plugin dir: $DEST"
cp "$NXR"/build/**/nxr_compute.mexmaca64 "$DEST"/nxr_compute.mexmaca64
mkdir -p "$DEST/+nxr/+manifold/+measure"
cp "$NXR"/bindings/mex/matlab/+nxr/+manifold/+measure/vertexFrame.m "$DEST/+nxr/+manifold/+measure/vertexFrame.m"
ls -la "$DEST/+nxr/+manifold/+measure/vertexFrame.m" "$DEST/nxr_compute.mexmaca64"
```

(If the `build/**` glob does not expand, locate the binary with `find "$NXR/build" -name nxr_compute.mexmaca64` and copy that path.)

- [ ] **Step 4: Reload the plugin in the live session, run the test to verify it PASSES**

Via the MATLAB MCP:

```matlab
clear nxr_compute
bst_plugin('Unload', 'nxr-compute');
bst_plugin('Load', 'nxr-compute');
test_vertex_frame_plugin
```

Expected: `ALL TESTS PASSED: test_vertex_frame_plugin`. (`clear nxr_compute` drops the stale MEX so the restaged binary loads — this clears only that function, NOT `clear` with no args, which would wipe `GlobalData`.) Paste the output.

- [ ] **Step 5: Commit the verification test (brainstorm3 repo)**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add dev/tests/test_vertex_frame_plugin.m
git commit -m "test(conn-phase): verify nxr vertexFrame export through the installed plugin"
```

---

## Notes for the implementer

- **Two repos:** Task 1 commits in `nxr-compute`; Task 2's test commits in `brainstorm3` (on branch `feature/connection-phase-readout`). Do not mix.
- **Never run bare `clear`** in the live MATLAB session (wipes `GlobalData`). `clear nxr_compute` (a single named MEX) is fine and necessary to swap the binary.
- **The build is slow** and re-configures CMake; that's expected for a new source file.
- **Deferred (operational publish step, not in this plan):** cutting a new `nxr-compute` GitHub release with the rebuilt macOS asset and bumping the `PlugDesc` version/URL in `toolbox/core/bst_plugin.m`. The local restage (Task 2) is sufficient for Plans B/C on this machine; the formal release is needed only for clean installs on other machines.
- **Gauge correctness** (that the exported frame is the *right rotation* for decoding eigenvectors, not just any orthonormal frame) is validated in Plan B via the decoded-Fiedler-vs-`smoothVertex` cross-check, not here.
```
