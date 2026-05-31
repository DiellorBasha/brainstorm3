# Design: nxr-compute as a Brainstorm Plugin (minimal slice)

- **Date:** 2026-05-31
- **Author:** Diellor Basha (with Claude)
- **Status:** Design — pending review
- **Branch:** feature/eigenmode-spectrum (fork of brainstorm3)

## 1. Goal

Integrate the `nxr-compute` C++ geometry engine into Brainstorm as a native,
on-demand plugin, and prove the integration end-to-end with a single, minimal
consumer: `tess_laplacian`.

This is the first vertical slice. The long-term intent is for nxr-compute to be a
general-purpose geometry **compute backend** that many Brainstorm functions call
on demand (vertices + faces in, arrays out). That broader build-out happens in a
later, separate effort. This document covers **only** the plumbing plus one proof
consumer.

## 2. Scope

**In scope**
- Register `nxr-compute` as a native Brainstorm plugin (`PlugDesc` entry in
  `toolbox/core/bst_plugin.m`).
- Ship a **macOS** binary (`nxr_compute.mexmaca64`) packaged with its `+nxr`
  MATLAB wrapper as a downloadable release zip.
- Wire `tess_laplacian` to be **nxr-first with a MATLAB fallback**.
- A parity test proving nxr and MATLAB produce equivalent operators, plus a
  plugin-lifecycle test.

**Out of scope (deferred to later efforts)**
- Linux/Windows binaries and CI to build them.
- A general `bst_nxr` bridge module (surface-file conversion, result mapping,
  scout helpers). Not needed for this slice — `tess_laplacian` already receives
  raw `Vertices`/`Faces`.
- All other nxr capabilities (eigenmodes, Hodge decomposition, trivial
  connections, heat-geodesics, vector transport, BFF, direction fields, …).
- Any GUI `process_*.m` wrapper.
- Exposing barycentric/galerkin mass via nxr (requires a C++/mex change — see §6).

## 3. Background

### 3.1 nxr-compute already ships a working MATLAB binding
The repo at `/Users/diellorbasha/workspace/research/code/nxr-compute`
(`github.com/neurodynamics-xr/nxr-compute`) already provides:
- A built MEX for this machine: `build/Release/nxr_compute.mexmaca64`.
- A `nxr_compute('<command>', ...)` dispatcher (35 commands).
- A clean functional wrapper package `+nxr` under `bindings/mex/matlab/`
  (e.g. `nxr.manifold.context`, `nxr.manifold.operator.stiffness`,
  `nxr.manifold.operator.mass`).

So this work is **packaging + surfacing an existing binding**, not authoring a
binding. nxr-compute is intentionally application-agnostic ("no file I/O, no
application knowledge"); we preserve that — no Brainstorm code goes into nxr.

### 3.2 Brainstorm plugins are optional and install-on-demand
`bst_plugin.m` registers 47 plugins; a stock Brainstorm boots with **none**
installed. SPM12 is the model to copy: `AutoLoad = 0`, used across ~10 core
functions (`bst_normalize_mni`, `mri_coregister`, …), each guarded by:

```matlab
[isInstalled, errMsg] = bst_plugin('Install', 'spm12');
if ~isInstalled, return; end
```

The rule is **never assume a plugin is present**; guard every use with an
on-demand `Install`. nxr is wired exactly like SPM. Whether a given consumer also
keeps a pure-MATLAB fallback is a per-function judgment call — warranted for
`tess_laplacian` because it (a) sits on a hot inner loop ("recomputed on every
project/filter call"), (b) already has an audited dependency-free implementation,
and (c) supports the eigenmode pipeline, which we want runnable on a stock install.

## 4. Architecture

Three seams; only the first two are built in this slice:

1. **`+nxr` library + mex** — shipped from nxr-compute, untouched. The Brainstorm
   side targets the stable functional wrapper (`nxr.manifold.*`), never the raw
   `nxr_compute('cmd', …)` dispatcher, so mex-internal churn never reaches us.
2. **Plugin registration** — one `PlugDesc` entry in `bst_plugin.m`.
3. **Consumer wiring** — a backend switch inside `tess_laplacian`. This is the
   only consumer in this slice; a general `bst_nxr` bridge is deferred.

## 5. Deliverable 1 — The plugin package (macOS)

A release zip attached to a GitHub release of `neurodynamics-xr/nxr-compute`,
laying the wrapper package and the binary side by side at the root so a single
folder on the MATLAB path resolves both:

```
nxr-compute-<ver>-mac.zip
├── +nxr/                     # from bindings/mex/matlab/+nxr/
├── nxr_compute.mexmaca64     # from build/Release/
└── README.md                 # from bindings/mex/matlab/README.md
```

The `+nxr/...context.m` wrapper invokes the bare mex `nxr_compute('…', …)`, so the
binary must be named `nxr_compute.mexmaca64` and live on the path next to `+nxr`.

**Packaging step:** add `scripts/package-plugin.sh` to the nxr-compute repo that
runs `scripts/build.sh Release`, then assembles the layout above into a versioned
zip. (Lives in the nxr repo because it packages nxr artifacts; it contains no
Brainstorm knowledge.)

## 6. Deliverable 2 — `PlugDesc` entry in `bst_plugin.m`

Modeled on OpenMEEG/DUNEuro (native code, per-OS, binary-verified `TestFile`):

```matlab
% === ANATOMY: NXR-COMPUTE ===
PlugDesc(end+1)              = GetStruct('nxr-compute');
PlugDesc(end).Version        = '<release tag>';
PlugDesc(end).Category       = 'Anatomy';
PlugDesc(end).AutoUpdate     = 0;
PlugDesc(end).AutoLoad       = 0;            % SPM-style install-on-demand
PlugDesc(end).URLinfo        = 'https://github.com/neurodynamics-xr/nxr-compute';
PlugDesc(end).ReadmeFile     = 'README.md';
PlugDesc(end).CompiledStatus = 1;            % native code, download-only
switch bst_get('OsType')
    case 'mac64arm'
        PlugDesc(end).URLzip   = '<github release asset URL, mac arm>';
        PlugDesc(end).TestFile = 'nxr_compute.mexmaca64';
    % 'linux64' / 'win64' / 'mac64' arms added when those binaries exist
end
```

The per-OS `switch` is laid down now so adding other platforms later is a
fill-in-the-blank, not a restructure. On unsupported platforms `URLzip`/`TestFile`
stay empty, and install fails cleanly — consumers fall back to MATLAB.

## 7. Deliverable 3 — `tess_laplacian` nxr-first with MATLAB fallback

A backend switch at the top of `tess_laplacian` (`toolbox/anatomy/tess_laplacian.m`).
nxr serves the stiffness `L` for **all** mass types (it is mass-independent); nxr
serves the mass `M` **only for Voronoi** (the only variant the mex exposes today).

| `MassType` request | `L` (stiffness) | `M` (mass) |
|---|---|---|
| `voronoi`     | nxr | **nxr** |
| `barycentric` | nxr | MATLAB fallback |
| `galerkin`    | nxr | MATLAB fallback |

If the plugin is unavailable, the entire existing MATLAB implementation runs
unchanged.

```matlab
% --- backend selection ---
% Use nxr only if it is already loaded. tess_laplacian is a hot-loop
% operator, so it never triggers a (possibly networked) install itself;
% installation is a one-time action (Plugins menu, or a higher-level
% pipeline calling bst_plugin('Install','nxr-compute') once up front).
useNxr = nxr_is_loaded();   % local helper, cheap in-memory check

if useNxr
    L = nxr.manifold.operator.stiffness(nxr.manifold.context(Vertices, Faces));
    if strcmp(MassType, 'voronoi')
        M = nxr.manifold.operator.mass(nxr.manifold.context(Vertices, Faces));
    else
        M = local_mass_matlab(Vertices, Faces, MassType);  % existing MATLAB mass
    end
    if Symmetrize, L = (L + L')/2; end   % match existing behavior
else
    % ... existing pure-MATLAB implementation, unchanged ...
end
```

where `nxr_is_loaded()` is a local helper:

```matlab
function tf = nxr_is_loaded()
    tf = false;
    try
        PlugDesc = bst_plugin('GetInstalled', 'nxr-compute');
        tf = ~isempty(PlugDesc) && isfield(PlugDesc,'isLoaded') && PlugDesc.isLoaded;
    catch
        tf = false;   % bst_plugin unavailable / any error → MATLAB path
    end
end
```

The existing MATLAB mass code is factored into a local helper (`local_mass_matlab`)
so both paths share it; the existing stiffness code remains the fallback. (Building
the context twice above is for clarity; the implementation builds it once and reads
both `K` and `M` from it.)

### Reconciliation work (where the real care goes)
These must be verified, not assumed — the parity test (§8) is the gate:
- **Faces convention.** `+nxr/...context.m`'s docstring says faces are *zero-based*,
  but it `int32`s and passes them through. Brainstorm `Faces` are **1-based**.
  Determine empirically whether the mex marshalling subtracts 1; pass faces in the
  convention nxr actually expects. A wrong convention yields a garbage/degenerate
  operator the parity test will flag.
- **Sign convention.** Confirm nxr's `K` is PSD with the same sign as the existing
  `L` (negative off-diagonals, `u'Lu ≥ 0`).
- **Symmetry.** Apply `(L+L')/2` to match the existing `Symmetrize` default.
- **Manifold semantics.** Preserve the existing `CheckManifold` behavior; nxr
  requires a clean 2-manifold (throws `nxr:nonManifold` otherwise). Decide whether
  a non-manifold input forces the MATLAB fallback or surfaces the nxr error.

## 8. Deliverable 4 — Tests

Run via the MATLAB MCP (`run_matlab_test_file`) against a live Brainstorm path.

1. **Plugin lifecycle test.** `bst_plugin('Install','nxr-compute')` succeeds; the
   `TestFile` is found; `nxr_compute('version')` returns a string; the `+nxr`
   package resolves; unload restores prior path state.
2. **Parity test (pass/fail gate).** On a real downsampled cortex surface
   (icosphere-downsampled, clean 2-manifold):
   - `max(abs(L_nxr - L_matlab))` within tolerance, for the stiffness.
   - `max(abs(M_nxr - M_matlab))` within tolerance, for **Voronoi** mass.
   - Sanity: `L` is symmetric, row sums ≈ 0, PSD (smallest eigenvalue ≥ -tol).
   - Tolerance is "numerically close," not bit-exact — the two implementations
     may differ in obtuse-triangle Voronoi handling and assembly order. Record the
     observed max delta and set `tol` accordingly; investigate if it is large.
3. **Fallback test.** With the plugin force-unloaded, `tess_laplacian` still
   produces correct operators (the original MATLAB path), confirming graceful
   degradation.

## 9. Multi-platform (documented, not built)

Currently nxr-compute has **only** `mexmaca64`; there is no CI. Future work:
- A GitHub Actions matrix building `nxr_compute.mexa64` (Linux) and
  `nxr_compute.mexw64` (Windows), each packaged by `scripts/package-plugin.sh`
  and attached to the same release.
- Add `linux64` / `win64` (and Intel `mac64` if needed) arms to the `PlugDesc`
  `switch` (§6). The per-OS structure already exists, so this is additive.

## 10. Risks

- **Faces 0/1-based mismatch** — most likely failure mode; caught by the parity
  test (§8) and resolved in §7 reconciliation.
- **Voronoi mass divergence** — nxr's Voronoi vs the existing mixed-Voronoi obtuse
  handling (Meyer et al.) may not be bit-identical; the parity tolerance accounts
  for this. If the delta is large, restrict the nxr `M` path or reconcile schemes.
- **Manifold strictness** — nxr throws on non-manifold input where the MATLAB path
  only warns; the backend switch must handle this so behavior is predictable.

## 11. Decisions log (from brainstorming)

- **Role:** nxr is a general-purpose geometry compute backend (this slice proves
  the pattern via one consumer).
- **Delivery:** pre-built binaries (no build-at-install); one zip per platform,
  hosted on GitHub releases.
- **Registration:** native `PlugDesc` entry in the fork's `bst_plugin.m`.
- **Integration depth:** library + one consumer (`tess_laplacian`); general
  `bst_nxr` bridge deferred.
- **Wiring model:** SPM-style guarded use. `tess_laplacian` (a hot-loop operator)
  uses nxr only when already loaded and never auto-installs on the operator path;
  installation is a one-time action. MATLAB fallback covers the not-loaded case.
- **Mass variants:** `L` always from nxr; `M` from nxr only for Voronoi (the only
  variant the mex exposes); barycentric/galerkin `M` via MATLAB. Extending the mex
  to all three variants is future work in the nxr repo.
- **First platform:** macOS only now; Linux/Windows packaging later (none exist in
  nxr-compute yet).
