# process_eigen orchestrator + process_eigenspectrum — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expose the eigen-domain analysis engine (`bst_eigen`) in the Brainstorm Process panel so eigenspectra can be batch-computed across many recordings/subjects and stored as `timefreq_` nodes — the spatial-spectral twin of `process_psd`.

**Architecture:** Three layers mirroring the timefreq family: `process_eigenspectrum` (registered, panel-visible) delegates its `Run` to a **pure** `process_eigen` orchestrator (no `GetDescription`, dispatches method by caller name), which calls the existing `bst_eigen` compute+save engine. The only engine change is closing one resolver gap so `bst_eigen` finds the `eigen_` basis implicitly from the input's `SurfaceFile`.

**Tech Stack:** MATLAB; Brainstorm process system (`eval(macro_method)` verb dispatch, `GetDescription`/`Run`); tests run via the MATLAB MCP against the live S01 protocol.

## Global Constraints

- Spec: `dev/2026-06-23-process-eigen-design.md`.
- v1 builds **only** the orchestrator + spectrum. `process_eigenfilter`/`process_eigenwavelet` are named in the dispatch map but NOT created.
- Do not change `bst_eigen`'s dispatch, compute, or `BuildSpectrumTimefreq` node shape — only `GetEigenBasis` + add a `Tau` default.
- `process_eigen` is pure: it must NOT define `GetDescription` (so the panel loader never registers it).
- Variant resolved implicitly from `SurfaceFile`; Variant chosen by dropdown, default `'Laplace-Beltrami'`; `nModes`/`Tau` NOT exposed.
- Welch-style time windowing is a first-class option, semantics copied from `process_psd` (`win_length`/`win_overlap`/`win_std`).
- New process Index = **484** (verified free in 478–495).
- Follow Brainstorm option-reading idioms exactly (see `process_timefreq.m:88–185`): `radio_linelabel` Value is the key string; `timewindow`/`win_length`/`win_overlap` Values are `{value, units, prec}` cells read via `.Value{1}`.
- MATLAB rule (memory): never `clear` in the live session; edited `.m` files auto-reload via `rehash` if needed.
- Tests are functions that `error(...)` on any failure and print `ALL TESTS PASSED`; place them in `dev/tests/`; they SKIP cleanly if the live protocol lacks a fixture.
- Each commit message ends with the two Co-Authored-By / Claude-Session trailers used in this repo.

---

### Task 1: Close the eigen-basis resolver gap in `bst_eigen`

Make `bst_eigen` resolve the `eigen_` node implicitly from a results file's `SurfaceFile` + Variant, instead of erroring. This is the foundational, independently valuable change: after it, `bst_eigen(resultsFile, OPTIONS)` with `OPTIONS.EigenFile=[]` works end-to-end.

**Files:**
- Modify: `toolbox/eigen/bst_eigen.m` (add `Def_OPTIONS.Tau`; rewrite the `GetEigenBasis` resolver block ~lines 300–319)
- Test: `dev/tests/test_bst_eigen_resolve.m` (create)

**Interfaces:**
- Consumes: `bst_get('EigenFileForSurface', SurfaceFile, Variant, nModes, Tau)` → `[sSubject, iSubject, iSurface, iEigen]`; returns `[]` for all when no match (bst_get.m:244–247, 1385–1411). Relative path is `sSubject.Surface(iSurface).Eigen(iEigen).FileName`.
- Produces: `bst_eigen(Data, OPTIONS)` with `OPTIONS.EigenFile=[]`, `OPTIONS.Variant`, `OPTIONS.nModes`, `OPTIONS.Tau` resolves the basis from `SurfaceFile`. New error id `bst_eigen:NoEigenForSurface` when none found. `OPTIONS.Tau` (default `[]`) now exists.

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_bst_eigen_resolve.m`:

```matlab
function test_bst_eigen_resolve()
% TEST_BST_EIGEN_RESOLVE: bst_eigen resolves the eigen_ basis implicitly from a results
% file's SurfaceFile (Variant=Laplace-Beltrami), and raises a clear error when none exists.
% SKIPs if the live protocol lacks an LBO eigen node or a constrained results file on its surface.
    nFail = 0; chk = @i_chk;
    [ef, surf] = i_find_eigen('Laplace-Beltrami');
    if isempty(ef); fprintf('SKIP test_bst_eigen_resolve: no LBO eigen node.\n'); return; end
    fprintf('LBO eigen: %s\n  surface: %s\n', ef, surf);

    % (1) resolver returns a node for the surface
    [sS, ~, iSurf, iE] = bst_get('EigenFileForSurface', surf, 'Laplace-Beltrami');
    nFail = nFail + chk('EigenFileForSurface resolves LBO', ~isempty(iE) && ~isempty(sS));
    if ~isempty(iE)
        nFail = nFail + chk('resolved entry has a FileName', ~isempty(sS.Surface(iSurf).Eigen(iE).FileName));
    end

    % (2) end-to-end implicit resolution on a real constrained results file
    rf = i_find_results(surf);
    if isempty(rf)
        fprintf('SKIP end-to-end: no constrained results on surface.\n');
    else
        fprintf('results: %s\n', rf);
        O = bst_eigen(); O.Method = 'spectrum'; O.Variant = 'Laplace-Beltrami';
        O.EigenFile = []; O.iTargetStudy = 'NoSave';
        [out, msg, err] = bst_eigen(rf, O);
        nFail = nFail + chk('implicit resolve: no error', err == 0 && isempty(msg));
        ok = iscell(out) && ~isempty(out) && isstruct(out{1});
        nFail = nFail + chk('implicit resolve: returns a struct', ok);
        if ok
            nFail = nFail + chk('spectrum has Freqs (sqrt-lambda)', isfield(out{1}, 'Freqs') && ~isempty(out{1}.Freqs));
            nFail = nFail + chk('spectrum Method == spectrum', isfield(out{1}, 'Method') && strcmp(out{1}.Method, 'spectrum'));
        end
    end

    % (3) negative: a variant absent on this surface -> clear typed error
    bad = 'Connection Laplacian';
    [~, ~, ~, iEb] = bst_get('EigenFileForSurface', surf, bad);
    if isempty(iEb) && ~isempty(rf)
        Ob = bst_eigen(); Ob.Method = 'spectrum'; Ob.Variant = bad;
        Ob.EigenFile = []; Ob.iTargetStudy = 'NoSave';
        gotErr = false;
        try, bst_eigen(rf, Ob); catch ME, gotErr = strcmp(ME.identifier, 'bst_eigen:NoEigenForSurface'); end
        nFail = nFail + chk('missing variant -> NoEigenForSurface error', gotErr);
    else
        fprintf('SKIP negative: a Connection basis exists or no results.\n');
    end

    fprintf('\n==== test_bst_eigen_resolve: %d failed ====\n', nFail);
    if nFail > 0; error('test_bst_eigen_resolve FAILED'); end
    disp('ALL TESTS PASSED');
end

function r = i_chk(nm, cond)
    if cond; r = 0; fprintf('  ok   %s\n', nm); else; r = 1; fprintf('  FAIL %s\n', nm); end
end
function [ef, surf] = i_find_eigen(want)
    ef = ''; surf = '';
    PI = bst_get('ProtocolInfo'); if isempty(PI); return; end
    d = dir(fullfile(PI.SUBJECTS, '**', 'eigen_*.mat')); [~, o] = sort([d.datenum], 'descend');
    for i = o(:)'
        rel = strrep(fullfile(d(i).folder, d(i).name), [PI.SUBJECTS filesep], '');
        try
            m = in_bst_eigen(rel, 'Variant', 'ParentSurface');
            if strcmpi(m.Variant, want); ef = rel; surf = m.ParentSurface; return; end
        catch
        end
    end
end
function rf = i_find_results(surf)
    % First constrained (nComponents==1) results file mapped to this surface.
    rf = '';
    sP = bst_get('ProtocolStudies'); if isempty(sP); return; end
    for s = 1:numel(sP.Study)
        R = sP.Study(s).Result;
        for r = 1:numel(R)
            try
                m = in_bst_results(R(r).FileName, 0, 'SurfaceFile', 'nComponents');
                if isfield(m, 'SurfaceFile') && ~isempty(m.SurfaceFile) ...
                        && file_compare(m.SurfaceFile, surf) ...
                        && (isempty(m.nComponents) || isequal(m.nComponents, 1))
                    rf = R(r).FileName; return;
                end
            catch
            end
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run (MATLAB MCP `run_matlab_file`): `dev/tests/test_bst_eigen_resolve.m`
Expected: FAIL — the end-to-end check errors with `bst_eigen:AutoResolveTODO` ("Auto-resolving the eigen_ basis … is not implemented yet").

- [ ] **Step 3: Add the `Tau` default option**

In `toolbox/eigen/bst_eigen.m`, in the `===== DEFAULT OPTIONS =====` block, add after the `nModes` line (currently bst_eigen.m:79):

```matlab
Def_OPTIONS.Tau           = [];          % Dirac-type smoothing for basis resolution; [] => any
```

- [ ] **Step 4: Implement the implicit resolution in `GetEigenBasis`**

In `toolbox/eigen/bst_eigen.m`, replace the resolver block in `GetEigenBasis` (currently bst_eigen.m:304–312):

```matlab
    EigenFile = OPTIONS.EigenFile;
    if isempty(EigenFile) && ~isempty(SurfaceFile)
        % TODO: EigenFile = bst_get('EigenFileForSurface', SurfaceFile, OPTIONS.Variant);
        error('bst_eigen:AutoResolveTODO', ...
            'Auto-resolving the eigen_ basis from the surface is not implemented yet; pass OPTIONS.EigenFile.');
    end
    if isempty(EigenFile)
        error('bst_eigen:NoEigenFile', 'No eigen_ basis specified (OPTIONS.EigenFile).');
    end
```

with:

```matlab
    EigenFile = OPTIONS.EigenFile;
    if isempty(EigenFile) && ~isempty(SurfaceFile)
        % Resolve the eigen_ node implicitly from the surface + operator family (the spatial
        % analogue of bst_timefreq deriving its axis from the recordings). Finding is bst_get's
        % job; a single surface can host several eigenbases, so Variant selects among them.
        Variant = OPTIONS.Variant;
        if isempty(Variant); Variant = 'Laplace-Beltrami'; end
        [sSubject, ~, iSurface, iEigen] = bst_get('EigenFileForSurface', SurfaceFile, Variant, OPTIONS.nModes, OPTIONS.Tau);
        if isempty(iEigen)
            error('bst_eigen:NoEigenForSurface', ...
                'No ''%s'' eigenbasis found for this surface — compute it first.', Variant);
        end
        EigenFile = sSubject.Surface(iSurface).Eigen(iEigen).FileName;
    end
    if isempty(EigenFile)
        error('bst_eigen:NoEigenFile', 'No eigen_ basis specified (OPTIONS.EigenFile).');
    end
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `dev/tests/test_bst_eigen_resolve.m`
Expected: `ALL TESTS PASSED` (or clean SKIP lines if no fixture, but the resolver check (1) must pass when an LBO node exists).

- [ ] **Step 6: Commit**

```bash
git add toolbox/eigen/bst_eigen.m dev/tests/test_bst_eigen_resolve.m
git commit -m "feat(eigen): bst_eigen resolves eigen_ basis implicitly from SurfaceFile

Close the GetEigenBasis TODO: call bst_get('EigenFileForSurface', Surface,
Variant, nModes, Tau) when OPTIONS.EigenFile is empty; clear NoEigenForSurface
error when none exists. Add OPTIONS.Tau default. Test: test_bst_eigen_resolve.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01C3Lf7KA3E5Y7maW6Knm2LM"
```

---

### Task 2: `process_eigen` — pure shared orchestrator

Create the pure orchestrator that maps the calling process name to a `bst_eigen` Method, translates panel options to `bst_eigen` OPTIONS, and runs the engine per input (so one bad input doesn't abort the batch). No `GetDescription` ⇒ never appears in the panel.

**Files:**
- Create: `toolbox/process/functions/process_eigen.m`
- Test: `dev/tests/test_process_eigen.m` (create)

**Interfaces:**
- Consumes: `bst_eigen(FileName, OPTIONS)` from Task 1; `bst_process('GetInputStruct', FileNames)`; `func2str(sProcess.Function)`.
- Produces: `process_eigen('Run', sProcess, sInputs)` → `OutputFiles` (cell of saved filenames). Dispatch map: `process_eigenspectrum`→`'spectrum'`, `process_eigenfilter`→`'filter'`, `process_eigenwavelet`→`'wavelet'`; unknown caller ⇒ `bst_report('Error',…)` + empty output (no throw). Reads options: `variant`/`measure` (`radio_linelabel`, Value=key string), `timewindow`/`win_length`/`win_overlap` (`.Value{1}`), `win_std` (checkbox → `WinFunc` `'mean'`/`'mean+std'`).

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_process_eigen.m`:

```matlab
function test_process_eigen()
% TEST_PROCESS_EIGEN: the pure orchestrator maps caller name -> method, translates options,
% runs bst_eigen per input, and rejects unknown callers without throwing.
% SKIPs if the live protocol lacks an LBO eigen node + a constrained results file on its surface.
    nFail = 0; chk = @i_chk;
    [ef, surf] = i_find_eigen('Laplace-Beltrami');
    if isempty(ef); fprintf('SKIP test_process_eigen: no LBO eigen node.\n'); return; end
    rf = i_find_results(surf);
    if isempty(rf); fprintf('SKIP test_process_eigen: no constrained results on surface.\n'); return; end

    % Synthetic sProcess mimicking process_eigenspectrum (str2func name need not exist on disk)
    sProcess = struct();
    sProcess.Function = str2func('process_eigenspectrum');
    sProcess.options.variant.Value     = 'Laplace-Beltrami';
    sProcess.options.measure.Value     = 'power';
    sProcess.options.timewindow.Value  = [];
    sProcess.options.win_length.Value  = {[], 's', []};
    sProcess.options.win_overlap.Value = {50, '%', 1};
    sProcess.options.win_std.Value     = 0;
    sInputs = bst_process('GetInputStruct', rf);

    Out = process_eigen('Run', sProcess, sInputs);
    ok = iscell(Out) && ~isempty(Out) && ischar(Out{1});
    nFail = nFail + chk('process_eigen returns a saved node', ok);
    if ok
        R = in_bst_timefreq(Out{1}, 'Method', 'Freqs');
        nFail = nFail + chk('node Method == spectrum', strcmp(R.Method, 'spectrum'));
        nFail = nFail + chk('node has Freqs', ~isempty(R.Freqs));
        try, file_delete(file_fullpath(Out{1}), 1); db_reload_studies(sInputs.iStudy); catch; end
    end

    % Unknown caller -> empty output, no throw
    sBad = sProcess; sBad.Function = str2func('process_eigenbogus');
    ob = process_eigen('Run', sBad, sInputs);
    nFail = nFail + chk('unknown caller -> empty output', isempty(ob));

    fprintf('\n==== test_process_eigen: %d failed ====\n', nFail);
    if nFail > 0; error('test_process_eigen FAILED'); end
    disp('ALL TESTS PASSED');
end

function r = i_chk(nm, cond)
    if cond; r = 0; fprintf('  ok   %s\n', nm); else; r = 1; fprintf('  FAIL %s\n', nm); end
end
function [ef, surf] = i_find_eigen(want)
    ef = ''; surf = '';
    PI = bst_get('ProtocolInfo'); if isempty(PI); return; end
    d = dir(fullfile(PI.SUBJECTS, '**', 'eigen_*.mat')); [~, o] = sort([d.datenum], 'descend');
    for i = o(:)'
        rel = strrep(fullfile(d(i).folder, d(i).name), [PI.SUBJECTS filesep], '');
        try
            m = in_bst_eigen(rel, 'Variant', 'ParentSurface');
            if strcmpi(m.Variant, want); ef = rel; surf = m.ParentSurface; return; end
        catch
        end
    end
end
function rf = i_find_results(surf)
    rf = '';
    sP = bst_get('ProtocolStudies'); if isempty(sP); return; end
    for s = 1:numel(sP.Study)
        R = sP.Study(s).Result;
        for r = 1:numel(R)
            try
                m = in_bst_results(R(r).FileName, 0, 'SurfaceFile', 'nComponents');
                if isfield(m, 'SurfaceFile') && ~isempty(m.SurfaceFile) ...
                        && file_compare(m.SurfaceFile, surf) ...
                        && (isempty(m.nComponents) || isequal(m.nComponents, 1))
                    rf = R(r).FileName; return;
                end
            catch
            end
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `dev/tests/test_process_eigen.m`
Expected: FAIL — `Undefined function or variable 'process_eigen'` (file not created yet).

- [ ] **Step 3: Implement `process_eigen.m`**

Create `toolbox/process/functions/process_eigen.m`:

```matlab
function varargout = process_eigen( varargin )
% PROCESS_EIGEN: Pure orchestrator for eigen-domain analysis (the spatial-spectral twin of
% PROCESS_TIMEFREQ). It is NOT a registered process (no GetDescription => never shown in the
% pipeline panel). Sibling processes (process_eigenspectrum, and future process_eigenfilter/
% process_eigenwavelet) delegate their Run to this function, exactly as process_psd/
% process_hilbert delegate to process_timefreq('Run', ...). The eigen METHOD is inferred from
% the calling process name; panel options are translated to bst_eigen OPTIONS; bst_eigen does
% the read -> resolve eigen_ basis -> compute -> save.

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Diellor Basha, 2026

eval(macro_method);
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    % Method from the calling process name (mirrors process_timefreq.m:91-100)
    switch func2str(sProcess.Function)
        case 'process_eigenspectrum', Method = 'spectrum';
        case 'process_eigenfilter',   Method = 'filter';
        case 'process_eigenwavelet',  Method = 'wavelet';
        otherwise
            bst_report('Error', sProcess, sInputs, 'Unsupported eigen process.');
            return;
    end
    % Build engine options: defaults, then overlay the panel options.
    OPTIONS = bst_eigen();
    OPTIONS.Method       = Method;
    OPTIONS.iTargetStudy = [];   % each output node lands in its own input's study
    % Variant (radio_linelabel: Value is the key string), default Laplace-Beltrami
    if isfield(sProcess.options, 'variant') && ~isempty(sProcess.options.variant.Value)
        OPTIONS.Variant = sProcess.options.variant.Value;
    end
    % Measure (power/magnitude)
    if isfield(sProcess.options, 'measure') && ~isempty(sProcess.options.measure.Value)
        OPTIONS.Measure = sProcess.options.measure.Value;
    end
    % Time window
    if isfield(sProcess.options, 'timewindow') && ~isempty(sProcess.options.timewindow.Value) ...
            && iscell(sProcess.options.timewindow.Value)
        OPTIONS.TimeWindow = sProcess.options.timewindow.Value{1};
    end
    % Welch-style windowing (window length empty => single full window)
    if isfield(sProcess.options, 'win_length') && ~isempty(sProcess.options.win_length.Value) ...
            && iscell(sProcess.options.win_length.Value)
        OPTIONS.WinLength  = sProcess.options.win_length.Value{1};
        OPTIONS.WinOverlap = sProcess.options.win_overlap.Value{1};
    end
    % Aggregate across windows
    if isfield(sProcess.options, 'win_std') && ~isempty(sProcess.options.win_std.Value)
        if sProcess.options.win_std.Value
            OPTIONS.WinFunc = 'mean+std';
        else
            OPTIONS.WinFunc = 'mean';
        end
    end
    % Run the engine per input so one failure does not abort the batch.
    for iIn = 1:numel(sInputs)
        try
            [out, Messages, isError] = bst_eigen(sInputs(iIn).FileName, OPTIONS);
            if isError
                bst_report('Error', sProcess, sInputs(iIn), Messages);
            else
                if ~isempty(Messages)
                    bst_report('Warning', sProcess, sInputs(iIn), Messages);
                end
                OutputFiles = [OutputFiles, out]; %#ok<AGROW>
            end
        catch ME
            bst_report('Error', sProcess, sInputs(iIn), ME.message);
        end
    end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `dev/tests/test_process_eigen.m`
Expected: `ALL TESTS PASSED`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/process/functions/process_eigen.m dev/tests/test_process_eigen.m
git commit -m "feat(process): process_eigen pure orchestrator for eigen-domain analysis

Spatial-spectral twin of process_timefreq: no GetDescription (not panel-
registered); maps caller name -> bst_eigen Method, translates panel options,
runs bst_eigen per input (batch-robust). Test: test_process_eigen.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01C3Lf7KA3E5Y7maW6Knm2LM"
```

---

### Task 3: `process_eigenspectrum` — registered panel process

Create the panel-visible sibling: a `GetDescription` defining the Variant dropdown, measure radio, and Welch windowing options; a one-line `Run` delegating to `process_eigen`.

**Files:**
- Create: `toolbox/process/functions/process_eigenspectrum.m`
- Test: `dev/tests/test_process_eigenspectrum.m` (create)

**Interfaces:**
- Consumes: `process_eigen('Run', sProcess, sInputs)` from Task 2.
- Produces: registered process — `Comment='Eigenspectrum (spatial PSD)'`, `Category='Custom'`, `SubGroup='Frequency'`, `Index=484`, `InputTypes={'results'}`, `OutputTypes={'timefreq'}`. Options: `variant` (radio_linelabel, default `'Laplace-Beltrami'`), `measure` (radio_linelabel, default `'power'`), `timewindow`, `win_length` (default `{1,'s',[]}`), `win_overlap` (default `{50,'%',1}`), `win_std` (checkbox, default 0).

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_process_eigenspectrum.m`:

```matlab
function test_process_eigenspectrum()
% TEST_PROCESS_EIGENSPECTRUM: the registered process advertises the right contract and runs
% end-to-end through process_eigen -> bst_eigen, producing a timefreq_eigenspectrum node.
% End-to-end portion SKIPs if the live protocol lacks the fixture.
    nFail = 0; chk = @i_chk;

    % ---- GetDescription contract ----
    sProcess = process_eigenspectrum('GetDescription');
    nFail = nFail + chk('Comment',     strcmp(sProcess.Comment, 'Eigenspectrum (spatial PSD)'));
    nFail = nFail + chk('SubGroup',    strcmp(sProcess.SubGroup, 'Frequency'));
    nFail = nFail + chk('Index 484',   isequal(sProcess.Index, 484));
    nFail = nFail + chk('InputTypes',  isequal(sProcess.InputTypes, {'results'}));
    nFail = nFail + chk('OutputTypes', isequal(sProcess.OutputTypes, {'timefreq'}));
    nFail = nFail + chk('variant default LBO', strcmp(sProcess.options.variant.Value, 'Laplace-Beltrami'));
    nFail = nFail + chk('measure default power', strcmp(sProcess.options.measure.Value, 'power'));
    nFail = nFail + chk('win_length present', isfield(sProcess.options, 'win_length'));
    nFail = nFail + chk('win_overlap present', isfield(sProcess.options, 'win_overlap'));
    nFail = nFail + chk('win_std present', isfield(sProcess.options, 'win_std'));

    % ---- end-to-end via the registered process ----
    [ef, surf] = i_find_eigen('Laplace-Beltrami');
    if isempty(ef)
        fprintf('SKIP e2e: no LBO eigen node.\n');
    else
        rf = i_find_results(surf);
        if isempty(rf)
            fprintf('SKIP e2e: no constrained results on surface.\n');
        else
            sProcess.Function = @process_eigenspectrum;
            sInputs = bst_process('GetInputStruct', rf);
            Out = process_eigenspectrum('Run', sProcess, sInputs);
            ok = iscell(Out) && ~isempty(Out) && ischar(Out{1});
            nFail = nFail + chk('Run produces a node', ok);
            if ok
                nFail = nFail + chk('node prefix timefreq_eigenspectrum', ...
                    ~isempty(strfind(Out{1}, 'timefreq_eigenspectrum'))); %#ok<STREMP>
                R = in_bst_timefreq(Out{1}, 'Method');
                nFail = nFail + chk('node Method == spectrum', strcmp(R.Method, 'spectrum'));
                try, file_delete(file_fullpath(Out{1}), 1); db_reload_studies(sInputs.iStudy); catch; end
            end
        end
    end

    fprintf('\n==== test_process_eigenspectrum: %d failed ====\n', nFail);
    if nFail > 0; error('test_process_eigenspectrum FAILED'); end
    disp('ALL TESTS PASSED');
end

function r = i_chk(nm, cond)
    if cond; r = 0; fprintf('  ok   %s\n', nm); else; r = 1; fprintf('  FAIL %s\n', nm); end
end
function [ef, surf] = i_find_eigen(want)
    ef = ''; surf = '';
    PI = bst_get('ProtocolInfo'); if isempty(PI); return; end
    d = dir(fullfile(PI.SUBJECTS, '**', 'eigen_*.mat')); [~, o] = sort([d.datenum], 'descend');
    for i = o(:)'
        rel = strrep(fullfile(d(i).folder, d(i).name), [PI.SUBJECTS filesep], '');
        try
            m = in_bst_eigen(rel, 'Variant', 'ParentSurface');
            if strcmpi(m.Variant, want); ef = rel; surf = m.ParentSurface; return; end
        catch
        end
    end
end
function rf = i_find_results(surf)
    rf = '';
    sP = bst_get('ProtocolStudies'); if isempty(sP); return; end
    for s = 1:numel(sP.Study)
        R = sP.Study(s).Result;
        for r = 1:numel(R)
            try
                m = in_bst_results(R(r).FileName, 0, 'SurfaceFile', 'nComponents');
                if isfield(m, 'SurfaceFile') && ~isempty(m.SurfaceFile) ...
                        && file_compare(m.SurfaceFile, surf) ...
                        && (isempty(m.nComponents) || isequal(m.nComponents, 1))
                    rf = R(r).FileName; return;
                end
            catch
            end
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `dev/tests/test_process_eigenspectrum.m`
Expected: FAIL — `Undefined function ... 'process_eigenspectrum'` (file not created yet).

- [ ] **Step 3: Verify the Index is free**

Run: `grep -rhoE "sProcess\.Index\s*=\s*484\b" toolbox/process/functions/*.m`
Expected: no output (484 unused). If any output appears, pick the next free value in 483–489 and update the test's `Index 484` assertion to match.

- [ ] **Step 4: Implement `process_eigenspectrum.m`**

Create `toolbox/process/functions/process_eigenspectrum.m`:

```matlab
function varargout = process_eigenspectrum( varargin )
% PROCESS_EIGENSPECTRUM: Windowed eigenspectrum of surface-mapped sources (the spatial-spectral
% analogue of PROCESS_PSD / Welch). Projects each source map onto an operator eigenbasis and
% returns mode power vs spatial frequency (Freqs = sqrt(lambda)) as a timefreq_ node. The eigen_
% basis is resolved implicitly from the input's SurfaceFile; the operator family is chosen by the
% Variant option. Thin wrapper: delegates Run to the pure orchestrator process_eigen, exactly as
% process_psd delegates to process_timefreq.

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Diellor Basha, 2026

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'Eigenspectrum (spatial PSD)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Frequency';
    sProcess.Index       = 484;
    sProcess.Description = '';
    % Source maps only (eigenspectrum needs a surface basis)
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Option: operator family (eigen_ basis resolved implicitly from the surface)
    sProcess.options.variant.Comment = {'Laplace-Beltrami', 'Connection Laplacian', 'Dirac', 'Dirac-Face', 'Hodge-Face', 'Operator:'; ...
        'Laplace-Beltrami', 'Connection Laplacian', 'Dirac', 'Dirac-Face', 'Hodge-Face', ''};
    sProcess.options.variant.Type    = 'radio_linelabel';
    sProcess.options.variant.Value   = 'Laplace-Beltrami';
    % Option: spectrum measure
    sProcess.options.measure.Comment = {'Power', 'Magnitude', 'Measure:'; 'power', 'magnitude', ''};
    sProcess.options.measure.Type    = 'radio_linelabel';
    sProcess.options.measure.Value   = 'power';
    % Option: time window
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    % Option: Welch window length (empty => single full window)
    sProcess.options.win_length.Comment = 'Window length: ';
    sProcess.options.win_length.Type    = 'value';
    sProcess.options.win_length.Value   = {1, 's', []};
    % Option: Welch window overlap
    sProcess.options.win_overlap.Comment = 'Window overlap ratio: ';
    sProcess.options.win_overlap.Type    = 'value';
    sProcess.options.win_overlap.Value   = {50, '%', 1};
    % Option: aggregate (mean vs mean+std across windows)
    sProcess.options.win_std.Comment = '<HTML><FONT color="#a0a0a0">Also save the std across windows</FONT>';
    sProcess.options.win_std.Type    = 'checkbox';
    sProcess.options.win_std.Value   = 0;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    % Delegate to the shared eigen orchestrator (method inferred from this process name).
    OutputFiles = process_eigen('Run', sProcess, sInputs);
end
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `dev/tests/test_process_eigenspectrum.m`
Expected: `ALL TESTS PASSED`.

- [ ] **Step 6: Verify the process registers in the panel loader**

Run (MATLAB MCP `evaluate_matlab_code`):

```matlab
sProcesses = panel_process_select('GetAvailableProcesses');
isSpec = any(arrayfun(@(p) ~isempty(p.Function) && strcmp(func2str(p.Function), 'process_eigenspectrum'), sProcesses));
isOrch = any(arrayfun(@(p) ~isempty(p.Function) && strcmp(func2str(p.Function), 'process_eigen'), sProcesses));
fprintf('process_eigenspectrum registered: %d\n', isSpec);
fprintf('process_eigen registered (should be 0): %d\n', isOrch);
```

Expected: `process_eigenspectrum registered: 1` and `process_eigen registered (should be 0): 0`. (If `panel_process_select('GetAvailableProcesses')` signature differs in this build, instead confirm by opening the Process panel in the GUI and checking **Frequency › Eigenspectrum (spatial PSD)** appears and **Eigen** does not.)

- [ ] **Step 7: Commit**

```bash
git add toolbox/process/functions/process_eigenspectrum.m dev/tests/test_process_eigenspectrum.m
git commit -m "feat(process): process_eigenspectrum panel process (spatial PSD)

Registered Frequency-group process (Index 484, results->timefreq): Variant
dropdown (default LBO), power/magnitude, Welch windowing. Run delegates to
process_eigen. Test: test_process_eigenspectrum.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01C3Lf7KA3E5Y7maW6Knm2LM"
```

---

## Final verification

- [ ] Run the full new suite, all green:
  - `dev/tests/test_bst_eigen_resolve.m`
  - `dev/tests/test_process_eigen.m`
  - `dev/tests/test_process_eigenspectrum.m`
- [ ] Run the existing eigen regression to confirm no breakage from the `bst_eigen` change:
  - `dev/tests/test_connection_wavelet_scalogram.m`
  - `dev/tests/test_eigen_wavelet_method.m` (if present)
- [ ] Manual smoke test in the GUI: drag a `results` file into the Process panel, select **Frequency › Eigenspectrum (spatial PSD)**, Variant = Laplace-Beltrami, Run; confirm a `timefreq_eigenspectrum` node appears under the input and opens as a spectrum (power vs sqrt-λ).

## Notes for the implementer

- `bst_eigen` already reads `results` with full load (kernel auto-applied), restricts the time window, computes per hemisphere, and saves the `timefreq_eigenspectrum` node — Tasks 2–3 add no compute or save logic, only option plumbing.
- The `radio_linelabel` widget stores its selected **key string** directly in `.Value` (see `process_psd.m:54-58`); that is what `process_eigen` reads — do not index into the label cell.
- Leaving `win_length` empty yields a single full-window spectrum; this is the simplest path and is exercised by the Task 2 test.
- Memory: this is dev code in `research/code/`; do not touch `library/software/`. Don't delete protocol `.mat` files directly except via `file_delete` + `db_reload_studies` (used only for test-node cleanup).
