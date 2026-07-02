# RESUME — SP2a (round out the atom designer) + program state (2026-07-02)

> **Purpose:** hand off to a FRESH session. Read top-to-bottom. The immediate task is **SP2a**: the design
> is DONE and spec'd; the next step is `writing-plans` → subagent-driven execution. Everything below the
> SP2a section is durable program context (SP1 shipped, SP2 design, concerns backlog, execution policy).

---

## 0. Where we are (one paragraph)

**SP1 (metadata-driven eigen stack + Dirac-Connectome factory operator) is SHIPPED + PUSHED** to
`origin/development` (HEAD `d4a3dd7d`, 2026-07-02). **SP2 (one design surface) is fully DESIGNED**; it
decomposes into **SP2a — round out the standalone designer** (spec written, ready to plan) and **SP2b —
excise the design half from the Dynamics panel + embed the shared designer** (next spec, not yet written).
The immediate job is to PLAN + BUILD SP2a.

## 1. Resume recipe for SP2a

1. Read the SP2a spec: **`dev/2026-07-02-designer-round-out-design.md`** (approved in brainstorm; do a spec
   review pass with the user first if they want changes).
2. Invoke **`superpowers:writing-plans`** to turn it into a task-by-task TDD plan.
3. Execute via **`superpowers:subagent-driven-development`** under the EXECUTION POLICY in §5.
4. Branch: create `feature/designer-round-out` off `development` (in-place, NOT a worktree — the live MATLAB
   session needs THIS checkout on its path for the controller live pass; see §5/§6).

## 2. SP2a — what it is (design locked)

**Goal:** make the standalone `view_atom_designer` (+ `panel_atom_designer`) support the FULL operator range,
data-free. Enabled by SP1 (`bst_eigen('Axes', variant)` now routes any operator by `field_type`;
Dirac-Connectome is a factory operator).

**The Atom/Fiber taxonomy reorg (the foundation — user was firm on this):**
- An **atom** is ONE supercategory; scalar-vs-vector is its **fiber** (a descriptor), NEVER a separate verb.
  There is **NO `AtomVec`** (rejected by name — it wrongly implies a vector atom is a different thing).
- **`bst_eigenfilter('Fiber', ax)` is RETAINED** — it answers "what fiber?" (`scalar`/`complex-tangent`/
  `quaternion`, driven by `field_type`).
- **`bst_eigenfilter('Atom', …)` is the SINGLE realiser** — extend it to return the atom already decoded for
  its fiber: `[W, gv, V3, isSigned]` (`V3=[]` for scalar; imag-3-vec for quaternion; `a·e1+b·e2` for tangent).
  It consults `Fiber` internally. The quaternion/tangent decode currently in `panel_bst_dynamics`
  `i_atom_realise_core` (`:679–697`) MOVES INTO `Atom`; `i_atom_realise_core` is then updated to CONSUME
  `Atom`'s `V3` (thin). Other `[W,gv]` callers unaffected (backward-compatible extra outputs).

**Round out the designer:** operator selector → 4 (Geometric/Connectomic/Dirac/Dirac-Connectome; Connection
Laplacian OUT — unpersisted tangent frame); fiber-driven render (scalar density, or Dirac magnitude+quivers);
app-side seed-direction control (Dirac presets, hidden for scalar); **default direction relocated** from the
panel's local `i_atom_default_dir` into a NEW shared `toolbox/gui/bst_atom_default_dir.m` (app-side, not the
library — user's rule: orchestrators carry no defaults). Seed-pick uses the designer's existing
`WaveletDesignerPick` (clean; not the geodesic Region tool — this also fixes concern C6).

**Components:** `bst_eigenfilter.m` (Atom extended, Fiber unchanged); `view_atom_designer.m` (i_eval_atom uses
Atom's V3; render vector); `panel_atom_designer.m` (4-operator selector + direction control); NEW
`bst_atom_default_dir.m`; `panel_bst_dynamics.m` (2 small byte-equivalent touches: call the shared default-dir
helper; i_atom_realise_core consumes Atom's V3). **Data-free — SP2a touches NOTHING about the inverse kernel /
sensor data.**

**Test surface:** `sub-MTL0005/tess_cortex_pial_low.mat` (the source `sub-MTL0002`'s R.SurfaceFile), protocol
`preventad`. Unit: `Atom` fiber-decoded output == the panel's old `i_atom_realise_core` for scalar + Dirac
(byte-equiv). Live: open the designer, each operator renders (scalar density / Dirac magnitude+quivers),
direction control reorients Dirac.

## 3. SP2b (the NEXT spec after SP2a — DESIGN already decided, not yet written)

Excise `panel_bst_dynamics`'s design half + **embed `panel_atom_designer` as a shared docked component in the
Dynamics session** (excise model **A**, user-chosen). Shared state = the filterbank/atom-spec table
(`db_template('atomgroup')`). The embedded designer paints its impulse on the Dynamics cortex figure and edits
`st.T.Groups(curAtom)`. Dynamics keeps only the atom LIST + real-source verbs (Apply/Analyze/Localize-bands/
Helmholtz/Measure). Subsumes concerns C6 (seed-pick), part of C7 (info pane). Write this spec via brainstorming
after SP2a lands.

## 4. Program artifacts (the map)

- **SP2a spec:** `dev/2026-07-02-designer-round-out-design.md` ← the thing to plan next.
- **SP1 (done):** spec `dev/2026-07-02-metadata-driven-eigen-stack-design.md`, plan
  `dev/2026-07-02-metadata-driven-eigen-stack-plan.md` (both committed, `8d3f8e43`/`0f814881`).
- **Concerns backlog + north-star:** `dev/2026-07-02-dynamics-panel-concerns.md` — N1 (operator-composition
  framework), N1.a (ik_/naming taxonomy), C1–C8 (C3 = the designer extraction SP2 realizes; C1/C2 operator
  lifecycle; C4/C5 differential-mode-kernel; C6 seed-pick; C7 info pane; C8 kernel menu). ⚑ UNCOMMITTED.
- **GUI workflow + tests:** `dev/2026-07-02-dynamics-gui-workflow-and-tests.md`. ⚑ UNCOMMITTED.
- **SDD ledger (durable progress):** `.superpowers/sdd/progress.md` (git-ignored scratch). Records SP1
  task-by-task + the historical ledger pointer (prior sub-projects A–D etc., all merged; full backup path in
  the ledger). ⚑ Check it after any compaction; trust it + `git log` over recollection.
- **Memories:** `operator-composition-northstar`, `atom-designer-architecture` (has the C3 direction +
  SP1-shipped note), `dirac-connectome-operator`, `dynamics-design-preview-modes`.

## 5. ⚑⚑ EXECUTION POLICY (hard-won across prior sessions — HONOR IT)

- **Implementer subagents do STATIC CHECKS ONLY:** `check_matlab_code` (MCP, read-only) + grep + write test
  `.m` files + commit. They **DO NOT run/eval MATLAB or launch Brainstorm** — it hung/WIPED the live session
  in prior sessions.
- **The CONTROLLER runs the consolidated LIVE/headless test pass** in the established session
  (`preventad`, `sub-MTL0002` — do NOT `clear`). `check_matlab_code` misses runtime errors; the controller
  pass is the gate. Headless unit tests via `matlab -batch "addpath('dev/tests'); test_X"` (separate process,
  safe). Live/GUI via the MATLAB-MCP `evaluate_matlab_code` (controller only).
- **Controller adjudicates** reviewer findings against "preserve prior behavior" before dispatching fixes.
- **Models:** transcription-with-complete-code = haiku; surgical edits in large files = sonnet; reviewers
  scaled to diff; **final whole-branch review = opus** (it earned its keep on SP1 — caught a stale-node
  regression the synthetic tests missed).
- **Tests:** `dev/tests/test_*.m`, assert-based; surface-dependent ones use `getenv('BST_TEST_SURF')` +
  assumeTrue skip; ⚑ `bst_eigen('Axes', …)` REQUIRES `SampleRate` + `TimeWindow` fields (BuildTimeAxis).
- ⚑ Common controller test-harness fixes: a test's second function must be CALLED by the primary (else it's an
  uncallable local → vacuous "PASS"); watch loop-var clobbering.

## 6. Environment / gotchas

- Mac Mini, ~17 GB RAM. MATLAB R2023b + MATLAB-MCP (`brainstorm-dev` plugin). Protocol `preventad`.
- ⚑ The MATLAB-MCP `run_matlab_file` tool RESET the path and KILLED the live Brainstorm session once this
  session — use `evaluate_matlab_code` (keeps the live session), not `run_matlab_file`, for live work.
- Never `clear` in the live session (wipes GlobalData); edited `.m` auto-reload; use `rehash`.
- Dev→feature branch model; `development` is the working branch; production is `library/software` (main only).
- ⚑ **Uncommitted docs** (in the working tree, not yet committed — commit-when-asked was the session rule):
  the SP2a spec, the concerns doc, the workflow doc, and this RESUME. Recommend committing them for
  durability before/at the start of the next session (`docs(...)` prefix on `development`).

## 7. Recommended first action next session
Confirm with the user, then: commit the uncommitted design docs → (optional spec-review pass on the SP2a
spec) → `superpowers:writing-plans` on `dev/2026-07-02-designer-round-out-design.md` → subagent-driven
execution under §5. SP2b follows once SP2a is merged.
