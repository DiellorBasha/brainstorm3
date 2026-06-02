# EigenModes Panel UI Refinement — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`).

**Goal:** Replace the EigenModes panel's two band-edge sliders + "Single" radio with an intuitive **center mode slider + width text field + shape radios**. The slider gains eigenmode-axis tick labels and `[ ]` window markers. The lever engine is untouched — only the panel's control surface changes.

**Architecture:** The panel translates Center + Width + Shape into existing lever verbs (`SetBand`, `SetWindowShape`). A pure `BandFromCenterWidth(center,width,K)` helper does the band math (unit-tested). `RefreshControls` reflects lever state back into the new controls and rebuilds the slider label table (axis labels + markers).

**Tech Stack:** MATLAB, Brainstorm Swing panels (`gui_river`/`gui_component`/`JSlider`/`JTextField`). Tests via MATLAB MCP `mcp__plugin_brainstorm-dev_MATLAB__evaluate_matlab_code`.

**Spec:** `docs/superpowers/specs/2026-06-01-eigenmode-panel-ui-design.md`

## Verified facts
- `panel_eigenmodes.m` uses `eval(macro_method)`, so any subfunction is callable as `panel_eigenmodes('Name', ...)`.
- Lever verbs (unchanged): `SetBand(lo,hi)` (clamps, sets `iCurrentMode=round((lo+hi)/2)`, recomputes weights, broadcasts), `SetWindowShape('single'|'box'|'tapered'|'gain')`, `GetWeights`, `GetCurrentMode`, `GetState`.
- State fields: `iCurrentMode`, `Band=[lo hi]`, `BandSpan`, `WindowShape`, `Weights`, `nModes` (=K_paired), `CacheEig` (has `.Values`, `.CompRank`).
- The merged `RefreshControls`/`UpdatePanel`/`SetSelectEnabled` reference `jSliderLo`,`jSliderHi`,`jLabelBand`,`jRadioSingle` — all to be replaced.
- The panel smoke test `dev/tests/test_eigenmode_lever_panel.m` asserts `jSliderLo`/`jSliderHi` exist — to be updated.
- Swing gotcha: `JSlider.setLabelTable` needs a `java.util.Hashtable` with **`java.lang.Integer` keys** mapping to `JComponent`s; `setPaintLabels(true)` to show them.

---

## Task 1: Pure band-from-center-width helper

**Files:** Modify `toolbox/gui/panel_eigenmodes.m`; create `dev/tests/test_eigenmode_panel_centerwidth.m`.

- [ ] **Step 1: failing test** — `dev/tests/test_eigenmode_panel_centerwidth.m`:
```matlab
function test_eigenmode_panel_centerwidth
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(fullfile(repoRoot, 'toolbox', 'gui'));
K = 600;
[lo,hi] = panel_eigenmodes('BandFromCenterWidth', 42, 13, K);
assert(lo==29 && hi==55, 'center 42 width 13 -> [29 55]');
[lo,hi] = panel_eigenmodes('BandFromCenterWidth', 42, 0, K);
assert(lo==42 && hi==42, 'width 0 -> single [42 42]');
[lo,hi] = panel_eigenmodes('BandFromCenterWidth', 3, 10, K);
assert(lo==1 && hi==13, 'clamp low edge to 1');
[lo,hi] = panel_eigenmodes('BandFromCenterWidth', 598, 10, K);
assert(lo==588 && hi==600, 'clamp high edge to K');
[lo,hi] = panel_eigenmodes('BandFromCenterWidth', 300, 9999, K);
assert(lo==1 && hi==600, 'width >= K -> full [1 K]');
fprintf('ALL TESTS PASSED: test_eigenmode_panel_centerwidth\n');
end
```

- [ ] **Step 2: run, verify FAIL** (`BandFromCenterWidth` missing).
  `run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_panel_centerwidth.m')`

- [ ] **Step 3: add the helper** to `panel_eigenmodes.m` (near `BuildWeights`):
```matlab
%% ===== PURE: clamped band from a center + half-width =====
function [lo, hi] = BandFromCenterWidth(center, width, K) %#ok<DEFNU>
    center = min(max(round(center), 1), K);
    width  = max(round(width), 0);
    lo = min(max(center - width, 1), K);
    hi = min(max(center + width, 1), K);
end
```

- [ ] **Step 4: run, verify PASS.**

- [ ] **Step 5: commit**
```bash
git add toolbox/gui/panel_eigenmodes.m dev/tests/test_eigenmode_panel_centerwidth.m
git commit -m "EigenModes panel: pure BandFromCenterWidth helper"
```

---

## Task 2: Rewrite the panel control surface

Replace the two band sliders + Single radio with a center slider + width field; add axis/marker label table; update callbacks, RefreshControls, UpdatePanel, SetSelectEnabled, and the smoke test.

**Files:** Modify `toolbox/gui/panel_eigenmodes.m`, `dev/tests/test_eigenmode_lever_panel.m`.

- [ ] **Step 1: update the smoke test** `dev/tests/test_eigenmode_lever_panel.m` to assert the NEW controls (this is the failing test for the GUI rewrite):
```matlab
function test_eigenmode_lever_panel
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end
bstPanel = panel_eigenmodes('CreatePanel');
assert(~isempty(bstPanel), 'CreatePanel must return a BstPanel');
ctrl = get(bstPanel, 'sControls');
assert(~isempty(ctrl), 'panel controls must be present');
assert(isfield(ctrl,'jSliderCenter') && isfield(ctrl,'jTextWidth') ...
    && isfield(ctrl,'jCheckActive') && isfield(ctrl,'jLabelReadout'), ...
    'expected center slider + width field + active + readout');
assert(isfield(ctrl,'jRadioBox') && isfield(ctrl,'jRadioTaper') && isfield(ctrl,'jRadioGauss'), ...
    'expected Box/Taper/Gauss shape radios');
fprintf('ALL TESTS PASSED: test_eigenmode_lever_panel\n');
end
```

- [ ] **Step 2: run it, verify FAIL** (old CreatePanel lacks `jSliderCenter`/`jTextWidth`/`jRadioGauss`).
  `run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_panel.m')`

- [ ] **Step 3: rewrite `CreatePanel`** in `panel_eigenmodes.m`. Replace the band-sliders/Single-radio block with a center slider + width field + Box/Taper/Gauss radios. Keep Active + readout. You MAY adjust Swing layout details so it builds, but keep these control names: `jSliderCenter, jTextWidth, jRadioBox, jRadioTaper, jRadioGauss, jCheckActive, jLabelReadout`, plus a `jLabelCenter`/`jLabelWidth`/`jLabelShape` for the text labels.
```matlab
function bstPanelNew = CreatePanel() %#ok<DEFNU>
    panelName = 'EigenModes';
    import java.awt.*;
    import javax.swing.*;
    jPanelNew = gui_river([2,2], [4,4,6,6], 'Spatial scale (eigenmodes)');

    % Active toggle + readout
    jCheckActive = gui_component('Checkbox', jPanelNew, '', 'Active', [], ...
        'Live-filter the displayed source map', @(h,ev)CheckActive_Callback());
    jLabelReadout = gui_component('Label', jPanelNew, 'hfill', '');
    jLabelReadout.setHorizontalAlignment(JLabel.RIGHT);

    % Center mode slider (thumb = center; axis labels + [ ] window markers via label table)
    jLabelCenter = gui_component('Label', jPanelNew, 'br', 'Center mode');
    jSliderCenter = JSlider(1, 100, 1);
    jSliderCenter.setPaintLabels(1);
    java_setcb(jSliderCenter, 'StateChangedCallback',   @(h,ev)CenterPreview_Callback());
    java_setcb(jSliderCenter, 'MouseReleasedCallback',  @(h,ev)Center_Callback());
    jPanelNew.add('br hfill', jSliderCenter);

    % Width text field
    jLabelWidth = gui_component('Label', jPanelNew, 'br', 'Width (+/- modes):');
    jTextWidth  = gui_component('Text', jPanelNew, '', '0', [], '', @(h,ev)Width_Callback());
    jTextWidth.setColumns(4);

    % Shape radios (greyed when width = 0)
    jGroup = ButtonGroup();
    jLabelShape = gui_component('Label', jPanelNew, 'br', 'Shape:');
    jRadioBox   = gui_component('Radio', jPanelNew, '',  'Box',   jGroup, '', @(h,ev)Shape_Callback('box'));
    jRadioTaper = gui_component('Radio', jPanelNew, '',  'Taper', jGroup, '', @(h,ev)Shape_Callback('tapered'));
    jRadioGauss = gui_component('Radio', jPanelNew, '',  'Gauss', jGroup, '', @(h,ev)Shape_Callback('gain'));
    jRadioBox.setSelected(1);

    ctrl = struct('jPanelTop',jPanelNew, 'jCheckActive',jCheckActive, 'jLabelReadout',jLabelReadout, ...
                  'jLabelCenter',jLabelCenter, 'jSliderCenter',jSliderCenter, ...
                  'jLabelWidth',jLabelWidth, 'jTextWidth',jTextWidth, ...
                  'jLabelShape',jLabelShape, 'jRadioBox',jRadioBox, 'jRadioTaper',jRadioTaper, 'jRadioGauss',jRadioGauss);
    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);
end
```
(Verify the `gui_component('Text', ...)` form against another panel that uses a text field — adapt the exact arg pattern; it must create an editable `JTextField` stored as `jTextWidth` whose callback fires `Width_Callback` on enter/focus-lost.)

- [ ] **Step 4: add the label-table builder** (axis labels + `[ ]` markers; Integer keys):
```matlab
function ApplySliderLabels(jSlider, K, lo, hi, width)
    import javax.swing.JLabel;
    tbl = java.util.Hashtable();
    pos = unique(round(linspace(1, K, 5)));
    for p = pos
        tbl.put(java.lang.Integer(p), JLabel(num2str(p)));
    end
    if (width > 0)
        tbl.put(java.lang.Integer(lo), JLabel('['));
        tbl.put(java.lang.Integer(hi), JLabel(']'));
    end
    jSlider.setLabelTable(tbl);
    jSlider.setPaintLabels(1);
end
```

- [ ] **Step 5: rewrite the callbacks.** Replace `Slider_Callback`/`Shape_Callback` and add `Center_Callback`/`CenterPreview_Callback`/`Width_Callback`:
```matlab
function CheckActive_Callback()
    ctrl = bst_get('PanelControls', 'EigenModes');
    SetActive(ctrl.jCheckActive.isSelected());
end

% Commit on slider release: SetBand(center +/- width)
function Center_Callback()
    ctrl = bst_get('PanelControls', 'EigenModes');
    st = GetState();
    c = ctrl.jSliderCenter.getValue();
    w = ReadWidth(ctrl, st.nModes);
    [lo, hi] = BandFromCenterWidth(c, w, st.nModes);
    SetBand(lo, hi);          % fires NotifyChanged -> RefreshControls
end

% Live preview while dragging (markers + readout only; no lever commit)
function CenterPreview_Callback()
    ctrl = bst_get('PanelControls', 'EigenModes');
    if ~ctrl.jSliderCenter.getValueIsAdjusting(), return; end
    st = GetState();
    c = ctrl.jSliderCenter.getValue();
    w = ReadWidth(ctrl, st.nModes);
    [lo, hi] = BandFromCenterWidth(c, w, st.nModes);
    ApplySliderLabels(ctrl.jSliderCenter, st.nModes, lo, hi, w);
end

function Width_Callback()
    ctrl = bst_get('PanelControls', 'EigenModes');
    st = GetState();
    w = ReadWidth(ctrl, st.nModes);
    c = ctrl.jSliderCenter.getValue();
    [lo, hi] = BandFromCenterWidth(c, w, st.nModes);
    if (w == 0)
        SetWindowShape('single');
    else
        SetWindowShape(CurrentShape(ctrl));
    end
    SetBand(lo, hi);
end

function Shape_Callback(shape)
    ctrl = bst_get('PanelControls', 'EigenModes');
    w = ReadWidth(ctrl, GetState().nModes);
    if (w == 0)
        SetWindowShape('single');
    else
        SetWindowShape(shape);
    end
end

% Read + clamp the width text field; falls back to current half-width on parse failure
function w = ReadWidth(ctrl, K)
    s = char(ctrl.jTextWidth.getText());
    w = str2double(s);
    if isnan(w) || ~isreal(w)
        st = GetState(); w = round((st.Band(2) - st.Band(1)) / 2);
    end
    w = min(max(round(w), 0), max(K-1,0));
end

function shape = CurrentShape(ctrl)
    if ctrl.jRadioTaper.isSelected(),     shape = 'tapered';
    elseif ctrl.jRadioGauss.isSelected(), shape = 'gain';
    else                                  shape = 'box';
    end
end
```

- [ ] **Step 6: rewrite `RefreshControls`** to populate the new controls from state:
```matlab
function RefreshControls()
    ctrl = bst_get('PanelControls', 'EigenModes');
    st = GetState();
    K  = max(st.nModes, 1);
    lo = st.Band(1); hi = st.Band(2);
    width = round((hi - lo) / 2);
    % Center slider
    ctrl.jSliderCenter.setMaximum(K);
    ctrl.jSliderCenter.setValue(st.iCurrentMode);
    ApplySliderLabels(ctrl.jSliderCenter, K, lo, hi, width);
    % Width field
    ctrl.jTextWidth.setText(num2str(width));
    % Shape radios: select current, grey out at width 0 (single)
    isBand = (width > 0) && ~strcmpi(st.WindowShape, 'single');
    ctrl.jRadioBox.setEnabled(isBand);
    ctrl.jRadioTaper.setEnabled(isBand);
    ctrl.jRadioGauss.setEnabled(isBand);
    switch st.WindowShape
        case 'tapered', ctrl.jRadioTaper.setSelected(true);
        case 'gain',    ctrl.jRadioGauss.setSelected(true);
        otherwise,      ctrl.jRadioBox.setSelected(true);
    end
    % Readout
    nKeep = nnz(st.Weights > 1e-6);
    lamStr = '';
    if ~isempty(st.CacheEig) && isfield(st.CacheEig,'Values') && ~isempty(st.CacheEig.Values) ...
            && isfield(st.CacheEig,'CompRank')
        cr = st.CacheEig.CompRank(:); vals = st.CacheEig.Values(:);
        iLo = find(cr == lo, 1); iHi = find(cr == hi, 1);
        if ~isempty(iLo) && ~isempty(iHi)
            lamStr = sprintf('    lambda %.3g - %.3g', vals(iLo), vals(iHi));
        end
    end
    if (width == 0)
        ctrl.jLabelReadout.setText(sprintf('Mode %d / %d%s', st.iCurrentMode, K, lamStr));
    else
        ctrl.jLabelReadout.setText(sprintf('Keeping %d modes%s', nKeep, lamStr));
    end
end
```

- [ ] **Step 7: update `UpdatePanel` + `SetSelectEnabled`** to the new control names. In `UpdatePanel`, change the slider-max line from `jSliderLo/jSliderHi.setMaximum(K)` to `ctrl.jSliderCenter.setMaximum(K);` (RefreshControls also sets it — fine). In `SetSelectEnabled`, replace the control-name list with:
```matlab
    sel = {'jSliderCenter','jTextWidth','jRadioBox','jRadioTaper','jRadioGauss'};
```
Search the file for any remaining reference to `jSliderLo`, `jSliderHi`, `jLabelBand`, `jRadioSingle` and remove/replace them (there should be none left after this task).

- [ ] **Step 8: verify**
```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_panel.m')
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_panel_centerwidth.m')
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_state.m')
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_lifecycle.m')
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_viewer_synth.m')
```
All `ALL TESTS PASSED`. `checkcode('.../toolbox/gui/panel_eigenmodes.m')` — no new errors. (If `gui_component('Text',...)` or the label-table API differs, fix by matching an existing panel that uses a `JTextField`/`setLabelTable`; keep the control names.)

- [ ] **Step 9: commit**
```bash
git add toolbox/gui/panel_eigenmodes.m dev/tests/test_eigenmode_lever_panel.m
git commit -m "EigenModes panel: center slider + width field + Box/Taper/Gauss; axis + window markers"
```

---

## Task 3: Live e2e (best-effort)

- [ ] **Step 1:** With a source map (or eigenmode viewer) open on a cortex with eigenmodes, via MCP: `gui_brainstorm('ShowToolTab','EigenModes')`, then drive the panel programmatically — set the center slider value + width field text and fire the callbacks — and confirm the lever `Band`/`Weights` update and the readout text reflects "Keeping N modes" / "Mode k". Confirm width 0 greys the shape radios. If GUI/data unavailable, document and rely on the headless suite.
- [ ] **Step 2: commit** any e2e helper if created; otherwise note completion.

## Done criteria
- Panel shows a center slider (with eigenmode axis labels + `[ ]` window markers), a width field, and Box/Taper/Gauss (greyed at width 0). No more lo/hi/Single confusion.
- Center + width map to `SetBand`/`SetWindowShape`; engine unchanged; all headless tests green.
