function varargout = panel_helmholtz(varargin)
% PANEL_HELMHOLTZ: Controls for the Helmholtz components view (view_helmholtz). A Smoothing
% section (shared Dirac-eigenmode kernel + scale, on/off) low-passes the active frame before
% the decomposition; a Component radio (Total |J| / Irrotational / Solenoidal / Harmonic)
% swaps the quiver + colormap; a Marker threshold slider prunes weak singular points. Plus
% Show vectors, Show singular points, a readout, and Close.
% Authors: Diellor Basha, 2026
    eval(macro_method);
end

function bstPanelNew = CreatePanel(hFig, Lambda) %#ok<DEFNU>
    import javax.swing.*;
    panelName = 'Helmholtz';
    jPanelNew = gui_component('Panel');
    jOpt = JPanel(); jOpt.setLayout(BoxLayout(jOpt, BoxLayout.Y_AXIS));
    jSec = gui_river([2 2], [2 8 3 6], 'Helmholtz / Hodge components');

    % --- Smoothing (Dirac eigenmodes) ---
    gui_component('label', jSec, 'br', 'Smoothing (Dirac eigenmodes):');
    [keys, displays] = bst_eigfilter_panel('Kernels');
    jKernel = gui_component('combobox', jSec, 'br hfill', [], {displays}, [], [], []);
    iHeat = find(strcmp(keys,'heat'),1); if ~isempty(iHeat); jKernel.setSelectedIndex(iHeat-1); end
    jParams = gui_river([2 2], [0 2 0 2]);  jSec.add('br hfill', jParams);
    jSmoothOn = gui_component('checkbox', jSec, 'br', 'Smoothing on');
    bst_eigfilter_panel('BuildSliders', jParams, bst_eigfilter_panel('CurrentKernel', jKernel, keys), Lambda, @() OnSmooth(panelName));
    java_setcb(jKernel,   'ActionPerformedCallback', @(h,e) OnKernel(panelName));
    java_setcb(jSmoothOn, 'ActionPerformedCallback', @(h,e) OnSmooth(panelName));

    % --- Component ---
    gui_component('label', jSec, 'br', 'Component:');
    names  = {'Total','Irrot','Solen','Harm'};
    labels = {'Total field |J|','Irrotational (grad phi)','Solenoidal (curl psi)','Harmonic (h)'};
    grp = ButtonGroup(); jRadio = javaArray('javax.swing.JRadioButton', numel(names));
    for i = 1:numel(names)
        jRadio(i) = gui_component('radio', jSec, 'br', labels{i});
        grp.add(jRadio(i));
        java_setcb(jRadio(i), 'ActionPerformedCallback', @(h,e) OnComponent(panelName, names{i}));
    end
    jRadio(1).setSelected(true);   % Total: native start
    jVec  = gui_component('checkbox', jSec, 'br', 'Show vectors');           jVec.setSelected(true);
    jMark = gui_component('checkbox', jSec, 'br', 'Show singular points');   jMark.setSelected(true);
    java_setcb(jVec,  'ActionPerformedCallback', @(h,e) OnVectors(panelName));
    java_setcb(jMark, 'ActionPerformedCallback', @(h,e) OnMarkers(panelName));

    % --- Marker threshold (magnitude gate) ---
    gui_component('label', jSec, 'br', 'Marker threshold:');
    jThresh = JSlider(0, 100, 0);  jThresh.setPreferredSize(java_scaled('dimension', 120, 22));
    jSec.add('br hfill', jThresh);
    java_setcb(jThresh, 'StateChangedCallback', @(h,e) OnGate(panelName));

    jReadout = gui_component('label', jSec, 'br', '');
    jClose  = gui_component('button', jSec, 'br', 'Close');
    java_setcb(jClose, 'ActionPerformedCallback', @(h,e) OnClose(panelName));

    jOpt.add(jSec); jPanelNew.add(jOpt, java.awt.BorderLayout.NORTH);
    ctrl = struct('hFig',hFig, 'jVec',jVec, 'jMark',jMark, 'jReadout',jReadout, ...
                  'jKernel',jKernel, 'KernelKeys',{keys}, 'jParams',jParams, ...
                  'jSmoothOn',jSmoothOn, 'Lambda',Lambda, 'jThresh',jThresh);
    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);
end

function OnComponent(panelName, name) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    view_helmholtz('SetComponent', ctrl.hFig, name);
end
function OnVectors(panelName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    view_helmholtz('SetVectors', ctrl.hFig, ctrl.jVec.isSelected());
end
function OnMarkers(panelName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    view_helmholtz('SetMarkers', ctrl.hFig, ctrl.jMark.isSelected());
end
function OnKernel(panelName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    key = bst_eigfilter_panel('CurrentKernel', ctrl.jKernel, ctrl.KernelKeys);
    bst_eigfilter_panel('BuildSliders', ctrl.jParams, key, ctrl.Lambda, @() OnSmooth(panelName));
    OnSmooth(panelName);
end
function OnSmooth(panelName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    name   = bst_eigfilter_panel('CurrentKernel', ctrl.jKernel, ctrl.KernelKeys);
    params = bst_eigfilter_panel('ReadParams', ctrl.jParams, ctrl.Lambda);
    view_helmholtz('SetSmoothing', ctrl.hFig, ctrl.jSmoothOn.isSelected(), name, params);
end
function OnGate(panelName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    if ctrl.jThresh.getValueIsAdjusting(); return; end
    view_helmholtz('SetGate', ctrl.hFig, double(ctrl.jThresh.getValue())/100);
end
function tf = i_valid(ctrl)
    tf = ~isempty(ctrl) && isfield(ctrl,'hFig') && ~isempty(ctrl.hFig) && all(ishandle(ctrl.hFig));
end
function SetReadout(text) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'Helmholtz'); if isempty(ctrl); return; end
    ctrl.jReadout.setText(text);
end
function OnClose(panelName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName);
    if i_valid(ctrl); view_helmholtz('Close', ctrl.hFig); else; gui_hide(panelName); end
end
