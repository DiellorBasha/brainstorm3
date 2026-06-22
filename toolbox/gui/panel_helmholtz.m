function varargout = panel_helmholtz(varargin)
% PANEL_HELMHOLTZ: Controls for the Helmholtz components view (view_helmholtz). A Smoothing
% section (shared Dirac-eigenmode kernel + scale, on/off) low-passes the active frame before
% the decomposition; a Component radio (Total |J| / Irrotational / Solenoidal / Harmonic)
% swaps the cortex colormap + markers (the quiver stays the total field); a persistence gate
% slider prunes low-persistence singular points. Plus Show vectors, Show singular points, a readout, and
% Close.
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
    [keys, displays] = panel_eigenfilter_design('Kernels');
    jKernel = gui_component('combobox', jSec, 'br hfill', [], {displays}, [], [], []);
    iHeat = find(strcmp(keys,'heat'),1); if ~isempty(iHeat); jKernel.setSelectedIndex(iHeat-1); end
    jParams = gui_river([2 2], [0 2 0 2]);  jSec.add('br hfill', jParams);
    jSmoothOn = gui_component('checkbox', jSec, 'br', 'Smoothing on');
    panel_eigenfilter_design('BuildSliders', jParams, panel_eigenfilter_design('CurrentKernel', jKernel, keys), Lambda, @() OnSmooth(panelName));
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

    % --- Derivative order ---
    gui_component('label', jSec, 'br', 'Derivative:');
    dnames = {'Field','Velocity','Acceleration'};
    grpD = ButtonGroup(); jDeriv = javaArray('javax.swing.JRadioButton', numel(dnames));
    for i = 1:numel(dnames)
        jDeriv(i) = gui_component('radio', jSec, 'br', dnames{i});
        grpD.add(jDeriv(i));
        java_setcb(jDeriv(i), 'ActionPerformedCallback', @(h,e) OnDeriv(panelName, i-1));
    end
    jDeriv(1).setSelected(true);   % Field (order 0)
    jVec  = gui_component('checkbox', jSec, 'br', 'Show vectors');           jVec.setSelected(true);
    jMark = gui_component('checkbox', jSec, 'br', 'Show singular points');   jMark.setSelected(true);
    java_setcb(jVec,  'ActionPerformedCallback', @(h,e) OnVectors(panelName));
    java_setcb(jMark, 'ActionPerformedCallback', @(h,e) OnMarkers(panelName));

    % --- Marker threshold (persistence gate) ---
    gui_component('label', jSec, 'br', 'Persistence gate:');
    jThresh = JSlider(0, 100, 0);  jThresh.setPreferredSize(java_scaled('dimension', 120, 22));
    jSec.add('br hfill', jThresh);
    java_setcb(jThresh, 'StateChangedCallback', @(h,e) OnGate(panelName));

    jTrack = gui_component('checkbox', jSec, 'br', 'Track trajectory');
    java_setcb(jTrack, 'ActionPerformedCallback', @(h,e) OnTrack(panelName));
    jReadout = gui_component('label', jSec, 'br', '');
    jClose  = gui_component('button', jSec, 'br', 'Close');
    java_setcb(jClose, 'ActionPerformedCallback', @(h,e) OnClose(panelName));

    jOpt.add(jSec); jPanelNew.add(jOpt, java.awt.BorderLayout.NORTH);
    ctrl = struct('hFig',hFig, 'jVec',jVec, 'jMark',jMark, 'jReadout',jReadout, ...
                  'jKernel',jKernel, 'KernelKeys',{keys}, 'jParams',jParams, ...
                  'jSmoothOn',jSmoothOn, 'Lambda',Lambda, 'jThresh',jThresh, 'jTrack',jTrack, 'jDeriv',jDeriv);
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
function OnTrack(panelName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    view_helmholtz('SetTrack', ctrl.hFig, ctrl.jTrack.isSelected());
end
function OnDeriv(panelName, order) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    view_helmholtz('SetDeriv', ctrl.hFig, order);
end
function OnKernel(panelName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    key = panel_eigenfilter_design('CurrentKernel', ctrl.jKernel, ctrl.KernelKeys);
    panel_eigenfilter_design('BuildSliders', ctrl.jParams, key, ctrl.Lambda, @() OnSmooth(panelName));
    OnSmooth(panelName);
end
function OnSmooth(panelName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    name   = panel_eigenfilter_design('CurrentKernel', ctrl.jKernel, ctrl.KernelKeys);
    params = panel_eigenfilter_design('ReadParams', ctrl.jParams, ctrl.Lambda);
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
