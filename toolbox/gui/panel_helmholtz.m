function varargout = panel_helmholtz(varargin)
% PANEL_HELMHOLTZ: Controls for the Helmholtz components view (view_helmholtz). Picks which
% Hodge component to display -- Total field |J| / Irrotational grad(phi) / Solenoidal
% curl(psi) / Harmonic h -- each swapping the quiver to that component's vector field and
% the cortex colormap to its scalar potential. Plus Show vectors, Show singular points
% (component-aware), a readout, and Close.
% Authors: Diellor Basha, 2026
    eval(macro_method);
end

function bstPanelNew = CreatePanel(hFig) %#ok<DEFNU>
    import javax.swing.*;
    panelName = 'Helmholtz';
    jPanelNew = gui_component('Panel');
    jOpt = JPanel(); jOpt.setLayout(BoxLayout(jOpt, BoxLayout.Y_AXIS));
    jSec = gui_river([2 2], [2 8 3 6], 'Helmholtz / Hodge components');

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
    jReadout = gui_component('label', jSec, 'br', '');
    jClose  = gui_component('button', jSec, 'br', 'Close');
    java_setcb(jClose, 'ActionPerformedCallback', @(h,e) OnClose(panelName));

    jOpt.add(jSec); jPanelNew.add(jOpt, java.awt.BorderLayout.NORTH);
    ctrl = struct('hFig',hFig, 'jVec',jVec, 'jMark',jMark, 'jReadout',jReadout);
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
