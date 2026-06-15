function varargout = panel_helmholtz(varargin)
% PANEL_HELMHOLTZ: Controls for the Helmholtz/vorticity view (view_helmholtz). Picks the
% scalar that colors the cortex -- Norm (|J|, native) / Curl.n (vorticity) / Divergence /
% Stream psi / Potential phi -- plus a Show vortex cores toggle, a core-count readout, and
% Close. The native source vector field is always shown.
% Authors: Diellor Basha, 2026
    eval(macro_method);
end

function bstPanelNew = CreatePanel(hFig) %#ok<DEFNU>
    import javax.swing.*;
    panelName = 'Helmholtz';
    jPanelNew = gui_component('Panel');
    jOpt = JPanel(); jOpt.setLayout(BoxLayout(jOpt, BoxLayout.Y_AXIS));
    jSec = gui_river([2 2], [2 8 3 6], 'Helmholtz / vorticity');

    gui_component('label', jSec, 'br', 'Colour the cortex by:');
    names    = {'Norm','Curl','Div','Psi','Phi'};
    labels   = {'Norm |J| (source)','Curl . n (vorticity)','Divergence','Stream function','Potential'};
    grp = ButtonGroup(); jRadio = javaArray('javax.swing.JRadioButton', numel(names));
    for i = 1:numel(names)
        jRadio(i) = gui_component('radio', jSec, 'br', labels{i});
        grp.add(jRadio(i));
        java_setcb(jRadio(i), 'ActionPerformedCallback', @(h,e) OnScalar(panelName, names{i}));
    end
    jRadio(1).setSelected(true);   % Norm: the native starting view
    jCores  = gui_component('checkbox', jSec, 'br', 'Show vortex cores');  jCores.setSelected(true);
    java_setcb(jCores, 'ActionPerformedCallback', @(h,e) OnCores(panelName));
    jReadout = gui_component('label', jSec, 'br', '');
    jClose  = gui_component('button', jSec, 'br', 'Close');
    java_setcb(jClose, 'ActionPerformedCallback', @(h,e) OnClose(panelName));

    jOpt.add(jSec); jPanelNew.add(jOpt, java.awt.BorderLayout.NORTH);
    ctrl = struct('hFig',hFig, 'jCores',jCores, 'jReadout',jReadout);
    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);
end

function OnScalar(panelName, scalarName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    view_helmholtz('SetScalar', ctrl.hFig, scalarName);
end
function OnCores(panelName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    view_helmholtz('SetCores', ctrl.hFig, ctrl.jCores.isSelected());
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
