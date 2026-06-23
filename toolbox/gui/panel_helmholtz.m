function varargout = panel_helmholtz(varargin)
% PANEL_HELMHOLTZ: Compact controls for the Helmholtz components view (view_helmholtz).
% Two sections: a Smooth toggle + Dirac-eigenmode kernel (low-passes the active frame
% before the decomposition), and a bordered Component section with state buttons --
% Potential (irrotational phi) / Stream (solenoidal psi). With both off, the default view
% is shown: the total field |J| coloring + the total-field vector quiver.
%
% NOTE: the Derivative, Show-singular-points, persistence gate, trajectory tracking and the
% (genus-0 trivial) Harmonic component were removed -- that machinery is being rebuilt on
% the atom system (bst_dynamics). This panel is designed to fold into a future
% panel_bst_dynamics.
% Authors: Diellor Basha, 2026
    eval(macro_method);
end

function bstPanelNew = CreatePanel(hFig, Lambda) %#ok<DEFNU>
    import java.awt.*;
    import javax.swing.*;
    panelName = 'Helmholtz';
    BTN_W = java_scaled('value', 64);
    BTN_H = java_scaled('value', 22);

    jPanelNew = gui_component('Panel');
    jOpt = JPanel();  jOpt.setLayout(BoxLayout(jOpt, BoxLayout.Y_AXIS));

    % ===== Smooth: checkbox to the LEFT of the kernel dropdown (no title) =====
    jSmoothSec = gui_river([2 2], [2 6 2 6]);
    jSmoothOn = gui_component('checkbox', jSmoothSec, '', 'Smooth', [], 'Low-pass the frame in the Dirac eigenbasis');
    [keys, displays] = panel_eigenfilter_design('Kernels');
    jKernel = gui_component('combobox', jSmoothSec, 'hfill', [], {displays}, [], [], []);
    iHeat = find(strcmp(keys,'heat'),1);  if ~isempty(iHeat); jKernel.setSelectedIndex(iHeat-1); end
    jParams = gui_river([2 2], [0 2 0 2]);  jSmoothSec.add('br hfill', jParams);
    panel_eigenfilter_design('BuildSliders', jParams, panel_eigenfilter_design('CurrentKernel', jKernel, keys), Lambda, @() OnSmooth(panelName));
    java_setcb(jKernel,   'ActionPerformedCallback', @(h,e) OnKernel(panelName));
    java_setcb(jSmoothOn, 'ActionPerformedCallback', @(h,e) OnSmooth(panelName));
    jOpt.add(jSmoothSec);

    % ===== Component: bordered, state toggle buttons (both off => total field + vectors) =====
    jCompSec = gui_river([0 4], [2 6 4 6], 'Component');
    jPot = gui_component('toggle', jCompSec, 'center', 'Potential', {Insets(0,0,0,0), Dimension(BTN_W, BTN_H)}, 'Irrotational potential \Phi (sources / sinks)');
    jStr = gui_component('toggle', jCompSec, '',       'Stream',    {Insets(0,0,0,0), Dimension(BTN_W, BTN_H)}, 'Solenoidal stream \Psi (vortices)');
    java_setcb(jPot, 'ActionPerformedCallback', @(h,e) OnComp(panelName, 'Irrot'));
    java_setcb(jStr, 'ActionPerformedCallback', @(h,e) OnComp(panelName, 'Solen'));
    jOpt.add(jCompSec);

    % ===== readout + close =====
    jFoot = gui_river([2 2], [2 6 2 6]);
    jReadout = gui_component('label', jFoot, 'hfill', '');
    jClose   = gui_component('button', jFoot, 'br right', 'Close', {Insets(2,8,2,8), Dimension(java_scaled('value',58), BTN_H)});
    java_setcb(jClose, 'ActionPerformedCallback', @(h,e) OnClose(panelName));
    jOpt.add(jFoot);

    jPanelNew.add(jOpt, BorderLayout.NORTH);
    ctrl = struct('hFig',hFig, 'jPot',jPot, 'jStr',jStr, 'jReadout',jReadout, ...
                  'jKernel',jKernel, 'KernelKeys',{keys}, 'jParams',jParams, ...
                  'jSmoothOn',jSmoothOn, 'Lambda',Lambda);
    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);
end

%% ===== Component state buttons (mutually exclusive; all off => Total) =====
function OnComp(panelName, which) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName); if ~i_valid(ctrl); return; end
    if strcmp(which, 'Irrot')
        if ctrl.jPot.isSelected(); ctrl.jStr.setSelected(false); name = 'Irrot'; else; name = 'Total'; end
    else  % Solen
        if ctrl.jStr.isSelected(); ctrl.jPot.setSelected(false); name = 'Solen'; else; name = 'Total'; end
    end
    view_helmholtz('SetComponent', ctrl.hFig, name);
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
