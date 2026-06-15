function varargout = bst_eigfilter_panel(varargin)
% BST_EIGFILTER_PANEL: Shared "Filter kernel" UI section (kernel dropdown + mode-index
% scale sliders), reused by panel_wavelet_designer and panel_spatial_filter.
%
% API (dispatched via macro_method):
%   [keys, displays] = bst_eigfilter_panel('Kernels')
%   key   = bst_eigfilter_panel('CurrentKernel', jKernel, keys)
%   bst_eigfilter_panel('BuildSliders', jParams, kernelKey, Lambda, onSettle)
%   names = bst_eigfilter_panel('ParamNames', jParams)
%   params = bst_eigfilter_panel('ReadParams', jParams, Lambda)
%
% The scale sliders live in mode-index space (1..K): mode k maps to eigenvalue
% Lambda(k), then to the kernel's scale parameter (t = 1/lambda_k for heat/mexhat/dog,
% beta = lambda_k for tikhonov). dog's t1<t2 is enforced in ReadParams. onSettle is a
% no-arg handle called when a slider drag settles (the owner's live recompute).
%
% Authors: Diellor Basha, 2026
    eval(macro_method);
end

function [keys, displays] = Kernels() %#ok<DEFNU>
    keys = {'mexhat','dog','heat','inverse_heat','tikhonov'};
    displays = cell(1, numel(keys));
    for i = 1:numel(keys)
        try
            m = bst_eigfilter_kernel('info', keys{i});  displays{i} = m.display;
        catch
            displays{i} = keys{i};
        end
    end
end

function key = CurrentKernel(jKernel, keys) %#ok<DEFNU>
    idx = max(1, min(numel(keys), jKernel.getSelectedIndex() + 1));
    key = keys{idx};
end

function BuildSliders(jParams, kernelKey, Lambda, onSettle) %#ok<DEFNU>
    meta = bst_eigfilter_kernel('info', kernelKey);
    K = numel(Lambda);
    jParams.removeAll();
    pf = fieldnames(meta.params);
    nP = numel(pf);
    for i = 1:nP
        nm = pf{i};
        % stagger defaults so multi-scale kernels (dog: t1,t2) start at distinct modes
        defMode = max(1, min(K, round(K * i/(nP+1))));
        [js, jTitle] = i_labeled_slider(jParams, ...
            sprintf('%s: mode %d', i_param_label(nm), defMode), 'coarse', 'fine', 1, K, defMode);
        jParams.putClientProperty(['slider_' nm], js);
        jParams.putClientProperty(['title_'  nm], jTitle);
        java_setcb(js, 'StateChangedCallback', @(h,e) i_slider_changed(jParams, nm, Lambda, h, onSettle));
    end
    jParams.putClientProperty('ParamNames', strjoin(pf(:).', ','));
    jParams.revalidate(); jParams.repaint();
end

function names = ParamNames(jParams) %#ok<DEFNU>
    s = jParams.getClientProperty('ParamNames');
    if isempty(s); names = {}; else; names = strsplit(char(s), ','); end
end

function params = ReadParams(jParams, Lambda) %#ok<DEFNU>
    params = struct();
    K = numel(Lambda);
    names = ParamNames(jParams);
    for i = 1:numel(names)
        nm = names{i};
        js = jParams.getClientProperty(['slider_' nm]);
        if isempty(js); continue; end
        k = max(1, min(K, double(js.getValue())));
        params.(nm) = i_param_value(nm, Lambda(k));
    end
    % dog requires t1 < t2 (the two sliders are independent): order + separate
    if isfield(params,'t1') && isfield(params,'t2')
        lo = min(params.t1, params.t2);  hi = max(params.t1, params.t2);
        if hi <= lo * (1 + 1e-3); hi = lo * 1.5; end
        params.t1 = lo;  params.t2 = hi;
    end
end

%% ===== internal =====
function i_slider_changed(jParams, name, Lambda, js, onSettle)
    jt = jParams.getClientProperty(['title_' name]);
    if ~isempty(jt)
        k = max(1, min(numel(Lambda), double(js.getValue())));
        jt.setText(sprintf('%s: mode %d', i_param_label(name), k));
    end
    if ~js.getValueIsAdjusting() && ~isempty(onSettle); onSettle(); end
end

function v = i_param_value(name, lamk)
    if strcmpi(name, 'beta'); v = max(lamk, eps); else; v = 1 ./ max(lamk, eps); end
end

function lab = i_param_label(name)
    switch lower(name)
        case 't',    lab = 'Scale';
        case 't1',   lab = 'Band edge 1';
        case 't2',   lab = 'Band edge 2';
        case 'beta', lab = 'Scale';
        otherwise,   lab = name;
    end
end

function js = i_slider(jParent, constraints, mn, mx, val)
    import javax.swing.*;
    js = JSlider(mn, mx, val);
    js.setPreferredSize(java_scaled('dimension', 40, 22));
    jParent.add(constraints, js);
end

function [js, jTitle] = i_labeled_slider(jParent, titleText, loLabel, hiLabel, mn, mx, val)
    jTitle = gui_component('label', jParent, 'br', titleText);
    gui_component('label', jParent, 'br', loLabel);
    js = i_slider(jParent, 'hfill', mn, mx, val);
    gui_component('label', jParent, '', hiLabel);
end
