function varargout = panel_eigenfilter_design(varargin)
% PANEL_EIGENFILTER_DESIGN: Shared eigenfilter "design" UI section (kernel dropdown +
% mode-index scale sliders), reused by panel_eigenfilter_options.
% Also renders the filter's spectral response h(lambda) (DrawResponse).
%
% API (dispatched via macro_method):
%   [keys, displays] = panel_eigenfilter_design('Kernels')
%   key   = panel_eigenfilter_design('CurrentKernel', jKernel, keys)
%   panel_eigenfilter_design('BuildSliders', jParams, kernelKey, Lambda, onSettle)
%   names = panel_eigenfilter_design('ParamNames', jParams)
%   params = panel_eigenfilter_design('ReadParams', jParams, Lambda)
%   panel_eigenfilter_design('DrawResponse', hAxes, kernelName, params, Lambda)
%
% The scale sliders live in mode-index space (1..K): mode k maps to eigenvalue
% Lambda(k), then to the kernel's scale parameter (t = 1/lambda_k for heat/mexhat/diffgauss,
% beta = lambda_k for tikhonov). diffgauss's t1<t2 is enforced in ReadParams. onSettle is a
% no-arg handle called when a slider drag settles (the owner's live recompute).
%
% Authors: Diellor Basha, 2026
    eval(macro_method);
end

function [keys, displays] = Kernels() %#ok<DEFNU>
    keys = {'mexhat','diffgauss','heat','inverse_heat','tikhonov'};
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
        % stagger defaults so multi-scale kernels (diffgauss: t1,t2) start at distinct modes
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

function DrawResponse(hAxes, kernelName, params, Lambda) %#ok<DEFNU>
% Plot the filter's spectral response h(lambda) (gain vs mode index) into hAxes.
    Lambda = double(Lambda(:));
    if isempty(Lambda)
        cla(hAxes);
        return;
    end
    g = bst_eigfilter_kernel(kernelName, params);
    if iscell(g)            % sliders yield single-scale params; guard the bank case
        g = g{1};
    end
    h = bst_eigfilter_evaluate(g, Lambda);
    k = (1:numel(Lambda))';
    plot(hAxes, k, h, 'LineWidth', 2);
    set(hAxes, 'XGrid', 'on', 'YGrid', 'on');
    xlabel(hAxes, 'Mode index k');
    ylabel(hAxes, 'Gain h(\lambda)');
    title(hAxes, sprintf('Spectral response: %s', kernelName), 'Interpreter', 'none');
    xlim(hAxes, [1, max(2, numel(Lambda))]);
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
    % diffgauss requires t1 < t2 (the two sliders are independent): order + separate
    if isfield(params,'t1') && isfield(params,'t2')
        lo = min(params.t1, params.t2);  hi = max(params.t1, params.t2);
        if hi <= lo * (1 + 1e-3); hi = lo * 1.5; end
        params.t1 = lo;  params.t2 = hi;
    end
end

%% ===== atom tool: physical-units sliders driven by bst_eigfilter_controls =====
% The atom designer + panel share one param spec (bst_eigfilter_controls). These render that spec
% (physical units: mm / Hz / m·s / nu — NOT mode index) as JSliders, for ALL kernels, and read the
% three slot values back as [s1 s2 s3] for bst_eigfilter_controls('ToKernel', kernel, vals, lmax).

function [keys, displays] = AtomKernels() %#ok<DEFNU>
% All eigfilter kernels, ordered (diffusion, then static, ts, js), with display names.
    keys = bst_eigfilter_kernel('list');
    rank = zeros(numel(keys), 1);  displays = cell(1, numel(keys));
    for i = 1:numel(keys)
        try, m = bst_eigfilter_kernel('info', keys{i}); catch, m = struct('display',keys{i},'domain','static'); end
        if isfield(m,'display') && ~isempty(m.display), displays{i} = m.display; else, displays{i} = keys{i}; end
        dom = ''; if isfield(m,'domain'), dom = m.domain; end
        switch lower(dom), case 'ts', rank(i) = 2; case 'js', rank(i) = 3; otherwise, rank(i) = 1; end
    end
    rank(strcmp(keys(:),'diffusion')) = 0;                 % diffusion first (the designer default)
    [~, ord] = sortrows([rank, (1:numel(keys))']);
    keys = keys(ord);  displays = displays(ord);
end

function BuildAtomSliders(jParams, kernel, bounds, onSettle) %#ok<DEFNU>
% Render kernel's bst_eigfilter_controls('Sliders') spec as physical-units JSliders (0..1000 -> [lo,hi]).
    S = bst_eigfilter_controls('Sliders', kernel, bounds);
    jParams.removeAll();
    for i = 1:3                                                 % clear stale slot handles (removeAll keeps client properties)
        jParams.putClientProperty(sprintf('atomslot_%d', i), []);
    end
    for i = 1:numel(S)
        if isempty(S(i).label), continue; end              % disabled slot -> no row
        rng = max(S(i).hi - S(i).lo, eps);
        pos = round(1000 * (S(i).def - S(i).lo) / rng);
        [js, jTitle] = i_labeled_slider(jParams, sprintf(['%s: ' S(i).fmt], S(i).label, S(i).def), ...
                                        'lo', 'hi', 0, 1000, max(0, min(1000, pos)));
        jParams.putClientProperty(sprintf('atomslot_%d', i), js);
        jParams.putClientProperty(sprintf('atomttl_%d', i),  jTitle);
        jParams.putClientProperty(sprintf('atomlo_%d', i),   S(i).lo);
        jParams.putClientProperty(sprintf('atomhi_%d', i),   S(i).hi);
        jParams.putClientProperty(sprintf('atomfmt_%d', i),  S(i).fmt);
        jParams.putClientProperty(sprintf('atomlbl_%d', i),  S(i).label);
        java_setcb(js, 'StateChangedCallback', @(h,e) i_atom_slider_changed(jParams, i, h, onSettle));
    end
    jParams.revalidate(); jParams.repaint();
end

function vals = ReadAtomVals(jParams) %#ok<DEFNU>
% Read the three slots as physical values [s1 s2 s3]; disabled slots read 0.
    vals = [0 0 0];
    for i = 1:3
        js = jParams.getClientProperty(sprintf('atomslot_%d', i));
        if isempty(js), continue; end
        lo = double(jParams.getClientProperty(sprintf('atomlo_%d', i)));
        hi = double(jParams.getClientProperty(sprintf('atomhi_%d', i)));
        vals(i) = lo + (hi - lo) * double(js.getValue()) / 1000;
    end
end

function SetAtomVals(jParams, vals) %#ok<DEFNU>
% Set the active slots to physical values (clamped); used when loading a stored atom's generator.
    for i = 1:min(3, numel(vals))
        js = jParams.getClientProperty(sprintf('atomslot_%d', i));
        if isempty(js), continue; end
        lo = double(jParams.getClientProperty(sprintf('atomlo_%d', i)));
        hi = double(jParams.getClientProperty(sprintf('atomhi_%d', i)));
        js.setValue(max(0, min(1000, round(1000 * (vals(i) - lo) / max(hi - lo, eps)))));
        i_atom_slider_changed(jParams, i, js, []);
    end
end

%% ===== internal =====
function i_atom_slider_changed(jParams, slot, js, onSettle)
    lo  = double(jParams.getClientProperty(sprintf('atomlo_%d', slot)));
    hi  = double(jParams.getClientProperty(sprintf('atomhi_%d', slot)));
    val = lo + (hi - lo) * double(js.getValue()) / 1000;
    jt  = jParams.getClientProperty(sprintf('atomttl_%d', slot));
    if ~isempty(jt)
        lbl = char(jParams.getClientProperty(sprintf('atomlbl_%d', slot)));
        fmt = char(jParams.getClientProperty(sprintf('atomfmt_%d', slot)));
        jt.setText(sprintf(['%s: ' fmt], lbl, val));
    end
    if (nargin >= 4) && ~isempty(onSettle) && ~js.getValueIsAdjusting(), onSettle(); end
end

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
