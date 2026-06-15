function varargout = panel_spatial_filter(varargin)
% PANEL_SPATIAL_FILTER: Live in-place spatial filter of a Dirac source map shown on a
% cortex figure. Launched from the figure popup; filters every time step's spatial
% field with bst_dirac_eigenmodes_filter and swaps the in-memory ImageGridAmp
% (non-destructive: the original is restored on toggle-off / close). Controls are the
% shared Filter-kernel section (bst_eigfilter_panel) + a Filter on/off toggle +
% Save filtered file + Close.
%
% Dispatched: Start(hFig), GetState(panelName), Apply/Restore/OnToggle/OnKernelChanged/
%             SaveFiltered/Close.
% Authors: Diellor Basha, 2026
    eval(macro_method);
end

%% ===== LAUNCH (from the source figure popup) =====
function Start(hFig) %#ok<DEFNU>
    global GlobalData;
    panelName = 'SpatialFilter';
    % resolve the displayed source results
    TessInfo = getappdata(hFig, 'Surface');
    iTess = find(arrayfun(@(t) ~isempty(t.DataSource) && strcmpi(t.DataSource.Type,'Source'), TessInfo), 1);
    if isempty(iTess)
        bst_error('No source map on this figure.', 'Spatial filter', 0); return;
    end
    [iDS, iResult] = bst_memory('GetDataSetResult', TessInfo(iTess).DataSource.FileName);
    if isempty(iResult)
        bst_error('Could not resolve the source results.', 'Spatial filter', 0); return;
    end
    R = GlobalData.DataSet(iDS).Results(iResult);
    if isempty(R.nComponents) || (R.nComponents ~= 3)
        bst_error('Spatial filter requires an unconstrained (3-component) source.', 'Spatial filter', 0); return;
    end
    SurfaceFile = R.SurfaceFile;

    % Dirac eigenbasis + operator mass (find-or-create)
    bst_progress('start', 'Spatial filter', 'Loading Dirac eigenbasis...');
    EigenMat = tess_eigen(SurfaceFile, 'Dirac');
    OpMat    = load(file_fullpath(EigenMat.OperatorFile));
    bst_progress('stop');
    nVert = double(max(cellfun(@(x) max(x(:)), EigenMat.GlobalVertices)));

    % materialize the full displayed field [3nVert x nT] and back it up
    J0 = i_full_field(iDS, iResult, nVert);
    if isempty(J0)
        bst_error(sprintf('Field size mismatch (expected 3*%d rows).', nVert), 'Spatial filter', 0); return;
    end

    % session state on the figure
    St = struct('iDS',iDS, 'iResult',iResult, 'hFig',hFig, 'iTess',iTess, ...
                'EigenMat',EigenMat, 'Mass',{OpMat.Mass}, 'Lambda',double(EigenMat.Lambda{1}(:)), ...
                'Orig',J0, 'OrigKernel',R.ImagingKernel, 'isOn',false);
    setappdata(hFig, 'SpatialFilterState', St);

    % build + dock the panel
    gui_hide(panelName);
    bstPanel = CreatePanel(St);
    gui_show(bstPanel, 'BrainstormTab', 'tools');
    try, gui_brainstorm('SetSelectedTab', panelName, 0); catch, end %#ok<CTCH>
    % restore the original if the figure closes
    set(hFig, 'DeleteFcn', @(h,e) Close(panelName));
end

%% ===== PANEL =====
function bstPanelNew = CreatePanel(St)
    import javax.swing.*;
    panelName = 'SpatialFilter';
    jPanelNew = gui_component('Panel');
    jOpt = JPanel(); jOpt.setLayout(BoxLayout(jOpt, BoxLayout.Y_AXIS));

    jSec = gui_river([2 2], [2 8 3 6], 'Spatial filter (Dirac)');
    [keys, displays] = bst_eigfilter_panel('Kernels');
    gui_component('label', jSec, 'br', 'Kernel:');
    jKernel = gui_component('combobox', jSec, 'br hfill', [], {displays}, [], [], []);
    % default to heat (low-pass)
    iHeat = find(strcmp(keys,'heat'),1); if ~isempty(iHeat); jKernel.setSelectedIndex(iHeat-1); end
    jParams = gui_river([2 2], [0 2 0 2]);
    jSec.add('br hfill', jParams);
    jFilterOn = gui_component('checkbox', jSec, 'br', 'Filter on');
    jSave   = gui_component('button', jSec, 'br', 'Save filtered file');
    jClose  = gui_component('button', jSec, '', 'Close');
    jOpt.add(jSec);
    jPanelNew.add(jOpt, java.awt.BorderLayout.NORTH);

    ctrl = struct('jKernel',jKernel, 'KernelKeys',{keys}, 'jParams',jParams, ...
                  'jFilterOn',jFilterOn, 'jSave',jSave, 'jClose',jClose, 'hFig',St.hFig);

    bst_eigfilter_panel('BuildSliders', jParams, bst_eigfilter_panel('CurrentKernel', jKernel, keys), St.Lambda, @() OnKernelOrScale(panelName));
    java_setcb(jKernel,   'ActionPerformedCallback', @(h,e) OnKernelChanged(panelName));
    java_setcb(jFilterOn, 'ActionPerformedCallback', @(h,e) OnToggle(panelName));
    java_setcb(jSave,     'ActionPerformedCallback', @(h,e) SaveFiltered(panelName));
    java_setcb(jClose,    'ActionPerformedCallback', @(h,e) Close(panelName));

    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);
end

function [St, ctrl] = GetState(panelName) %#ok<DEFNU>
    if (nargin < 1) || isempty(panelName); panelName = 'SpatialFilter'; end
    ctrl = bst_get('PanelControls', panelName);
    St = [];
    if isempty(ctrl) || ~isfield(ctrl,'hFig') || ~ishandle(ctrl.hFig); return; end
    St = getappdata(ctrl.hFig, 'SpatialFilterState');
end

%% ===== PURE: filter the whole series with a kernel =====
function Jf = ComputeFiltered(St, kernelName, params) %#ok<DEFNU>
    g  = bst_eigfilter_kernel(kernelName, params);
    Jf = real(bst_dirac_eigenmodes_filter(St.EigenMat, St.Mass, St.Orig, 'custom', 'TransferFn', g));
end

%% ===== APPLY (swap in the filtered series) =====
function Apply(panelName)
    global GlobalData;
    [St, ctrl] = GetState(panelName);
    if isempty(St); return; end
    name   = bst_eigfilter_panel('CurrentKernel', ctrl.jKernel, ctrl.KernelKeys);
    params = bst_eigfilter_panel('ReadParams', ctrl.jParams, St.Lambda);
    Jf = ComputeFiltered(St, name, params);
    GlobalData.DataSet(St.iDS).Results(St.iResult).ImageGridAmp  = Jf;
    GlobalData.DataSet(St.iDS).Results(St.iResult).ImagingKernel = [];
    St.isOn = true; setappdata(ctrl.hFig, 'SpatialFilterState', St);
    i_refresh(ctrl.hFig, St.iTess);
end

%% ===== RESTORE (put the original back) =====
function Restore(panelName)
    global GlobalData;
    [St, ctrl] = GetState(panelName);
    if isempty(St) || ~ishandle(ctrl.hFig); return; end
    GlobalData.DataSet(St.iDS).Results(St.iResult).ImageGridAmp  = St.Orig;
    GlobalData.DataSet(St.iDS).Results(St.iResult).ImagingKernel = St.OrigKernel;
    St.isOn = false; setappdata(ctrl.hFig, 'SpatialFilterState', St);
    i_refresh(ctrl.hFig, St.iTess);
end

function OnToggle(panelName) %#ok<DEFNU>
    [St, ctrl] = GetState(panelName);
    if isempty(St); return; end
    if ctrl.jFilterOn.isSelected(); Apply(panelName); else; Restore(panelName); end
end

function OnKernelChanged(panelName) %#ok<DEFNU>
    [St, ctrl] = GetState(panelName); %#ok<ASGLU>
    if isempty(ctrl); return; end
    key = bst_eigfilter_panel('CurrentKernel', ctrl.jKernel, ctrl.KernelKeys);
    bst_eigfilter_panel('BuildSliders', ctrl.jParams, key, St.Lambda, @() OnKernelOrScale(panelName));
    OnKernelOrScale(panelName);
end

% kernel/scale changed: re-apply only if the filter is currently on
function OnKernelOrScale(panelName)
    [St, ctrl] = GetState(panelName); %#ok<ASGLU>
    if ~isempty(ctrl) && ctrl.jFilterOn.isSelected(); Apply(panelName); end
end

%% ===== refresh the figure display =====
function i_refresh(hFig, iTess)
    TessInfo = getappdata(hFig, 'Surface');
    for k = 1:numel(TessInfo); TessInfo(k).DataMinMax = []; end
    setappdata(hFig, 'Surface', TessInfo);
    panel_surface('UpdateSurfaceData', hFig);
    panel_surface('UpdateSurfaceColormap', hFig);
    try, figure_3d('SetShowSourceVectors', hFig, iTess, 1); catch, end %#ok<CTCH>
end

%% ===== CLOSE (restore original, undock) =====
function Close(panelName) %#ok<DEFNU>
    [St, ctrl] = GetState(panelName);
    if ~isempty(St) && ~isempty(ctrl) && ishandle(ctrl.hFig)
        if St.isOn; Restore(panelName); end
        try, set(ctrl.hFig, 'DeleteFcn', ''); catch, end %#ok<CTCH>
        try, rmappdata(ctrl.hFig, 'SpatialFilterState'); catch, end %#ok<CTCH>
    end
    gui_hide(panelName);
end

%% ===== materialize the full displayed field as [3nVert x nT] =====
function J = i_full_field(iDS, iResult, nVert)
    global GlobalData;
    R = GlobalData.DataSet(iDS).Results(iResult);
    if ~isempty(R.ImageGridAmp)
        J = double(R.ImageGridAmp);
    else
        J = double(bst_memory('GetResultsValues', iDS, iResult, [], [], 0));   % [3nVert x nT], no orient
    end
    if isempty(J) || (size(J,1) ~= 3*nVert); J = []; end
end
