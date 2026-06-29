function hFig = view_atom_designer(SurfaceFile, variant, nModes, seed0)
% VIEW_ATOM_DESIGNER: Interactive joint time-cortex WAVELET ATOM designer. Opens a figure_3d cortex view
% with a top-right filter-design panel. LEFT-CLICK a vertex (on the cortex) to drop a delta and propagate
% it as a dynamic eigenwavelet (atom) through the cortex eigenbasis (geometry with Laplace-Beltrami, or
% geometry+connectome with LB-Connectome). DRAG / click off the cortex rotates the camera.
%
% The atom is displayed through Brainstorm's MANAGED SOURCE OVERLAY (a working results file + in-place
% updates), so the native colormap context menu (right-click), the colorbar, and the native time stepper
% (LEFT/RIGHT arrows over the 1 s @ 100 Hz atom time axis) all work as for any source map. The working
% results file is removed when the designer window is closed.
%
% Filter-design panel: a flat kernel dropdown grouped into 'dynamic' (atoms: dampedwave/wave/kleingordon/
% diffusion) and 'static' (eigenwavelets: heat/mexhat/...); sliders for the kernel's PHYSICAL parameters
% (spatial Scale in mm calibrated to the cortex via 2*pi/sqrt(lambda); wave Speed c in m/s, non-separable
% kernels only; Decay in s, damped only); and Save (level-sets the atom into a Scout + Event).
%
% USAGE:  hFig = view_atom_designer(SurfaceFile)                       % default LB-Connectome
%         hFig = view_atom_designer(SurfaceFile, 'Laplace-Beltrami')  % geometry-only
%         hFig = view_atom_designer(SurfaceFile, variant, nModes, seed0)
%
% SEE ALSO: bst_atom, bst_eigfilter_kernel, tess_eigen, view_surface_data, view_manifold

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% Copyright (c) University of Southern California & McGill University
% Distributed under the terms of the GNU General Public License (GPLv3).
% =============================================================================@
%
% Authors: Diellor Basha, 2026

    global GlobalData;
    if (nargin < 2) || isempty(variant), variant = 'LB-Connectome'; end
    if (nargin < 3) || isempty(nModes),  nModes  = 200; end

    % --- surface geometry ---
    sCx = in_tess_bst(SurfaceFile);  V = sCx.Vertices;  nV = size(V,1);
    meanEdge = i_mean_edge(V, sCx.Faces);

    % --- eigenbasis (cached) + physical-scale calibration ---
    bst_progress('start', 'Atom designer', sprintf('Building/loading %s eigenbasis...', variant));
    E  = tess_eigen(SurfaceFile, variant, 'nModes', nModes);   % cache + reuse the eigen_ file
    Op = in_bst_operator(E.OperatorFile);
    bst_progress('stop');
    ax = struct('Phi',{E.Phi}, 'Lambda',{E.Lambda}, 'Mass',{Op.Mass}, 'GlobalVertices',{E.GlobalVertices});
    lamAll = E.Lambda{1}(:); if numel(E.Lambda) > 1 && ~isempty(E.Lambda{2}), lamAll = [lamAll; E.Lambda{2}(:)]; end
    lmax = max(lamAll);  lminPos = min(lamAll(lamAll > 1e-9));
    mm   = @(l) 2*pi./sqrt(max(l,eps))*1000;                   % eigenvalue -> physical scale (mm)
    scaleMinMM = mm(lmax);  scaleMaxMM = mm(lminPos);          % finest .. coarsest cortical scale

    % --- registry: one FLAT kernel list, dynamic + static group headers ---
    allK = bst_eigfilter_kernel('list');  spatialK = {}; dynK = {};
    for ii = 1:numel(allK)
        try m = bst_eigfilter_kernel('info', allK{ii}); catch, continue; end
        if isfield(m,'domain') && ~isempty(m.domain), dynK{end+1}=allK{ii}; else, spatialK{end+1}=allK{ii}; end %#ok<AGROW>
    end
    pref = {'dampedwave','wave','kleingordon','diffusion'};
    dynK = [pref(ismember(pref,dynK)), setdiff(dynK, pref, 'stable')];
    DH = '──── dynamic ────';  SH = '──── static ────';
    kList = [{DH}, dynK, {SH}, spatialK];  iDH = 1;  iSH = numel(dynK) + 2;

    % --- state (closure-shared) ---
    state    = 'design';
    kernel   = 'dampedwave';  lastKIdx = 2;
    pScaleMM = round((scaleMinMM+scaleMaxMM)/2);   % spatial scale (mm)
    pSpeed   = 1.0;                                 % wave speed c (m/s), non-separable only
    pDecay   = 0.5;                                 % decay time (s), damped only
    nFrames  = 100;  ax.nT = nFrames;  ax.tlag = (0:nFrames-1)/100;   % standard 1 s @ 100 Hz time vector
    % default seed: vertex nearest the centroid of the first eigenbasis support
    if (nargin >= 4) && ~isempty(seed0)
        seedVtx = seed0;
    else
        gv1 = ax.GlobalVertices{1};  [~,j] = min(sum((V(gv1,:) - mean(V(gv1,:),1)).^2, 2));  seedVtx = gv1(j);
    end
    W = i_eval_atom(seedVtx, ax, kernel, i_phys2kernel(), V, nV);

    % --- working results file -> display through the managed SOURCE overlay (native colormap/bar/stepper) ---
    [sSubj] = bst_get('SurfaceFile', SurfaceFile);
    [ss, iSs] = bst_get('StudyWithSubject', sSubj.FileName);
    kk = find(arrayfun(@(s)~isempty(s.Result), ss), 1);  if isempty(kk), iStudyG = iSs(1); else, iStudyG = iSs(kk); end
    R = db_template('resultsmat');  R.ImageGridAmp = W;  R.Time = ax.tlag;  R.nComponents = 1;
    R.SurfaceFile = file_short(SurfaceFile);  R.HeadModelType = 'surface';  R.Comment = 'Atom designer (working)';
    resFile = db_add(iStudyG, R);
    [hFig, iDS] = view_surface_data(SurfaceFile, resFile);
    if isempty(hFig), bst_error('Could not open the source display.', 'Atom designer', 0); return; end
    set(hFig, 'Name', ['Atom designer: ' SurfaceFile]);
    hAxes  = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    iRes   = bst_memory('GetResultInDataSet', iDS, resFile);
    TI     = getappdata(hFig, 'Surface');  iTess = find(~cellfun('isempty', {TI.SurfaceFile}), 1);

    % --- top-right filter-design panel: kernel + sliders, compactly spaced ---
    hP = uipanel('Parent',hFig, 'Title','Filter design', 'Units','normalized', 'Position',[0.70 0.68 0.295 0.30], ...
                 'FontUnits','points','FontSize',bst_get('FigFont'));
    yk = 0.80; ys = [0.60 0.43 0.26];  lw=0.30; sw=0.42; vw=0.16;
    uicontrol(hP,'Style','text','String','Kernel','Units','normalized','Position',[0.04 yk lw 0.13],'HorizontalAlignment','left');
    hKern = uicontrol(hP,'Style','popupmenu','String',kList,'Value',2,'Units','normalized','Position',[0.35 yk 0.61 0.14],'Callback',@KernelChanged);
    [hScale,hScaleV] = i_slider(hP,'Scale', ys(1), lw,sw,vw, scaleMinMM, scaleMaxMM, pScaleMM, @ParamChanged);
    [hSpeed,hSpeedV] = i_slider(hP,'Speed', ys(2), lw,sw,vw, 0.1, 10,  pSpeed, @ParamChanged);
    [hDecay,hDecayV] = i_slider(hP,'Decay', ys(3), lw,sw,vw, 0.05, 2,  pDecay, @ParamChanged);
    hSave = uicontrol(hP,'Style','pushbutton','String','Save','Units','normalized','Position',[0.80 0.04 0.16 0.14], ...
                      'Callback',@SaveAtom, 'TooltipString','Save atom -> Scout + Event');
    hLabel = uicontrol(hFig,'Style','text','String','', 'Units','Pixels','Position',[8 6 820 18], ...
        'HorizontalAlignment','left','FontUnits','points','FontSize',bst_get('FigFont'));
    SyncControls();

    % --- chain figure_3d handlers: keep camera, context menu AND the native time stepper (no KeyPress override) ---
    origDown = get(hFig,'WindowButtonDownFcn');  origUp = get(hFig,'WindowButtonUpFcn');  origClose = get(hFig,'CloseRequestFcn');
    downXY = [];
    set(hFig, 'WindowButtonDownFcn', @OnDown, 'WindowButtonUpFcn', @OnUp, 'CloseRequestFcn', @OnClose);

    % ===== nested callbacks =====
    function KernelChanged(src,~)
        idx = get(src,'Value');
        if idx==iDH || idx==iSH, set(src,'Value',lastKIdx); return; end
        lastKIdx = idx;  opt = get(src,'String');  kernel = opt{idx};
        SyncControls();  Regen();
    end
    function ParamChanged(~,~)
        pScaleMM=get(hScale,'Value'); pSpeed=get(hSpeed,'Value'); pDecay=get(hDecay,'Value');
        set(hScaleV,'String',num2str(round(pScaleMM))); set(hSpeedV,'String',sprintf('%.2g',pSpeed)); set(hDecayV,'String',sprintf('%.2g',pDecay));
        Regen();
    end
    function SyncControls()
        m = bst_eigfilter_kernel('info', kernel);
        isDyn = isfield(m,'domain') && ~isempty(m.domain);
        isSep = ~isfield(m,'separable') || m.separable;
        set([hSpeed hSpeedV], 'Enable', i_en(isDyn && ~isSep));
        set([hDecay hDecayV], 'Enable', i_en(strcmpi(kernel,'dampedwave')));
        Status();
    end
    function Status()
        m = bst_eigfilter_kernel('info', kernel);  isSep = ~isfield(m,'separable') || m.separable;
        set(hLabel,'String',sprintf('[%s] %s @ vtx %d%s | scale %.0f mm (cortex %.0f-%.0f) | right-click=colormap, arrows=time', ...
            state, kernel, seedVtx, i_en2b(~isSep,pSpeed), pScaleMM, scaleMinMM, scaleMaxMM));
    end
    function Regen(), if ~isempty(seedVtx), Generate(); end, end
    function Generate()
        try
            W = i_eval_atom(seedVtx, ax, kernel, i_phys2kernel(), V, nV);
            GlobalData.DataSet(iDS).Results(iRes).ImageGridAmp = W;       % in-place update, no file I/O
            T2 = getappdata(hFig,'Surface');  T2(iTess).DataMinMax = [min(W(:)) max(W(:))];  setappdata(hFig,'Surface',T2);
            panel_surface('UpdateSurfaceData', hFig, iTess);  panel_surface('UpdateSurfaceColormap', hFig);
            Status();
        catch ME
            set(hLabel,'String',['Could not propagate: ' regexprep(ME.message,'\s+',' ')]);
        end
    end
    function OnDown(h,ev), downXY = get(hFig,'CurrentPoint'); i_call(origDown,h,ev); end
    function OnUp(h,ev)
        i_call(origUp,h,ev);
        if strcmpi(state,'design') && strcmpi(get(hFig,'SelectionType'),'normal') && ~isempty(downXY)
            if norm(get(hFig,'CurrentPoint')-downXY) < 4               % click, not a drag
                [v,hit] = PickVertex();  if hit, seedVtx=v; Generate(); end
            end
        end
    end
    function [v,hit] = PickVertex()
        cp=get(hAxes,'CurrentPoint'); o=cp(1,:); d=cp(2,:)-cp(1,:); d=d/max(norm(d),eps);
        w=V-o; t=w*d'; proj=o+t.*d; dist=sqrt(sum((V-proj).^2,2));
        near=find(dist < 2*meanEdge);
        if isempty(near), v=[]; hit=false; return; end
        [~,j]=min(t(near)); v=near(j); hit=true;                        % front-most grazed vertex
    end
    function kp = i_phys2kernel()
        kp = struct('lmax', lmax);
        switch lower(kernel)
            case {'wave','dampedwave','kleingordon'}
                kp.alpha = pSpeed * sqrt(lmax) / 2;                    % c (m/s) -> alpha
                if strcmpi(kernel,'dampedwave'), kp.beta = 1/max(pDecay,eps); end
                if strcmpi(kernel,'kleingordon'), kp.mu = 0.1*lmax; end
            case 'diffusion'
                kp.tau = max((pScaleMM/1000)^2 * lmax, eps);
            case 'heat'
                lamS = (2*pi/(pScaleMM/1000))^2; kp.t = log(2)/max(lamS,eps);
            case 'mexhat'
                lamS = (2*pi/(pScaleMM/1000))^2; kp.t = 1/max(lamS,eps);
        end
    end
    function SaveAtom(~,~)
        if isempty(W), set(hLabel,'String','[design] nothing to save - drop an atom first.'); return; end
        state='save';  Status();
        T=bst_dynamics('New', sprintf('atom %s @vtx%d', kernel, seedVtx)); T.SurfaceFile=SurfaceFile;
        G=bst_dynamics('NewGroup','atom'); G.vertices=seedVtx; G.pos=V(seedVtx,:);
        LS=bst_atom('Levelset', W(ax.GlobalVertices{1},:), ax.GlobalVertices{1}, 0.5);
        G.region={LS.scoutVertices(:)'}; G.times=[ax.tlag(LS.eventSamples(1)); ax.tlag(LS.eventSamples(end))];
        T=bst_dynamics('AddGroup',T,G);
        outDir=bst_fileparts(file_fullpath(sSubj.FileName));
        of=file_unique(bst_fullfile(outDir,'dynamics_atom.mat')); bst_dynamics('Save', of, T);
        state='design';
        set(hLabel,'String',sprintf('[saved] atom -> %s (Scout %d vtx, Event %.2f-%.2fs)', file_short(of), numel(LS.scoutVertices), G.times(1), G.times(2)));
    end
    function OnClose(h,ev)
        try, i_call(origClose,h,ev); catch, if ishandle(hFig), delete(hFig); end, end
        try, file_delete(file_fullpath(resFile), 1); db_reload_studies(iStudyG); catch, end   % remove the working file
    end
end

% ===== helpers (local) =====
function W = i_eval_atom(s, ax, kernel, kp, V, nV)
    G = bst_dynamics('NewGroup','atom');  G.vertices = s;  G.pos = V(s,:);
    [Wloc, gv] = bst_atom('Evaluate', G, 1, ax, kernel, kp);
    W = zeros(nV, ax.nT);  W(gv,:) = Wloc;
end
function i_call(fcn,h,ev)
    if isempty(fcn), return; end
    if isa(fcn,'function_handle'), fcn(h,ev);
    elseif iscell(fcn), feval(fcn{1},h,ev,fcn{2:end});
    elseif ischar(fcn), eval(fcn); end %#ok<EVLDR>
end
function s = i_en(tf),  if tf, s='on'; else, s='off'; end, end
function s = i_en2b(show,c), if show, s=sprintf(' | speed %.2g m/s', c); else, s=''; end, end
function [hS,hV] = i_slider(p, name, y, lw, sw, vw, mn, mx, v0, cb)
    uicontrol(p,'Style','text','String',name,'Units','normalized','Position',[0.04 y lw 0.12],'HorizontalAlignment','left');
    hS = uicontrol(p,'Style','slider','Min',mn,'Max',mx,'Value',max(min(v0,mx),mn),'Units','normalized','Position',[0.35 y+0.01 sw 0.10],'Callback',cb);
    hV = uicontrol(p,'Style','text','String',num2str(round(v0,2)),'Units','normalized','Position',[0.80 y vw 0.12],'HorizontalAlignment','left');
end
function L = i_mean_edge(V,F)
    e=[F(:,[1 2]);F(:,[2 3]);F(:,[3 1])]; L=mean(sqrt(sum((V(e(:,1),:)-V(e(:,2),:)).^2,2)));
end
