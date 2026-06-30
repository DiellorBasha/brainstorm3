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
% SEE ALSO: bst_eigen('Axes'), bst_eigenfilter('Atom'), bst_dynamics, bst_eigfilter_kernel, view_surface_data

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

    % --- eigenbasis + canonical axes (bst_eigen('Axes')) + physical-scale calibration; (re)built by i_build_basis ---
    mm = @(l) 2*pi./sqrt(max(l,eps))*1000;                     % eigenvalue -> physical scale (mm)
    nFrames = 100;                                             % standard 1 s @ 100 Hz time axis (ax.nT/tlag from bst_eigen)
    [lmax, lminPos, scaleMinMM, scaleMaxMM, rateMinMM2, rateMaxMM2] = deal(0);   % shared with i_build_basis
    ax = struct();  i_build_basis(variant);                    % ax = bst_eigen('Axes', ...) + scale/rate bounds

    % --- state (closure-shared) --- (the kernel list lives in the docked panel: panel_eigenfilter_design('AtomKernels'))
    state    = 'design';
    kernel   = 'diffusion';                          % default kernel (the panel's Filter combobox selects it)
    normMode = '';                                  % active colormap normalization (set by i_normalize)
    showFib  = false;  hFibers = [];  fibPts = [];  fibEndV = [];   % connectome fiber overlay (lazy-loaded)
    % default seed: vertex nearest the centroid of the first eigenbasis support
    if (nargin >= 4) && ~isempty(seed0)
        seedVtx = seed0;
    else
        gv1 = ax.GlobalVertices{1};  [~,j] = min(sum((V(gv1,:) - mean(V(gv1,:),1)).^2, 2));  seedVtx = gv1(j);
    end
    W = i_eval_atom(seedVtx, ax, kernel, struct('lmax',lmax), V, nV);   % init frame from kernel defaults (sliders not built yet); Generate() recomputes with the slider params at the end

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
    iRes   = bst_memory('GetResultInDataSet', iDS, resFile);
    TI     = getappdata(hFig, 'Surface');  iTess = find(~cellfun('isempty', {TI.SurfaceFile}), 1);

    % --- design controls live in a DOCKED panel (panel_atom_designer), NOT on the figure ---
    cb = struct('Kernel',   @()bst_call(@KernelChanged), ...
                'Operator', @()bst_call(@OperatorChanged), ...
                'Param',    @()bst_call(@ParamChanged), ...
                'Fibers',   @(s)bst_call(@()OnToggleConnectome(s)), ...
                'Save',     @()bst_call(@SaveAtom));
    try, gui_hide('AtomDesigner'); catch, end %#ok<CTCH>
    bstPanelDes = panel_atom_designer('CreatePanel');
    gui_show(bstPanelDes, 'BrainstormTab', 'tools');
    panel_atom_designer('Configure', cb, i_bounds(), kernel, variant);
    Status();

    % --- vertex picking: use figure_3d's native WaveletDesignerPick hook (select3d-based via
    %     panel_coordinates('SelectPoint'): occlusion-correct, with native click/drag separation) rather
    %     than overriding the mouse callbacks. Camera, colormap context menu and the native time stepper
    %     are all preserved untouched. Only the close handler is chained (to remove the working file). ---
    setappdata(hFig, 'WaveletDesignerPick', @i_seed);
    origClose = get(hFig,'CloseRequestFcn');
    set(hFig, 'CloseRequestFcn', @OnClose);
    Generate();                                  % normalize the initial frame + set the kernel-matched colormap

    % ===== nested callbacks =====
    function i_build_basis(newVar)
        % (Re)build the canonical axes (eigenbasis x time x temporal-frequency) for the chosen operator via
        % bst_eigen('Axes'), then refresh the physical-scale bounds. Works at init and on live operator switches.
        variant = newVar;
        bst_progress('start', 'Atom designer', sprintf('Building/loading %s eigenbasis...', variant));
        ax = bst_eigen('Axes', struct('SurfaceFile',SurfaceFile, 'Variant',variant, ...
                       'nModes',nModes, 'TimeWindow',[0 (nFrames-1)/100], 'SampleRate',100));
        bst_progress('stop');
        lamAll = ax.Lambda{1}(:); if numel(ax.Lambda) > 1 && ~isempty(ax.Lambda{2}), lamAll = [lamAll; ax.Lambda{2}(:)]; end
        lmax = max(lamAll);  lminPos = min(lamAll(lamAll > 1e-9));
        scaleMinMM = mm(lmax);  scaleMaxMM = mm(lminPos);          % finest .. coarsest cortical scale
        rateMinMM2 = scaleMinMM^2;  rateMaxMM2 = scaleMaxMM^2;     % diffusion rate bounds (mm^2/s) = scale^2
    end
    function b = i_bounds()
        b = struct('scaleMinMM',scaleMinMM, 'scaleMaxMM',scaleMaxMM, 'rateMinMM2',rateMinMM2, 'rateMaxMM2',rateMaxMM2);
    end
    function OperatorChanged()
        newVar = panel_atom_designer('CurrentOperator');
        if strcmpi(newVar, variant), return; end                   % unchanged (e.g. the configure-time callback)
        variant = newVar;
        i_build_basis(variant);
        gv = ax.GlobalVertices{1};                                 % keep the seed if still supported, else recenter
        if ~ismember(seedVtx, gv), [~,j] = min(sum((V(gv,:)-mean(V(gv,:),1)).^2,2)); seedVtx = gv(j); end
        panel_atom_designer('RebuildSliders', kernel, i_bounds());  % new spectrum -> refresh slider ranges
        Regen();  i_reset_time();                                  % new basis -> diffuse anew from t=0
    end
    function KernelChanged()
        newK = panel_atom_designer('CurrentKernel');
        if strcmp(newK, kernel), return; end                       % unchanged (e.g. the configure-time callback)
        kernel = newK;
        panel_atom_designer('RebuildSliders', kernel, i_bounds());  % contextual params for the new kernel
        Regen();  i_reset_time();                                  % new kernel -> diffuse anew from t=0
    end
    function ParamChanged()
        Regen();                                                   % panel sliders update their own readouts
    end
    function Status()
        panel_atom_designer('SetStatus', sprintf('[%s] %s @ vtx %d | norm: %s | arrows=time', ...
            state, kernel, seedVtx, normMode));
    end
    function [W, isSigned] = i_normalize(W)
        % Pick the sign-class -> normalization + colormap. Honour the kernel's DECLARED class first so
        % the colormap is STABLE per kernel (a band-pass WAVELET like gabor/mexhat is always signed ->
        % diverging; a low-pass / prior like heat/matern is always one-signed -> density). Only fall back
        % to the empirical mass test when the kernel declares neither (the wave family). This prevents
        % e.g. gabor flipping to the sequential 'hot' map when its band slides onto near-DC modes at
        % large scale (where the field happens to be mass-positive but it is still a wavelet).
        cls  = i_kernel_signclass(kernel);               % 'signed' | 'onesigned' | '' (indeterminate)
        gv   = ax.GlobalVertices{1};
        mvec = full(sum(ax.Mass{1}, 2));                 % lumped vertex areas on the eigenbasis support
        Wg   = W(gv, :);
        s    = mvec.' * Wg;                              % [1 x nT] signed mass integral per frame
        l1   = mvec.' * abs(Wg);                         % [1 x nT] total absolute mass per frame
        empOneSigned = isempty(cls) && (min(abs(s)./max(l1,eps)) > 0.1);
        if strcmp(cls,'onesigned') || empOneSigned       % mass-conserving -> probability density
            sden = s;  sden(abs(sden) < eps) = 1;        % guard frames with ~zero mass
            W = W ./ sden;  normMode = 'density (unit mass)';  isSigned = false;
        else                                             % wavelet / signed -> relative amplitude
            pk = max(abs(W(:)));  if pk > 0, W = W / pk; end
            normMode = 'peak';  isSigned = true;
        end
    end
    function c = i_kernel_signclass(k)
        % A band-pass kernel is a wavelet (zero-mean -> signed); a priorAdmissible low-pass is one-signed.
        c = '';
        try
            m = bst_eigfilter_kernel('info', k);
            if isfield(m,'bandpass') && m.bandpass,                c = 'signed';
            elseif isfield(m,'priorAdmissible') && m.priorAdmissible, c = 'onesigned';
            end
        catch %#ok<CTCH>
        end
    end
    function Regen(), if ~isempty(seedVtx), Generate(); end, end
    function Generate()
        try
            W = i_eval_atom(seedVtx, ax, kernel, i_phys2kernel(), V, nV);
            [W, isSigned] = i_normalize(W);                              % one-signed (density) vs signed (peak)
            if isSigned, i_set_cmap('stat2'); else, i_set_cmap('source'); end   % diverging-symmetric vs sequential
            GlobalData.DataSet(iDS).Results(iRes).ImageGridAmp = W;       % in-place update, no file I/O
            T2 = getappdata(hFig,'Surface');  T2(iTess).DataMinMax = [min(W(:)) max(W(:))];  setappdata(hFig,'Surface',T2);
            panel_surface('UpdateSurfaceData', hFig, iTess);  panel_surface('UpdateSurfaceColormap', hFig);
            if showFib, i_recolor_fibers(); end                          % keep the fiber overlay in sync with the atom
            Status();
        catch ME
            panel_atom_designer('SetStatus',['Could not propagate: ' regexprep(ME.message,'\s+',' ')]);
        end
    end
    function OnToggleConnectome(stateOn)
        showFib = logical(stateOn);
        if showFib
            if isempty(fibPts), i_load_fibers(); end
            if isempty(fibPts), showFib = false; return; end   % nothing to show
            [~, hFibers] = figure_3d('PlotFibers', hFig, fibPts, repmat([.5 .5 .5], size(fibPts,1), 1));
            setappdata(hFig, 'AtomFiberRecolor', @i_recolor_fibers);    % recolor on every time frame (bst_figures hook)
            i_recolor_fibers();
        else
            if ~isempty(hFibers), try, delete(hFibers); catch, end, end %#ok<CTCH>
            hFibers = [];
            if isappdata(hFig,'AtomFiberRecolor'), rmappdata(hFig,'AtomFiberRecolor'); end
        end
    end
    function i_load_fibers()
        % Locate a Fibers surface (subject's own, else the default-anatomy HCP-1065) and load the full
        % streamlines. The subject cortex and the template fibers live in DIFFERENT SCS (different brain),
        % so the raw template cloud is mis-scaled/-placed; we similarity-align it to this cortex (below).
        ff = '';  isDef = false;  k = find(strcmpi({sSubj.Surface.SurfaceType},'Fibers'), 1);
        if ~isempty(k), ff = sSubj.Surface(k).FileName;
        else
            sDef = bst_get('Subject', 0);  kd = [];
            if ~isempty(sDef) && isfield(sDef,'Surface') && ~isempty(sDef.Surface)
                kd = find(strcmpi({sDef.Surface.SurfaceType},'Fibers'), 1);
            end
            if ~isempty(kd), ff = sDef.Surface(kd).FileName;  isDef = true; end
        end
        if isempty(ff), panel_atom_designer('SetStatus','[connectome] no fibers on this subject or the default anatomy.'); return; end
        FibMat = load(file_fullpath(ff), 'Points');  P = double(FibMat.Points);   % [nF x nP x 3]
        cap = 5000;  nF = size(P,1);
        if nF > cap, P = P(sort(randperm(nF, cap)),:,:); end                       % stay under PlotFibers' cap
        if isDef
            % The HCP-1065 fibers live in the default-anatomy (ICBM152) SCS. Bring them to THIS subject's SCS
            % THROUGH MNI: default-SCS -> MNI -> subject-SCS via cs_convert (a true 3D transform using each
            % MRI's MNI normalisation, so white-matter interiors map correctly, not just the endpoints).
            sz = size(P);  Pf = reshape(P, [], 3);  Psub = [];
            try
                sMriDef  = in_mri_bst(sDef.Anatomy(sDef.iAnatomy).FileName);
                sMriSubj = in_mri_bst(sSubj.Anatomy(sSubj.iAnatomy).FileName);
                Psub = cs_convert(sMriSubj, 'mni', 'scs', cs_convert(sMriDef, 'scs', 'mni', Pf, 1), 1);
            catch, end %#ok<CTCH>
            if ~isempty(Psub) && (mean(isnan(Psub(:))) < 0.5)
                P = reshape(Psub, sz);  P = P(~any(any(isnan(P),3),2), :, :);      % drop any out-of-FOV fibers
            else
                % no MNI normalisation available -> coarse centroid+scale alignment of the endpoint cloud
                e  = [squeeze(P(:,1,:)); squeeze(P(:,end,:))];
                cT = mean(e,1);  sT = sqrt(mean(sum((e-cT).^2,2)));
                cS = mean(V,1);  sS = sqrt(mean(sum((V-cS).^2,2)));
                P  = reshape((Pf - cT) * (sS/max(sT,eps)) + cS, sz);
                panel_atom_designer('SetStatus','[connectome] no MNI transform; coarse similarity alignment used.');
            end
        end
        fibPts  = P;
        fibEndV = [dsearchn(V, squeeze(P(:,1,:))), dsearchn(V, squeeze(P(:,end,:)))];   % endpoint -> nearest vertex
    end
    function i_recolor_fibers()
        % Colour each fiber by the atom: kernel value at its two endpoint vertices, interpolated along the
        % streamline, mapped through the active overlay colormap + CLim. Recolours in place (no re-plot).
        if ~showFib || isempty(hFibers) || isempty(fibPts) || ~ishandle(hFig), return; end
        try
            iT   = i_current_frame();
            Wcur = GlobalData.DataSet(iDS).Results(iRes).ImageGridAmp(:, iT);   % [nV x 1] atom at this frame
            nP   = size(fibPts, 2);  tt = linspace(0, 1, nP);
            S    = Wcur(fibEndV(:,1)) * (1-tt) + Wcur(fibEndV(:,2)) * tt;       % [nF x nP] endpoint interpolation
            figure_3d('ColorFibers', hFibers, i_scalar2rgb(S));
        catch, end %#ok<CTCH>
    end
    function iT = i_current_frame()
        ct = GlobalData.UserTimeWindow.CurrentTime;  if isempty(ct), ct = ax.tlag(1); end
        [~, iT] = min(abs(ax.tlag - ct));  iT = max(1, min(iT, ax.nT));
    end
    function RGB = i_scalar2rgb(S)
        TIc   = getappdata(hFig, 'Surface');  sCmap = bst_colormaps('GetColormap', TIc(iTess).ColormapType);
        CMap  = sCmap.CMap;  N = size(CMap, 1);
        hAx   = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');  cl = get(hAx, 'CLim');
        if sCmap.isAbsoluteValues, S = abs(S); end                          % source: |value|; stat2: signed
        idx   = round((S - cl(1)) / max(cl(2)-cl(1), eps) * (N-1)) + 1;  idx = max(1, min(idx, N));
        RGB   = reshape(CMap(idx(:), :), size(S,1), size(S,2), 3);
    end
    function i_set_cmap(cmapType)
        % Match the colorbar to the kernel's sign-class: 'source' (sequential, absolute) for one-signed
        % fields (diffusion/heat), 'stat2' (diverging cmap_mandrill, signed -> auto 0-centred symmetric range)
        % for sign-changing fields (waves/mexhat). Switches the overlay's type only; the user's global
        % 'source' colormap is never modified. No-op when already on the right type (so it never flickers).
        TIc = getappdata(hFig, 'Surface');
        if isempty(TIc) || (iTess > numel(TIc)) || strcmpi(TIc(iTess).ColormapType, cmapType), return; end
        TIc(iTess).ColormapType = cmapType;  setappdata(hFig, 'Surface', TIc);
        other = 'source'; if strcmpi(cmapType,'source'), other = 'stat2'; end
        try, bst_colormaps('AddColormapToFigure', hFig, cmapType);  bst_colormaps('RemoveColormapFromFigure', hFig, other); catch, end
    end
    function i_seed(vi)                          % figure_3d native pick callback: clicked cortical vertex -> seed
        if ~strcmpi(state,'design') || isempty(vi), return; end
        seedVtx = vi;  Generate();  i_reset_time();       % new seed -> diffuse anew from t=0
    end
    function i_reset_time()                       % rewind the global time cursor to the start of the atom (t=0)
        try, panel_time('SetCurrentTime', ax.tlag(1)); catch, end
    end
    function kp = i_phys2kernel()
        kp = bst_eigfilter_controls('ToKernel', kernel, panel_atom_designer('ReadVals'), lmax);
    end
    function SaveAtom(~,~)
        if isempty(W), panel_atom_designer('SetStatus','[design] nothing to save - drop an atom first.'); return; end
        state='save';  Status();
        T=bst_dynamics('New', sprintf('atom %s @vtx%d', kernel, seedVtx)); T.SurfaceFile=SurfaceFile;
        G=bst_dynamics('NewGroup','atom'); G.vertices=seedVtx; G.pos=V(seedVtx,:);
        LS=bst_dynamics('Levelset', W(ax.GlobalVertices{1},:), ax.GlobalVertices{1}, 0.5);
        G.region={LS.scoutVertices(:)'}; G.times=[ax.tlag(LS.eventSamples(1)); ax.tlag(LS.eventSamples(end))];
        T=bst_dynamics('AddGroup',T,G);
        outDir=bst_fileparts(file_fullpath(sSubj.FileName));
        of=file_unique(bst_fullfile(outDir,'dynamics_atom.mat')); bst_dynamics('Save', of, T);
        state='design';
        panel_atom_designer('SetStatus',sprintf('[saved] atom -> %s (Scout %d vtx, Event %.2f-%.2fs)', file_short(of), numel(LS.scoutVertices), G.times(1), G.times(2)));
    end
    function OnClose(h,ev)
        try, gui_hide('AtomDesigner'); catch, end %#ok<CTCH>          % undock the design panel
        try, i_call(origClose,h,ev); catch, if ishandle(hFig), delete(hFig); end, end
        try, file_delete(file_fullpath(resFile), 1); db_reload_studies(iStudyG); catch, end   % remove the working file
    end
end

% ===== helpers (local) =====
function W = i_eval_atom(s, ax, kernel, kp, V, nV) %#ok<INUSL>
    [Wloc, gv] = bst_eigenfilter('Atom', ax, kernel, kp, s);   % V unused now; signature kept for callers
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
