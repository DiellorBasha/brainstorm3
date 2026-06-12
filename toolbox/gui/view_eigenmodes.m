function varargout = view_eigenmodes(varargin)
% VIEW_EIGENMODES: Standalone vector-field viewer for Dirac eigenmodes.
%
% USAGE:  hFig = view_eigenmodes(EigenFile)
%         V3   = view_eigenmodes('ReconstructModeField', EigenMat, k, nVert)
%
% Loads an eigen_ DB node (db_template('eigenmat')) and renders each Dirac
% eigenvector as an ambient 3D quiver field on the parent cortex. Cycle modes
% with the keyboard (Left/Right = -/+1, PgUp/PgDn = +/-10), like
% view_leadfield_vectors cycles sensor channels.
%
% Only the Dirac variant is implemented; LBO / Connection Laplacian raise a
% bst_error (added later as variant branches). Dirac eigenvectors are real
% [4*nVh x K] quaternion fields; the per-vertex 3-vector is the quaternion
% vector part (rows 2:4 of each 4-block; the w slot 1:4:end is dropped) --
% exactly as bst_dirac/local_reconstruct extracts the ambient field.
%
% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Diellor Basha, 2026

if (nargin >= 1) && ischar(varargin{1}) && any(strcmp(varargin{1}, {'ReconstructModeField'}))
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
[varargout{1:nargout}] = ViewFigure(varargin{:});
end


%% ===== PURE: per-vertex ambient 3-vector for Dirac eigenmode k =====
function V3 = ReconstructModeField(EigenMat, k, nVert)
% V3 [nVert x 3]: zeros off-support; quaternion vector part (i,j,k)->(x,y,z),
% w slot dropped; scattered to global vertices via EigenMat.GlobalVertices.
    if ~isfield(EigenMat,'Phi') || isempty(EigenMat.Phi) || numel(EigenMat.Phi) ~= 2
        error('view_eigenmodes:badEigen', 'EigenMat.Phi must be a 1x2 per-hemisphere cell.');
    end
    V3 = zeros(nVert, 3);
    for hh = 1:2
        vH  = EigenMat.GlobalVertices{hh}(:);
        Phi = EigenMat.Phi{hh};
        if (k < 1) || (k > size(Phi,2))
            error('view_eigenmodes:badMode', 'Mode index %d out of range 1..%d.', k, size(Phi,2));
        end
        if size(Phi,1) ~= 4*numel(vH)
            error('view_eigenmodes:shapeMismatch', ...
                'Hemisphere %d: Phi has %d rows, expected 4*nV=%d.', hh, size(Phi,1), 4*numel(vH));
        end
        col = double(Phi(:, k));
        V3(vH,1) = col(2:4:end);   % quaternion i -> x
        V3(vH,2) = col(3:4:end);   % quaternion j -> y
        V3(vH,3) = col(4:4:end);   % quaternion k -> z
    end
end


%% ===== GUI: standalone quiver viewer =====
function hFig = ViewFigure(EigenFile)
    hFig = [];
    % --- load + validate eigen node ---
    EigenFull = file_fullpath(EigenFile);
    if ~file_exist(EigenFull)
        bst_error('Eigen file not found.', 'View eigenmodes', 0);
        return;
    end
    EigenMat = load(EigenFull);
    if ~isfield(EigenMat,'Variant') || isempty(EigenMat.Variant)
        bst_error('Eigen file has no Variant field.', 'View eigenmodes', 0);
        return;
    end
    % --- variant dispatch (Dirac implemented; others deferred) ---
    switch lower(EigenMat.Variant)
        case 'dirac'
            % implemented below
        case {'laplace-beltrami','connection laplacian'}
            bst_error(sprintf(['Vector viewer currently supports Dirac eigenmodes only.' 10 ...
                'This is a "%s" node.'], EigenMat.Variant), 'View eigenmodes', 0);
            return;
        otherwise
            bst_error(sprintf('Unknown eigen variant: %s', EigenMat.Variant), 'View eigenmodes', 0);
            return;
    end
    if ~isfield(EigenMat,'Phi') || isempty(EigenMat.Phi) || numel(EigenMat.Phi) ~= 2 ...
            || isempty(EigenMat.Phi{1}) || isempty(EigenMat.Phi{2})
        bst_error('Dirac eigen node has empty Phi.', 'View eigenmodes', 0);
        return;
    end

    Surface  = EigenMat.ParentSurface;
    TessMat  = in_tess_bst(Surface);
    Vertices = TessMat.Vertices;
    nVert    = size(Vertices, 1);
    K        = min(size(EigenMat.Phi{1},2), size(EigenMat.Phi{2},2));
    Tau      = NaN;
    if isfield(EigenMat,'Provenance') && isstruct(EigenMat.Provenance) ...
            && isfield(EigenMat.Provenance,'Tau') && ~isempty(EigenMat.Provenance.Tau)
        Tau = EigenMat.Provenance.Tau;
    end

    % --- display cortex (translucent gray) ---
    hFig = view_surface(Surface, 0.5, [0.5 0.5 0.5], 'NewFigure');
    if isempty(hFig)
        bst_error('Could not open the surface figure.', 'View eigenmodes', 0);
        return;
    end
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    % CRITICAL: hold the axes so quiver3 (high-level) does not run newplot and
    % reset the 'Axes3D' axes (which would delete the cortex patch).
    hold(hAxes, 'on');
    figure_3d('SetStandardView', hFig, 'left');
    set(hFig, 'Name', ['Eigenmodes: ' EigenMat.Variant ' | ' Surface]);

    % --- state (closure vars) ---
    iMode              = 1;
    quiverSize         = 1;
    quiverWidth        = 1;
    thresholdAmplitude = 1;     % fraction of cumulative norm kept (1 = all)
    thresholdBalance   = 0;     % 0 = keep small (<=), 1 = keep large (>)
    useNormalize       = false;

    % --- legend ---
    hLabel = uicontrol('Style','text','String','...','Units','Pixels', ...
        'Position',[6 1 1600 35],'HorizontalAlignment','left', ...
        'FontUnits','points','FontSize',bst_get('FigFont'), ...
        'ForegroundColor',[.3 1 .3],'BackgroundColor',[0 0 0],'Parent',hFig);

    % --- keyboard ---
    KeyPressFcn_bak = get(hFig, 'KeyPressFcn');
    set(hFig, 'KeyPressFcn', @KeyPress_Callback);

    DrawArrows();

    % ===== NESTED: draw the current mode's field =====
    function DrawArrows()
        delete(findobj(hAxes, '-depth', 1, 'Tag', 'eigArrows'));
        V3 = ReconstructModeField(EigenMat, iMode, nVert);
        if useNormalize
            nv = sqrt(sum(V3.^2, 2));
            nz = nv > eps;
            V3(nz,:) = V3(nz,:) ./ nv(nz);
        end
        % cumulative-norm amplitude gate (like view_leadfield_vectors)
        normV = sqrt(sum(V3.^2, 2));
        [sv, ind] = sort(normV, 'ascend');
        cdf = cumsum(sv);
        if cdf(end) > 0, cdf = cdf / cdf(end); end
        if thresholdBalance == 0
            keep = find(cdf <= thresholdAmplitude);
        else
            keep = find(cdf > thresholdAmplitude);
        end
        Vre = Vertices(ind, :);
        Dre = zeros(numel(ind), 3);
        Dre(keep, :) = V3(ind(keep), :);
        quiver3(Vre(:,1), Vre(:,2), Vre(:,3), Dre(:,1), Dre(:,2), Dre(:,3), quiverSize, ...
            'Parent', hAxes, 'LineWidth', quiverWidth, 'Color', [.3 1 .3], 'Tag', 'eigArrows');
        % legend
        lamL = EigenMat.Lambda{1}(iMode);
        lamR = EigenMat.Lambda{2}(iMode);
        tauStr = ''; if ~isnan(Tau), tauStr = sprintf(' | tau=%.3g', Tau); end
        normStr = ''; if useNormalize, normStr = ' | unit'; end
        set(hLabel, 'String', sprintf(['Mode %d / %d   |   lambdaL=%.4g, lambdaR=%.4g%s%s   ' ...
            '[arrows: %d | H for help]'], iMode, K, lamL, lamR, tauStr, normStr, numel(keep)));
    end

    % ===== NESTED: keyboard navigation =====
    function KeyPress_Callback(h, keyEvent)
        switch (keyEvent.Key)
            case 'leftarrow'
                if     ismember('shift',   keyEvent.Modifier), quiverSize  = quiverSize  / 1.2;
                elseif ismember('control', keyEvent.Modifier), quiverWidth = quiverWidth / 1.2;
                elseif ismember('alt',     keyEvent.Modifier), thresholdAmplitude = thresholdAmplitude - 0.01;
                else,  iMode = iMode - 1; end
            case 'rightarrow'
                if     ismember('shift',   keyEvent.Modifier), quiverSize  = quiverSize  * 1.2;
                elseif ismember('control', keyEvent.Modifier), quiverWidth = quiverWidth * 1.2;
                elseif ismember('alt',     keyEvent.Modifier), thresholdAmplitude = thresholdAmplitude + 0.01;
                else,  iMode = iMode + 1; end
            case 'uparrow'
                if     ismember('shift',   keyEvent.Modifier), quiverSize  = quiverSize  * 1.2;
                elseif ismember('control', keyEvent.Modifier), quiverWidth = quiverWidth * 1.2;
                elseif ismember('alt',     keyEvent.Modifier), thresholdAmplitude = thresholdAmplitude + 0.01;
                else,  return; end
            case 'downarrow'
                if     ismember('shift',   keyEvent.Modifier), quiverSize  = quiverSize  / 1.2;
                elseif ismember('control', keyEvent.Modifier), quiverWidth = quiverWidth / 1.2;
                elseif ismember('alt',     keyEvent.Modifier), thresholdAmplitude = thresholdAmplitude - 0.01;
                else,  return; end
            case 'pageup',   iMode = iMode + 10;
            case 'pagedown', iMode = iMode - 10;
            case 'n',        useNormalize = ~useNormalize;
            case 'return'
                if ismember('alt', keyEvent.Modifier), thresholdBalance = ~thresholdBalance; else, return; end
            case 'h'
                java_dialog('msgbox', ['<HTML><TABLE>' ...
                    '<TR><TD><B>Left/Right</B></TD><TD>Previous/next mode</TD></TR>' ...
                    '<TR><TD><B>PgUp/PgDn</B></TD><TD>+/- 10 modes</TD></TR>' ...
                    '<TR><TD><B>Shift+Left/Right</B></TD><TD>Arrow length -/+</TD></TR>' ...
                    '<TR><TD><B>Control+Left/Right</B></TD><TD>Arrow width -/+</TD></TR>' ...
                    '<TR><TD><B>Alt+Left/Right</B></TD><TD>Amplitude threshold -/+</TD></TR>' ...
                    '<TR><TD><B>Alt+Enter</B></TD><TD>Toggle threshold direction</TD></TR>' ...
                    '<TR><TD><B>N</B></TD><TD>Toggle unit-normalized arrows</TD></TR>' ...
                    '<TR><TD><B>0-9</B></TD><TD>Change view</TD></TR>' ...
                    '</TABLE>'], 'Keyboard shortcuts', [], 0);
                return;
            otherwise
                if ~isempty(KeyPressFcn_bak), KeyPressFcn_bak(h, keyEvent); end
                return;
        end
        if iMode < 1, iMode = K; end
        if iMode > K, iMode = 1; end
        if thresholdAmplitude < 0, thresholdAmplitude = 0; end
        if thresholdAmplitude > 1, thresholdAmplitude = 1; end
        DrawArrows();
    end
end
