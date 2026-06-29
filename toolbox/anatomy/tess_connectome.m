function R = tess_connectome(SurfaceFile, varargin)
% TESS_CONNECTOME: Build a structural connectome on a cortical surface from tractography fibers.
%
% USAGE:  R = tess_connectome(SurfaceFile)
%         R = tess_connectome(SurfaceFile, 'Resolution','vertex', 'SmoothHops',3, ...
%                             'Atlas','Desikan-Killiany', 'EdgeThr',1e-4)
%
% DESCRIPTION:
%     Assembles the structural connectome of a cortex from fiber endpoints. The fibers are the
%     subject's OWN (a Fibers surface in the subject's anatomy) when available; otherwise the protocol
%     default HCP-1065 template fibers are co-registered to this surface via the FreeSurfer registration
%     sphere (tess_interp_tess2tess) - i.e. the subject is "assumed" to carry the HCP connectome,
%     projected onto its own cortex. This reuses Brainstorm's fiber machinery rather than duplicating it:
%       - fibers_helper('AssignToScouts')  for the ROI endpoint->scout assignment,
%       - tess_interp_tess2tess            for the template->subject registration.
%
%     Two resolutions:
%       'roi'    : region x region connectome over an atlas (default Desikan-Killiany); C(a,b) = number
%                  of fibers whose endpoints fall in scouts a and b (via AssignToScouts centroids).
%       'vertex' : vertex x vertex connectome; endpoints -> nearest vertex (dsearchn), then a
%                  continuous-kernel mesh smoothing W = S^p Wraw S^p' bridges the sulcal disconnection
%                  (raw endpoints leave ~46% of vertices isolated); the largest connected component is
%                  returned in .keep.
%
% OPTIONS:
%     'Resolution'  : 'vertex' (default) | 'roi'
%     'Atlas'       : atlas name for 'roi' (default 'Desikan-Killiany')
%     'SmoothHops'  : mesh-smoothing hops for 'vertex' (default 3, ~10mm)
%     'EdgeThr'     : drop edges below this fraction of max weight, 'vertex' (default 1e-4)
%
% OUTPUT:
%     R : struct with fields
%         W           : sparse connectome (region x region for 'roi'; vertex x vertex for 'vertex')
%         Resolution  : 'roi' | 'vertex'
%         keep        : largest-CC vertex indices ('vertex' only; [] for 'roi')
%         Labels      : scout labels ('roi' only; [] for 'vertex')
%         Provenance  : struct(fibers, resolution, atlas, nFibers, registration)
%
% SEE ALSO: fibers_helper, tess_interp_tess2tess, tess_operators, proto_fiber_connectome

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
% WARRANTY, EXPRESS OR IMPLIED.
% =============================================================================@
%
% Author: Diellor Basha, 2026

    % ===== OPTIONS =====
    Opts = struct('Resolution','vertex','Atlas','Desikan-Killiany','SmoothHops',3,'EdgeThr',1e-4);
    for i = 1:2:numel(varargin), Opts.(varargin{i}) = varargin{i+1}; end

    sTess = in_tess_bst(SurfaceFile); V = sTess.Vertices; F = sTess.Faces; nV = size(V,1);
    sSubj = bst_get('SurfaceFile', SurfaceFile);

    % ===== RESOLVE FIBER ENDPOINTS ON THIS SURFACE =====
    [endPts, regStr, nFib] = local_get_endpoints(sSubj, SurfaceFile, V);

    % ===== ASSEMBLE CONNECTOME =====
    R = struct('W',[],'Resolution',lower(Opts.Resolution),'keep',[],'Labels',[],'Provenance',[]);
    switch lower(Opts.Resolution)
        case 'roi'
            [R.W, R.Labels] = local_roi_connectome(sTess, Opts.Atlas, endPts);
        case 'vertex'
            [R.W, R.keep]   = local_vertex_connectome(V, F, endPts, Opts);
        otherwise
            error('Unknown Resolution "%s" (use ''roi'' or ''vertex'').', Opts.Resolution);
    end
    R.Provenance = struct('fibers',regStr,'resolution',R.Resolution,'atlas',Opts.Atlas, ...
                          'nFibers',nFib,'registration',regStr);
    fprintf('tess_connectome [%s]: %d fibers, %s | %d edges%s\n', R.Resolution, nFib, regStr, ...
        nnz(triu(R.W)), local_ifelse(strcmpi(R.Resolution,'vertex'), sprintf(' | %d-vtx CC', numel(R.keep)), ''));
end

% ---- fiber endpoints on THIS surface (own fibers, or default HCP-1065 via registration) ----
function [endPts, regStr, nFib] = local_get_endpoints(sSubj, SurfaceFile, V)
    iFib = find(strcmpi({sSubj.Surface.SurfaceType},'Fibers'), 1);
    if ~isempty(iFib)
        FibMat = load(file_fullpath(sSubj.Surface(iFib).FileName),'Points');
        endPts = double(FibMat.Points(:,[1 end],:)); regStr = ['own:' sSubj.Surface(iFib).Comment];
    else
        % default HCP-1065 template fibers
        sDef = bst_get('Subject',0); iFibD = find(strcmpi({sDef.Surface.SurfaceType},'Fibers'),1);
        if isempty(iFibD), error('No fibers on this subject and none on the default anatomy.'); end
        FibMat = load(file_fullpath(sDef.Surface(iFibD).FileName),'Points');
        ep = double(FibMat.Points(:,[1 end],:));
        % default cortex (mid, low res) as the registration source
        icD = find(~cellfun('isempty',regexp({sDef.Surface.FileName},'cortex_mid_low\.mat$','once')),1);
        if isempty(icD), icD = sDef.iCortex; end
        defCortex = sDef.Surface(icD).FileName;
        if strcmp(file_fullpath(defCortex), file_fullpath(SurfaceFile))
            endPts = ep; regStr = 'default:HCP-1065';                       % same surface, no registration
        else
            endPts = local_project_endpoints(ep, defCortex, SurfaceFile, V);
            regStr = 'default:HCP-1065 (registered)';
        end
    end
    nFib = size(endPts,1);
end

% ---- project template endpoints onto the subject cortex via the registration sphere ----
function endPts = local_project_endpoints(ep, srcCortex, destSurfFile, destV)
    Wmat = tess_interp_tess2tess(srcCortex, destSurfFile, 0, 0);     % [nDest x nSrc] sphere interpolation
    srcV = in_tess_bst(srcCortex).Vertices;
    [~, src2dest] = max(Wmat, [], 1);                                % best dest vertex per src vertex
    e1 = src2dest(dsearchn(srcV, squeeze(ep(:,1,:))));
    e2 = src2dest(dsearchn(srcV, squeeze(ep(:,2,:))));
    endPts = cat(2, reshape(destV(e1,:),[],1,3), reshape(destV(e2,:),[],1,3));
end

% ---- ROI connectome via fibers_helper('AssignToScouts') ----
function [C, labels] = local_roi_connectome(sTess, atlasName, endPts)
    ai = find(strcmp({sTess.Atlas.Name}, atlasName), 1);
    if isempty(ai), error('Atlas "%s" not found on the surface.', atlasName); end
    Sc = sTess.Atlas(ai).Scouts; nReg = numel(Sc);
    centroids = zeros(nReg,3);
    for r = 1:nReg, centroids(r,:) = mean(sTess.Vertices(Sc(r).Vertices,:),1); end
    FibMat = struct('Points',endPts,'Scouts',struct('ConnectFile',[],'Assignment',[]));
    FibMat = fibers_helper('AssignToScouts', FibMat, '', centroids);
    asg = FibMat.Scouts(end).Assignment;
    ok = asg(:,1)>0 & asg(:,2)>0 & asg(:,1)~=asg(:,2);
    C = accumarray(asg(ok,:), 1, [nReg nReg]); C = C + C';
    labels = {Sc.Label};
end

% ---- vertex connectome: nearest-vertex + continuous-kernel smoothing ----
function [W, keep] = local_vertex_connectome(V, F, endPts, Opts)
    nV = size(V,1);
    v1 = dsearchn(V, squeeze(endPts(:,1,:))); v2 = dsearchn(V, squeeze(endPts(:,2,:)));
    Wraw = sparse([v1;v2],[v2;v1],1,nV,nV); Wraw = Wraw - spdiags(diag(Wraw),0,nV,nV);
    ii = [F(:,1);F(:,2);F(:,3)]; jj = [F(:,2);F(:,3);F(:,1)];
    A = double(sparse([ii;jj],[jj;ii],1,nV,nV) > 0);
    Sm = A + speye(nV); Sm = spdiags(1./sum(Sm,2),0,nV,nV)*Sm; Sp = Sm^Opts.SmoothHops;
    W = Sp*Wraw*Sp'; W = (W+W')/2; W = W - spdiags(diag(W),0,nV,nV);
    W(W < Opts.EdgeThr*max(W(:))) = 0;
    cc = conncomp(graph(W>0)); tb = tabulate(cc); [~,big] = max(tb(:,2)); keep = find(cc==big);
end

function s = local_ifelse(c,a,b), if c, s=a; else, s=b; end, end
