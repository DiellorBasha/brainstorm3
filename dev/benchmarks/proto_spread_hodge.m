function R = proto_spread_hodge(SurfaceFile, A, Opts)
% PROTO_SPREAD_HODGE: Layer-2 differential-geometry flow analysis of a longitudinal surface field
% series A=[nVert x nT]. Uses the DEC operators on a manifold node (tess_manifold):
%   gradient (bst_gradient)  -> the per-face deposition FLOW direction
%   divergence (bst_divergence) -> SOURCES/sinks (the deposition foci; gradient converges at maxima)
%   curl (bst_curl)          -> circulation (Hodge: a pure gradient is irrotational ~0)
% Reported per timepoint, plus the gradient of the ACCUMULATION field (a(end)-a(1)) = where amyloid
% grows + in what direction. Also projects the flow onto the CONNECTION-LAPLACIAN eigenmodes (the
% vector-field basis) for a spectral flow readout.
%
% USAGE: R = proto_spread_hodge(SurfaceFile, A, Opts)   Opts: .nVecModes=60
%
% Author: Diellor Basha, 2026 (prototype)
    if nargin<3, Opts=struct(); end
    if ~isfield(Opts,'nVecModes'), Opts.nVecModes=60; end
    Mani=tess_manifold(SurfaceFile);
    sW=in_tess_bst(SurfaceFile); nV=size(sW.Vertices,1); nT=size(A,2);
    src=zeros(nV,nT); curlFrac=zeros(1,nT); epiSrc=zeros(1,nT);
    for t=1:nT
        a=A(:,t); a(~isfinite(a))=0;
        G=bst_gradient(a,Mani); Dv=bst_divergence(G,Mani); Cu=bst_curl(G,Mani);
        src(:,t)=Dv; curlFrac(t)=norm(Cu)/max(norm(Dv),eps);   % rotational vs source content
        [~,epiSrc(t)]=min(Dv);                                  % gradient converges (div<0) at the focus
    end
    % accumulation field flow (a(end) - a(1))
    dA=A(:,end)-A(:,1); dA(~isfinite(dA))=0;
    Gacc=bst_gradient(dA,Mani); accDiv=bst_divergence(Gacc,Mani); accCurl=bst_curl(Gacc,Mani);
    [~,accFocus]=min(accDiv);                                   % focus of accumulation
    accCurlFrac=norm(accCurl)/max(norm(accDiv),eps);
    % connection-Laplacian eigenmodes: spectral content of the accumulation flow (energy by mode)
    vecEnergy=[];
    try
        Evec=tess_eigen(SurfaceFile,'Connection Laplacian','nModes',Opts.nVecModes);
        vecEnergy=local_vecspectrum(Gacc, Evec, Mani);
    catch ME, fprintf('  connection-Laplacian projection skipped: %s\n', regexprep(ME.message,'\s+',' ')); end
    R=struct('src',src,'curlFrac',curlFrac,'epiSrc',epiSrc,'accDiv',accDiv,'accFocus',accFocus, ...
             'accCurlFrac',accCurlFrac,'Gacc',Gacc,'vecEnergy',vecEnergy,'V',sW.Vertices);
end

function e=local_vecspectrum(G, Evec, Mani) %#ok<INUSD>
    % crude spectral energy of the face vector field by per-mode inner products if the eigenbasis
    % exposes face vector modes; otherwise returns [] (placeholder for the full connection-Laplacian path).
    e=[];
    if isfield(Evec,'Phi') && ~isempty(Evec.Phi)
        % connection-Laplacian Phi are complex vertex fields; map is non-trivial -> defer, report norm only
        e=struct('note','connection-Laplacian eigenbasis available; full vector projection deferred', ...
                 'nModes', sum(cellfun(@(p)size(p,2),Evec.Phi)));
    end
end
