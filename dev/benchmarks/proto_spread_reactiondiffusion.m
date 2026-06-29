function R = proto_spread_reactiondiffusion(SurfaceFile, A, dt)
% PROTO_SPREAD_REACTIONDIFFUSION: Layer-3 single-tracer reaction-diffusion fit of a longitudinal
% surface field series. Recovers growth r_a and diffusion D_a (via pet_spread_invert's single-tracer
% path), the Fisher-KPP front width sqrt(D_a/r_a), and the marginal eigenvalue lambda*=r_a/D_a.
% Cross-checks lambda* against the field's gradient spectral scale (the Layer-1<->Layer-3 link).
%
% USAGE: R = proto_spread_reactiondiffusion(SurfaceFile, A, dt)   % A = [nLH x nT] or [nVert x nT]
%
% Author: Diellor Basha, 2026 (prototype)
    Op=tess_operators(SurfaceFile,'Laplace-Beltrami'); gv=double(Op.GlobalVertices{1}(:));
    if size(A,1)~=numel(gv), A=A(gv,:); end                 % reduce a full-surface field to LH
    A(~isfinite(A))=0;
    est=pet_spread_invert(SurfaceFile, A, A, dt);            % single-tracer: ra,Da from the Abeta path
    if est.Da>0 && est.ra>0, frontMM=sqrt(est.Da/est.ra)*1000; lstar=est.ra/est.Da;
    else, frontMM=NaN; lstar=NaN; end
    % spectral cross-check: lambda* vs the gradient-spectral scale of the mid-spread field
    M=Op.Mass{1}; K=Op.Operator{1};
    Eig=tess_eigen(SurfaceFile,'Laplace-Beltrami','nModes',300); Phi=Eig.Phi{1}; Lam=Eig.Lambda{1}(:);
    fmid=A(:,max(1,round(size(A,2)/2))); g=M\(K*fmid); c=Phi'*(M*g); e=c.^2;
    lpeak=sum(Lam.*e)/max(sum(e),eps);
    reactionDominated = isfinite(frontMM) && frontMM < 3;   % front < ~3mm => negligible diffusion
    R=struct('ra',est.ra,'Da',est.Da,'frontMM',frontMM,'lstar',lstar,'lpeak',lpeak, ...
             'reactionDominated',reactionDominated);
    fprintf('reaction-diffusion fit: ra=%.3f Da=%.2e | front width=%.1fmm | lambda*=%.0f lpeak=%.0f | %s\n', ...
        est.ra, est.Da, frontMM, lstar, lpeak, local_verdict(reactionDominated, est.Da));
end

function s=local_verdict(rd, Da)
    if Da<=0, s='UNRELIABLE (Da<=0; sparse-time regime)';
    elseif rd, s='REACTION-DOMINATED (diffusion ~0; local growth, no spread)';
    else, s='reaction-diffusion (genuine spatial spread)'; end
end
