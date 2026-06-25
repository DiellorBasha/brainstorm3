function varargout = process_helmholtz( varargin )
% PROCESS_HELMHOLTZ: Helmholtz-Hodge decomposition of a 3-D cortical source vector field.
%
% Compute (pure, stateless, I/O-free) is the flat-covariant algorithm: it sequences
% bst_divergence + bst_curl (ambient div/curl from the Covariant operator) and recovers the
% scalar potentials phi/psi via the weak Hodge solve (bst_poisson, factor cached by
% tess_cholesky), plus the irrotational/solenoidal/harmonic reconstruction. Run loops Compute
% over a source series and saves results maps. Same algorithm for the GUI ephemeral feedback
% (view_helmholtz / panel_bst_dynamics) and the on-file save.
%
% USAGE:  Ht = process_helmholtz('Compute', J, Cov)   % J [3nV x nT], Cov = 'Covariant' node
%         OutputFiles = process_helmholtz('Run', sProcess, sInputs)
%
% Authors: Diellor Basha, 2026

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

eval(macro_method);
end

%% ===== EXTERNAL CALL: pure flat-covariant Helmholtz of a 3-D source frame =====
% Ported verbatim from the old bst_helmholtz i_prepare_vertex/i_frame_vertex (vertex domain).
function Ht = Compute(J, Cov)
    s = +1;
    nVtot = max(cellfun(@(c) max(double(c(:))), Cov.GlobalVertices));
    nT = size(J, 2);
    z1 = zeros(nVtot, nT);  z3 = zeros(nVtot, 3);
    Ht = struct('Div',z1, 'Curl',z1, 'Phi',z1, 'Psi',z1, 'Fmag',z1, 'Hmag',z1, ...
                'Vtot',z3, 'Virr',z3, 'Vsol',z3, 'Hresid',z3, 'HarmFrac',0);
    % Strong fields from the dedicated engines (single home of div/curl).
    Ht.Div  = bst_divergence(J, [], 'Ambient', [], Cov);
    Ht.Curl = bst_curl(J, [], 'Ambient', [], Cov);
    harmNum = 0;  harmDen = 0;
    for hh = 1:numel(Cov.Covariant)
        C = Cov.Covariant{hh};  vH = double(Cov.GlobalVertices{hh}(:));
        nFh = size(C.Faces, 1);  nVh = numel(vH);
        Gx = C.ScalarGrad(1:nFh,:);  Gy = C.ScalarGrad(nFh+1:2*nFh,:);  Gz = C.ScalarGrad(2*nFh+1:3*nFh,:);
        Nf = C.FaceNormal;  Af = C.FaceArea;  nv = C.Frame.normal;
        W   = spdiags(Af, 0, nFh, nFh);
        Wfv = bst_face2vertex(C.Faces, Af);           % shared math helper
        Fvf = sparse([(1:nFh)';(1:nFh)';(1:nFh)'], [C.Faces(:,1);C.Faces(:,2);C.Faces(:,3)], 1/3, nFh, nVh);
        nx=Nf(:,1); ny=Nf(:,2); nz=Nf(:,3);
        Sx = spdiags(ny,0,nFh,nFh)*Gz - spdiags(nz,0,nFh,nFh)*Gy;
        Sy = spdiags(nz,0,nFh,nFh)*Gx - spdiags(nx,0,nFh,nFh)*Gz;
        Sz = spdiags(nx,0,nFh,nFh)*Gy - spdiags(ny,0,nFh,nFh)*Gx;
        dF = tess_cholesky(Cov, hh, 1);                 % cached pinned factor (pin vertex 1 => free=2:nVh; matches bst_poisson)
        for t = 1:nT
            Jx = J(3*(vH-1)+1, t);  Jy = J(3*(vH-1)+2, t);  Jz = J(3*(vH-1)+3, t);  Jv = [Jx Jy Jz];
            Jf = [Fvf*Jx, Fvf*Jy, Fvf*Jz];
            divw  = s * (Gx'*W*Jf(:,1) + Gy'*W*Jf(:,2) + Gz'*W*Jf(:,3));   % weak divergence source
            vortw = s * (Sx'*W*Jf(:,1) + Sy'*W*Jf(:,2) + Sz'*W*Jf(:,3));   % weak vorticity source
            phi = tess_cholesky('solve', dF, divw);   phi = phi - mean(phi);
            psi = tess_cholesky('solve', dF, vortw);  psi = psi - mean(psi);
            Virr = Wfv * (s * [Gx*phi, Gy*phi, Gz*phi]);
            Vsol = Wfv * (s * cross(Nf, [Gx*psi, Gy*psi, Gz*psi], 2));
            Jn   = Jx.*nv(:,1) + Jy.*nv(:,2) + Jz.*nv(:,3);
            Hres = Jv - Virr - Vsol - Jn.*nv;
            Ht.Phi(vH,t)=phi;  Ht.Psi(vH,t)=psi;  Ht.Fmag(vH,t)=sqrt(Jx.^2+Jy.^2+Jz.^2);
            Ht.Hmag(vH,t)=sqrt(sum(Hres.^2,2));
            if t == nT     % last-frame vector fields (single-frame contract)
                Ht.Vtot(vH,:)=Jv;  Ht.Virr(vH,:)=Virr;  Ht.Vsol(vH,:)=Vsol;  Ht.Hresid(vH,:)=Hres;
                av = full(sum(Cov.Mass{hh}, 2));
                harmNum = harmNum + sum(av .* sum(Hres.^2,2));
                harmDen = harmDen + sum(av .* sum(Jv.^2,2));
            end
        end
    end
    if harmDen > 0, Ht.HarmFrac = harmNum / harmDen; else, Ht.HarmFrac = 0; end
end
