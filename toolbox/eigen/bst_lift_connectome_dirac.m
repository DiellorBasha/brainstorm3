function [Phiq, Lamq, Mq] = bst_lift_connectome_dirac(Phis, Lams, Ms)
% BST_LIFT_CONNECTOME_DIRAC: Lift a scalar eigenbasis to an ambient quaternion (Dirac) basis.
%   [Phiq,Lamq,Mq] = bst_lift_connectome_dirac(Phis,Lams,Ms)
%
% USAGE:
%   [Phiq, Lamq, Mq] = bst_lift_connectome_dirac(Phis, Lams, Ms)
%
% DESCRIPTION:
%     Each scalar mode phi_k becomes three quaternion modes with phi_k in the imaginary x/y/z slots
%     (w=0), eigenvalue lam_k each; mass Mq = kron(Ms,I4). Realizes L_conn (x) I3 (ambient,
%     component-wise).
%
% INPUTS:
%     Phis [n x K]       - Scalar eigenbasis matrix (n vertices, K modes)
%     Lams [K x 1]       - Scalar eigenvalues
%     Ms   [n x n]       - Scalar mass matrix (e.g. vertex areas)
%
% OUTPUTS:
%     Phiq [4n x 3K]     - Lifted quaternion eigenbasis (interleaved [w;x;y;z] per vertex)
%     Lamq [3K x 1]      - Lifted eigenvalues (each Lams value repeated 3 times)
%     Mq   [4n x 4n]     - Lifted mass matrix (kron(Ms, I4))
%
% SEE ALSO: bst_eigen, bst_eigenfilter, manifold_quat_imag

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

    Phis = double(Phis);  [n, K] = size(Phis);
    Phiq = zeros(4*n, 3*K);
    Lamq = zeros(3*K, 1);
    for k = 1:K
        for d = 1:3                                  % d=1,2,3 -> imag x,y,z at rows (v-1)*4 + (d+1)
            c = (k-1)*3 + d;
            Phiq((0:n-1)*4 + (d+1), c) = Phis(:, k);
            Lamq(c) = Lams(k);
        end
    end
    Mq = kron(Ms, speye(4));
end
