function J = bst_dirac_eigenmode_field(DiracEigen, GainRows, CompHM)
% BST_DIRAC_EIGENMODE_FIELD: Reconstruct cortical 3-vector field(s) from Dirac
% eigenmode leadfield coefficients (forward change-of-basis, back to vertex space).
%
% USAGE:  J = bst_dirac_eigenmode_field(DiracEigen, GainRows, CompHM)
%
% DESCRIPTION:
%     For each row of GainRows (one channel's Dirac eigenmode leadfield, a [1 x 2K]
%     vector of per-mode coefficients), map it back through the Dirac eigenvectors
%     to a per-vertex 3-vector field on the cortex:
%         psi_h = Phi_D,h * gainRow(modes of hemisphere h)'   % [4Vh x 1] real
%         J(v, :) = imag-part = [psi_h(4(v-1)+2), +3, +4]     % drop the w slot
%     The two hemispheres are scattered to the full cortex. The result is that
%     channel's leadfield band-limited to the K Dirac eigenmodes (-> the raw
%     unconstrained leadfield as K grows). Pure forward operation (no inverse).
%
% INPUTS:
%   DiracEigen : 1x2 per-hemisphere struct (TessMat.DiracEigen): .Vectors [4Vh x K]
%   GainRows   : [m x 2K] eigenmode-leadfield rows (m channels; columns = modes,
%                stacked L then R as in the composed head model Gain)
%   CompHM     : composed head model with .nModes (=2K), .ModeHemisphere [2K x 1],
%                .HemiGlobalVertices {Lverts, Rverts}
% OUTPUT:
%   J : [m x 3*nVert] reconstructed cortical leadfield (x,y,z per vertex):
%       J(c, 3*(v-1)+(1:3)) is channel c's reconstructed 3-vector at vertex v.

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

    Kh    = CompHM.nModes / 2;
    nVert = sum(cellfun(@numel, CompHM.HemiGlobalVertices));
    m     = size(GainRows, 1);
    J     = zeros(m, 3*nVert);
    for hh = 1:2
        vH   = CompHM.HemiGlobalVertices{hh}(:);
        Phi  = double(DiracEigen(hh).Vectors(:, 1:Kh));     % [4Vh x Kh]
        cols = (CompHM.ModeHemisphere(:) == hh);            % [2K x 1] logical
        R    = Phi * GainRows(:, cols).';                   % [4Vh x m]
        J(:, (vH-1)*3 + 1) = R(2:4:end, :).';               % x  (drop w rows 1:4:end)
        J(:, (vH-1)*3 + 2) = R(3:4:end, :).';               % y
        J(:, (vH-1)*3 + 3) = R(4:4:end, :).';               % z
    end
end
