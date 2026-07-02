function V = manifold_quat_imag(Q)
% MANIFOLD_QUAT_IMAG: Imaginary 3-vector of an interleaved quaternion field.
%
% USAGE:  V = manifold_quat_imag(Q)
%
% INPUT:
%    - Q : [4*n x T] interleaved quaternion rows [w1,x1,y1,z1, w2,x2,y2,z2, ...]
% OUTPUT:
%    - V : [3*n x T] imaginary 3-vector rows      [x1,y1,z1, x2,y2,z2, ...]
%
% This is the inverse of embedding a per-vertex 3-vector field into the imaginary
% quaternion slots (w=0): it extracts the physical current from a Dirac quaternion
% field. Shared by the Dirac source-mapping / atom / scalogram paths that all use
% the same interleaved [w,x,y,z] convention (bst_dirac, bst_eigen, bst_eigenwavelet,
% panel_bst_dynamics).
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

if mod(size(Q,1), 4) ~= 0
    error('manifold_quat_imag:size', 'Q must have 4*n rows (got %d).', size(Q,1));
end
n  = size(Q,1) / 4;
Qr = reshape(Q, 4, n, []);            % [4 x n x T]
V  = reshape(Qr(2:4,:,:), 3*n, []);   % [3n x T]  rows [x1,y1,z1, x2,...]
end
