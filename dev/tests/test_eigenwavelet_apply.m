function test_eigenwavelet_apply(EigenFile)
% Analysis/Synthesis: shapes, tight-frame round-trip, heat scale-space, localized atom.

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

if nargin < 1 || isempty(EigenFile); EigenFile = 'Subject01/eigen_260622_1920.mat'; end
E = in_bst_eigen(EigenFile); O = in_bst_operator(E.OperatorFile);
% Frame range must span ALL hemispheres (a frame applied to a hemi whose lambda exceeds
% the design lmax would zero those modes). Use the global [min, max].
lams = E.Lambda(~cellfun(@isempty, E.Lambda));
lr = [min(cellfun(@min, lams)), max(cellfun(@max, lams))];
nSrc = max(cellfun(@(v) max([v(:);0]), E.GlobalVertices));
% in-subspace signal so a tight frame reconstructs exactly (no truncation residual)
rng(3);
Fin = zeros(nSrc, 4);
for hh=1:numel(E.Phi)
    Phi=E.Phi{hh}; if isempty(Phi); continue; end
    Fin(E.GlobalVertices{hh}(:),:) = Phi*randn(size(Phi,2),4);
end
% --- Analysis shape ---
fr = bst_eigenwavelet('Design','itersine',8,lr);
W = bst_eigenwavelet('Analysis', Fin, E, O, fr);
assert(isequal(size(W),[nSrc,4,numel(fr.g)]), 'Analysis shape wrong.');
% --- Tight-frame round-trip: synthesis(analysis(f)) ~ f for in-subspace f ---
Frec = bst_eigenwavelet('Synthesis', W, E, O, fr);
relerr = norm(Frec(:)-Fin(:))/norm(Fin(:));
assert(relerr < 1e-8, sprintf('itersine round-trip not exact: relerr=%.2e', relerr));
% --- Heat scale-space: coarsest member smoother than finest (less high-mode energy) ---
fh = bst_eigenwavelet('Design','heat',5,lr);
Wh = bst_eigenwavelet('Analysis', Fin, E, O, fh);
hi_coarse = i_highmode_energy(Wh(:,:,1), E, O);
hi_fine   = i_highmode_energy(Wh(:,:,end), E, O);
assert(hi_coarse < hi_fine, 'heat coarsest member should hold less high-mode energy.');
% --- Atom: synthesis of a single-member delta is spatially localized at the seed ---
v0 = E.GlobalVertices{1}(round(numel(E.GlobalVertices{1})/2));
Wd = zeros(nSrc,1,numel(fr.g)); Wd(v0,1,4)=1;   % delta in mid-band member 4
atom = bst_eigenwavelet('Synthesis', Wd, E, O, fr);
[sv,si] = sort(abs(atom),'descend');  %#ok<ASGLU>
assert(ismember(v0, si(1:5)), 'atom not localized near the seed vertex.');
disp('test_eigenwavelet_apply PASSED');
end

function e = i_highmode_energy(F, E, O)
    e = 0;
    for hh=1:numel(E.Phi)
        Phi=E.Phi{hh}; if isempty(Phi); continue; end
        idx=E.GlobalVertices{hh}(:); C = Phi'*(O.Mass{hh}*F(idx,:));
        K=size(C,1); hi=round(K/2):K; e = e + sum(sum(C(hi,:).^2));
    end
end
