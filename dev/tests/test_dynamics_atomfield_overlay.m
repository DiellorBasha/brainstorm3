% test_dynamics_atomfield_overlay - the atom-field scatter/frame-pick used by the overlay (pure verb)
nV = 50; k = 8; nT = 10;
gv = (3:3:3*k)';                       % 8 support vertices scattered into a 50-vertex surface
W  = zeros(k, nT);  W(4, :) = 1:nT;    % support row 4 ramps 1..nT over time
s5 = view_dynamics('AtomFrameScalar', W, gv, 5, nV);
assert(numel(s5)==nV, 'scal is full-surface length');
assert(s5(gv(4))==5, 'frame 5 value scattered to its global vertex');
assert(nnz(s5)==1 && all(s5(setdiff((1:nV)', gv(4)))==0), 'only the support vertex is nonzero');
sHi = view_dynamics('AtomFrameScalar', W, gv, 999, nV);
assert(sHi(gv(4))==nT, 'out-of-range frame index clamps to the last frame');
sLo = view_dynamics('AtomFrameScalar', W, gv, 0, nV);
assert(sLo(gv(4))==1, 'frame index < 1 clamps to the first frame');
disp('OK');
