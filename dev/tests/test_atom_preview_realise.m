% test_atom_preview_realise - the atom-preview normalization (density vs peak), pure verb
nGv = 40; nT = 12;
M = speye(nGv);                                   % unit vertex areas -> mass-weighted sum = plain sum
% (1) one-signed mass-conserving field -> DENSITY: each frame integrates (mass-weighted) to 1
Wpos = abs(reshape(sin(1:(nGv*nT)), nGv, nT)) + 0.5;   % strictly positive
[Wd, sgnD] = panel_bst_dynamics('i_atom_normalize', Wpos, M);
assert(sgnD==false, 'one-signed field -> density (not signed)');
massPerFrame = full(sum(M,2)).' * Wd;
assert(all(abs(massPerFrame - 1) < 1e-9), 'each frame integrates to unit mass');
% (2) zero-mean oscillatory field -> PEAK: max|.| == 1, flagged signed
Wsig = reshape(sin(1:(nGv*nT)), nGv, nT);  Wsig = Wsig - mean(Wsig,1);   % per-frame zero mean
[Wp, sgnP] = panel_bst_dynamics('i_atom_normalize', Wsig, M);
assert(sgnP==true, 'zero-mean field -> peak (signed)');
assert(abs(max(abs(Wp(:))) - 1) < 1e-12, 'peak-normalized to unit max amplitude');
disp('OK');
