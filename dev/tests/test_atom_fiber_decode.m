function test_atom_fiber_decode
% Byte-equivalence: bst_eigenfilter('Atom') V3 output == the panel's former i_atom_realise_core decode,
% for the scalar (Laplace-Beltrami) and quaternion (Dirac) fibers. V3=[] for scalar; imag-3-vec for Dirac.
    surf = getenv('BST_TEST_SURF');
    assert(~isempty(surf), 'Set BST_TEST_SURF to a cortex surface file (skips if unset).');

    kp = struct('lmax', []);
    % ---- scalar fiber ----
    axS  = i_axes(surf, 'Laplace-Beltrami');
    seedS = axS.GlobalVertices{1}(1);
    [Ws, gvs, V3s, sgnS] = bst_eigenfilter('Atom', axS, 'heat', kp, seedS, []);
    assert(isempty(V3s), 'scalar fiber must yield V3=[]');
    assert(isempty(sgnS), 'scalar fiber must yield isSigned=[]');
    [Ws2, gvs2] = bst_eigenfilter('Atom', axS, 'heat', kp, seedS, []);   % backward-compat 2-output
    assert(isequal(Ws, Ws2) && isequal(gvs, gvs2), 'scalar [W,gv] must be unchanged');

    % ---- quaternion fiber ----
    axQ  = i_axes(surf, 'Dirac');
    seedQ = axQ.GlobalVertices{1}(1);
    dir   = [1 0 0];
    [Wq, gvq, V3q, sgnQ] = bst_eigenfilter('Atom', axQ, 'heat', kp, seedQ, dir);
    V3ref = i_ref_quat_decode(Wq, gvq, axQ);                             % former panel decode, reproduced
    assert(isequal(V3q, V3ref), 'quaternion V3 must byte-match the former decode');
    assert(size(V3q,1) == i_nsrc(axQ) && size(V3q,2) == 3, 'V3 must be [nSrc x 3]');
    assert(sgnQ == false, 'vector fiber must yield isSigned=false');
    disp('test_atom_fiber_decode PASSED');
end

function ax = i_axes(surf, variant)
    ax = bst_eigen('Axes', struct('SurfaceFile',surf, 'Variant',variant, ...
                   'nModes',60, 'TimeWindow',[0 0.5], 'SampleRate',100));
end

function nV = i_nsrc(ax)
    nV = 0; for h=1:numel(ax.GlobalVertices)
        if ~isempty(ax.GlobalVertices{h}), nV = max(nV, max(ax.GlobalVertices{h}(:))); end
    end
end

function V3 = i_ref_quat_decode(W, gv, ax)
    % Verbatim reproduction of the former panel_bst_dynamics i_atom_realise_core quaternion branch.
    nV = i_nsrc(ax);  n = numel(gv);
    im = reshape(manifold_quat_imag(W(:,1)), 3, n).';
    V3 = zeros(nV,3);  V3(gv,:) = im;
end
