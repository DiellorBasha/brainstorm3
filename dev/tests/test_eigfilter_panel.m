function test_eigfilter_panel()
% Headless test of the shared kernel-section helper: build sliders into a JPanel,
% read back the kernel name + params for a synthetic Lambda.
% Authors: Diellor Basha, 2026
    nFail = 0;
    Lambda = sort(rand(400,1) * 3e-5);     % synthetic Dirac eigenvalues
    [keys, displays] = bst_eigfilter_panel('Kernels');
    nFail = nFail + chk('lists curated kernels', numel(keys) >= 4 && numel(keys)==numel(displays));
    nFail = nFail + chk('heat is in the list', any(strcmp(keys,'heat')));

    % build a heat section into a river panel, read params back
    jP = gui_river([2 2], [0 2 0 2]);
    bst_eigfilter_panel('BuildSliders', jP, 'heat', Lambda, @() []);
    names = bst_eigfilter_panel('ParamNames', jP);
    nFail = nFail + chk('heat has one param (t)', numel(names)==1 && strcmp(names{1},'t'));
    p = bst_eigfilter_panel('ReadParams', jP, Lambda);
    js = jP.getClientProperty('slider_t'); k = double(js.getValue());
    nFail = nFail + chk('t = 1/lambda_k', abs(p.t - 1/max(Lambda(k),eps)) < 1e-9);

    % dog: two params, ordered t1<t2 regardless of slider positions
    bst_eigfilter_panel('BuildSliders', jP, 'dog', Lambda, @() []);
    names = bst_eigfilter_panel('ParamNames', jP);
    nFail = nFail + chk('dog has t1,t2', numel(names)==2);
    jP.getClientProperty('slider_t1').setValue(300);   % force t1 mode > t2 mode
    jP.getClientProperty('slider_t2').setValue(50);
    p = bst_eigfilter_panel('ReadParams', jP, Lambda);
    nFail = nFail + chk('dog t1 < t2 enforced', p.t1 < p.t2);

    fprintf('\n==== test_eigfilter_panel: %d failed ====\n', nFail);
    if nFail > 0, error('test_eigfilter_panel FAILED'); end
end

function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
