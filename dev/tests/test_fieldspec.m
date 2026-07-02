function test_fieldspec()
% FieldSpec derives (field_type,domain,width,nComponents,C,kind) from operator metadata,
% and falls back to Phi-layout inference when metadata is absent.
addpath(genpath(fullfile(pwd,'toolbox')));
fprintf('--- bst_eigen FieldSpec ---\n');

% helper: fake ax with a given registry field_type/domain + a Phi of the right width
mk = @(ft, dom, width, nGv) struct( ...
    'Operator', struct('Registry', struct('Primary', i_prim(ft,dom))), ...
    'Phi', {{randn(width*nGv, 5), []}}, 'GlobalVertices', {{(1:nGv)', []}});

% metadata-driven cases
s = bst_eigen('FieldSpec', mk('real','vertex',1,50));
assert(strcmp(s.field_type,'real') && strcmp(s.domain,'vertex') && s.width==1 && s.nComponents==1 && s.C==1 && strcmp(s.kind,'scalar'));
s = bst_eigen('FieldSpec', mk('quaternion','vertex',4,50));
assert(strcmp(s.field_type,'quaternion') && s.width==4 && s.nComponents==3 && s.C==4 && strcmp(s.kind,'quaternion'));
s = bst_eigen('FieldSpec', mk('complex','vertex',1,50));
assert(strcmp(s.field_type,'complex') && s.C==2 && s.nComponents==1 && strcmp(s.kind,'tangent'));

% domain default = vertex when metadata omits it
p = i_prim('quaternion',''); p = rmfield(p,'domain');
axNoDom = struct('Operator',struct('Registry',struct('Primary',p)), 'Phi',{{randn(4*50,5),[]}}, 'GlobalVertices',{{(1:50)',[]}});
s = bst_eigen('FieldSpec', axNoDom);  assert(strcmp(s.domain,'vertex'));

% pre-registry fallback: no Registry -> infer C from Phi rows / nV
axFb = struct('Phi',{{randn(4*50,5),[]}}, 'GlobalVertices',{{(1:50)',[]}});
s = bst_eigen('FieldSpec', axFb);  assert(s.C==4 && strcmp(s.kind,'quaternion'));
fprintf('PASS\n');
end
function p = i_prim(ft,dom), p = struct('field_type',ft,'domain',dom); end
