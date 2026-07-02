function test_fieldspec()
% FieldSpec derives (field_type,domain,width,nComponents,C,kind) from operator metadata (falling back to
% Phi-layout inference when metadata is absent); and bst_nxr_registry('fieldspec',Variant) declares the
% field_type/domain of every operator variant. One callable entry runs both checks.
addpath(genpath(fullfile(pwd,'toolbox')));
fprintf('--- bst_eigen FieldSpec ---\n');

% fake ax with a given registry field_type/domain + a Phi of the right width
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

i_registry_fieldspec();     % every operator variant declares field_type/domain
fprintf('PASS\n');
end

% bst_nxr_registry('fieldspec', Variant) declares (field_type, domain) for every operator variant,
% including the connectome family that has no nxr binary id.
function i_registry_fieldspec()
fprintf('--- bst_nxr_registry fieldspec ---\n');
want = {'Laplace-Beltrami','real','vertex'; 'LB-Connectome','real','vertex'; ...
        'Dirac','quaternion','vertex'; 'Dirac-Connectome','quaternion','vertex'; ...
        'Dirac-Face','quaternion','face'; 'Connection Laplacian','complex','vertex'};
for i = 1:size(want,1)
    fs = bst_nxr_registry('fieldspec', want{i,1});
    assert(~isempty(fs), 'no fieldspec for %s', want{i,1});
    assert(strcmp(fs.field_type, want{i,2}) && strcmp(fs.domain, want{i,3}), 'wrong for %s', want{i,1});
end
end

function p = i_prim(ft,dom), p = struct('field_type',ft,'domain',dom); end
