function test_eigen_metadata_routing()
% The FieldSpec-driven embedding must reproduce the old switch-Variant U_h byte-for-byte.
% Synthetic per-variant: build a fake EigenMat + source map F, compare embed to a reference.
addpath(genpath(fullfile(pwd,'toolbox')));
fprintf('--- metadata routing == old switch ---\n');
rng(3);  nGv = 40;  nT = 6;  gv = (1:nGv)';
% reference embeddings (copied from the pre-refactor switch bodies)
ref.real = @(F) F(gv,:);
ref.quat = @(F) i_embed4(F, gv, nT);              % [0;x;y;z]
for tc = {{'real',1},{'quaternion',3}}
    ft = tc{1}{1};  ncomp = tc{1}{2};
    F = randn(ncomp*nGv, nT);
    ax = struct('Variant',i_variant(ft), 'Operator',struct('Registry',struct('Primary',struct('field_type',ft,'domain','vertex'))), ...
                'Phi',{{randn(i_width(ft)*nGv,5),[]}}, 'GlobalVertices',{{gv,[]}}, 'Lambda',{{(1:5)',[]}});
    Op = struct('Mass',{{speye(i_width(ft)*nGv),[]}});
    [U_h,~,~,~,msg] = bst_eigen('ExtractHemiFieldTest', F, ax, Op, 1);  % thin test hook, see step 3
    assert(isempty(msg), 'msg: %s', msg);
    if strcmp(ft,'real'), R = ref.real(F); else, R = ref.quat(F); end
    assert(isequal(size(U_h),size(R)) && max(abs(U_h(:)-R(:)))<1e-12, '%s embed mismatch', ft);
end
fprintf('PASS\n');
end
function U=i_embed4(F,gv,nT), n=numel(gv); U=zeros(4*n,nT); U(2:4:end,:)=F((gv-1)*3+1,:); U(3:4:end,:)=F((gv-1)*3+2,:); U(4:4:end,:)=F((gv-1)*3+3,:); end
function w=i_width(ft), if strcmp(ft,'quaternion'), w=4; else, w=1; end, end
function v=i_variant(ft), if strcmp(ft,'quaternion'), v='Dirac'; else, v='Laplace-Beltrami'; end, end
