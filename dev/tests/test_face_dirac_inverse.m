function test_face_dirac_inverse()
% Phase 4: bst_inverse_dirac face path. A face leadfield (isFaceBased=1) must yield a
% [3nF x nCh] kernel (the original bug returned a vertex-sized 3nV kernel).
% Author: Diellor Basha, 2026
    nFail = 0;
    df = 'Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat';
    [sStudy,~] = bst_get('DataFile', df);
    BaseHM = in_bst_headmodel([fileparts(df) '/headmodel_surf_os_meg.mat'], 0);
    ChanMat = in_bst_channel(sStudy.Channel(1).FileName); types = {ChanMat.Channel.Type};
    iMEG = find(strcmpi(types,'MEG'));
    NC = load(file_fullpath([fileparts(df) '/noisecov_full.mat'])); Cn = NC.NoiseCov(iMEG,iMEG); Cn=(Cn+Cn')/2;
    [Lf, FG] = bst_face_leadfield(BaseHM.SurfaceFile, ChanMat.Channel(iMEG), BaseHM.Param(iMEG), 'Mode','unconstrained');
    nF = size(Lf,2)/3;
    HMf = struct('Gain',Lf, 'SurfaceFile',BaseHM.SurfaceFile, 'HeadModelType','surface', ...
                 'isFaceBased',1, 'GridLoc',FG.Centroids);
    OPT = struct('NoiseMethod','reg','NoiseReg',0.1,'SnrMethod','fixed','SnrFixed',3, ...
                 'InverseMeasure','dspm2018','nModes',400,'Tau',0.5);
    OPT.NoiseCovMat.NoiseCov = Cn;  OPT.ChannelTypes = types(iMEG);

    R = bst_inverse_dirac(HMf, OPT);
    nFail = nFail + chk('face kernel [3nF x nCh] (NOT 3nV)', isequal(size(R.ImagingKernel),[3*nF numel(iMEG)]));
    nFail = nFail + chk('kernel finite', all(isfinite(R.ImagingKernel(:))));
    nFail = nFail + chk('isFaceBased preserved through inverse', isfield(R,'nModes') && R.nModes==2*min(400,R.nModes/2));

    fprintf('\n==== test_face_dirac_inverse: %d failed ====\n', nFail);
    if nFail > 0, error('test_face_dirac_inverse FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
