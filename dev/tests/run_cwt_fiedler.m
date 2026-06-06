% RUN_CWT_FIEDLER  Run the full CWT/Fiedler pipeline and save a compact result.
if ~brainstorm('status'), brainstorm nogui; end
HM='Subject01/S01_AEF_20131218_01_notch/headmodel_face_eigenmode_260606_1131.mat';
DATA='Subject01/S01_AEF_20131218_02_notch/data_block001_band.mat';
SURF='Subject01/tess_cortex_pial_low.mat';
outFile=fullfile(fileparts(mfilename('fullpath')),'cwt_fiedler_result.mat');

t0=tic;
R=bst_cwt_fiedler_pipeline(HM,DATA,SURF,'TargetFreq',10,'FreqRange',[5 20],...
    'VoicesPerOctave',8,'SNR',3,'Verbose',true);
fprintf('TOTAL: %.1fs\n',toc(t0));

% Save only the compact summary (not the giant arrays)
Summary.DispInfo  = R.DispInfo;
Summary.peak_info = R.peak_info;
Summary.Wave      = R.Wave;
Summary.S_joint   = R.S_joint;
Summary.f_ax      = R.f_ax;
Summary.k_ax      = R.k_ax;
Summary.phi_F     = R.phi_F;
save(outFile,'Summary','-v7');
fprintf('Saved compact summary → %s\n',outFile);

fprintf('\n=== CWT/FIEDLER SUMMARY ===\n');
fprintf('Dispersion: v=%.2f m/s  R²=%.3f\n', R.DispInfo.v_linear, R.DispInfo.R2);
fprintf('Source: face %d  [%.1f %.1f %.1f] mm\n', R.Wave.SourceFace, R.Wave.SourcePos*1000);
fprintf('Fiedler wave speed: %.2f ± %.2f m/s\n', R.Wave.v_fiedler, R.Wave.v_fiedler_std);
fprintf('PLV: max=%.3f  mean=%.3f\n', max(R.Wave.PLV,[],'omitnan'), mean(R.Wave.PLV,'omitnan'));
