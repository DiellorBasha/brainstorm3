function C = bench_config(preset)
% BENCH_CONFIG: Configuration for the eigenmode accuracy benchmark.
% USAGE:  C = bench_config()         % full benchmark
%         C = bench_config('smoke')  % fast pipeline-validation preset
if nargin < 1 || isempty(preset); preset = 'full'; end

% Each anatomy: a protocol name + the subject whose cortex/headmodel to use.
audi = struct('label','auditory','protocol','TutorialAuditory','subject','Subject01');
neur = struct('label','neuromag','protocol','TutorialNeuromag','subject','Subject01');

C = struct();
C.methods      = {'wmne','dspm','sloreta','eig_mne_log','eig_dspm_log'};
C.eigMethods   = {'eig_mne_log','eig_dspm_log'};
C.stdMethods   = {'wmne','dspm','sloreta'};
C.regimes      = {'focal','patch','distributed'};
C.regimeOpts   = struct('focal',{{}}, 'patch',{{'Radius',2}}, 'distributed',{{'Sigma',0.01}});
C.snr_db       = [2 4 6 10 20];
C.k_total      = [600 1200 2000];     % total modes ~ {300,600,1000}/hemisphere
C.nModes_eig   = 1000;                % per hemisphere, at compute time
C.nReplicates  = 15;
C.nTime        = 20;
C.seed         = 20260602;            % master seed (date-derived, fixed for reproducibility)
C.plateauTol_mm= 1.0;                 % focal LocError improvement below this => plateau

switch lower(preset)
    case 'smoke'
        C.anatomies   = {audi};
        C.regimes     = {'focal'};
        C.regimeOpts  = struct('focal',{{}});
        C.snr_db      = [4 10];
        C.k_total     = [600 2000];
        C.nReplicates = 2;
        C.outDir      = fullfile(fileparts(mfilename('fullpath')), 'smoke_run');
    otherwise
        C.anatomies   = {audi, neur};
        C.outDir      = fullfile(fileparts(mfilename('fullpath')), 'eigenmode_accuracy_run');
end
end
