function R = tutorial_benchmark_esi(varargin)
% TUTORIAL_BENCHMARK_ESI: ESI benchmark driver - simulate, reconstruct, score, report.
%
% USAGE:  R = tutorial_benchmark_esi('Regimes',{...},'SNRs',[...],'nLoc',N,'nNoise',M,'OutDir',dir)
%
% Sweeps regimes x SNRs x source locations x noise draws over an OMEGA MEG study,
% runs the comparator panel, scores each estimate vs ground truth, and aggregates
% into descriptive statistics. Returns [] (skips) if no suitable study is found.
%
% Authors: Diellor Basha, 2026
Regimes = {'focal','patch','distributed'}; SNRs = [0 3 6 10];
nLoc = 20; nNoise = 5;
OutDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), '..', 'dev','benchmarks','run');
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'regimes', Regimes = varargin{i+1};
        case 'snrs',    SNRs    = varargin{i+1};
        case 'nloc',    nLoc    = varargin{i+1};
        case 'nnoise',  nNoise  = varargin{i+1};
        case 'outdir',  OutDir  = varargin{i+1};
    end
end
if ~brainstorm('status'); brainstorm nogui; end

target = find_benchmark_study();
if isempty(target); R = []; return; end

baseHM = in_bst_headmodel(target.baseHmFile, 1);    % constrained [nCh x nVert]
goodMask = all(isfinite(double(baseHM.Gain)), 2);
L = double(baseHM.Gain(goodMask, :));
Surf = in_tess_bst(baseHM.SurfaceFile);
if ~isfield(Surf,'VertConn') || isempty(Surf.VertConn)
    Surf.VertConn = tess_vertconn(Surf.Vertices, Surf.Faces);
end
GridLoc = Surf.Vertices;                             % per-vertex positions (surface model)
NC = load(file_fullpath(target.ncFile)); C = NC.NoiseCov(goodMask, goodMask);
SurfStruct = struct('Vertices', Surf.Vertices, 'VertConn', Surf.VertConn);

rows = struct('regime',{},'snr',{},'method',{},'metric',{},'value',{},'realization',{});
real = 0;
for ir = 1:numel(Regimes)
  for il = 1:nLoc
    S = bst_benchmark_sources(SurfStruct, Regimes{ir}, 'Seed', il);
    for is = 1:numel(SNRs)
      for inum = 1:nNoise
        real = real + 1;
        Sim = bst_benchmark_simulate(L, S.Sources, C, 'SNR', SNRs(is), 'Seed', 1000*il+inum);
        Est = bst_benchmark_inverse(Sim.F, target.baseHmFile, target.ncFile, target.chFile, goodMask, SNRs(is));
        tEval = round(size(S.Sources,2)/2);
        fn = fieldnames(Est);
        for k = 1:numel(fn)
            M = bst_benchmark_metrics(S.GT, Est.(fn{k})(:,tEval), GridLoc, S.SeedVertex);
            mn = {'LocError','AUC','NRMSE','Correlation','SpatialDispersion'};
            for q = 1:numel(mn)
                rows(end+1) = struct('regime',Regimes{ir},'snr',SNRs(is),'method',fn{k}, ...
                    'metric',mn{q},'value',M.(mn{q}),'realization',real); %#ok<AGROW>
            end
        end
      end
    end
  end
end

R = bst_benchmark_report(rows, 'RefMethod','eig_mne_log', 'Seed',1, 'OutDir',OutDir);
fprintf('ESI benchmark complete: %d rows -> %s\n', numel(rows), OutDir);
end

function target = find_benchmark_study()
target = [];
sStudies = bst_get('ProtocolStudies');
if isempty(sStudies) || ~isfield(sStudies,'Study'); return; end
for iS = 1:numel(sStudies.Study)
    s = sStudies.Study(iS);
    if isempty(s.iHeadModel) || numel(s.HeadModel)<1; continue; end
    if isempty(s.NoiseCov) || isempty(s.NoiseCov(1).FileName); continue; end
    if isempty(s.Channel) || ~any(strcmpi(s.Channel.DisplayableSensorTypes,'MEG')); continue; end
    iBase = [];
    for ih = 1:numel(s.HeadModel)
        try hm = in_bst_headmodel(s.HeadModel(ih).FileName,0); catch; continue; end
        if (~isfield(hm,'isEigenmode')||~hm.isEigenmode) && strcmpi(hm.HeadModelType,'surface'); iBase=ih; break; end
    end
    if isempty(iBase); continue; end
    hmB = in_bst_headmodel(s.HeadModel(iBase).FileName,0,'SurfaceFile');
    [~, isE] = in_tess_eigenmodes(hmB.SurfaceFile);
    if ~isE; continue; end
    target = struct('baseHmFile', s.HeadModel(iBase).FileName, ...
        'ncFile', s.NoiseCov(1).FileName, ...
        'chFile', bst_get('ChannelFileForStudy', s.FileName));
    return;
end
end
