function test_eigfilter_design_smoke
% Smoke: the extended EigenModes panel builds, and kernel-mode state verbs drive
% weights without error. Skips if no eigenmode surface is loaded.
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3');
if ~brainstorm('status'); brainstorm nogui; end

p = panel_eigenmodes('CreatePanel');
assert(~isempty(p), 'panel must build.');

sSubjects = bst_get('ProtocolSubjects');
SurfaceFile = '';
if ~isempty(sSubjects)
    for iS = 1:numel(sSubjects.Subject)
        for iC = 1:numel(sSubjects.Subject(iS).Surface)
            sf = sSubjects.Subject(iS).Surface(iC).FileName;
            [~, ok] = in_tess_eigenmodes(sf);
            if ok; SurfaceFile = sf; break; end
        end
        if ~isempty(SurfaceFile); break; end
    end
end
if isempty(SurfaceFile); disp('SKIP: no eigenmode surface.'); return; end

[Eig,~] = in_tess_eigenmodes(SurfaceFile);
Kp = double(max(Eig.CompRank));
panel_eigenmodes('ResetState', SurfaceFile, Kp);
global GlobalData;
GlobalData.UserModes.CacheSurfaceFile = SurfaceFile;
GlobalData.UserModes.CacheEig  = Eig;
GlobalData.UserModes.CacheMass = Eig.MassMatrix;
panel_eigenmodes('SetKernelName', 'heat');
W = panel_eigenmodes('GetWeights');
assert(numel(W) == Kp && all(isfinite(W)), 'kernel weights must be finite, length Kpaired.');
panel_eigenmodes('SetKernelParam', 1, 0.2);
assert(all(isfinite(panel_eigenmodes('GetWeights'))), 'weights finite after param change.');
PairedGrid = zeros(size(Eig.Vectors,1), Kp);
for k=1:Kp; PairedGrid(:,k)=sum(Eig.Vectors(:,Eig.CompRank==k),2); end
panel_eigenmodes('SetDisplayMode','delta');
panel_eigenmodes('SetDeltaVertex', 1);
col = panel_eigenmodes('GetDisplayColumn', SurfaceFile, PairedGrid);
assert(isequal(size(col),[size(Eig.Vectors,1) 1]) && all(isfinite(col)), 'delta column ok.');

disp('ALL TESTS PASSED');
end
