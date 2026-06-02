function test_harmonic_inverse_e2e
% Smoke: compute a Harmonic inverse and verify its consistency contract:
%   - bst_inverse_eigenmodes('harmonic') returns a finite [K x nGoodCh] kernel
%   - the vertex map is Phi * EigenKernel
% Requires a protocol with a surface head model + eigenmodes + noise cov.
% Skips cleanly otherwise.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sProtocol = bst_get('ProtocolStudies');
if isempty(sProtocol) || ~isfield(sProtocol, 'Study') || isempty(sProtocol.Study)
    disp('SKIP: no protocol loaded.');
    return;
end
% Find a study with surface head model + eigenmodes + noise cov
target = [];
for iS = 1:numel(sProtocol.Study)
    s = sProtocol.Study(iS);
    if isfield(s,'iHeadModel') && ~isempty(s.iHeadModel) && (s.iHeadModel >= 1) ...
            && (length(s.HeadModel) >= s.iHeadModel) ...
            && isfield(s,'NoiseCov') && ~isempty(s.NoiseCov) && ~isempty(s.NoiseCov(1).FileName)
        try
            hm = in_bst_headmodel(s.HeadModel(s.iHeadModel).FileName, 0, 'HeadModelType', 'SurfaceFile');
            if strcmpi(hm.HeadModelType, 'surface')
                [~, isEig] = in_tess_eigenmodes(hm.SurfaceFile);
                if isEig
                    target = struct('surf', hm.SurfaceFile, ...
                        'hmFile', s.HeadModel(s.iHeadModel).FileName, ...
                        'ncFile', s.NoiseCov(1).FileName);
                    break;
                end
            end
        catch
        end
    end
end
if isempty(target)
    disp('SKIP: no study with surface head model + eigenmodes + noise cov.');
    return;
end

% Build a good-channel mask: drop non-sensor / NaN-gain rows so the smoke test
% mirrors the realistic compute paths (which pass good_channel-filtered sets).
HM = in_bst_headmodel(target.hmFile, 1);   % ApplyOrient=1 -> constrained [nch x nVert]
goodMask = all(isfinite(double(HM.Gain)), 2);
assert(any(goodMask), 'No finite-gain channels found.');

% Compute the harmonic kernel via bst_inverse_eigenmodes (engine under both compute surfaces)
[InvE, errE] = bst_inverse_eigenmodes(target.hmFile, target.surf, target.ncFile, ...
    'Method', 'harmonic', 'nModes', 0, 'GoodChannel', goodMask);
assert(isempty(errE), ['bst_inverse_eigenmodes harmonic failed: ' errE]);
assert(~isempty(InvE.ImagingKernel), 'Harmonic kernel must be non-empty.');
K = InvE.nModes;
assert(size(InvE.ImagingKernel, 1) == K, 'EigenKernel must have K rows.');
assert(all(isfinite(InvE.ImagingKernel(:))), 'Harmonic kernel must be finite (rank-safe).');

% Vertex map = Phi * M̃
[Eig, ~] = in_tess_eigenmodes(target.surf);
Phi  = double(Eig.Vectors(:, 1:K));
Vmap = Phi * InvE.ImagingKernel;
assert(size(Vmap, 1) == size(Eig.Vectors, 1), 'Vertex kernel must have nVert rows.');
assert(size(Vmap, 2) == size(InvE.ImagingKernel, 2), 'Vertex kernel column count must match channels.');

disp('ALL TESTS PASSED');
end
