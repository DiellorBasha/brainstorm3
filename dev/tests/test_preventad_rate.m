function test_preventad_rate
% Verify the bst-meg worker brings mixed-rate noise and rest runs to one rate (preventad_rate):
%   the nine OMEGA subjects with a noise/rest rate mismatch (request a31963ec) each get a
%   target rate every recording can reach, process_resample('Compute') lands every recording
%   on it, and process_notch('Compute') accepts the notch list at that rate.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot, fileparts(thisDir));
if ~brainstorm('status')
    if usejava('awt'), brainstorm nogui; else, brainstorm server; end   % -nodisplay: no GUI at all
end

% subject, rates of its recordings (noise, rest), expected target  -- excluded_sfreq_mismatch.csv
cases = { ...
    'sub-0002',   [6000 2400],           1200;
    'sub-0010',   [6000 2400],           1200;
    'sub-0011',   [6000 2400],           1200;
    'sub-0019',   [6000 2400],           1200;
    'sub-0043',   [6000 2400],           1200;
    'sub-0047',   [6000 2400],           1200;
    'sub-0020',   [1200 2400 1200 2400], 1200;   % ses-01 and ses-02 both have 1200 Hz noise
    'sub-0062',   [2400 6000],           1200;
    'sub-PD1306', [600 2400 2400 2400 2400], 600;
    'sub-0001',   [2400 2400],           1200};  % the common case is unchanged
for i = 1:size(cases, 1)
    [fsT, notch] = preventad_rate(cases{i,2});
    assert(fsT == cases{i,3}, '%s: target %g, expected %g', cases{i,1}, fsT, cases{i,3});
    assert(all(notch < fsT/2) && all(mod(notch, 60) == 0), '%s: notch above Nyquist', cases{i,1});
    for fs = unique(cases{i,2})
        t = (0:round(2*fs)-1) / fs;                       % 2 s, 3 channels
        x = randn(3, numel(t)) + sin(2*pi*60*t);
        if abs(fs - fsT) < 0.05
            y = x;  tOut = t;                             % the worker passes it through
        else
            [y, tOut] = process_resample('Compute', x, t, fsT);
        end
        assert(abs(1/(tOut(2)-tOut(1)) - fsT) < 1e-6, '%s: %g Hz did not land on %g Hz', cases{i,1}, fs, fsT);
        [~, ~, msg] = process_notch('Compute', y, fsT, notch);
        assert(isempty(msg), '%s: notch refused at %g Hz: %s', cases{i,1}, fsT, msg);
    end
    fprintf('PASS %-11s rates %-26s -> %4d Hz, notch %s\n', cases{i,1}, mat2str(cases{i,2}), fsT, mat2str(notch));
end
assert(isequal(nargout_notch_trap(), true), 'the old fixed notch list must be refused at 600 Hz');
fprintf('PASS old fixed notch [60..300] is refused at 600 Hz (why the list now depends on the rate)\n');
end

function refused = nargout_notch_trap()
[~, ~, msg] = process_notch('Compute', randn(2, 1200), 600, [60 120 180 240 300]);
refused = ~isempty(msg);
end
