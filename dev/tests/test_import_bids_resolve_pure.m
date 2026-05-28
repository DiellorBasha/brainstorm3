function test_import_bids_resolve_pure
% Verify ResolveAnatDownsample decides (nVertices, Method, warning) per anatomy format:
%   FreeSurfer + icosphere   -> ico total count + 'icosphere', no warning
%   FreeSurfer + reducepatch -> passthrough nVertices + 'reducepatch', no warning
%   non-FreeSurfer + icosphere -> reducepatch fallback + non-empty warning
%   non-FreeSurfer + reducepatch -> passthrough, no warning
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% FreeSurfer + icosphere/ico5 -> 20484, icosphere, no warning
[nv, m, w] = process_import_bids('ResolveAnatDownsample', 'FreeSurfer', 'icosphere', 'ico5', 15000);
assert(nv == 20484 && strcmp(m, 'icosphere') && isempty(w), 'FreeSurfer icosphere ico5 must give 20484/icosphere/no-warning.');

% FreeSurfer + reducepatch -> passthrough 15000, reducepatch, no warning
[nv, m, w] = process_import_bids('ResolveAnatDownsample', 'FreeSurfer', 'reducepatch', 'ico5', 15000);
assert(nv == 15000 && strcmp(m, 'reducepatch') && isempty(w), 'FreeSurfer reducepatch must pass through 15000/reducepatch.');

% CAT12 + icosphere -> reducepatch fallback + warning (warning names the format and vertex count)
[nv, m, w] = process_import_bids('ResolveAnatDownsample', 'CAT12', 'icosphere', 'ico5', 15000);
assert(nv == 15000 && strcmp(m, 'reducepatch'), 'Non-FreeSurfer icosphere must fall back to reducepatch.');
assert(~isempty(w) && contains(w, 'CAT12') && contains(w, '15000'), 'Warning must name the format and the vertex count.');

% CAT12 + reducepatch -> passthrough, no warning
[nv, m, w] = process_import_bids('ResolveAnatDownsample', 'CAT12', 'reducepatch', 'ico5', 12000);
assert(nv == 12000 && strcmp(m, 'reducepatch') && isempty(w), 'Non-FreeSurfer reducepatch must pass through, no warning.');

fprintf('ALL TESTS PASSED: test_import_bids_resolve_pure\n');
end
