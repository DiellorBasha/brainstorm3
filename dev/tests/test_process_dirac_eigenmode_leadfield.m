function test_process_dirac_eigenmode_leadfield
% Static dispatch sanity (no DB): GetDescription + FormatComment.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);

sProcess = process_dirac_eigenmode_leadfield('GetDescription');
assert(strcmp(sProcess.Comment, 'Compute Dirac eigenmode leadfield'), 'Comment');
assert(isfield(sProcess.options,'nmodes') && isfield(sProcess.options,'tau'), 'options nmodes+tau');
assert(strcmpi(sProcess.SubGroup, 'Sources'), 'SubGroup');

sProcess.options.nmodes.Value = {0, '', 0};
c0 = process_dirac_eigenmode_leadfield('FormatComment', sProcess);
assert(contains(lower(c0), 'default'), 'FormatComment default');

sProcess.options.nmodes.Value = {50, '', 0};
c1 = process_dirac_eigenmode_leadfield('FormatComment', sProcess);
assert(contains(c1, '50'), 'FormatComment count');

disp('ALL TESTS PASSED');
end
