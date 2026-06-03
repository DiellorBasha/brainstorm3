function test_benchmark_inverse_nmodes_e2e
% Verify the optional nModes arg is accepted (test that the function signature changed).
% This is a signature test - it confirms the 7-arg variant is callable.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);

% Call with 7 arguments; before the change, this should fail with "Too many input arguments".
% After the change, this will hit protocol-not-loaded and skip, but the signature will be accepted.
try
    % Dummy data - just to test that 7 args is accepted
    Est = bst_benchmark_inverse([], '', '', '', [], 6, 600);
    % If we got here, the signature accepts 7 args - good!
    disp('PASS: bst_benchmark_inverse accepts nModes arg (7 args).');
catch ME
    if contains(ME.message, 'Too many input arguments')
        disp('FAIL: bst_benchmark_inverse still only accepts 6 args.');
        rethrow(ME);
    else
        % Other error (expected, since we passed empty data) - still means signature is OK
        disp('PASS: bst_benchmark_inverse accepts nModes arg (7 args) - called with dummy data.');
    end
end
end
