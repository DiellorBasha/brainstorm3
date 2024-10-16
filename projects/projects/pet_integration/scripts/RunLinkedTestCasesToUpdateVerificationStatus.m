%  scripts/RunLinkedTestsToUpdateVerificationStatus.m
% 
% Invoke the command line API to execute M-Unit test cases linked as
% verifying requirements, that will update the verification status
% afterwards.

% Copyright 2021 The MathWorks, Inc.

% Get API objects for both requirement sets
functionalReqs = slreq.open('shortest_path_func_reqs.slreqx');
testReqs = slreq.open('shortest_path_tests_reqs.slreqx');

% Run the linked tests to update the verification status.
Hdlg = helpdlg('Running tests to update verification status.  Please wait.','Running Tests');
functionalReqs.runTests();
testReqs.runTests();
delete(Hdlg);
