%  scripts/GenerateReqsReport.m
% 
% Invoke the command line API to update verification status that will 
% cause the M-Unit test cases linked as verifying requirements to be 
% executed.

% Copyright 2021 The MathWorks, Inc.

% Get API objects for both requirement sets
functionalReqs = slreq.open('shortest_path_func_reqs.slreqx');
testReqs = slreq.open('shortest_path_tests_reqs.slreqx');

% Report options
reqReportOptions = slreq.getReportOptions;
reqReportOptions.includes.toc = true;
reqReportOptions.includes.publishedDate = true;
reqReportOptions.includes.revision = true;
reqReportOptions.includes.properties = false;
reqReportOptions.includes.links = true;
reqReportOptions.includes.changeInformation = false;
reqReportOptions.includes.groupLinksBy = 'Artifact';
reqReportOptions.includes.keywords = false;
reqReportOptions.includes.comments = false;
reqReportOptions.includes.implementationStatus = true;
reqReportOptions.includes.verificationStatus = true;
reqReportOptions.includes.emptySections = false;
reqReportOptions.includes.rationale = false;
reqReportOptions.includes.customAttributes = false;


% Generate the report
reportPath = slreq.generateReport([functionalReqs testReqs], reqReportOptions);

% Open the report
uiopen(reportPath,1);