%  scripts/RunTestsWithCoverage.m
% 
%  Invoke MATLAB unit tests for this project while logging code coverage
%  Generate a coverage report after completing all tests

% Copyright 2021 The MathWorks, Inc.

import matlab.unittest.TestRunner
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoberturaFormat

testFile = 'graph_unit_tests.m';
sourceCodeFile = which('shortest_path.m');
suite = testsuite(testFile);

runner = TestRunner.withTextOutput;
runner.addPlugin(CodeCoveragePlugin.forFile(sourceCodeFile))
result = runner.run(suite);


