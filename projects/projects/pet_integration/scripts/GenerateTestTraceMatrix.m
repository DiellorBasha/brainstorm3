% scripts/GenerateTestTraceMatrix.m
%
% Create a traceability matrix that shows how unit tests in graph_unit_tests.m
% link to functional equirements shortest_path_func_reqs.slreqx and testing
% requirements shortest_path_tests_reqs.slreqx

% Copyright 2021 The MathWorks, Inc.

mtxOpts = slreq.getTraceabilityMatrixOptions;
mtxOpts.leftArtifacts = {which('graph_unit_tests')};
mtxOpts.topArtifacts = {which('shortest_path_func_reqs.slreqx'), ...
                        which('shortest_path_tests_reqs.slreqx')};

slreq.generateTraceabilityMatrix(mtxOpts);