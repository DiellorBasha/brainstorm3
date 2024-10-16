% scripts/GenerateTraceDiagram.m
%
% Create a traceability diagram for a requirement that has tests and
% implementation

% Copyright 2021 The MathWorks, Inc.

% Find the requirement "Returns -1 if no path from startIdx to endIdx exists"
functionalReqs = slreq.open('shortest_path_func_reqs.slreqx');
targetReq = functionalReqs.find('Type', 'Requirement', 'Summary', ...
    'Returns -1 if no path from startIdx to endIdx exists');

% Produce the trace diagram for it
slreq.generateTraceabilityDiagram(targetReq);