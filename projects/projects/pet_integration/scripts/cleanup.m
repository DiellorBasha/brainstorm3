%  scripts/cleanup.m
% 
%  Cleanup actions for ShortestPath project.  

% Copyright 2021 The MathWorks, Inc.

% Close functional requirements
functionalReqs = slreq.find('Type','ReqSet','Name','shortest_path_func_reqs');
if ~isempty(functionalReqs)
    functionalReqs.close();
end

% Close test requirements
testReqs = slreq.find('Type','ReqSet','Name','shortest_path_tests_reqs');
if ~isempty(testReqs)
    testReqs.close();
end


