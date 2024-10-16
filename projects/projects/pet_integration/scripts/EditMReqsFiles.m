%  scripts/EditMReqsFiles.m
% 
%  Open the requirement files and edit the test and design M files.  

% Copyright 2021 The MathWorks, Inc.

% Load and open requirements
slreq.open('shortest_path_func_reqs.slreqx');
slreq.open('shortest_path_tests_reqs.slreqx');

% Edit the design and test MATLAB code
edit shortest_path.m;
edit graph_unit_tests.m;




